"""Fresh-batch adaptive amplitude intervals with an anytime union-bound budget.

The monotonic-branch principle follows Grinko et al. (2021). This implementation
uses independent fresh batches and delta/(r*(r+1)), not a fixed budget reused
after adaptive looks. The scheduler never receives the unknown amplitude.
"""
import math
import numpy as np
from scipy.stats import beta


def cp_interval(successes, shots, delta):
    successes = np.asarray(successes)
    shots = np.asarray(shots)
    lower = np.where(successes == 0, 0., beta.ppf(delta / 2, successes, shots - successes + 1))
    upper = np.where(successes == shots, 1., beta.ppf(1 - delta / 2, successes + 1, shots - successes))
    return lower, upper


def choose_power(lo, hi, old_m=1, maximum_power=1048576):
    """Largest admissible odd multiplier; theta is measured in radians."""
    maximum = min(2 * maximum_power + 1, max(1, int(math.pi / (2 * (hi - lo)))))
    if maximum % 2 == 0:
        maximum -= 1
    for m in range(maximum, max(1, old_m) - 1, -2):
        branch_l = math.floor(2 * m * lo / math.pi + 1e-13)
        branch_u = math.floor(2 * m * hi / math.pi - 1e-13)
        if branch_l == branch_u:
            return m, branch_l
    # The previous interval was admitted on this multiplier, including endpoints.
    branch = math.floor(2 * old_m * lo / math.pi + 1e-13)
    assert 2 * old_m * hi / math.pi <= branch + 1 + 1e-11
    return old_m, branch


def estimate_batch(sampler, number, epsilon, delta_output, maximum_rounds=64):
    """sampler(indices, multipliers, shots) returns fresh Bernoulli successes."""
    lo = np.zeros(number)
    hi = np.full(number, math.pi / 2)
    old = np.ones(number, dtype=np.int64)
    repeated = np.full(number, -1, dtype=np.int64)
    status = np.zeros(number, dtype=np.int8)  # 0 active, 1 complete, -1 contradiction.
    A_queries = np.zeros(number, dtype=np.int64)
    Q_queries = np.zeros(number, dtype=np.int64)
    readouts = np.zeros(number, dtype=np.int64)
    maximum_k = np.zeros(number, dtype=np.int64)
    rounds = []
    for r in range(1, maximum_rounds + 1):
        indices = np.flatnonzero(status == 0)
        if not len(indices):
            break
        picked = [choose_power(lo[i], hi[i], old[i]) for i in indices]
        m = np.array([x[0] for x in picked], dtype=np.int64)
        branch = np.array([x[1] for x in picked], dtype=np.int64)
        same = m == old[indices]
        repeated[indices] = np.where(same, repeated[indices] + 1, 0)
        # First round uses256; repeated fresh rounds double, bounded at65536.
        exponent = np.minimum(np.maximum(repeated[indices], 0), 8)
        shots = 256 * np.power(2, exponent)
        successes = sampler(indices, m, shots)
        delta_round = delta_output / (r * (r + 1))
        p_lo, p_hi = cp_interval(successes, shots, delta_round)
        u = np.arcsin(np.sqrt(p_lo))
        v = np.arcsin(np.sqrt(p_hi))
        even = branch % 2 == 0
        next_lo = np.where(even, branch * math.pi / 2 + u, (branch + 1) * math.pi / 2 - v) / m
        next_hi = np.where(even, branch * math.pi / 2 + v, (branch + 1) * math.pi / 2 - u) / m
        next_lo = np.maximum(lo[indices], next_lo)
        next_hi = np.minimum(hi[indices], next_hi)
        bad = next_lo > next_hi
        status[indices[bad]] = -1
        lo[indices] = next_lo
        hi[indices] = next_hi
        a_lo = np.sin(next_lo) ** 2
        a_hi = np.sin(next_hi) ** 2
        done = (~bad) & ((a_hi - a_lo) <= 2 * epsilon)
        status[indices[done]] = 1
        k = (m - 1) // 2
        A_queries[indices] += shots * m
        Q_queries[indices] += shots * k
        readouts[indices] += shots
        maximum_k[indices] = np.maximum(maximum_k[indices], k)
        old[indices] = m
        rounds.append(dict(round=r,indices=indices,multipliers=m,branches=branch,shots=shots,successes=successes,
                           theta_lower=next_lo.copy(),theta_upper=next_hi.copy(),delta=delta_round))
    return dict(status=status,a_lower=np.sin(lo)**2,a_upper=np.sin(hi)**2,
                A_queries=A_queries,Q_queries=Q_queries,readouts=readouts,maximum_k=maximum_k,rounds=rounds)


def simulated_estimates(amplitudes, epsilon, delta_output, seed):
    """Ideal classical likelihood simulation; does not execute a quantum circuit."""
    theta = np.arcsin(np.sqrt(np.clip(amplitudes, 0, 1)))
    rng = np.random.default_rng(seed)
    def sampler(indices, multipliers, shots):
        return rng.binomial(shots, np.sin(multipliers * theta[indices]) ** 2)
    return estimate_batch(sampler, len(theta), epsilon, delta_output)
