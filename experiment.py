"""Exact-current-law finite-time verification; synthetic, not biological evidence.

No fitting, target labels, molecular trajectories or quantum-device execution.
Five contact-network nodes; two fixed shared harmonic displacement coordinates.
This script records every setting and measured result without selecting for effect.
"""
import os
from contractions import contract
from runtime_info import usage
from pathlib import Path
import csv, hashlib, json, math, platform, sys, time
import numpy as np
import scipy
from scipy import sparse
from scipy.special import expit, ive, logsumexp
from scipy.sparse.linalg import expm_multiply
ROOT = Path(os.environ.get('PULSAR_OUTPUT_DIR', Path(__file__).resolve().parent / 'results'))
R0 = np.array([[0, 0, 0], [1.1, 0, 0], [0.45, 0.95, 0], [0.35, 0.3, 0.95], [0.65, 0.28, -0.85]], float)
EDGES = [(i, j) for i in range(5) for j in range(i + 1, 5) if (i, j) != (3, 4)]
STIFFNESS = 10.0
BETA = 1.0
MOBILITY = 1.0
D = MOBILITY / BETA
SENDER = 3
RECEIVER = 4
POWERS = [(a, b) for degree in (2, 3, 4) for a in range(degree + 1) for b in [degree - a]]

def primitive(x):
    if isinstance(x, np.ndarray):
        return x.tolist()
    if isinstance(x, np.generic):
        return x.item()
    if isinstance(x, dict):
        return {str(k): primitive(v) for (k, v) in x.items()}
    if isinstance(x, (list, tuple)):
        return [primitive(v) for v in x]
    return x

def maxabs(x):
    return float(np.max(np.abs(x)))

def dense_taylor_expm(A):
    """Independent small-matrix scaling/squaring Taylor reference.

    Scale the infinity norm to at most 0.5, take 40 Taylor terms, then square.
    This avoids the dense library matmul path that emitted local runtime warnings.
    Roundoff is still assessed by comparison; no interval certificate is claimed.
    """
    A = np.asarray(A, float)
    bound = np.linalg.norm(A, np.inf)
    scale = max(0, int(math.ceil(math.log2(bound / 0.5)))) if bound else 0
    scaled = A / 2 ** scale
    result = np.eye(len(A))
    term = result.copy()
    for k in range(1, 41):
        term = np.einsum('ij,jk->ik', term, scaled, optimize=False) / k
        result += term
    for _ in range(scale):
        result = np.einsum('ij,jk->ik', result, result, optimize=False)
    return result

def geometry():
    K = np.zeros((15, 15))
    for (i, j) in EDGES:
        a = R0[i] - R0[j]
        v = np.zeros(15)
        v[3 * i:3 * i + 3] = a / np.linalg.norm(a)
        v[3 * j:3 * j + 3] = -a / np.linalg.norm(a)
        K += STIFFNESS * np.outer(v, v)
    (vals, vecs) = np.linalg.eigh(K)
    positive = np.flatnonzero(vals > 1e-09)
    assert len(positive) == 9 and np.sum(np.abs(vals) < 1e-09) == 6
    B = vecs[:, positive[:2]].copy()
    for col in range(2):
        if B[np.argmax(np.abs(B[:, col])), col] < 0:
            B[:, col] *= -1
    lam = vals[positive[:2]]
    rigid = []
    centered = R0 - R0.mean(axis=0)
    for axis in np.eye(3):
        rigid.extend([np.tile(axis, (5, 1)).ravel(), np.cross(np.tile(axis, (5, 1)), centered).ravel()])
    rigid = np.array(rigid).T
    return (K, B, lam, {'hessian': K, 'eigenvalues': vals, 'positive_mode_count': len(positive), 'rigid_null_residual': maxabs(contract(K, rigid)), 'B_orthonormal_residual': maxabs(contract(B.T, B) - np.eye(2)), 'B_rigid_overlap': maxabs(contract(B.T, rigid)), 'retained_eigenvalues': lam, 'B': B})
(K, B, LAM, GEO) = geometry()
TAU = 1 / (MOBILITY * LAM[0])
DEG = np.bincount(np.array(EDGES).ravel(), minlength=5)

def observables(q):
    r = R0[None, :, :] + contract(q, B.T).reshape(-1, 5, 3)
    E = np.zeros((len(q), 5))
    U = np.zeros(len(q))
    stretch = []
    for (i, j) in EDGES:
        l2 = np.sum((R0[i] - R0[j]) ** 2)
        current2 = np.sum((r[:, i] - r[:, j]) ** 2, axis=1)
        u = STIFFNESS / (8 * l2) * (current2 - l2) ** 2
        U += u
        E[:, i] += u / DEG[i]
        E[:, j] += u / DEG[j]
        stretch.append(np.sqrt(current2 / l2))
    return (U, E, r, np.array(stretch).T)

def polynomial_coefficients():
    A = np.zeros((5, len(POWERS)))
    for (i, j) in EDGES:
        a = R0[i] - R0[j]
        T = B[3 * i:3 * i + 3] - B[3 * j:3 * j + 3]
        l2 = contract(a, a)
        v = contract(2 * a, T)
        Q = contract(T.T, T)
        g = {(1, 0): v[0], (0, 1): v[1], (2, 0): Q[0, 0], (1, 1): 2 * Q[0, 1], (0, 2): Q[1, 1]}
        co = {power: 0.0 for power in POWERS}
        for (p, c) in g.items():
            for (z, w) in g.items():
                co[p[0] + z[0], p[1] + z[1]] += c * w * STIFFNESS / (8 * l2)
        for (k, p) in enumerate(POWERS):
            A[i, k] += co[p] / DEG[i]
            A[j, k] += co[p] / DEG[j]
    return A
A = polynomial_coefficients()

def monomials(q):
    return np.array([q[:, 0] ** a * q[:, 1] ** b for (a, b) in POWERS]).T

def gaussian_moment(power):
    val = 1.0
    for (k, lam) in zip(power, LAM):
        if k % 2:
            return 0.0
        val *= math.prod(range(1, k, 2)) * (1 / (BETA * lam)) ** (k / 2)
    return val
MUZ = np.array([gaussian_moment(p) for p in POWERS])
GAUSS = np.array([[gaussian_moment((p[0] + q[0], p[1] + q[1])) for q in POWERS] for p in POWERS]) - np.outer(MUZ, MUZ)
HARM_COV = contract(contract(A, GAUSS), A.T)
SH = np.sqrt(np.diag(HARM_COV))
NORMAL = np.outer(SH, SH)

def lattice(n, extent):
    half = extent / np.sqrt(BETA * LAM[0])
    axis = np.linspace(-half, half, n)
    delta = axis[1] - axis[0]
    (xx, yy) = np.meshgrid(axis, axis, indexing='ij')
    q = np.column_stack([xx.ravel(), yy.ravel()])
    ids = np.arange(n * n).reshape(n, n)
    x = np.concatenate([ids[:-1, :].ravel(), ids[:, :-1].ravel()])
    y = np.concatenate([ids[1:, :].ravel(), ids[:, 1:].ravel()])
    return (q, delta, x, y, half)

def generator(U, delta, x, y):
    G = len(U)
    scale = 2 * D / delta ** 2
    qxy = scale * expit(-BETA * (U[y] - U[x]))
    qyx = scale * expit(BETA * (U[y] - U[x]))
    exitrate = np.bincount(x, weights=qxy, minlength=G) + np.bincount(y, weights=qyx, minlength=G)
    row = np.concatenate([x, y, np.arange(G)])
    col = np.concatenate([y, x, np.arange(G)])
    L = sparse.csr_matrix((np.concatenate([qxy, qyx, -exitrate]), (row, col)), shape=(G, G))
    off = -np.sqrt(qxy * qyx)
    H = sparse.csr_matrix((np.concatenate([off, off, exitrate]), (row, col)), shape=(G, G))
    pi = np.exp(-BETA * U - logsumexp(-BETA * U))
    checks = {'generator_conservation': maxabs(contract(L, np.ones(G))), 'stationarity': maxabs(contract(pi, L)), 'detailed_balance': maxabs(pi[x] * qxy - pi[y] * qyx), 'H_sqrtpi_null': maxabs(contract(H, np.sqrt(pi))), 'H_symmetry': maxabs((H - H.T).data) if (H - H.T).nnz else 0.0}
    return (L, H, pi, checks)

def run_case(n, extent, model, finite_fields=False, polycheck=False, densecheck=False):
    start = time.perf_counter()
    (q, delta, x, y, half) = lattice(n, extent)
    (UA, E, r, stretch) = observables(q)
    U = UA if model == 'biquadratic' else 0.5 * np.sum(q * q * LAM, axis=1)
    (L, H, pi, checks) = generator(U, delta, x, y)
    means = contract(pi, E)
    centered = E - means
    e = np.sqrt(pi)[:, None] * centered
    evolved = expm_multiply(-TAU * H, e, traceA=-TAU * H.diagonal().sum())
    cov = contract(e.T, e)
    delayed = contract(e.T, evolved)
    R = -BETA * (cov - delayed)
    C = R / (BETA * NORMAL)
    Z = monomials(q)
    ez = np.sqrt(pi)[:, None] * (Z - contract(pi, Z))
    zevolved = expm_multiply(-TAU * H, ez, traceA=-TAU * H.diagonal().sum())
    Rrec = contract(contract(-BETA * A, contract(ez.T, ez - zevolved)), A.T)
    checks.update({'fixed_observable_polynomial_max_error': maxabs(E - contract(Z, A.T)), 'reconstruction_R_max_error': maxabs(R - Rrec), 'response_symmetry': maxabs(R - R.T), 'largest_R_eigenvalue': float(np.linalg.eigvalsh((R + R.T) / 2).max()), 'maximum_R_diagonal': float(np.diag(R).max()), 'R_at_zero_max_error': maxabs(-BETA * (cov - contract(e.T, expm_multiply(0.0 * H, e, traceA=0.0)))), 'harmonic_normalization_covariance_max_grid_error': maxabs(cov - HARM_COV) if model == 'harmonic' else None})
    idx = np.arange(n * n).reshape(n, n)
    boundary = np.unique(np.concatenate([idx[0], idx[-1], idx[:, 0], idx[:, -1]]))
    apices = np.linalg.norm(r[:, 3] - r[:, 4], axis=1)
    distort = np.max(np.abs(stretch - 1), axis=1)
    result = {'model': model, 'n': n, 'G': n * n, 'extent_slowest_harmonic_sigma': extent, 'half_width': half, 'spacing': delta, 'horizon': TAU, 'U_min': float(U.min()), 'U_max': float(U.max()), 'means': means, 'variances': np.diag(cov), 'R': R, 'C': C, 'apex_response_R': float(R[RECEIVER, SENDER]), 'apex_response_C': float(C[RECEIVER, SENDER]), 'checks': checks, 'validity_diagnostics': {'outer_grid_shell_probability': float(pi[boundary].sum()), 'probability_any_contact_relative_length_change_gt_0p2': float(pi[distort > 0.2].sum()), 'probability_any_contact_below_sqrt_half_reference': float(pi[np.min(stretch, axis=1) < np.sqrt(0.5)].sum()), 'probability_noncontact_apices_distance_lt_0p5': float(pi[apices < 0.5].sum()), 'mean_maximum_absolute_contact_strain': float(contract(pi, distort))}, 'finite_field_checks': []}
    if finite_fields:
        for h in [0.01, 0.001, 0.0001, 1e-05]:
            (Lh, _, _, _) = generator(U + h * E[:, SENDER], delta, x, y)
            fh = expm_multiply(TAU * Lh, E[:, RECEIVER], traceA=TAU * Lh.diagonal().sum())
            slope = (contract(pi, fh) - means[RECEIVER]) / h
            result['finite_field_checks'].append({'h': h, 'forward_slope': float(slope), 'difference_from_covariance_response': float(slope - R[RECEIVER, SENDER])})
    if densecheck:
        hd = H.toarray()
        gd = L.toarray()
        rd = -BETA * (cov - contract(contract(e.T, dense_taylor_expm(-TAU * hd)), e))
        direct = -BETA * (cov - contract(centered.T, pi[:, None] * contract(dense_taylor_expm(TAU * gd), centered)))
        checks['dense_expm_R_max_difference'] = maxabs(rd - R)
        checks['backward_generator_expm_R_max_difference'] = maxabs(direct - R)
        eps = 1e-05
        Hess = np.zeros((2, 2))
        zero = np.zeros((1, 2))
        f0 = observables(zero)[0][0]
        for a in range(2):
            for b in range(2):
                ea = np.eye(2)[a] * eps
                eb = np.eye(2)[b] * eps
                Hess[a, b] = (observables((ea + eb)[None, :])[0][0] - observables((ea - eb)[None, :])[0][0] - observables((-ea + eb)[None, :])[0][0] + observables((-ea - eb)[None, :])[0][0]) / (4 * eps * eps)
        checks['finite_difference_reduced_Hessian_error'] = maxabs(Hess - np.diag(LAM))
    if polycheck:
        nu = 8 * D / delta ** 2
        z = nu * TAU
        coeff = [float(ive(0, z))]
        tail = 1 - coeff[0]
        order = 0
        while tail > 1e-11:
            order += 1
            coeff.append(float(2 * ive(order, z)))
            tail = 1 - sum(coeff)
            if order > 10000:
                raise RuntimeError('Unexpected polynomial size')
        M = sparse.eye(n * n, format='csr') - H / nu
        T0 = e.copy()
        total = coeff[0] * T0
        if order:
            T1 = contract(M, e)
            total += coeff[1] * T1
            for k in range(2, order + 1):
                T2 = 2 * contract(M, T1) - T0
                total += coeff[k] * T2
                (T0, T1) = (T1, T2)
        RP = -BETA * (cov - contract(e.T, total))
        result['polynomial_check'] = {'nu': nu, 'z': z, 'order_K': order, 'coefficient_sum': sum(coeff), 'omitted_coefficient_sum': tail, 'evolved_vector_max_difference': maxabs(total - evolved), 'response_R_max_difference': maxabs(RP - R), 'normalized_C_max_difference': maxabs((RP - R) / (BETA * NORMAL)), 'raw_entry_error_bound_max': float(BETA * np.outer(np.linalg.norm(e, axis=0), np.linalg.norm(e, axis=0)).max() * tail), 'configuration_qubits': int(math.ceil(math.log2(n * n))), 'two_configuration_register_qubits': int(2 * math.ceil(math.log2(n * n))), 'coefficient_qubits': int(math.ceil(math.log2(order + 1))), 'overlap_readout_qubits': 1, 'status': 'classical Chebyshev recurrence; not a compiled or executed quantum circuit'}
    result['wall_time_seconds'] = time.perf_counter() - start
    return result

def main():
    ROOT.mkdir(parents=True, exist_ok=True)
    begun = time.time()
    results = []
    designs = [(n, 4.0) for n in [4, 8, 16, 32, 64]] + [(48, e) for e in [3.0, 4.0, 5.0, 6.0]]
    for (n, extent) in designs:
        for model in ['harmonic', 'biquadratic']:
            res = run_case(n, extent, model, finite_fields=n in [16, 32] and extent == 4, polycheck=n in [4, 8, 16, 32] and extent == 4, densecheck=n in [4, 8] and extent == 4)
            results.append(res)
            print(json.dumps({'n': n, 'extent': extent, 'model': model, 'R43': res['apex_response_R'], 'C43': res['apex_response_C'], 'seconds': res['wall_time_seconds']}), flush=True)
            (ROOT / 'partial-results.json').write_text(json.dumps(primitive(results), indent=2) + '\n')
    out = {'status': 'completed deterministic synthetic contact-network calculation; no protein, device or biological-performance claim', 'timestamp_unix': begun, 'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), 'software': {'python': sys.version, 'numpy': np.__version__, 'scipy': scipy.__version__, 'platform': platform.platform()}, 'parameters': {'positions': R0, 'edges_zero_based': EDGES, 'stiffness': STIFFNESS, 'beta': BETA, 'mobility': MOBILITY, 'time_horizon': TAU, 'sender_zero_based': SENDER, 'receiver_zero_based': RECEIVER, 'reference_domain_extent': 4, 'field_steps': [0.01, 0.001, 0.0001, 1e-05], 'grids': designs, 'degree': DEG}, 'geometry': GEO, 'polynomial': {'monomial_powers': POWERS, 'residue_coefficients': A, 'analytic_reduced_harmonic_covariance': HARM_COV, 'harmonic_standard_deviations': SH, 'normalization': 'exact unbounded two-coordinate Gaussian moments of the fixed quartic Ei; common to every case'}, 'cases': results, 'resources': {'total_wall_seconds': time.time() - begun, 'max_rss_bytes': usage()['max_rss_bytes'], 'rss_unit': 'bytes; process high-water mark, not per-case incremental memory', 'process_user_cpu_seconds': usage()['user_cpu_seconds'], 'process_system_cpu_seconds': usage()['system_cpu_seconds']}, 'limitations': ['two of nine internal coordinates retained; omitted-mode error not measured', 'dimensionless synthetic geometry and uniform contact stiffness are not calibrated protein parameters', 'finite reflecting domains and grid refinement do not certify continuum accuracy', 'five nodes cannot test a top-five biological prioritization endpoint', 'polynomial evaluation is classical and does not establish circuit execution or quantum advantage']}
    (ROOT / 'results.json').write_text(json.dumps(primitive(out), indent=2) + '\n')
    fields = ['model', 'n', 'G', 'extent_slowest_harmonic_sigma', 'spacing', 'apex_response_R', 'apex_response_C', 'wall_time_seconds']
    with (ROOT / 'figure-data.csv').open('w', newline='') as f:
        wr = csv.DictWriter(f, fieldnames=fields)
        wr.writeheader()
        for r in results:
            wr.writerow({k: r[k] for k in fields})
    with (ROOT / 'finite-field-data.csv').open('w', newline='') as f:
        wr = csv.writer(f)
        wr.writerow(['model', 'n', 'h', 'covariance_R', 'forward_field_slope', 'absolute_difference'])
        for r in results:
            for row in r['finite_field_checks']:
                wr.writerow([r['model'], r['n'], row['h'], r['apex_response_R'], row['forward_slope'], abs(row['difference_from_covariance_response'])])
    print(json.dumps({'complete': True, 'wall_seconds': out['resources']['total_wall_seconds']}), flush=True)
if __name__ == '__main__':
    main()
