"""Cost the frozen KRAS d2/33² model. This does not execute its quantum circuit.

All matrix exponentials here are explicitly classical validation/planning references.
No protein coordinate, receiver, energy, normalization or candidate rule is changed.
"""
import argparse
import hashlib
import json
import math
import platform
import resource
import time
from pathlib import Path

import numpy as np
import scipy
from scipy import sparse, special
from scipy.sparse.linalg import expm_multiply

CAPS = dict(total_local_shots=2_000_000, statevector_bytes=512*1024**2,
            primitive_wall_seconds=20, primitive_rss_bytes=1024**3,
            primitive_native_operations=200_000)
ETA = .05
POLYNOMIAL_C_CAP = .0005  # Half the existing .001 polynomial/other deterministic allocation.
SAMPLING_C_CAP = .004


def mm(a, b):
    # Avoid platform-BLAS warning artefacts in small dense products; same contraction.
    return np.einsum('ik,kj->ij', a, b, optimize=False)


def maxabs(a):
    return float(np.max(np.abs(a)))


def tail_bound(z, degree):
    """Chernoff bound for sum_{k>K} (2 exp(-z) I_k(z)).

    The Bessel generating function is the MGF of a symmetric Skellam variable:
    E exp(theta D) = exp(z(cosh(theta)-1)). Apply Markov to each tail and
    choose theta=asinh((K+1)/z). Mathematical bound, evaluated in float64;
    it is not a directed-rounding numerical certificate.
    """
    k = degree+1
    return min(1., 2*math.exp(-k*math.asinh(k/z)+z*(math.hypot(1, k/z)-1)))


def coefficients(z, degree):
    c = special.ive(np.arange(degree+1), z)
    c[1:] *= 2
    return c


def choose_degree(z, rho2, cap):
    for degree in range(2049):
        if rho2*tail_bound(z, degree) <= cap:
            return degree
    raise RuntimeError('Degree exceeds declared analysis ceiling of 2048')


def shot_count(epsilon, pairs, failure):
    """Sufficient, not optimal: ±1 Hoeffding and a union bound over all pairs."""
    if epsilon <= 0 or failure <= 0:
        raise ValueError('A positive error and failure probability are required')
    return math.ceil(2*math.log(2*pairs/failure)/(min(epsilon, 2.)**2))


def basis_from(matrix):
    u, singular, _ = np.linalg.svd(matrix, full_matrices=False)
    tolerance = max(matrix.shape)*np.finfo(float).eps*singular[0]
    rank = int(np.sum(singular > tolerance))
    u = u[:, :rank].copy()
    for column in range(rank):
        if u[np.argmax(np.abs(u[:, column])), column] < 0:
            u[:, column] *= -1
    return u, singular, tolerance


def scores(C, receiver):
    return np.sqrt(np.mean(C[:, receiver]**2, axis=1))


def order_scores(value, candidates, canonical):
    return np.array(sorted(np.flatnonzero(candidates), key=lambda i:(-value[i], canonical[i])))


def interval_decision(means, D, equilibrium, scale, retained_sum, tail, epsilon,
                      gamma, rho, receiver, candidates, canonical, device_score_bound=0.):
    """A production-valid decision interface: it has no exact target scores/ranks.

    Inputs means are measured shared-basis overlap means. In this artifact ONLY,
    calls use exact polynomial means and are labelled oracle-centred planning.
    This criterion accounts for sampling, polynomial tail and declared device bias
    on this fixed finite model. It cannot certify errors from omitted protein modes.
    """
    delayed = retained_sum*mm(mm(D, means), D.T)
    C_hat = -(equilibrium-delayed)/np.outer(scale, scale)
    score_hat = scores(C_hat, receiver)
    gamma_f = float(np.sqrt(np.mean(gamma[receiver]**2)))
    rho_f = float(np.sqrt(np.mean(rho[receiver]**2)))
    halfwidth = tail*rho*rho_f + retained_sum*epsilon*gamma*gamma_f + device_score_bound
    ranked = order_scores(score_hat, candidates, canonical)
    separation = float(np.min(score_hat[ranked[:5]]-halfwidth[ranked[:5]])-
                       np.max(score_hat[ranked[5:]]+halfwidth[ranked[5:]]))
    return dict(certified=bool(separation > 0), separation=separation,
                top5_canonical=canonical[ranked[:5]].astype(int).tolist(),
                max_candidate_halfwidth=float(np.max(halfwidth[candidates])))


def polynomial_overlap(M, u, coeff):
    old = u.copy()
    acc = coeff[0]*old
    if len(coeff) > 1:
        now = M@u
        acc += coeff[1]*now
        for c in coeff[2:]:
            new = 2*(M@now)-old
            acc += c*new
            old, now = now, new
    return mm(u.T, acc)/coeff.sum()


def resources(states, degree):
    n = math.ceil(math.log2(states))
    a = max(1, math.ceil(math.log2(degree+1)))
    walks = 2**a-1  # Existing binary SELECT loops over all powers, including padding.
    q = 2*n+a+1
    return dict(configuration_qubits_each=n, coefficient_qubits=a,
                readout_qubits=1, workspace_qubits_in_existing_abstract_encoding=0,
                abstract_total_qubits_before_any_compiler_workspace=q,
                full_complex128_statevector_bytes=16*2**q,
                binary_select_walk_invocations=walks,
                V_or_inverse_invocations=2*walks,
                row_controlled_PREP_invocations=2*walks*2**n,
                register_pair_swaps=walks*n,
                zero_workspace_reflections=walks,
                observable_PREP_and_inverse=2, coefficient_PREP_and_inverse=2,
                full_statevector_within_local_cap=bool(16*2**q <= CAPS['statevector_bytes']))


def run(inputs, output):
    start = time.perf_counter()
    output.mkdir(parents=True, exist_ok=True)
    for filename, sha in json.loads((inputs/'manifest.json').read_text()).items():
        assert hashlib.sha256((inputs/filename).read_bytes()).hexdigest() == sha, filename
    r = np.load(inputs/'primary.npz')
    static = np.load(inputs/'static-model.npz')
    meta = json.loads((inputs/'primary.json').read_text())
    scale = np.load(inputs/'all-mode-harmonic-scales.npz')['harmonic_sd']
    receiver = static['receiver'].astype(int)
    candidates = static['candidate'].astype(bool)
    canonical = static['canonical']
    G = len(r['pi'])
    H = sparse.csr_matrix((r['H_data'], r['H_indices'], r['H_indptr']), shape=(G,G))
    tau, delta = meta['tau_model_time'], meta['delta_A']
    nu = 4*meta['d']/delta**2
    z = nu*tau
    M = sparse.eye(G, format='csr')-H/nu
    pi = r['pi']
    # P = diag(sqrt(pi))^{-1} M diag(sqrt(pi)). All pi are positive here.
    assert np.all(np.isfinite(pi)) and np.all(pi > 0) and np.all(scale > 0)
    P = sparse.diags(1/np.sqrt(pi))@M@sparse.diags(np.sqrt(pi))
    F = np.sqrt(pi)[:,None]*r['monomial_centered']
    E = mm(F, r['A'].T)
    u, singular, tolerance = basis_from(F)
    rank = u.shape[1]
    pairs = rank*(rank+1)//2
    D = mm(E.T, u)
    gamma = np.sum(np.abs(D), axis=1)/scale
    rho = np.linalg.norm(D, axis=1)/scale
    rho2, gamma2 = float(np.max(rho)**2), float(np.max(gamma)**2)
    gamma_f = float(np.sqrt(np.mean(gamma[receiver]**2)))
    rho_f = float(np.sqrt(np.mean(rho[receiver]**2)))
    physical_u, physical_singular, physical_tol = basis_from(E)
    physical_D = mm(E.T, physical_u)
    physical_gamma = np.sum(abs(physical_D), axis=1)/scale
    norms = np.linalg.norm(F, axis=0)
    raw_gamma = np.sum(abs(r['A']*norms[None,:]), axis=1)/scale
    equilibrium = mm(E.T, E)
    t = time.perf_counter()
    evolved = expm_multiply(-tau*H, u, traceA=-tau*H.diagonal().sum())
    shared_exact = mm(u.T, evolved)
    C_exact = -(equilibrium-mm(mm(D, shared_exact), D.T))/np.outer(scale, scale)
    shared_solver_seconds = time.perf_counter()-t
    exact_score = scores(C_exact, receiver)
    ranking = order_scores(exact_score, candidates, canonical)
    gap = float(exact_score[ranking[4]]-exact_score[ranking[5]])
    check = dict(max_C_difference_from_frozen=maxabs(C_exact-r['C']),
                 max_score_difference_from_frozen=maxabs(exact_score-r['score']),
                 max_weighted_observable_reconstruction=maxabs(mm(u,D.T)-E),
                 max_orthonormality_error=maxabs(mm(u.T,u)-np.eye(rank)),
                 max_covariance_difference_from_frozen=maxabs(mm(F.T,F)-r['monomial_covariance']),
                 max_normalized_equilibrium_difference_from_frozen=maxabs(equilibrium/np.outer(scale,scale)-r['equilibrium']),
                 maximum_stochastic_row_error=maxabs(np.asarray(P.sum(axis=1)).ravel()-1),
                 minimum_stored_transition=float(P.data.min()),
                 max_detailed_balance_error=maxabs((sparse.diags(pi)@P-P.T@sparse.diags(pi)).data),
                 max_exit_rate=float(H.diagonal().max()),
                 max_exit_rate_over_uniformization=float(H.diagonal().max()/nu))
    assert check['max_C_difference_from_frozen'] < 1e-10
    assert check['max_score_difference_from_frozen'] < 1e-10
    assert check['max_normalized_equilibrium_difference_from_frozen'] < 1e-10
    assert check['maximum_stochastic_row_error'] < 1e-12
    assert check['minimum_stored_transition'] >= 0
    assert np.array_equal(canonical[ranking[:5]], meta['top5_canonical'])
    fixed = []
    for budget in [.001, POLYNOMIAL_C_CAP]:
        K = choose_degree(z, rho2, budget)
        coeff = coefficients(z, K)
        s, tail = float(coeff.sum()), tail_bound(z,K)
        means = polynomial_overlap(M,u,coeff)
        C_poly = -(equilibrium-s*mm(mm(D,means),D.T))/np.outer(scale,scale)
        epsilon = min(2.,SAMPLING_C_CAP/(s*gamma2))
        shots = shot_count(epsilon,pairs,ETA)
        decision = interval_decision(means,D,equilibrium,scale,s,tail,epsilon,gamma,rho,
                                     receiver,candidates,canonical)
        fixed.append(dict(polynomial_C_cap=budget,degree=K,retained_sum=s,
                          tail_bound=tail,numerical_tail=1-s,tail_C_bound=tail*rho2,
                          actual_polynomial_C_error=maxabs(C_poly-C_exact),
                          overlap_sampling_halfwidth=epsilon,shots_per_overlap=shots,
                          total_shots=shots*pairs,simultaneous_failure_probability=ETA,
                          oracle_centred_interval_diagnostic=decision,resources=resources(G,K)))
        assert 1-s <= tail+1e-13
        assert maxabs(C_poly-C_exact) <= tail*rho2+1e-13
    # Oracle-calibrated estimate: exact reference shortlist and exact pairwise gaps
    # are used ONLY here, never inside interval_decision or the adaptive schedule.
    top, other = ranking[:5], ranking[5:]
    K = 0
    while True:
        tail = tail_bound(z,K)
        deterministic_pair = tail*rho_f*(rho[top,None]+rho[other][None,:])
        margins = exact_score[top,None]-exact_score[other][None,:]
        if np.all(deterministic_pair <= .1*margins):
            break
        K += 1
        assert K <= 2048
    coeff = coefficients(z,K)
    s = float(coeff.sum())
    epsilon_rank = float(np.min((margins-deterministic_pair)/
                           (s*gamma_f*(gamma[top,None]+gamma[other][None,:]))))
    epsilon = min(SAMPLING_C_CAP/(s*gamma2),.5*epsilon_rank)
    shots = shot_count(epsilon,pairs,ETA)
    oracle = dict(status='Oracle-calibrated sufficient planning estimate; no executed primary quantum circuit',
                  degree=K,retained_sum=s,tail_bound=tail,tail_C_bound=tail*rho2,
                  maximum_admissible_overlap_error_by_exact_ranking=epsilon_rank,
                  planned_overlap_halfwidth=epsilon,shots_per_overlap=shots,
                  total_shots=pairs*shots,resources=resources(G,K))
    family_comparisons=[]
    for family,g in [('weighted_monomial_SVD',gamma),('weighted_physical_observable_SVD',physical_gamma),
                     ('separately_normalized_monomials',raw_gamma)]:
        gf=float(np.sqrt(np.mean(g[receiver]**2)))
        admissible=float(np.min((margins-deterministic_pair)/(s*gf*(g[top,None]+g[other][None,:]))))
        precision=min(SAMPLING_C_CAP/(s*float(np.max(g)**2)),.5*admissible)
        count=shot_count(precision,pairs,ETA)
        family_comparisons.append(dict(family=family,shared_pairs=pairs,
            gamma_max_squared=float(np.max(g)**2),receiver_rms_gamma=gf,
            oracle_degree=K,oracle_halfwidth=precision,oracle_sufficient_total_shots=count*pairs,
            scope='Algebraically exact known-observable representation; oracle-gap-calibrated planning only'))
    # Every stage's degree and precision use known reconstruction constants only.
    # Alpha_l telescopes to ETA. New K means a changed experiment: charge fresh shots.
    stages, cumulative = [], 0
    for stage in range(16):
        pcap = POLYNOMIAL_C_CAP/2**stage
        scap = SAMPLING_C_CAP/2**stage
        alpha = ETA/((stage+1)*(stage+2))
        K = choose_degree(z,rho2,pcap)
        coeff = coefficients(z,K)
        s, tail = float(coeff.sum()), tail_bound(z,K)
        epsilon = min(2.,scap/(s*gamma2))
        shots = shot_count(epsilon,pairs,alpha)
        cumulative += shots*pairs
        means = polynomial_overlap(M,u,coeff)
        ideal = interval_decision(means,D,equilibrium,scale,s,tail,epsilon,gamma,rho,
                                 receiver,candidates,canonical)
        device = interval_decision(means,D,equilibrium,scale,s,tail,epsilon,gamma,rho,
                                  receiver,candidates,canonical,device_score_bound=.001)
        stages.append(dict(stage=stage,degree=K,polynomial_C_cap=pcap,sampling_C_cap=scap,
                           stage_failure_probability=alpha,tail_bound=tail,
                           retained_sum=s,overlap_halfwidth=epsilon,
                           shots_per_overlap=shots,total_fresh_shots=pairs*shots,
                           cumulative_shots=cumulative,
                           within_local_shot_cap=cumulative <= CAPS['total_local_shots'],
                           ideal_oracle_centred_diagnostic=ideal,
                           device_point001_oracle_centred_diagnostic=device,
                           resources=resources(G,K)))
        if ideal['certified']:
            break
    else:
        raise AssertionError('No ideal oracle-centred separation by declared stage limit')
    padded = 2**math.ceil(math.log2(G))
    pad_states = padded-G
    cov_normalized = mm((F/norms).T,F/norms)
    np.savez_compressed(output/'analysis-arrays.npz',basis=u,reconstruction=D,
                        physical_basis=physical_u,physical_reconstruction=physical_D,
                        gamma=gamma,rho=rho,scale=scale,shared_exact=shared_exact,
                        C_exact=C_exact,exact_score=exact_score,ranking=ranking,
                        P_data=P.data,P_indices=P.indices,P_indptr=P.indptr,
                        monomial_correlation=cov_normalized,
                        fixed_coefficients=coefficients(z,fixed[1]['degree']),
                        oracle_coefficients=coefficients(z,oracle['degree']))
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if platform.system() != 'Darwin': rss *= 1024
    results = dict(scope='Frozen KRAS 4OBE, biquadratic d2/33²/e4; cost analysis, no primary quantum execution',
                   input_hashes=json.loads((inputs/'manifest.json').read_text()),
                   versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),
                   model=dict(states=G,grid=[33,33],residues=len(scale),receiver_count=int(len(receiver)),
                              candidate_count=int(candidates.sum()),tau=tau,delta=delta,nu=nu,z=z,
                              same_full_harmonic_dimensions=492,physical_modes_retained=2,
                              common_normalization='frozen all-mode harmonic SD of fixed quartic Ei'),
                   observables=dict(monomial_rank=rank,monomial_rank_tolerance=tolerance,
                         monomial_singular_values=singular.tolist(),physical_observable_rank=physical_u.shape[1],
                         physical_rank_tolerance=physical_tol,physical_singular_values=physical_singular.tolist(),
                         unique_symmetric_overlaps=pairs,primary_basis='weighted monomial SVD; canonical column signs',
                         reconstruction_residual=check['max_weighted_observable_reconstruction'],
                         primary_gamma_max_squared=gamma2,rho_max_squared=rho2,
                         receiver_rms_gamma=gamma_f,receiver_rms_rho=rho_f,
                         physical_SVD_gamma_max_squared=float(np.max(physical_gamma)**2),
                         separately_normalized_monomial_gamma_max_squared=float(np.max(raw_gamma)**2)),
                   checks=check,
                   ranking=dict(top6_canonical=canonical[ranking[:6]].astype(int).tolist(),
                                top6_scores=exact_score[ranking[:6]].tolist(),fifth_sixth_gap=gap,
                                uniform_C_error_sufficient_for_separation=gap/2,
                                maximum_absolute_C=maxabs(C_exact),
                                zero_response_passes_loose_sampling_C_cap=maxabs(C_exact)<SAMPLING_C_CAP),
                   fixed_budget=fixed,oracle_calibrated_rank_planning=oracle,
                   exact_readout_basis_comparison=family_comparisons,
                   adaptive_schedule=dict(status='Oracle-centred planning only, no sampled measurements',
                       policy='Degree and precision halve per stage without target ranking; stage alpha=eta/((l+1)(l+2)); fresh samples when K changes',
                       covers='Simultaneous sampling, polynomial tail and declared device bias on fixed finite model only',
                       stages=stages),
                   loading=dict(padded_states=padded,padding_states=pad_states,
                       valid_stochastic_nonzeros=int(P.nnz),padded_stochastic_nonzeros=int(P.nnz+pad_states),
                       maximum_valid_row_support=int(np.max(np.diff(P.indptr))),
                       padded_dense_float64_transition_bytes=padded*padded*8,
                       padded_dense_float64_amplitude_table_bytes=padded*padded*8,
                       actual_sparse_P_array_bytes=int(P.data.nbytes+P.indices.nbytes+P.indptr.nbytes),
                       valid_monomial_table_bytes=int(F.nbytes),valid_physical_observable_table_bytes=int(E.nbytes),
                       primary_basis_float64_table_bytes=int(u.nbytes),
                       padded_basis_float64_table_bytes=int(padded*rank*8),
                       note='Logical row-PREP counts describe existing exhaustive constructor, not synthesized CX bounds or an optimal loading algorithm'),
                   caps=CAPS,
                   timings=dict(analysis_total_seconds=time.perf_counter()-start,
                       loaded_same_model_shared_solver_and_reconstruction_seconds=shared_solver_seconds,
                       analysis_process_peak_RSS_bytes=int(rss),
                       original_primary_grid_seconds=meta['seconds'],
                       timing_scope='Analysis includes loaded-data SVD, classical references and degree/shot planning; excludes upstream PDB/Hessian/normalizer, installation, quantum compilation'),
                   interpretation=['This fixed physical slice previously failed the full-mode harmonic compression gate.',
                       'An exact rank-12 readout factorization does not recover the 490 omitted physical modes.',
                       'No primary quantum gate count, hardware execution, quantum speedup or biological advantage is demonstrated.',
                       'Conservative sufficient shot budgets are not lower bounds on optimal estimators.'])
    (output/'results.json').write_text(json.dumps(results,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(rank=rank,gamma2=gamma2,gap=gap,fixed=fixed[1],oracle=oracle,
                          last_stage=stages[-1],checks=check),indent=2))


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs',type=Path,default=Path(__file__).parent/'inputs')
    parser.add_argument('--output',type=Path,default=Path(__file__).parent/'results')
    args=parser.parse_args()
    run(args.inputs,args.output)
