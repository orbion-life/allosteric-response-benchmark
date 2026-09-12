"""Frozen, coarse KRAS fixture and response-unit conditioning calculations."""
from pathlib import Path
import hashlib
import json
import math
import time
import numpy as np
from scipy.linalg import expm
from scipy.special import expit, logsumexp, ive


def mm(a, b):
    a, b = np.asarray(a), np.asarray(b)
    return np.einsum({(1,1): 'i,i->', (1,2): 'i,ij->j',
                      (2,1): 'ij,j->i', (2,2): 'ij,jk->ik'}[(a.ndim,b.ndim)],
                     a, b, optimize=False)


def coefficient_tail_bound(K, z):
    """Analytic upper bound on 2 exp(-z) sum_{k>K} I_k(z), evaluated in float.

    I_k(z) <= (z/2)^k exp[z^2/(4(k+1))]/k! from its positive series.
    The remaining factorial series is bounded geometrically. No directed-rounding
    certification is claimed for this floating-point evaluation of the bound.
    """
    if z == 0:
        return 0.0
    ratio = z / (2*(K+2))
    if ratio >= 1:
        return float('inf')
    return float(math.exp(math.log(2)-z+z*z/(4*(K+2))
                          +(K+1)*math.log(z/2)-math.lgamma(K+2)-math.log1p(-ratio)))


def shots_for_uniform_response(gamma, unique_pairs, s, halfwidth=.004, failure=.05):
    amplification = float(np.max(gamma)**2)
    overlap_error = min(2.0, halfwidth/(s*amplification))
    per = math.ceil(2*math.log(2*unique_pairs/failure)/overlap_error**2)
    return {"unique_observable_pairs": int(unique_pairs), "gamma_max_squared": amplification,
            "response_halfwidth": halfwidth, "family_failure_probability": failure,
            "required_uniform_overlap_halfwidth": overlap_error,
            "sufficient_shots_per_pair": per, "sufficient_total_shots": int(per*unique_pairs),
            "bound_scope": "uniform independent +/-1 sampling, Hoeffding and union bound; sufficient, not necessary"}


def prepare_fixture(inputs, output, grid_side=4, extent=4.0, algorithm_budget=.0005):
    start = time.perf_counter()
    inputs, output = Path(inputs), Path(output)
    output.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((inputs/'input-manifest.json').read_text())
    for row in manifest:
        assert hashlib.sha256((inputs/row['file']).read_bytes()).hexdigest() == row['sha256']
    protocol = json.loads((inputs/'preanalysis-protocol.json').read_text())
    meta = json.loads((inputs/'static-model.json').read_text())
    coordinate = np.load(inputs/'coordinate-d1.npz')
    static = np.load(inputs/'static-model.npz')
    r0, edges = coordinate['r0'], coordinate['edges']
    B, A, powers = coordinate['B'], coordinate['A'], coordinate['powers']
    own_sh = coordinate['harmonic_sd']
    all_mode = np.load(inputs/'all-mode-harmonic-scales.npz')
    sh = all_mode['harmonic_sd']
    beta, mu, kappa = (float(protocol[k]) for k in ('beta','mu','kappa'))
    lam = float(coordinate['eigenvalues'][0])
    assert B.shape[1] == 1 and grid_side == 4
    assert np.max(abs(B-static['B'][:,:1])) < 1e-12
    tau = 1/(mu*lam)
    sigma = 1/math.sqrt(beta*lam)
    q = np.linspace(-extent*sigma, extent*sigma, grid_side)
    delta = float(q[1]-q[0])
    positions = r0[None,:,:] + (q[:,None]*B[:,0]).reshape(grid_side, len(r0), 3)
    degree = np.bincount(edges.ravel(), minlength=len(r0))
    E = np.zeros((grid_side, len(r0)))
    energy = np.zeros(grid_side)
    for i,j in edges:
        length2 = float(np.sum((r0[i]-r0[j])**2))
        contact = kappa/(8*length2)*(np.sum((positions[:,i]-positions[:,j])**2,axis=1)-length2)**2
        energy += contact
        E[:,i] += contact/degree[i]
        E[:,j] += contact/degree[j]
    monomials = np.stack([q**int(p[0]) for p in powers], axis=1)
    polynomial_error = float(np.max(abs(mm(monomials,A.T)-E)))
    assert polynomial_error < 1e-10
    pi = np.exp(-beta*energy-logsumexp(-beta*energy))
    assert pi.min() > 0 and sh.min() > 0
    Ddiff = mu/beta
    L = np.zeros((grid_side,grid_side))
    for x in range(grid_side):
        for y in (x-1,x+1):
            if 0<=y<grid_side:
                L[x,y] = 2*Ddiff/delta**2*expit(-beta*(energy[y]-energy[x]))
        L[x,x] = -L[x].sum()
    nu = 4*Ddiff/delta**2
    H = -(np.sqrt(pi)[:,None]*L)/np.sqrt(pi)[None,:]
    P = np.eye(grid_side)+L/nu
    M = np.eye(grid_side)-H/nu
    assert np.max(abs(H-H.T))<1e-12
    weighted = np.sqrt(pi)[:,None]*(E-mm(pi,E))
    phi = monomials-mm(pi,monomials)
    raw_basis = np.sqrt(pi)[:,None]*phi
    # A basis change is exact on this fixed grid: it is not a physical-motion
    # reduction. Four symmetric points make centered q^4 proportional to q^2.
    U, singular, Vh = np.linalg.svd(raw_basis, full_matrices=False)
    threshold = max(raw_basis.shape)*np.finfo(float).eps*singular[0]
    rank = int(np.sum(singular>threshold))
    basis = U[:,:rank].copy()
    reconstruction = mm(A, Vh[:rank,:].T)*singular[:rank]
    for col in range(rank):
        if basis[np.argmax(abs(basis[:,col])),col] < 0:
            basis[:,col] *= -1
            reconstruction[:,col] *= -1
    residual = float(np.max(abs(mm(basis,reconstruction.T)-weighted)))
    assert residual < 1e-11
    norms = np.sqrt(np.sum(weighted**2,axis=0))
    rho = norms/sh
    propagator = expm(-tau*H)
    covariance = mm(weighted.T,weighted)
    delayed_exact = mm(mm(weighted.T,propagator),weighted)
    R_exact = -beta*(covariance-delayed_exact)
    C_exact = R_exact/(beta*np.outer(sh,sh))
    receiver = coordinate['receiver']
    candidates = np.flatnonzero(coordinate['candidate'])
    canonical = coordinate['canonical']
    scores = np.sqrt(np.mean(C_exact[receiver,:]**2,axis=0))
    ranking = np.array(sorted(candidates,key=lambda i:(-scores[i],canonical[i])))
    top, outside = ranking[:5], ranking[5:]
    gap = float(scores[ranking[4]]-scores[ranking[5]])
    receiver_rho = float(np.sqrt(np.mean(rho[receiver]**2)))
    z = nu*tau
    for K in range(1,128):
        tail = coefficient_tail_bound(K,z)
        score_tail = tail*rho*receiver_rho
        retained_gap = float(np.min(scores[top]-score_tail[top])-np.max(scores[outside]+score_tail[outside]))
        if tail*np.max(rho)**2 <= algorithm_budget and retained_gap >= .9*gap:
            break
    else:
        raise RuntimeError('No bounded truncation found')
    coefficients = ive(np.arange(K+1), z)
    coefficients[1:] *= 2
    s = float(coefficients.sum())
    Ts=[np.eye(grid_side),M]
    for k in range(2,K+1):
        Ts.append(2*mm(M,Ts[-1])-Ts[-2])
    polynomial = sum(float(a)*T for a,T in zip(coefficients,Ts))
    response_polynomial = -beta*(covariance-mm(mm(weighted.T,polynomial),weighted))
    C_polynomial = response_polynomial/(beta*np.outer(sh,sh))
    overlap_polynomial = mm(mm(basis.T,polynomial),basis)/s
    assert float(np.max(abs(C_polynomial-C_exact))) <= algorithm_budget+1e-10
    raw_norms = np.sqrt(np.sum(raw_basis**2,axis=0))
    raw_gamma = mm(abs(A),raw_norms)/sh
    ortho_gamma = np.sum(abs(reconstruction),axis=1)/sh
    comparisons = {
        'direct_all_residue_pairs': shots_for_uniform_response(rho,len(r0)*(len(r0)+1)//2,s),
        'raw_three_monomial_basis': shots_for_uniform_response(raw_gamma,len(powers)*(len(powers)+1)//2,s),
        'orthonormal_grid_basis': shots_for_uniform_response(ortho_gamma,rank*(rank+1)//2,s),
    }
    own_comparisons = {
        'direct_all_residue_pairs': shots_for_uniform_response(norms/own_sh,len(r0)*(len(r0)+1)//2,s),
        'raw_three_monomial_basis': shots_for_uniform_response(mm(abs(A),raw_norms)/own_sh,len(powers)*(len(powers)+1)//2,s),
        'orthonormal_grid_basis': shots_for_uniform_response(np.sum(abs(reconstruction),axis=1)/own_sh,rank*(rank+1)//2,s),
    }
    receiver_gamma = float(np.sqrt(np.mean(ortho_gamma[receiver]**2)))
    rank_overlap_cap = min((scores[i]-scores[j]-score_tail[i]-score_tail[j]) /
                           (s*receiver_gamma*(ortho_gamma[i]+ortho_gamma[j]))
                           for i in top for j in outside)
    matrix_overlap_cap = comparisons['orthonormal_grid_basis']['required_uniform_overlap_halfwidth']
    selected_overlap_error = min(matrix_overlap_cap, .5*rank_overlap_cap)
    m = rank*(rank+1)//2
    selected_shots = math.ceil(2*math.log(2*m/.05)/selected_overlap_error**2)
    sampling_score_halfwidth = s*selected_overlap_error*ortho_gamma*receiver_gamma
    total_score_halfwidth = score_tail+sampling_score_halfwidth
    rank_lower = float(np.min(scores[top]-total_score_halfwidth[top]))
    rank_upper = float(np.max(scores[outside]+total_score_halfwidth[outside]))
    assert rank_lower>rank_upper
    uniform_C_gap_budget = .25*gap
    uniform_gap_cost = shots_for_uniform_response(ortho_gamma,m,s,halfwidth=uniform_C_gap_budget)
    rank_family_comparison = {}
    for label,gamma,count in [('direct_all_residue_pairs',rho,len(r0)*(len(r0)+1)//2),
                              ('raw_three_monomial_basis',raw_gamma,len(powers)*(len(powers)+1)//2),
                              ('orthonormal_grid_basis',ortho_gamma,m)]:
        gamma_receiver = float(np.sqrt(np.mean(gamma[receiver]**2)))
        cap = min((scores[i]-scores[j]-score_tail[i]-score_tail[j]) /
                  (s*gamma_receiver*(gamma[i]+gamma[j])) for i in top for j in outside)
        eps = min(comparisons[label]['required_uniform_overlap_halfwidth'], .5*cap)
        per = math.ceil(2*math.log(2*count/.05)/eps**2)
        rank_family_comparison[label] = {'unique_pairs':count,'chosen_overlap_halfwidth':float(eps),
                                        'sufficient_shots_per_pair':per,'sufficient_total_shots':count*per,
                                        'scope':'same fixed-grid ranking intervals, common full-mode normalization, same global failure probability'}
    data = dict(q=q, E=E, U=energy, pi=pi, L=L, H=H, P=P, M=M, coefficients=coefficients,
                e=weighted, raw_basis=raw_basis, orthonormal_basis=basis, reconstruction=reconstruction,
                raw_A=A, powers=powers, harmonic_sd=sh, own_d1_harmonic_sd=own_sh,
                covariance=covariance, C_exact=C_exact,
                C_polynomial=C_polynomial, overlap_polynomial=overlap_polynomial,
                gamma_orthonormal=ortho_gamma, gamma_raw=raw_gamma, rho_direct=rho,
                exact_scores=scores, ranking=ranking, score_tail_bound=score_tail,
                sampling_score_halfwidth=sampling_score_halfwidth,
                canonical=coordinate['canonical'], receiver=coordinate['receiver'],
                candidate=coordinate['candidate'])
    np.savez(output/'fixture.npz',**data)
    summary = {'input_structure':meta['input'],'chain':meta['chain'],'retained_residues':len(r0),
               'contact_count':len(edges),'positive_harmonic_modes':len(static['eigenvalues']),
               'retained_physical_modes':1,'grid_values':grid_side,'grid_extent_harmonic_sigma':extent,
               'grid_spacing':delta,'kappa':kappa,'beta':beta,'mu':mu,'tau':tau,'lambda1':lam,
               'normalization':'pilot analytic all-492-mode harmonic SD of fixed quartic Ei, common to primary model comparisons',
               'stationary_probabilities':pi.tolist(),'boundary_probability':float(pi[0]+pi[-1]),
               'energy_polynomial_max_error':polynomial_error,'weighted_basis_reconstruction_max_error':residual,
               'raw_monomial_count':len(powers),'retained_algebraic_rank':rank,'singular_values':singular.tolist(),
               'singular_value_cutoff':float(threshold),'rank_reduction_scope':'exact observable-span reduction on this coarse grid; does not restore omitted protein motion',
               'nu':nu,'z':z,'K':K,'retained_coefficient_sum':s,
               'omitted_coefficient_sum_float':float(1-s),'analytic_tail_upper_bound_float':tail,
               'requested_algorithm_C_error_budget':algorithm_budget,
               'truncation_selection':'smallest K>=1 meeting C budget and retaining at least 90% of the exact same-grid fifth/sixth ranking gap under conservative score-tail bounds',
               'all_entry_truncation_C_upper_bound':float(tail*np.max(rho)**2),
               'observed_all_entry_truncation_C_error':float(np.max(abs(C_polynomial-C_exact))),
               'maximum_absolute_C':float(np.max(abs(C_exact))),
               'zero_response_meets_loose_sampling_cap':bool(np.max(abs(C_exact))<.004),
               'exact_same_grid_top5_canonical':canonical[top].tolist(),
               'exact_same_grid_fifth_sixth_score_gap':gap,
               'shot_family_comparison':comparisons,'own_d1_normalization_comparison':own_comparisons,
               'rank_precision_family_comparison':rank_family_comparison,
               'conservative_uniform_C_gap_cost':uniform_gap_cost,
               'selected_sampling':{'family_size':m,'failure_probability':.05,
                   'ranking_uniform_overlap_cap':float(rank_overlap_cap),'matrix_uniform_overlap_cap':matrix_overlap_cap,
                   'chosen_overlap_halfwidth':float(selected_overlap_error),
                   'selection_rule':'minimum of matrix cap and half of the conservative pairwise rank-separation cap',
                   'shots_per_overlap':selected_shots,'total_shots':selected_shots*m,
                   'achieved_uniform_C_halfwidth_bound':float(s*selected_overlap_error*np.max(ortho_gamma)**2),
                   'top5_lower_score_bound':rank_lower,'outside_top5_upper_score_bound':rank_upper,
                   'scope':'simultaneous sampling and truncation bounds for this fixed coarse model only; excludes device bias, grid error and omitted physical modes'},
               'preprocessing_seconds':time.perf_counter()-start,
               'input_manifest':manifest,'scope':'complete-circuit engineering fixture; unvalidated d1/four-point physical representation'}
    (output/'fixture-summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    return data,summary
