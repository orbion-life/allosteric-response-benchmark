"""Follow-up diagnostics declared after the initial fixed-parameter run.

These exploratory checks retain all tested settings. They are not selection for
biological accuracy. κ10 initial results were seen before adding higher stiffness.
"""
import os
from contractions import contract
from runtime_info import usage
from pathlib import Path
import csv, itertools, json, math, time, warnings
import numpy as np
from scipy import sparse
from scipy.special import expit, logsumexp, ive
from scipy.sparse.linalg import expm_multiply
from scipy.linalg import block_diag
import experiment as base
ROOT = Path(os.environ.get('PULSAR_OUTPUT_DIR', Path(__file__).resolve().parent / 'results'))
DECLARED = {'stiffnesses': [10, 30, 100, 300, 1000], 'beta_at_kappa10': [0.5, 1, 2], 'strength_grid': 32, 'strength_domain': 4, 'fixed_spacing_domains': [(49, 3), (65, 4), (81, 5), (97, 6)], 'dimension_grids': {2: [16, 32], 3: [8, 12, 16, 24], 4: [4, 6, 8, 12, 16]}, 'dimension_stiffness': 100, 'display_grid_check_stiffness': 100, 'display_grid_sizes': [4, 8, 16, 32, 64], 'time_multiples': [0, 0.1, 0.25, 0.5, 1, 2, 5, 10], 'explicit_walk_powers': [1, 2, 3, 5, 8]}

def make_model(kappa, beta, dim):
    K = base.K * kappa / base.STIFFNESS
    (vals, vecs) = np.linalg.eigh(K)
    ix = np.flatnonzero(vals > 1e-08)[:dim]
    B = vecs[:, ix].copy()
    lam = vals[ix]
    for j in range(dim):
        if B[np.argmax(np.abs(B[:, j])), j] < 0:
            B[:, j] *= -1
    powers = [p for degree in [2, 3, 4] for p in itertools.product(range(degree + 1), repeat=dim) if sum(p) == degree]
    A = np.zeros((5, len(powers)))
    for (i, j) in base.EDGES:
        a = base.R0[i] - base.R0[j]
        T = B[3 * i:3 * i + 3] - B[3 * j:3 * j + 3]
        l2 = contract(a, a)
        v = contract(2 * a, T)
        Q = contract(T.T, T)
        co = {}
        for u in range(dim):
            p = [0] * dim
            p[u] = 1
            co[tuple(p)] = v[u]
            for w in range(u, dim):
                p = [0] * dim
                p[u] += 1
                p[w] += 1
                co[tuple(p)] = Q[u, w] * (1 if u == w else 2)
        full = {p: 0.0 for p in powers}
        for (p, c) in co.items():
            for (q, z) in co.items():
                full[tuple((a + b for (a, b) in zip(p, q)))] += c * z * kappa / (8 * l2)
        for (n, p) in enumerate(powers):
            A[i, n] += full[p] / base.DEG[i]
            A[j, n] += full[p] / base.DEG[j]

    def gm(power):
        if any((p % 2 for p in power)):
            return 0.0
        return math.prod((math.prod(range(1, p, 2)) * (1 / (beta * l)) ** (p / 2) for (p, l) in zip(power, lam)))
    means = np.array([gm(p) for p in powers])
    cov = np.array([[gm(tuple((a + b for (a, b) in zip(p, q)))) for q in powers] for p in powers]) - np.outer(means, means)
    sh = np.sqrt(np.diag(contract(contract(A, cov), A.T)))
    return (B, lam, A, powers, sh)

def setup(kappa, beta, dim, n, extent, model, field=0):
    (B, lam, A, powers, sh) = make_model(kappa, beta, dim)
    half = extent / np.sqrt(beta * lam[0])
    axis = np.linspace(-half, half, n)
    delta = axis[1] - axis[0]
    q = np.stack(np.meshgrid(*[axis] * dim, indexing='ij'), axis=-1).reshape(-1, dim)
    ids = np.arange(n ** dim).reshape([n] * dim)
    xs = []
    ys = []
    for d in range(dim):
        sl = [slice(None)] * dim
        sr = sl.copy()
        sl[d] = slice(None, -1)
        sr[d] = slice(1, None)
        xs.append(ids[tuple(sl)].ravel())
        ys.append(ids[tuple(sr)].ravel())
    x = np.concatenate(xs)
    y = np.concatenate(ys)
    r = base.R0[None, :, :] + contract(q, B.T).reshape(-1, 5, 3)
    E = np.zeros((len(q), 5))
    UA = np.zeros(len(q))
    st = []
    for (i, j) in base.EDGES:
        l2 = np.sum((base.R0[i] - base.R0[j]) ** 2)
        r2 = np.sum((r[:, i] - r[:, j]) ** 2, axis=1)
        u = kappa / (8 * l2) * (r2 - l2) ** 2
        UA += u
        E[:, i] += u / base.DEG[i]
        E[:, j] += u / base.DEG[j]
        st.append(np.sqrt(r2 / l2))
    U = UA if model == 'biquadratic' else 0.5 * np.sum(q * q * lam, axis=1)
    U = U + field * E[:, 3]
    diffusivity = 1 / beta
    rate = 2 * diffusivity / delta ** 2
    xy = rate * expit(-beta * (U[y] - U[x]))
    yx = rate * expit(beta * (U[y] - U[x]))
    ex = np.bincount(x, weights=xy, minlength=len(q)) + np.bincount(y, weights=yx, minlength=len(q))
    rows = np.concatenate([x, y, np.arange(len(q))])
    cols = np.concatenate([y, x, np.arange(len(q))])
    L = sparse.csr_matrix((np.r_[xy, yx, -ex], (rows, cols)), shape=(len(q), len(q)))
    off = -np.sqrt(xy * yx)
    H = sparse.csr_matrix((np.r_[off, off, ex], (rows, cols)), shape=L.shape)
    pi = np.exp(-beta * U - logsumexp(-beta * U))
    ec = E - contract(pi, E)
    e = np.sqrt(pi)[:, None] * ec
    return {'B': B, 'lam': lam, 'q': q, 'E': E, 'H': H, 'L': L, 'pi': pi, 'e': e, 'sh': sh, 'tau': 1 / lam[0], 'delta': delta, 'half': half, 'A': A, 'powers': powers, 'stretch': np.array(st).T, 'r': r}

def calculate(kappa=10, beta=1, dim=2, n=32, extent=4, model='biquadratic', curve=False):
    (ROOT / 'follow-up-design.json').write_text(json.dumps(DECLARED, indent=2) + '\n')
    start = time.perf_counter()
    o = setup(kappa, beta, dim, n, extent, model)
    H = o['H']
    e = o['e']
    cov = contract(e.T, e)
    f = expm_multiply(-o['tau'] * H, e, traceA=-o['tau'] * H.diagonal().sum())
    R = -beta * (cov - contract(e.T, f))
    C = R / (beta * np.outer(o['sh'], o['sh']))
    distort = np.max(np.abs(o['stretch'] - 1), axis=1)
    pi = o['pi']
    r = o['r']
    boundary = np.any(np.isclose(np.abs(o['q']), o['half']), axis=1)
    result = {'kappa': kappa, 'beta': beta, 'dimensions': dim, 'n': n, 'G': n ** dim, 'extent': extent, 'spacing': o['delta'], 'model': model, 'horizon': o['tau'], 'R': R, 'C': C, 'C43': float(C[4, 3]), 'R43': float(R[4, 3]), 'retained_eigenvalues': o['lam'], 'B': o['B'], 'harmonic_sd': o['sh'], 'checks': {'L_conservation': base.maxabs(contract(o['L'], np.ones(len(pi)))), 'stationarity': base.maxabs(contract(pi, o['L'])), 'response_symmetry': base.maxabs(R - R.T)}, 'diagnostics': {'outer_grid_shell_probability': float(pi[boundary].sum()), 'mean_maximum_contact_strain': float(contract(pi, distort)), 'probability_any_strain_gt_0p2': float(pi[distort > 0.2].sum()), 'probability_any_strain_gt_0p1': float(pi[distort > 0.1].sum()), 'probability_strong_compression': float(pi[np.min(o['stretch'], axis=1) < np.sqrt(0.5)].sum()), 'probability_apex_overlap_lt_0p5': float(pi[np.linalg.norm(r[:, 3] - r[:, 4], axis=1) < 0.5].sum())}, 'time_curve': []}
    if curve:
        for fac in DECLARED['time_multiples']:
            ef = expm_multiply(-fac * o['tau'] * H, e, traceA=-fac * o['tau'] * H.diagonal().sum())
            rr = -beta * (cov - contract(e.T, ef))
            cc = rr / (beta * np.outer(o['sh'], o['sh']))
            result['time_curve'].append({'t_over_tau': fac, 'R43': float(rr[4, 3]), 'C43': float(cc[4, 3])})
    result['seconds'] = time.perf_counter() - start
    return result

def explicit_walk(kappa=10):
    (ROOT / 'follow-up-design.json').write_text(json.dumps(DECLARED, indent=2) + '\n')
    start = time.perf_counter()
    o = setup(kappa, 1, 2, 4, 4, 'biquadratic')
    G = 16
    nu = 8 / o['delta'] ** 2
    P = np.eye(G) + o['L'].toarray() / nu
    blocks = []
    for row in P:
        v = np.sqrt(row)
        u = np.eye(G)[0] - v
        unit = np.eye(G) if np.linalg.norm(u) < 1e-14 else np.eye(G) - 2 * np.outer(u, u) / contract(u, u)
        blocks.append(unit)
    V = block_diag(*blocks)
    swap = np.zeros((G * G, G * G))
    for i in range(G):
        for j in range(G):
            swap[j * G + i, i * G + j] = 1
    UM = contract(contract(V.T, swap), V)
    inds = np.arange(G) * G
    diag = np.full(G * G, -1.0)
    diag[inds] = 1
    W = diag[:, None] * UM
    M = np.eye(G) - o['H'].toarray() / nu
    result = {'kappa': kappa, 'G': G, 'unitary_dimension': G * G, 'V_unitarity_error': base.maxabs(contract(V.T, V) - np.eye(G * G)), 'UM_involution_error': base.maxabs(contract(UM, UM) - np.eye(G * G)), 'block_M_error': base.maxabs(UM[np.ix_(inds, inds)] - M), 'walk_power_checks': []}
    maxk = max(DECLARED['explicit_walk_powers'])
    WP = np.eye(G * G)
    T0 = np.eye(G)
    T1 = M
    for k in range(1, maxk + 1):
        WP = contract(WP, W)
        if k == 1:
            T = T1
        else:
            T = contract(2 * M, T1) - T0
            (T0, T1) = (T1, T)
        if k in DECLARED['explicit_walk_powers']:
            result['walk_power_checks'].append({'k': k, 'projected_power_max_error': base.maxabs(WP[np.ix_(inds, inds)] - T)})
    result['seconds'] = time.perf_counter() - start
    result['status'] = 'explicit dense unitary algebra; no gate transpilation or quantum device execution'
    return result

def validate_display(model):
    o = setup(100, 1, 2, 32, 4, model)
    H = o['H']
    e = o['e']
    cov = contract(e.T, e)
    f = expm_multiply(-o['tau'] * H, e, traceA=-o['tau'] * H.diagonal().sum())
    R = -(cov - contract(e.T, f))
    normal = np.outer(o['sh'], o['sh'])
    checks = []
    for field in [0.01, 0.001, 0.0001, 1e-05]:
        oh = setup(100, 1, 2, 32, 4, model, field=field)
        fh = expm_multiply(o['tau'] * oh['L'], o['E'][:, 4], traceA=o['tau'] * oh['L'].diagonal().sum())
        slope = (contract(o['pi'], fh) - contract(o['pi'], o['E'][:, 4])) / field
        checks.append({'h': field, 'forward_slope': float(slope), 'difference_from_R': float(slope - R[4, 3])})
    Z = np.array([np.prod(o['q'] ** p, axis=1) for p in o['powers']]).T
    ez = np.sqrt(o['pi'])[:, None] * (Z - contract(o['pi'], Z))
    zf = expm_multiply(-o['tau'] * H, ez, traceA=-o['tau'] * H.diagonal().sum())
    Rrec = contract(contract(-o['A'], contract(ez.T, ez - zf)), o['A'].T)
    nu = 8 / o['delta'] ** 2
    z = nu * o['tau']
    coeff = [float(ive(0, z))]
    k = 0
    while 1 - sum(coeff) > 1e-11:
        k += 1
        coeff.append(float(2 * ive(k, z)))
    M = sparse.eye(len(o['pi']), format='csr') - H / nu
    T0 = e
    T1 = contract(M, e)
    total = coeff[0] * T0 + coeff[1] * T1
    for j in range(2, k + 1):
        T2 = contract(2 * M, T1) - T0
        total += coeff[j] * T2
        (T0, T1) = (T1, T2)
    RP = -(cov - contract(e.T, total))
    return {'model': model, 'kappa': 100, 'n': 32, 'R43': float(R[4, 3]), 'C43': float(R[4, 3] / normal[4, 3]), 'finite_field': checks, 'reconstruction_R_error': base.maxabs(Rrec - R), 'fixed_observable_polynomial_error': base.maxabs(o['E'] - contract(Z, o['A'].T)), 'chebyshev_order': k, 'coefficient_tail': 1 - sum(coeff), 'chebyshev_R_error': base.maxabs(RP - R), 'chebyshev_C_error': base.maxabs((RP - R) / normal), 'scope': 'classical numerical/operator validation of the proposed discrete response, not quantum execution'}

def main():
    ROOT.mkdir(parents=True, exist_ok=True)
    (ROOT / 'follow-up-design.json').write_text(json.dumps(DECLARED, indent=2) + '\n')
    start = time.perf_counter()
    out = {'design': DECLARED, 'status': 'exploratory follow-up diagnostics to initial kappa10 fixture', 'strength': [], 'beta': [], 'fixed_spacing_domain': [], 'dimension': [], 'grid_kappa100': [], 'walk': [], 'fixed_spacing_domain_kappa100': [], 'display_validation': []}
    tasks = []
    for k in DECLARED['stiffnesses']:
        for m in ['harmonic', 'biquadratic']:
            tasks.append(('strength', dict(kappa=k, model=m, curve=k in [10, 100, 1000])))
    for beta in [0.5, 2]:
        for m in ['harmonic', 'biquadratic']:
            tasks.append(('beta', dict(beta=beta, model=m)))
    for (n, extent) in DECLARED['fixed_spacing_domains']:
        for m in ['harmonic', 'biquadratic']:
            tasks.append(('fixed_spacing_domain', dict(n=n, extent=extent, model=m)))
    for (n, extent) in DECLARED['fixed_spacing_domains']:
        for m in ['harmonic', 'biquadratic']:
            tasks.append(('fixed_spacing_domain_kappa100', dict(kappa=100, n=n, extent=extent, model=m)))
    for (dim, ns) in DECLARED['dimension_grids'].items():
        for n in ns:
            for m in ['harmonic', 'biquadratic']:
                tasks.append(('dimension', dict(kappa=100, dim=dim, n=n, model=m)))
    for n in DECLARED['display_grid_sizes']:
        for m in ['harmonic', 'biquadratic']:
            tasks.append(('grid_kappa100', dict(kappa=100, n=n, model=m)))
    for (group, kwargs) in tasks:
        r = calculate(**kwargs)
        if group == 'dimension':
            sh4 = make_model(100, 1, 4)[-1]
            r['C_common_d4'] = r['R'] / np.outer(sh4, sh4)
            r['C43_common_d4'] = float(r['C_common_d4'][4, 3])
            r['normalization_comparison'] = 'C_common_d4 uses the same analytic four-coordinate harmonic SD for d2,d3,d4; ordinary C uses its declared within-dimension SD'
        out[group].append(r)
        print(json.dumps({'group': group, **{k: r[k] for k in ['kappa', 'beta', 'dimensions', 'n', 'model', 'C43', 'seconds']}}), flush=True)
        (ROOT / 'sensitivity-partial.json').write_text(json.dumps(base.primitive(out), indent=2) + '\n')
    for k in [10, 100]:
        out['walk'].append(explicit_walk(k))
    for model in ['harmonic', 'biquadratic']:
        out['display_validation'].append(validate_display(model))
    out['total_seconds'] = time.perf_counter() - start
    out['max_rss_bytes'] = usage()['max_rss_bytes']
    (ROOT / 'sensitivity-results.json').write_text(json.dumps(base.primitive(out), indent=2) + '\n')
    with (ROOT / 'sensitivity-figure-data.csv').open('w', newline='') as f:
        wr = csv.writer(f)
        wr.writerow(['family', 'kappa', 'beta', 'dimensions', 'n', 'G', 'extent', 'model', 'R43', 'C43', 'mean_max_contact_strain', 'prob_strain_gt_0p2', 'seconds'])
        for group in ['strength', 'beta', 'fixed_spacing_domain', 'fixed_spacing_domain_kappa100', 'dimension', 'grid_kappa100']:
            for r in out[group]:
                wr.writerow([group] + [r[k] for k in ['kappa', 'beta', 'dimensions', 'n', 'G', 'extent', 'model', 'R43', 'C43']] + [r['diagnostics']['mean_maximum_contact_strain'], r['diagnostics']['probability_any_strain_gt_0p2'], r['seconds']])
    with (ROOT / 'time-curve-data.csv').open('w', newline='') as f:
        wr = csv.writer(f)
        wr.writerow(['kappa', 'model', 't_over_tau', 'R43', 'C43'])
        for r in out['strength']:
            for row in r['time_curve']:
                wr.writerow([r['kappa'], r['model'], row['t_over_tau'], row['R43'], row['C43']])
    print(json.dumps({'complete': True, 'seconds': out['total_seconds']}), flush=True)
if __name__ == '__main__':
    main()
