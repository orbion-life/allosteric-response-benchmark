"""Harmonic response by exact polynomial moments and separable finite dynamics."""
import argparse
import hashlib
from itertools import product
import json
import math
from pathlib import Path
import resource
import sys
import time
import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.special import eval_hermitenorm, expit, logsumexp, roots_hermitenorm
from model import (ROOT, analytic_gaussian, build_model, contract_moments,
                   energies, grid_shape, polynomial_values, protocol, runtime_metadata)


def one_dimensional(A, h, diffusion, times):
    n = grid_shape([A], h)[0]
    x = -A+(np.arange(n)+.5)*h
    V = .5*x*x
    pi = np.exp(-V-logsumexp(-V))
    jump = np.diff(V)
    forward = 2*diffusion/h**2*expit(-jump)
    reverse = 2*diffusion/h**2*expit(jump)
    diag = np.zeros(n)
    diag[:-1] += forward
    diag[1:] += reverse
    off = -np.sqrt(forward*reverse)
    vals, Q = eigh_tridiagonal(diag, off)
    P = x[:, None]**np.arange(5)[None, :]
    weighted = np.sqrt(pi)[:, None]*P
    spectral = Q.T @ weighted
    moments = []
    for t in [0., *times]:
        moments.append(spectral.T @ (np.exp(-t*vals)[:, None]*spectral))
    first = pi @ P
    stat = diag*np.sqrt(pi)
    stat[:-1] += off*np.sqrt(pi[1:])
    stat[1:] += off*np.sqrt(pi[:-1])
    return moments, first, dict(n=n, matrix_entries=2*n-1,
                               dense_eigenvectors_bytes=Q.nbytes,
                               minimum_eigenvalue=float(vals[0]),
                               stationarity_residual=float(np.max(np.abs(stat))))


def tensorized(model, faces, h, sd):
    moments, first, diagnostics = [], [], []
    for k in range(3):
        m, f, d = one_dimensional(faces[k], h, model['mu']*model['eigenvalues'][k], model['times'])
        moments.append(m)
        first.append(f)
        diagnostics.append(d)
    covariances = []
    for at in range(4):
        cov, mean = contract_moments(model, [m[at] for m in moments], first)
        covariances.append(cov)
    norm = np.outer(sd, sd)
    G0, K = covariances[0]/norm, np.array(covariances[1:])/norm
    return dict(G0=G0, K=K, C=K-G0, sd=sd, mean=mean), diagnostics


def independent_hermite(model):
    """Five-node normal quadrature exactly integrates degree-eight products."""
    nodes, weights = roots_hermitenorm(5)
    weights /= np.sqrt(2*np.pi)
    inds = np.array(list(product(range(5), repeat=3)))
    x = nodes[inds]
    w = np.prod(weights[inds], axis=1)
    # Original geometric contact law, not the monomial coefficient evaluator.
    _, E = energies(model, x)
    powers = np.array([a for a in product(range(5), repeat=3) if sum(a) <= 4])
    He = np.ones((len(x), len(powers)))
    factorial = np.ones(len(powers))
    for k in range(3):
        He *= eval_hermitenorm(powers[None, :, k], x[:, None, k])
        factorial *= np.array([math.factorial(int(a)) for a in powers[:, k]])
    coefficients = E.T @ (w[:, None]*He)/factorial[None, :]
    centered = np.sum(powers, axis=1) != 0
    coeff = coefficients[:, centered]
    factor = factorial[centered]
    powers = powers[centered]
    cov = (coeff*factor) @ coeff.T
    sd = np.sqrt(np.diag(cov))
    delayed = []
    for t in model['times']:
        decay = np.exp(-t*model['mu']*(powers @ model['eigenvalues']))
        delayed.append((coeff*(factor*decay)) @ coeff.T)
    norm = np.outer(sd, sd)
    G0, K = cov/norm, np.array(delayed)/norm
    return dict(G0=G0, K=K, C=K-G0, sd=sd,
                mean=coefficients[:, np.sum(np.array([a for a in product(range(5), repeat=3) if sum(a) <= 4]), axis=1) == 0].ravel())


def max_errors(actual, expected):
    return {key: float(np.max(np.abs(actual[key]-expected[key]))) for key in ('G0', 'K', 'C')}


def small_full_grid(model, faces, h, sd):
    """Independent full 3D assembly and sparse exponential action for parity."""
    from scipy.sparse import coo_matrix
    from scipy.sparse.linalg import expm_multiply
    shape = grid_shape(faces, h)
    ix = np.indices(shape).reshape(3, -1).T
    x = -np.array(faces)+(ix+.5)*h
    _, E = energies(model, x)
    V = .5*np.sum(x*x, axis=1)
    pi = np.exp(-V-logsumexp(-V))
    F = np.sqrt(pi)[:, None]*(E-pi @ E)/sd[None, :]
    G = len(x)
    flat = np.arange(G).reshape(shape)
    rows, cols, vals = [], [], []
    diag = np.zeros(G)
    for k in range(3):
        lo, hi = [slice(None)]*3, [slice(None)]*3
        lo[k], hi[k] = slice(0, -1), slice(1, None)
        a, b = flat[tuple(lo)].ravel(), flat[tuple(hi)].ravel()
        diff = V[b]-V[a]
        ab = 2*model['mu']*model['eigenvalues'][k]/h**2*expit(-diff)
        ba = 2*model['mu']*model['eigenvalues'][k]/h**2*expit(diff)
        off = -np.sqrt(ab*ba)
        np.add.at(diag, a, ab)
        np.add.at(diag, b, ba)
        rows.extend([a, b]); cols.extend([b, a]); vals.extend([off, off])
    rows.append(np.arange(G)); cols.append(np.arange(G)); vals.append(diag)
    H = coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(G, G)).tocsr()
    G0 = F.T @ F
    K = np.array([F.T @ expm_multiply(-t*H, F, traceA=-t*diag.sum()) for t in model['times']])
    return dict(G0=G0, K=K, C=K-G0), dict(states=G, csr_bytes=H.data.nbytes+H.indices.nbytes+H.indptr.nbytes)


def run(output):
    started = time.perf_counter()
    output.mkdir(parents=True, exist_ok=False)
    p = protocol()
    np.random.seed(p['execution']['numpy_random_seed_per_worker'])
    runtime = runtime_metadata()
    model = build_model()
    analytic = analytic_gaussian(model)
    hermite = independent_hermite(model)
    exact_error = max_errors(hermite, analytic)
    scale_error = float(np.max(np.abs(hermite['sd']/analytic['sd']-1)))
    small = p['harmonic_milestone']
    full, full_meta = small_full_grid(model, small['small_full_grid_faces'], small['small_full_grid_spacing'], analytic['sd'])
    tensor_small, _ = tensorized(model, small['small_full_grid_faces'], small['small_full_grid_spacing'], analytic['sd'])
    parity_error = max_errors(tensor_small, full)
    np.savez_compressed(output/'analytic.npz', **analytic, eigenvalues=model['eigenvalues'], basis=model['basis'], times=model['times'])
    np.savez_compressed(output/'parity.npz', **{f'full_{k}': v for k, v in full.items()},
                        **{f'tensor_{k}': tensor_small[k] for k in ('G0', 'K', 'C')})
    cases, arrays = [], []
    for A, h in small['finite_calibration_cases']:
        at = time.perf_counter()
        result, dims = tensorized(model, [A]*3, h, analytic['sd'])
        errors = max_errors(result, analytic)
        name = f'harmonic-A{A}-h{h:g}'
        np.savez_compressed(output/(name+'.npz'), **result)
        arrays.append(result)
        cases.append(dict(name=name, faces=[A]*3, spacing=h,
                          equivalent_full_states=math.prod(grid_shape([A]*3, h)),
                          full_grid_allocated=False, one_dimensional= dims,
                          seconds=time.perf_counter()-at, errors_against_analytic=errors,
                          all_components_within_0_001=all(v <= .001 for v in errors.values())))
    comparisons = []
    for i, j, kind in [(0, 1, 'old_grid'), (1, 2, 'old_grid'), (2, 3, 'old_domain'), (3, 4, 'old_domain'), (4, 5, 'fine_grid'), (5, 6, 'fine_grid'), (6, 7, 'fine_grid'), (7, 8, 'fine_domain'), (8, 9, 'fine_domain')]:
        errors = max_errors(arrays[j], arrays[i])
        comparisons.append(dict(first=cases[i]['name'], second=cases[j]['name'], kind=kind,
                                errors=errors, all_components_within_0_001=all(v <= .001 for v in errors.values())))
    elapsed = time.perf_counter()-started
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform == 'darwin' else 1024)
    accuracy = max([*exact_error.values(), scale_error, *parity_error.values()]) <= 1e-10
    receipt = dict(status='COMPLETE', runtime=runtime, numpy_random_seed=p['execution']['numpy_random_seed_per_worker'], protocol_sha256=hashlib.sha256((ROOT/'protocol.json').read_bytes()).hexdigest(),
                   source_hashes={name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in ['model.py', 'calibration.py']},
                   exact_method_difference=exact_error, relative_scale_difference=scale_error,
                   full_small_grid=full_meta, tensor_full_difference=parity_error,
                   cases=cases, comparisons=comparisons, wall_seconds=elapsed,
                   worker_peak_rss_bytes=peak,
                   harmonic_implementation_milestone_pass=accuracy and elapsed <= 60 and peak <= 512_000_000,
                   scientific_scope='Analytic harmonic implementation and finite harmonic calibration only; nonlinear reference remains unresolved.')
    (output/'receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps({k: receipt[k] for k in ['exact_method_difference', 'tensor_full_difference', 'wall_seconds', 'worker_peak_rss_bytes', 'harmonic_implementation_milestone_pass']}, indent=2))
    if not accuracy:
        raise AssertionError('Harmonic algebra or full-grid parity failed; full receipt preserved.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    run(args.output)
