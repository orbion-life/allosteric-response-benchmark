"""Independent polynomial Wick, quadrature and derivative checks; fixed seeds."""
from collections import defaultdict
from functools import lru_cache
from itertools import product
from pathlib import Path
import argparse
import json
import time
import numpy as np
from scipy.linalg import eigh
from numpy.polynomial.hermite import hermgauss
from gaussian_protein import CovarianceObjective, newton_cg, contracted_covariance, host, moments_for_covariance, backend


def fixture():
    r0 = np.array([[0., 0., 0.], [4., 0., 0.], [1., 3.5, 0.], [.5, 1.2, 4.2]])
    edges = np.array([(i, j) for i in range(4) for j in range(i+1, 4)])
    H = np.zeros((12, 12))
    for i, j in edges:
        a = r0[i]-r0[j]
        T = np.zeros((3, 12)); T[:, 3*i:3*i+3] = np.eye(3); T[:, 3*j:3*j+3] = -np.eye(3)
        H += T.T @ np.outer(a, a) @ T/(a @ a)
    lam, B = eigh(H); positive = lam > 1e-8
    return {'r0': r0, 'edges': edges, 'B': B[:, positive], 'eigenvalues': lam[positive]}


def polynomial(a, offset):
    g = defaultdict(float)
    for i in range(3):
        p = [0]*6; p[offset+i] = 1; g[tuple(p)] += 2*a[i]
        p = [0]*6; p[offset+i] = 2; g[tuple(p)] += 1
    result = defaultdict(float)
    for p, c in g.items():
        for q, d in g.items():
            result[tuple(x+y for x, y in zip(p, q))] += c*d
    return result


def wick_covariance(a, b, joint):
    @lru_cache(None)
    def moment(ids):
        if not ids: return 1.
        if len(ids) % 2: return 0.
        return sum(joint[ids[0], ids[j]]*moment(ids[1:j]+ids[j+1:]) for j in range(1, len(ids)))
    def evaluate(p):
        return moment(tuple(i for i, n in enumerate(p) for _ in range(n)))
    P, Q = polynomial(a, 0), polynomial(b, 3)
    meanp = sum(c*evaluate(p) for p, c in P.items())
    meanq = sum(c*evaluate(p) for p, c in Q.items())
    product_mean = sum(c*d*evaluate(tuple(x+y for x, y in zip(p, q))) for p, c in P.items() for q, d in Q.items())
    return product_mean-meanp*meanq


def independent_energy(Y, obj):
    # Five nodes per six-dimensional fixture coordinate integrates degree 4 exactly.
    nodes, weights = hermgauss(5)
    indices = np.array(list(product(range(5), repeat=obj.d)))
    q = np.einsum('si,ji->sj', np.sqrt(2)*nodes[indices], np.linalg.cholesky(Y), optimize=False)
    probability = np.prod(weights[indices], axis=1)/np.pi**(obj.d/2)
    rr = obj.r0[None]+np.einsum('si,ji->sj', q, obj.W, optimize=False).reshape(len(q), obj.n, 3)
    U = np.zeros(len(q))
    for i, j in obj.edges:
        l2 = np.sum((obj.r0[i]-obj.r0[j])**2)
        dr2 = np.sum((rr[:, i]-rr[:, j])**2, axis=1)
        U += (dr2-l2)**2/(8*l2)
    return float(np.sum(probability*U))


def checks(xp=np):
    rng = np.random.default_rng(1729)
    model = fixture()
    obj = CovarianceObjective(**{k: model[k] for k in ['r0', 'edges', 'B']}, lam=model['eigenvalues'])
    A = rng.normal(size=(6, 6)); Y = .2*np.eye(6)+.05*(A @ A.T)
    value, G, inverse, _, _ = obj.evaluate(Y)
    energy = independent_energy(Y, obj)
    energy_error = abs(energy-(value+.5*np.linalg.slogdet(Y)[1]))
    gradient_errors, hessian_errors = [], []
    epsilon = 1e-5
    for _ in range(8):
        D = rng.normal(size=(6, 6)); D = (D+D.T)/2; D /= np.linalg.norm(D)
        plus = obj.evaluate(Y+epsilon*D); minus = obj.evaluate(Y-epsilon*D)
        gradient_errors.append(abs((plus[0]-minus[0])/(2*epsilon)-np.sum(G*D)))
        hessian_errors.append(float(np.max(np.abs((plus[1]-minus[1])/(2*epsilon)-obj.hessian_action(D, inverse)))))
    errors = []
    rotation_errors = []
    for _ in range(12):
        A = rng.normal(size=(6, 6)); joint = .02*(A @ A.T)+.05*np.eye(6)
        a, b = rng.normal(size=(2, 3))
        exact = wick_covariance(a, b, joint)
        result = contracted_covariance(a, b, joint[:3, :3], joint[3:, 3:], joint[:3, 3:])
        errors.append(abs(float(result)-exact))
        R, _ = np.linalg.qr(rng.normal(size=(3, 3)))
        rotated = contracted_covariance(R @ a, R @ b, R @ joint[:3, :3] @ R.T,
                                         R @ joint[3:, 3:] @ R.T, R @ joint[:3, 3:] @ R.T)
        rotation_errors.append(abs(float(result-rotated)))
    harmonic = CovarianceObjective(model['r0'], model['edges'], model['B'], model['eigenvalues'], gamma=0.)
    HY, HF = newton_cg(harmonic, initial_scale=.7)
    Yfit, fit = newton_cg(obj)
    otherY, otherfit = newton_cg(obj, initial_scale=.25)
    sigma = Yfit/np.sqrt(model['eigenvalues'][:, None]*model['eigenvalues'][None, :])
    times = np.array([.1, 1., 10.])/model['eigenvalues'][0]
    covs, means, _ = moments_for_covariance(model, sigma, times)
    record = {'quadrature_expected_energy_max_abs': energy_error,
              'gradient_finite_difference_max_abs': max(gradient_errors),
              'hessian_action_finite_difference_max_abs': max(hessian_errors),
              'arbitrary_joint_gaussian_wick_max_abs': max(errors),
              'joint_coordinate_rotation_max_abs': max(rotation_errors),
              'harmonic_limit_covariance_max_abs': float(np.max(np.abs(HY-np.eye(6)))),
              'two_initializations_covariance_max_abs': float(np.max(np.abs(Yfit-otherY))),
              'final_gradient_frobenius': fit['history'][-1]['whitened_gradient_frobenius']}
    record['passed'] = (energy_error <= 1e-10 and max(gradient_errors) <= 1e-8 and max(hessian_errors) <= 1e-8
                        and max(errors) <= 1e-10 and max(rotation_errors) <= 1e-10
                        and record['harmonic_limit_covariance_max_abs'] <= 1e-10
                        and record['two_initializations_covariance_max_abs'] <= 1e-9)
    arrays = {'Y': Yfit, 'Sigma': sigma, 'covariances': covs, 'means': means}
    if xp is not np:
        gpuobj = CovarianceObjective(model['r0'], model['edges'], model['B'], model['eigenvalues'], xp=xp)
        GY, GF = newton_cg(gpuobj)
        Gl = xp.asarray(model['eigenvalues']); GS = GY/xp.sqrt(Gl[:, None]*Gl[None, :])
        Gcov, Gmeans, _ = moments_for_covariance(model, GS, times, xp)
        gpuarrays = {'Y': host(GY), 'Sigma': host(GS), 'covariances': host(Gcov), 'means': host(Gmeans)}
        record['cpu_gpu_max_abs'] = {key: float(np.max(np.abs(arrays[key]-gpuarrays[key]))) for key in arrays}
        record['cpu_gpu_passed'] = max(record['cpu_gpu_max_abs'].values()) <= 1e-9
        record['gpu_fit_gradient_frobenius'] = GF['history'][-1]['whitened_gradient_frobenius']
        arrays.update({'gpu_'+key: value for key, value in gpuarrays.items()})
    return record, arrays


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--output', required=True); parser.add_argument('--backend', default='numpy')
    args = parser.parse_args(); out = Path(args.output); out.mkdir(parents=True, exist_ok=True)
    start = time.perf_counter(); result, arrays = checks(backend(args.backend)); result['seconds'] = time.perf_counter()-start
    (out/'independent-verification.json').write_text(json.dumps(result, indent=2))
    np.savez_compressed(out/'independent-check-arrays.npz', **arrays)
    print(json.dumps(result, indent=2))
    if not result['passed'] or result.get('cpu_gpu_passed') is False: raise SystemExit(1)
