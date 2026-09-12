"""Exploratory 3x3 convex-closure algebra check, not a protein or response fit.

This checks the objective, derivatives, positive definiteness, harmonic limit and
agreement of deterministic Newton solves. It evaluates no nonlinear reference,
contact label, receiver score, top-five ranking or quantum algorithm.
"""
from pathlib import Path
from itertools import product
import hashlib
import json
import platform
import time

import numpy as np
import scipy
from scipy.linalg import eigh
from numpy.polynomial.hermite import hermgauss

ROOT = Path(__file__).resolve().parent


def setup():
    points = np.array([[0., 0., 0.], [4., 0., 0.], [1., 3.5, 0.]])
    edges = [(0, 1), (0, 2), (1, 2)]
    hessian = np.zeros((9, 9))
    for i, j in edges:
        a = points[i] - points[j]
        t = np.zeros((3, 9))
        t[:, 3*i:3*i+3] = np.eye(3)
        t[:, 3*j:3*j+3] = -np.eye(3)
        hessian += t.T @ np.outer(a, a) @ t / np.dot(a, a)
    values, vectors = eigh(hessian)
    basis = vectors[:, values > 1e-8]
    assert basis.shape == (9, 3)
    K = basis.T @ hessian @ basis
    terms = []
    for i, j in edges:
        a = points[i] - points[j]
        T = basis[3*i:3*i+3] - basis[3*j:3*j+3]
        terms.append((a, T, 1 / (8 * np.dot(a, a))))
    return K, terms


def objective(C, K, terms, gamma=1., beta=1.):
    sign, logdet = np.linalg.slogdet(C)
    if sign <= 0 or np.linalg.eigvalsh(C)[0] <= 0:
        return np.inf
    value = .5 * np.trace(K @ C) - logdet / (2 * beta)
    for _, T, coefficient in terms:
        S = T @ C @ T.T
        value += gamma * coefficient * (2 * np.trace(S @ S) + np.trace(S)**2)
    return float(value)


def gradient(C, K, terms, gamma=1., beta=1.):
    G = .5 * K - np.linalg.inv(C) / (2 * beta)
    for _, T, coefficient in terms:
        S = T @ C @ T.T
        G += gamma * coefficient * T.T @ (4*S + 2*np.trace(S)*np.eye(3)) @ T
    return G


def hessian_action(C, D, terms, gamma=1., beta=1.):
    inverse = np.linalg.inv(C)
    answer = inverse @ D @ inverse / (2 * beta)
    for _, T, coefficient in terms:
        change = T @ D @ T.T
        answer += gamma * coefficient * T.T @ (4*change + 2*np.trace(change)*np.eye(3)) @ T
    return answer


def symmetric_basis(d):
    matrices = []
    for i in range(d):
        for j in range(i, d):
            value = np.zeros((d, d))
            value[i, j] = value[j, i] = 1. if i == j else 1 / np.sqrt(2)
            matrices.append(value)
    return np.array(matrices)


def optimize(initial, K, terms, gamma):
    C = initial.copy()
    basis = symmetric_basis(len(C))
    history = []
    for iteration in range(100):
        G = gradient(C, K, terms, gamma)
        g = np.einsum('aij,ij->a', basis, G)
        value = objective(C, K, terms, gamma)
        history.append({"iteration": iteration, "F_without_constant": value, "gradient_frobenius": float(np.linalg.norm(G))})
        if np.linalg.norm(G) < 1e-11:
            return C, history
        H = np.array([[np.sum(a * hessian_action(C, b, terms, gamma)) for b in basis] for a in basis])
        step = -np.linalg.solve(H, g)
        D = np.einsum('a,aij->ij', step, basis)
        slope = float(np.sum(G * D))
        alpha = 1.
        while (np.linalg.eigvalsh(C + alpha*D)[0] <= 0 or
               objective(C + alpha*D, K, terms, gamma) > value + 1e-4*alpha*slope):
            alpha *= .5
            if alpha < 2**-40:
                raise RuntimeError("Backtracking failed.")
        C = (C + alpha*D + (C + alpha*D).T) / 2
    raise RuntimeError("Stationarity tolerance not reached.")


def gaussian_energy_quadrature(C, terms):
    # Five-point Gauss-Hermite per axis integrates this quartic exactly up to rounding.
    nodes, weights = hermgauss(5)
    L = np.linalg.cholesky(C)
    value = 0.
    for indices in product(range(5), repeat=3):
        q = np.sqrt(2) * L @ nodes[list(indices)]
        weight = np.prod(weights[list(indices)]) / np.pi**1.5
        U = sum(coefficient * (2*np.dot(a, T @ q) + np.dot(T @ q, T @ q))**2
                for a, T, coefficient in terms)
        value += weight * U
    return float(value)


def main():
    start = time.monotonic()
    K, terms = setup()
    harmonic = np.linalg.inv(K)
    records = []
    optima = {}
    for gamma in [0., 1.]:
        solutions = []
        for scale in [.25, 1., 4.]:
            C, history = optimize(scale*harmonic, K, terms, gamma)
            solutions.append(C)
            records.append({"quartic_multiplier": gamma, "initial_covariance_scale": scale,
                            "covariance": C.tolist(), "history": history,
                            "minimum_covariance_eigenvalue": float(np.linalg.eigvalsh(C)[0])})
        difference = float(max(np.max(np.abs(C-solutions[0])) for C in solutions))
        assert difference < 1e-10
        optima[str(gamma)] = solutions[0]
    C = optima['1.0']
    harmonic_error = float(np.max(np.abs(optima['0.0']-harmonic)))
    quadrature = gaussian_energy_quadrature(C, terms)
    analytic_energy = objective(C, K, terms) + .5*np.linalg.slogdet(C)[1]
    gradient_matrix = gradient(C, K, terms)
    stationarity = float(np.linalg.norm(2*gradient_matrix))
    basis = symmetric_basis(3)
    epsilon = 1e-5
    gradient_errors, hessian_errors = [], []
    for D in basis:
        derivative = (objective(C+epsilon*D, K, terms)-objective(C-epsilon*D, K, terms))/(2*epsilon)
        gradient_errors.append(abs(derivative-np.sum(gradient_matrix*D)))
        numerical = (gradient(C+epsilon*D, K, terms)-gradient(C-epsilon*D, K, terms))/(2*epsilon)
        hessian_errors.append(float(np.max(np.abs(numerical-hessian_action(C, D, terms)))))
    H = np.array([[np.sum(a*hessian_action(C,b,terms)) for b in basis] for a in basis])
    centroid_gradient = sum(4*coefficient*T.T @ (np.trace(T@C@T.T)*a + 2*(T@C@T.T)@a)
                            for a,T,coefficient in terms)
    result = {"status": "exploratory algebra check only", "date": "2026-09-12",
              "scope": "One three-coordinate triangle fixture; no nonlinear reference or delayed-response accuracy was evaluated.",
              "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "geometry_A": [[0,0,0],[4,0,0],[1,3.5,0]], "beta": 1, "kappa": 1,
              "K": K.tolist(), "harmonic_covariance": harmonic.tolist(), "variational_covariance": C.tolist(),
              "records": records, "harmonic_limit_max_covariance_error": harmonic_error,
              "analytic_expected_U": analytic_energy, "Gauss_Hermite_expected_U": quadrature,
              "quadrature_absolute_error": abs(analytic_energy-quadrature),
              "stationarity_Frobenius": stationarity,
              "minimum_symmetric_parameter_Hessian_eigenvalue": float(np.linalg.eigvalsh(H)[0]),
              "minimum_eigenvalue_harmonic_minus_variational_covariance": float(np.linalg.eigvalsh(harmonic-C)[0]),
              "gradient_finite_difference_max_error": max(gradient_errors),
              "Hessian_finite_difference_max_error": max(hessian_errors),
              "unrelaxed_centroid_gradient": centroid_gradient.tolist(),
              "unrelaxed_centroid_gradient_norm": float(np.linalg.norm(centroid_gradient)),
              "free_energy_difference_to_harmonic_trial_same_U": objective(C,K,terms)-objective(harmonic,K,terms),
              "runtime_seconds": time.monotonic()-start,
              "software": {"python": platform.python_version(),"numpy":np.__version__,"scipy":scipy.__version__}}
    assert harmonic_error < 1e-10 and stationarity < 1e-10
    assert result['quadrature_absolute_error'] < 1e-12
    assert result['gradient_finite_difference_max_error'] < 1e-8
    assert result['Hessian_finite_difference_max_error'] < 1e-8
    assert result['minimum_symmetric_parameter_Hessian_eigenvalue'] > 0
    assert result['minimum_eigenvalue_harmonic_minus_variational_covariance'] >= -1e-10
    (ROOT/'algebra-check.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k.endswith('error') or k in ['status','stationarity_Frobenius','unrelaxed_centroid_gradient_norm','runtime_seconds']}))


if __name__ == '__main__':
    main()
