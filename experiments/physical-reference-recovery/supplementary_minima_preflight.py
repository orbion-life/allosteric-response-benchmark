"""Additional independently reconstructed congruent minima; no grid allocation."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.linalg import eigh
from model import ROOT, analytic_gaussian, build_model, contract, energies, protocol
from preflight import grid_record


def minima(model):
    centered = model['r0']-model['r0'].mean(0)
    covariance2 = contract(centered[:, :2].T, centered[:, :2])
    _, principal = eigh(covariance2)
    axes = [np.array([0., 0., 1.]), *[np.r_[u, 0.] for u in principal.T]]
    out = []
    for u in axes:
        rotation = 2*np.outer(u, u)-np.eye(3)
        displacement = contract(centered, rotation.T)-centered
        q = contract(model['basis'].T, displacement.ravel())
        reconstructed = contract(model['basis'], q).reshape(3, 3)
        x = np.sqrt(model['beta']*model['eigenvalues'])*q
        U, _ = energies(model, x[None])
        error = float(np.max(np.abs(reconstructed-displacement)))
        if error > 1e-12 or U[0] > 1e-20:
            raise AssertionError('Constructed congruent minimum fails direct reconstruction.')
        out.append(dict(axis=u.tolist(), q=q.tolist(), x=x.tolist(),
                        projection_residual=error, contact_energy=float(U[0])))
    return out


def run():
    p = protocol()
    cases = []
    for name, geometry in dict(development=p['geometry_A'], **p['other_geometries_A']).items():
        for kappa in (1, 10, 100):
            model = build_model(geometry, kappa)
            records = minima(model)
            faces = np.ceil(np.maximum(4, np.max(np.abs([r['x'] for r in records]), axis=0)+2)).astype(int)
            grids = [grid_record(faces.tolist(), h) for h in (.5, .25, .125)]
            grids += [grid_record((faces+a).tolist(), .125) for a in (1, 2)]
            sd = analytic_gaussian(model)['sd']
            for grid in grids:
                A = np.array(grid['faces'], dtype=float)
                maxE = contract(np.abs(model['coefficients']), np.prod(A[None, :]**model['powers'], axis=1))
                grid['observable_absolute_upper_bound_by_polynomial'] = maxE.tolist()
                grid['variance_normalization_amplification_upper_bound'] = float(np.max((maxE/(2*sd))**2))
                grid['maximum_rate_sum_bound'] = float(4*model['mu']*model['eigenvalues'].sum()/grid['spacing']**2)
                grid['spectral_upper_bound_by_gershgorin_exit_bound'] = 2*grid['maximum_rate_sum_bound']
                grid['rate_bound_note'] = 'Exit <=4 sum(mu lambda)/h^2, so spectral upper <=8 sum(mu lambda)/h^2. Polynomial nonnegative E bound gives Var(E)<=max(E)^2/4; no unknown equilibrium variance is assumed.'
            cases.append(dict(geometry=name, kappa=kappa, all_three_constructed_minima=records,
                              initial_faces=faces.tolist(), nonlinear_grids=grids))
    return dict(status='SUPPLEMENTARY STATIC PREFLIGHT; NO NEW RESPONSE CASES',
                interpretation='These extra congruent minima were identified after the original cost-probe protocol froze. The original two cost probes remain unchanged; all later coverage-aware allocation must use this expanded inventory. This is still not proof that all minima have been found.',
                protocol_sha256=hashlib.sha256((ROOT/'protocol.json').read_bytes()).hexdigest(),
                cases=cases)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.write_text(json.dumps(run(), indent=2)+'\n')
