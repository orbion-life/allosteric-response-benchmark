"""Static resource arithmetic only; never constructs a configuration grid."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import numpy as np
from model import ROOT, build_model, grid_shape, known_minimum, protocol


def grid_record(faces, spacing):
    shape = grid_shape(faces, spacing)
    G = math.prod(shape)
    links = sum((shape[k]-1)*math.prod(shape[j] for j in range(3) if j != k) for k in range(3))
    nnz = G+2*links
    return dict(faces=list(faces), spacing=spacing, shape=list(shape), states=G,
                undirected_links=links, csr_nonzeros=nnz,
                old_one_million_state_cap_pass=G <= 1_000_000,
                storage_bytes={
                    'csr_float64_int64_exact': 16*nnz+8*(G+1),
                    'matrix_free_diag_plus_three_link_fields_float64': 8*(G+links),
                    'three_observable_columns_float64': 24*G,
                    'one_state_vector_float64': 8*G,
                    'full_coordinates_plus_positions_float64_if_allocated': 96*G,
                    'forty_state_vector_planning_workspace_float64': 320*G,
                    'csr_plus_forty_vectors_planning_only': 16*nnz+8*(G+1)+320*G,
                    'matrix_free_plus_forty_vectors_planning_only': 8*(G+links)+320*G
                },
                maximum_rate_sum_bound=4*sum([1., 1., 1.])/spacing**2,
                rate_bound_note='Placeholder isotropic unit-diffusion value; actual sum(mu*lambda) multiplier is attached per model. Storage formulas exclude allocator, Python and library overhead; forty vectors are a declared planning scenario, not an observed peak.')


def run():
    p = protocol()
    geometries = dict(development=p['geometry_A'], **p['other_geometries_A'])
    cases = []
    for name, geom in geometries.items():
        for kappa in (1, 10, 100):
            model = build_model(geom, kappa)
            minimum = known_minimum(model)
            A = np.array(minimum['faces'])
            harmonic = [grid_record([4]*3, h) for h in (.5, .25, .125)]
            harmonic += [grid_record([a]*3, .125) for a in (5, 6)]
            nonlinear = [grid_record(A.tolist(), h) for h in (.5, .25, .125)]
            nonlinear += [grid_record((A+a).tolist(), .125) for a in (1, 2)]
            for r in harmonic+nonlinear:
                r['maximum_rate_sum_bound'] = float(4*model['mu']*model['eigenvalues'].sum()/r['spacing']**2)
                r['rate_bound_note'] = 'Actual per-model upper bound 4*sum(mu*lambda)/h^2 on exit-rate sum; a Gershgorin spectral bound is twice this. Arrays exclude Python/allocator/library overhead; forty-vector workspace is a planning scenario, not measured peak.'
            cases.append(dict(geometry=name, kappa=kappa, eigenvalues=model['eigenvalues'].tolist(),
                              minimum=minimum, harmonic=harmonic, nonlinear=nonlinear))
    return dict(status='STATIC PREFLIGHT ONLY; NO GRID OR PROPAGATION ALLOCATED',
                protocol_sha256=hashlib.sha256((ROOT/'protocol.json').read_bytes()).hexdigest(),
                total_harmonic_requirements=45, total_nonlinear_requirements=45,
                cases=cases)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    receipt = run()
    args.output.write_text(json.dumps(receipt, indent=2)+'\n')
    dev = receipt['cases'][0]
    print(json.dumps({'development_faces': dev['minimum']['faces'], 'nonlinear_states': [x['states'] for x in dev['nonlinear']]}, indent=2))
