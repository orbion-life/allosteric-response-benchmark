"""Largest-grid construction and action-only capacity; no propagation."""
import argparse
import hashlib
import json
from pathlib import Path
import resource
import sys
import time
import numpy as np
from model import ROOT, analytic_gaussian, build_model, runtime_metadata
from nonlinear_profile import assemble, build_grid
from preflight import grid_record


def run(output):
    output.mkdir(parents=True, exist_ok=False)
    p = json.loads((ROOT/'capacity-protocol.json').read_text())
    np.random.seed(p['numpy_random_seed'])
    started = time.perf_counter()
    runtime = runtime_metadata()
    model = build_model()
    exact = analytic_gaussian(model)
    at = time.perf_counter()
    shape, V, pi, F, boundary = build_grid(model, p['faces'], p['spacing'], exact['sd'])
    grid_seconds = time.perf_counter()-at
    if list(shape) != p['shape'] or len(V) != p['states']:
        raise AssertionError('Grid differs from predeclared capacity probe.')
    at = time.perf_counter()
    H, stencil, diag = assemble(model, shape, V, p['spacing'])
    assembly_seconds = time.perf_counter()-at
    preflight = grid_record(p['faces'], p['spacing'])
    if H.nnz != preflight['csr_nonzeros']:
        raise AssertionError('Sparse graph count differs from preflight.')
    at = time.perf_counter()
    workspace = np.empty((len(V), p['workspace_state_vectors']), dtype=float)
    workspace.fill(0.)  # Touch every page so RSS measures capacity, not reservation.
    workspace_seconds = time.perf_counter()-at
    csr = H @ F
    free = stencil.action(F)
    action_error = float(np.max(np.abs(csr-free))/max(1., np.max(np.abs(csr))))
    stationarity = float(np.max(np.abs(H @ np.sqrt(pi))))
    measurements = {}
    for name, action in [('csr', lambda x: H @ x), ('stencil', stencil.action)]:
        samples = []
        for _ in range(p['microbenchmark_actions_per_representation']):
            at = time.perf_counter()
            result = action(F)
            samples.append(time.perf_counter()-at)
        if not np.isfinite(result).all():
            raise FloatingPointError('Nonfinite capacity action.')
        measurements[name] = dict(three_column_action_seconds=samples, minimum=min(samples),
                                  median=float(np.median(samples)), maximum=max(samples))
    spectral_upper = float(np.max(np.asarray(abs(H).sum(axis=1)).ravel()))
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform == 'darwin' else 1024)
    receipt = dict(status='COMPLETE', capacity_protocol_sha256=hashlib.sha256((ROOT/'capacity-protocol.json').read_bytes()).hexdigest(),
                   source_hashes={name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in ['model.py','nonlinear_profile.py','capacity_probe.py']},
                   runtime=runtime, faces=p['faces'], spacing=p['spacing'], shape=list(shape), states=len(V),
                   grid_seconds=grid_seconds, assembly_seconds=assembly_seconds,
                   workspace_allocation_and_touch_seconds=workspace_seconds,
                   workspace_bytes=workspace.nbytes, workspace_state_vectors=p['workspace_state_vectors'],
                   actual_csr_entries=H.nnz, actual_csr_bytes=H.data.nbytes+H.indices.nbytes+H.indptr.nbytes,
                   stored_stencil_bytes=diag.nbytes+sum(x.nbytes for x in stencil.links),
                   three_column_bytes=F.nbytes, boundary_probability_diagnostic=boundary,
                   probability_underflow_nodes=int(np.sum(pi == 0)), maximum_energy=float(V.max()),
                   matrix_action_measurements=measurements, scaled_action_error=action_error,
                   stationarity_residual=stationarity, spectral_gershgorin_upper_bound=spectral_upper,
                   observable_norms=np.linalg.norm(F, axis=0).tolist(),
                   worker_peak_rss_bytes=peak, worker_timed_seconds=time.perf_counter()-started,
                   identities_pass=stationarity <= 1e-9 and action_error <= 1e-12,
                   scope='Construction, explicitly touched forty-vector workspace, identity checks, and10actions per representation only. No matrix exponential or nonlinear response was calculated. The spectral bound and norm are planning inputs on this finite grid.')
    (output/'receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps(receipt, indent=2))
    if not receipt['identities_pass']:
        raise AssertionError('Capacity identity check failed; receipt retained.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();run(args.output)
