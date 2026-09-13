"""Full nonlinear coarse/mid-grid costs; no continuum-accuracy claim."""
import argparse
import hashlib
import json
from pathlib import Path
import resource
import sys
import time
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import LinearOperator, expm_multiply
from scipy.special import expit, logsumexp
from model import ROOT, analytic_gaussian, build_model, energies, grid_shape, protocol, runtime_metadata
from chebyshev import action as chebyshev_action


class CountedOperator(LinearOperator):
    def __init__(self, H):
        self.H = H
        self.calls = dict(forward_vectors=0, forward_blocks=0, adjoint_vectors=0, adjoint_blocks=0, vector_equivalents=0)
        super().__init__(dtype=np.dtype(float), shape=H.shape)

    def action(self, x, label):
        self.calls[label] += 1
        self.calls['vector_equivalents'] += 1 if x.ndim == 1 else x.shape[1]
        return self.H @ x

    def _matvec(self, x):
        return self.action(x, 'forward_vectors')

    def _matmat(self, x):
        return self.action(x, 'forward_blocks')

    def _rmatvec(self, x):
        return self.action(x, 'adjoint_vectors')

    def _rmatmat(self, x):
        return self.action(x, 'adjoint_blocks')


class Stencil:
    def __init__(self, diag, links, shape):
        self.diag = diag.reshape(shape)
        self.links = links
        self.shape = shape

    def action(self, x):
        X = x.reshape((*self.shape, -1))
        out = self.diag[..., None]*X
        for k, off in enumerate(self.links):
            lo, hi = [slice(None)]*3, [slice(None)]*3
            lo[k], hi[k] = slice(0, -1), slice(1, None)
            out[tuple(lo)] += off[..., None]*X[tuple(hi)]
            out[tuple(hi)] += off[..., None]*X[tuple(lo)]
        return out.reshape(x.shape)


def build_grid(model, faces, h, sd):
    shape = grid_shape(faces, h)
    G = int(np.prod(shape))
    idx = np.indices(shape).reshape(3, -1).T
    x = -np.array(faces)+(idx+.5)*h
    U, E = energies(model, x)
    V = model['beta']*U
    pi = np.exp(-V-logsumexp(-V))
    F = np.sqrt(pi)[:, None]*(E-pi @ E)/sd[None, :]
    boundary = np.any((idx == 0) | (idx == np.array(shape)-1), axis=1)
    boundary_mass = float(pi[boundary].sum())
    del x, idx, E
    return shape, V, pi, F, boundary_mass


def assemble(model, shape, V, h):
    G = len(V)
    flat = np.arange(G).reshape(shape)
    diag = np.zeros(G)
    rows, cols, values, links = [], [], [], []
    for k in range(3):
        lo, hi = [slice(None)]*3, [slice(None)]*3
        lo[k], hi[k] = slice(0, -1), slice(1, None)
        a, b = flat[tuple(lo)].ravel(), flat[tuple(hi)].ravel()
        delta = V[b]-V[a]
        ab = 2*model['mu']*model['eigenvalues'][k]/h**2*expit(-delta)
        ba = 2*model['mu']*model['eigenvalues'][k]/h**2*expit(delta)
        off = -np.sqrt(ab*ba)
        np.add.at(diag, a, ab)
        np.add.at(diag, b, ba)
        linkshape = list(shape); linkshape[k] -= 1
        links.append(off.reshape(linkshape))
        rows.extend([a, b]); cols.extend([b, a]); values.extend([off, off])
    rows.append(np.arange(G)); cols.append(np.arange(G)); values.append(diag)
    H = coo_matrix((np.concatenate(values), (np.concatenate(rows), np.concatenate(cols))), shape=(G, G)).tocsr()
    return H, Stencil(diag, links, shape), diag


def run(output, faces, h):
    started = time.perf_counter()
    output.mkdir(parents=True, exist_ok=False)
    np.random.seed(protocol()['execution']['numpy_random_seed_per_worker'])
    runtime = runtime_metadata()
    timestamps = {}
    model = build_model()
    exact = analytic_gaussian(model)
    timestamps['native_model_and_exact_normalizer_seconds'] = time.perf_counter()-started
    at = time.perf_counter()
    shape, V, pi, F, boundary_mass = build_grid(model, faces, h, exact['sd'])
    timestamps['grid_equilibrium_observable_seconds'] = time.perf_counter()-at
    at = time.perf_counter()
    H, stencil, diag = assemble(model, shape, V, h)
    timestamps['sparse_and_stencil_assembly_seconds'] = time.perf_counter()-at
    # Separate action parity and timing. Both representations remain live here;
    # the measured profile peak is not a pure matrix-free peak.
    csr_action = H @ F
    stencil_action = stencil.action(F)
    relative_error = float(np.max(np.abs(csr_action-stencil_action))/max(1., np.max(np.abs(csr_action))))
    bench = {}
    for name, action in [('csr', lambda x: H @ x), ('matrix_free_stencil', stencil.action)]:
        at = time.perf_counter()
        for _ in range(10):
            result = action(F)
        bench[name+'_ten_three_column_actions_seconds'] = time.perf_counter()-at
        if not np.isfinite(result).all():
            raise FloatingPointError('Nonfinite matrix action.')
    stationarity = float(np.max(np.abs(H @ np.sqrt(pi))))
    if stationarity > 1e-9 or relative_error > 1e-12:
        raise AssertionError(f'Operator identity failure: stationary={stationarity}; action={relative_error}.')
    G0 = F.T @ F
    op = CountedOperator(H)
    propagation = []
    K = []
    for t in model['times']:
        before = dict(op.calls)
        at = time.perf_counter()
        evolved = expm_multiply(-t*op, F, traceA=-t*diag.sum())
        K.append(F.T @ evolved)
        propagation.append(dict(time=float(t), time_over_tau=float(t*model['mu']*model['eigenvalues'][0]),
                                seconds=time.perf_counter()-at,
                                counts={key: op.calls[key]-before[key] for key in before}))
        # A checkpoint preserves completed query diagnostics if the watchdog stops.
        (output/'progress.json').write_text(json.dumps(dict(status='RUNNING', completed_queries=propagation), indent=2)+'\n')
    K = np.array(K)
    upper = float(np.max(np.asarray(abs(H).sum(axis=1)).ravel()))
    independent = []
    K_cheb = []
    for t in model['times']:
        evolved, meta = chebyshev_action(H, F, t, upper)
        K_cheb.append(F.T @ evolved)
        meta['time'] = float(t)
        independent.append(meta)
    K_cheb = np.array(K_cheb)
    independent_error = float(np.max(np.abs(K_cheb-K)))
    at = time.perf_counter()
    np.savez_compressed(output/'response.npz', G0=G0, K=K, C=K-G0,
                        K_chebyshev=K_cheb, C_chebyshev=K_cheb-G0, sd=exact['sd'],
                        times=model['times'], eigenvalues=model['eigenvalues'], basis=model['basis'])
    timestamps['array_serialization_seconds'] = time.perf_counter()-at
    elapsed = time.perf_counter()-started
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform == 'darwin' else 1024)
    receipt = dict(status='COMPLETE', runtime=runtime, numpy_random_seed=protocol()['execution']['numpy_random_seed_per_worker'], protocol_sha256=hashlib.sha256((ROOT/'protocol.json').read_bytes()).hexdigest(),
                   source_hashes={name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in ['model.py', 'nonlinear_profile.py', 'chebyshev.py']},
                   faces=faces, spacing=h, shape=list(shape), states=len(V),
                   full_positive_coordinates=3, maximum_energy=float(V.max()),
                   probability_underflow_nodes=int(np.sum(pi == 0)), boundary_mass=boundary_mass,
                   stationarity_residual=stationarity, csr_stencil_scaled_action_error=relative_error,
                   csr_stencil_raw_max_action_difference=float(np.max(np.abs(csr_action-stencil_action))),
                   action_error_definition='maxabs(CSR(F)-stencil(F))/max(1,maxabs(CSR(F))); mixed absolute/relative, not a pure relative error.',
                   assembly_and_preprocessing= timestamps, action_microbenchmark=bench,
                   propagation=propagation, propagation_counts_total=op.calls,
                   independent_chebyshev=independent,
                   independent_max_normalized_delayed_response_difference=independent_error,
                   independent_comparison_pass=independent_error <= 1e-8,
                   actual_csr_nonzeros=H.nnz,
                   actual_csr_bytes=H.data.nbytes+H.indices.nbytes+H.indptr.nbytes,
                   stored_matrix_free_stencil_bytes=stencil.diag.nbytes+sum(x.nbytes for x in stencil.links),
                   three_observable_columns_bytes=F.nbytes,
                   wall_seconds=elapsed, worker_peak_rss_bytes=peak,
                   wall_seconds_scope='Worker timing from run entry through scientific array serialization; excludes interpreter/import startup, final JSON receipt formatting/write and shutdown. Use watchdog outer_worker_wall_seconds for the complete worker invocation.',
                   measurement_scope='Fresh process; input/model/Hessian/normalizers, full equilibrium and generator construction, CSR and stencil action check/benchmark, three full CSR Taylor queries, three independent Chebyshev queries and array serialization. Both operators are resident; this is not a measured matrix-free-only peak. Excludes dependency installation and response convergence checks.',
                   scientific_scope='Full nonlinear finite model on this stated box; no whole-plane convergence or domain-coverage acceptance.')
    (output/'receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps({k: receipt[k] for k in ['states', 'wall_seconds', 'worker_peak_rss_bytes', 'propagation_counts_total']}, indent=2))
    if independent_error > 1e-8:
        raise AssertionError('Independent propagation comparison failed; full receipt retained.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--faces', type=int, nargs=3, default=[4, 4, 17])
    parser.add_argument('--spacing', type=float, choices=[.5, .25], required=True)
    args = parser.parse_args()
    run(args.output, args.faces, args.spacing)
