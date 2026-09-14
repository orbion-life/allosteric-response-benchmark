"""Full-coordinate fixed-centroid Gaussian variational fit and exact OU moments.

This is a classical Gaussian surrogate. It is not nonlinear-equilibrium,
nonlinear-response, biological, or quantum-advantage validation.
"""
from pathlib import Path
import argparse
import csv
import datetime
import hashlib
import json
import platform
import resource
import time
import numpy as np


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def backend(name):
    if name == 'numpy':
        return np
    if name == 'cupy':
        import cupy
        return cupy
    raise ValueError(name)


def host(x):
    return x.get() if hasattr(x, 'get') else np.asarray(x)


def scalar(x):
    return float(host(x))


def sync(xp):
    if xp.__name__ == 'cupy':
        xp.cuda.Stream.null.synchronize()


def mm(a, b, xp=np):
    # Direct BLAS calls avoid NumPy/Apple Accelerate stale floating-point flags.
    # No floating-point warnings are suppressed. Batched 3x3 products use an
    # explicit contraction, checked independently against Wick moments.
    if xp is np:
        if a.ndim == b.ndim == 2:
            from scipy.linalg.blas import dgemm, zgemm
            return (zgemm if np.iscomplexobj(a) or np.iscomplexobj(b) else dgemm)(1., a, b)
        return np.einsum('...ij,...jk->...ik', a, b, optimize=False)
    return xp.matmul(a, b)


class CovarianceObjective:
    """beta*F in harmonic-whitened covariance Y; no parameter Hessian formed."""
    def __init__(self, r0, edges, B, lam, xp=np, beta=1., gamma=1.):
        self.xp, self.beta, self.gamma = xp, beta, gamma
        self.n, self.d = len(r0), len(lam)
        self.r0, self.edges = xp.asarray(r0), xp.asarray(edges)
        self.B, self.lam = xp.asarray(B), xp.asarray(lam)
        self.W = self.B / xp.sqrt(beta*self.lam)[None, :]
        self.a = self.r0[self.edges[:, 0]] - self.r0[self.edges[:, 1]]
        self.c = 1/(8*xp.sum(self.a*self.a, axis=1))
        self.eye = xp.eye(self.d)
        self.eye3 = xp.eye(3)
        self.actions = 0

    def edge_covariances(self, Y):
        X = mm(mm(self.W, Y, self.xp), self.W.T, self.xp).reshape(self.n, 3, self.n, 3).transpose(0, 2, 1, 3)
        i, j = self.edges[:, 0], self.edges[:, 1]
        S = X[i, i]-X[i, j]-X[j, i]+X[j, j]
        return (S+S.transpose(0, 2, 1))/2

    def quartic_adjoint(self, S):
        xp = self.xp
        traces = xp.trace(S, axis1=1, axis2=2)
        V = self.beta*self.gamma*self.c[:, None, None]*(4*S+2*traces[:, None, None]*self.eye3)
        L = xp.zeros((self.n, self.n, 3, 3))
        i, j = self.edges[:, 0], self.edges[:, 1]
        xp.add.at(L, (i, i), V)
        xp.add.at(L, (j, j), V)
        xp.add.at(L, (i, j), -V)
        xp.add.at(L, (j, i), -V)
        cart = L.transpose(0, 2, 1, 3).reshape(3*self.n, 3*self.n)
        out = mm(self.W.T, mm(cart, self.W, xp), xp)
        return (out+out.T)/2

    def evaluate(self, Y):
        xp = self.xp
        # Cholesky verifies positive definiteness, unlike determinant sign alone.
        L = xp.linalg.cholesky(Y)
        inverse = xp.linalg.solve(Y, self.eye)
        S = self.edge_covariances(Y)
        traces = xp.trace(S, axis1=1, axis2=2)
        quartic = self.beta*self.gamma*xp.sum(self.c*(2*xp.sum(S*S, axis=(1, 2))+traces**2))
        value = .5*xp.trace(Y)+quartic-xp.sum(xp.log(xp.diag(L)))
        G = .5*self.eye-.5*inverse+self.quartic_adjoint(S)
        return scalar(value), (G+G.T)/2, inverse, S, L

    def hessian_action(self, D, inverse):
        self.actions += 1
        out = .5*mm(mm(inverse, D, self.xp), inverse, self.xp)+self.quartic_adjoint(self.edge_covariances(D))
        return (out+out.T)/2

    def stable_difference(self, Y, D, G, chol):
        xp = self.xp
        relative = xp.linalg.solve(chol, D)
        relative = xp.linalg.solve(chol, relative.T).T
        values = xp.linalg.eigvalsh((relative+relative.T)/2)
        if scalar(xp.min(values)) <= -1:
            return float('inf')
        small = xp.abs(values) < 1e-4
        y = values
        series = y*y*(.5+y*(-1/3+y*(.25+y*(-.2+y/6))))
        remainder = xp.where(small, series, values-xp.log1p(values))
        DS = self.edge_covariances(D)
        tr = xp.trace(DS, axis1=1, axis2=2)
        quadratic = self.beta*self.gamma*xp.sum(self.c*(2*xp.sum(DS*DS, axis=(1, 2))+tr*tr))
        return scalar(xp.sum(G*D)+.5*xp.sum(remainder)+quadratic)


def newton_cg(obj, initial_scale=1., tolerance=1e-9, max_iterations=50, cg_limit=100, progress=None):
    xp = obj.xp
    Y = initial_scale*xp.eye(obj.d)
    history = []
    start = time.perf_counter()
    for iteration in range(max_iterations):
        value, G, inverse, S, chol = obj.evaluate(Y)
        norm = scalar(xp.linalg.norm(G))
        row = {'iteration': iteration, 'beta_F_minus_constant': value,
               'whitened_gradient_frobenius': norm, 'elapsed_seconds': time.perf_counter()-start}
        history.append(row)
        if progress:
            progress(row)
        if norm <= tolerance:
            return Y, {'converged': True, 'history': history, 'hessian_actions': obj.actions,
                       'seconds': time.perf_counter()-start, 'gradient_tolerance': tolerance,
                       'matrix_free_covariance_parameter_hessian': True}
        # Preconditioned CG in the D-by-D symmetric matrix space.
        rhs = -G
        D = xp.zeros_like(Y)
        residual = rhs.copy()
        z = 2*mm(mm(Y, residual, xp), Y, xp)
        direction = z.copy()
        rz = scalar(xp.sum(residual*z))
        relative_tolerance = min(.1, max(1e-7, .1*norm**.5))
        cg_converged = False
        for inner in range(cg_limit):
            action = obj.hessian_action(direction, inverse)
            curvature = scalar(xp.sum(direction*action))
            if not curvature > 0:
                raise ArithmeticError('Nonpositive Newton-CG curvature.')
            alpha = rz/curvature
            D += alpha*direction
            residual -= alpha*action
            if scalar(xp.linalg.norm(residual)) <= relative_tolerance*norm:
                cg_converged = True
                break
            z = 2*mm(mm(Y, residual, xp), Y, xp)
            rz_new = scalar(xp.sum(residual*z))
            direction = z+(rz_new/rz)*direction
            rz = rz_new
        if not cg_converged:
            raise ArithmeticError('Newton-CG inner solve reached its declared limit.')
        D = (D+D.T)/2
        slope = scalar(xp.sum(G*D))
        if not slope < 0:
            raise ArithmeticError('Newton direction was not a descent direction.')
        step = 1.
        for backtrack in range(41):
            decrease = obj.stable_difference(Y, step*D, G, chol)
            if decrease <= 1e-4*step*slope:
                break
            step *= .5
        else:
            raise ArithmeticError('Positive-definite Armijo line search failed.')
        row.update(cg_iterations=inner+1, cg_relative_tolerance=relative_tolerance,
                   step_length=step, stable_objective_decrease=decrease)
        Y = (Y+step*D+(Y+step*D).T)/2
    return Y, {'converged': False, 'history': history, 'hessian_actions': obj.actions,
               'seconds': time.perf_counter()-start, 'gradient_tolerance': tolerance,
               'matrix_free_covariance_parameter_hessian': True}


def contracted_covariance(a, b, S, T, C, xp=np):
    """Exact Gaussian covariance of (2a.x+x.x)^2 and (2b.y+y.y)^2.

    S=Cov(x), T=Cov(y), C=Cov(x,y); leading dimensions broadcast.
    """
    tr = lambda M: xp.trace(M, axis1=-2, axis2=-1)
    bil = lambda x, M, y: xp.einsum('...i,...ij,...j->...', x, M, y)
    h = 4*(tr(S)[..., None]*a+2*xp.einsum('...ij,...j->...i', S, a))
    k = 4*(tr(T)[..., None]*b+2*xp.einsum('...ij,...j->...i', T, b))
    F = 4*a[..., :, None]*a[..., None, :]+2*tr(S)[..., None, None]*xp.eye(3)+4*S
    G = 4*b[..., :, None]*b[..., None, :]+2*tr(T)[..., None, None]*xp.eye(3)+4*T
    Ct = C.swapaxes(-1, -2)
    CC = mm(C, Ct, xp)
    return (bil(h, C, k)+2*tr(mm(mm(mm(F, C, xp), G, xp), Ct, xp))+32*bil(a, C, b)*tr(CC)
            +64*bil(a, mm(CC, C, xp), b)+8*tr(CC)**2+16*tr(mm(CC, CC, xp)))


def residue_covariance(r0, edges, cart_equal, cart_lag, xp=np, block=32):
    """Stream edge-pair blocks immediately into an M-by-N accumulator."""
    n, m = len(r0), len(edges)
    edges_host = host(edges).astype(int)
    deg = np.bincount(edges_host.ravel(), minlength=n)
    if np.any(deg == 0):
        raise ValueError('Isolated residues have undefined degree-averaged observables.')
    # Use sparse incidence for both contractions, avoiding repeated N residue masks.
    if xp is np:
        from scipy.sparse import coo_matrix
    else:
        from cupyx.scipy.sparse import coo_matrix
    eh = xp.asarray(edges_host)
    W = coo_matrix((xp.asarray(1/deg[edges_host.ravel()]),
                    (xp.asarray(edges_host.ravel()), xp.repeat(xp.arange(m), 2))), shape=(n, m)).tocsr()
    r0 = xp.asarray(r0)
    a = r0[eh[:, 0]]-r0[eh[:, 1]]
    coef = 1/(8*xp.sum(a*a, axis=1))
    V = cart_equal.reshape(n, 3, n, 3).transpose(0, 2, 1, 3)
    Q = cart_lag.reshape(n, 3, n, 3).transpose(0, 2, 1, 3)
    i, j = eh[:, 0], eh[:, 1]
    S = V[i, i]-V[i, j]-V[j, i]+V[j, j]
    trS = xp.trace(S, axis1=1, axis2=2)
    edge_mean = coef*(4*xp.einsum('mi,mij,mj->m', a, S, a)+2*xp.sum(S*S, axis=(1, 2))+trS*trS)
    means = W @ edge_mean
    KW = xp.empty((m, n), dtype=cart_lag.dtype)
    u, v = eh[None, :, 0], eh[None, :, 1]
    for start in range(0, m, block):
        ix = xp.arange(start, min(start+block, m))
        x, y = eh[ix, 0, None], eh[ix, 1, None]
        C = Q[x, u]-Q[x, v]-Q[y, u]+Q[y, v]
        K = contracted_covariance(a[ix, None], a[None], S[ix, None], S[None], C, xp)*coef[ix, None]*coef[None]
        KW[ix] = (W @ K.T).T
    covariance = W @ KW
    return covariance, means


def moments_for_covariance(model, Sigma, times, xp=np, beta=1., mobility=1., block=32):
    B = xp.asarray(model['B'])
    Sigma = xp.asarray(Sigma)
    values, vectors = xp.linalg.eigh((Sigma+Sigma.T)/2)
    if scalar(xp.min(values)) <= 0:
        raise ArithmeticError('Nonpositive fitted covariance eigenvalue.')
    modes = mm(B, vectors, xp)
    cart = mm(modes*values[None, :], modes.T, xp)
    covariances, timings = [], []
    means = None
    for t in [0.]+list(times):
        sync(xp)
        start = time.perf_counter()
        lag = cart if t == 0 else mm(modes*(values*xp.exp(-mobility*t/(beta*values)))[None, :], modes.T, xp)
        cov, means = residue_covariance(model['r0'], model['edges'], cart, lag, xp, block)
        sync(xp)
        covariances.append(cov)
        timings.append({'time': float(t), 'seconds': time.perf_counter()-start})
    return xp.stack(covariances), means, timings


def run_model(model_path, scales_path, output_dir, backend_name='numpy', label='protein', block=32):
    begin = time.perf_counter()
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    model = dict(np.load(model_path))
    sd = np.load(scales_path)['harmonic_sd']
    xp = backend(backend_name)
    tau = 1/model['eigenvalues'][0]
    times = tau*np.array([.1, 1., 10.])
    manifest = {'label': label, 'model_path': str(model_path), 'model_sha256': digest(model_path),
                'scales_path': str(scales_path), 'scales_sha256': digest(scales_path),
                'residues': len(model['r0']), 'contacts': len(model['edges']),
                'internal_coordinates': len(model['eigenvalues']), 'beta': 1., 'mobility': 1.,
                'times_over_tau': [.1, 1., 10.], 'tau': float(tau), 'times': times.tolist(),
                'centroid': 'fixed at native structure', 'domain': 'all retained internal-coordinate space R^D for Gaussian surrogate',
                'potential': 'original biquadratic contact potential expectation in Gaussian family',
                'observables': 'original degree-averaged quartic contact energies',
                'normalization': 'original all-positive-mode harmonic SD of same quartic observables',
                'source_sha256': digest(__file__), 'created_utc': datetime.datetime.now(datetime.timezone.utc).isoformat()}
    (out/'input-and-model-manifest.json').write_text(json.dumps(manifest, indent=2))
    obj = CovarianceObjective(model['r0'], model['edges'], model['B'], model['eigenvalues'], xp)
    def progress(row):
        (out/'fit-progress.json').write_text(json.dumps(row, indent=2))
        print(label, backend_name, 'fit', row['iteration'], row['whitened_gradient_frobenius'], flush=True)
    Y, fit = newton_cg(obj, progress=progress)
    if not fit['converged']:
        (out/'fit-receipt.json').write_text(json.dumps(fit, indent=2))
        raise ArithmeticError('Covariance solver failed stationarity criterion.')
    lam = xp.asarray(model['eigenvalues'])
    Sigma = Y/xp.sqrt(lam[:, None]*lam[None, :])
    sigma_host = host(Sigma)
    # Release fitted covariance promptly for the original-measure integration lane.
    np.savez_compressed(out/'fitted-covariance.npz', covariance=sigma_host,
                        whitened_covariance=host(Y), model_sha256=manifest['model_sha256'])
    fit.update(minimum_covariance_eigenvalue=float(np.linalg.eigvalsh(sigma_host)[0]),
               normalized_stationarity_frobenius=fit['history'][-1]['whitened_gradient_frobenius']/np.sqrt(len(lam)))
    (out/'fit-receipt.json').write_text(json.dumps(fit, indent=2))
    covs, means, timings = moments_for_covariance(model, Sigma, times, xp, block=block)
    covs, means = host(covs), host(means)
    scales = np.outer(sd, sd)
    C = (covs[1:]-covs[0])/scales[None]
    receivers = model['receiver']
    scores = np.sqrt(np.mean(C[:, :, receivers]**2, axis=-1))
    candidates = np.flatnonzero(model['candidate'])
    orders = np.array([candidates[np.lexsort((model['canonical'][candidates], -score[candidates]))] for score in scores])
    np.savez_compressed(out/'raw-moments.npz', covariance=covs[0], delayed_covariance=covs[1:],
                        residue_means=means, harmonic_sd=sd, times=times)
    np.savez_compressed(out/'response-matrices.npz', C=C, covariance=covs[0]/scales,
                        delayed_covariance=covs[1:]/scales, scores=scores, orders=orders,
                        canonical=model['canonical'], receiver=receivers, candidate=model['candidate'], times=times)
    with (out/'residue-shortlists.csv').open('w') as f:
        writer = csv.writer(f); writer.writerow(['time_over_tau', 'rank', 'residue_canonical', 'row_index', 'score'])
        for time_ratio, order, score in zip([.1, 1., 10.], orders, scores):
            for rank, index in enumerate(order, 1):
                writer.writerow([time_ratio, rank, int(model['canonical'][index]), int(index), float(score[index])])
    symmetry = float(np.max(np.abs(covs-covs.swapaxes(-1, -2))))
    minimum_covariance_eigenvalue = float(np.linalg.eigvalsh((covs[0]+covs[0].T)/2)[0])
    check = {'engineering_status': 'COMPLETE', 'scientific_status': 'PASS' if symmetry <= 1e-9 and minimum_covariance_eigenvalue >= -1e-10 and np.isfinite(C).all() else 'FAIL',
             'scope': 'Gaussian-surrogate fit and complete exact OU moment computation only',
             'nonlinear_response_fidelity': 'NOT_EVALUATED', 'biological_benefit': 'NOT_EVALUATED',
             'quantum_advantage': 'NOT_EVALUATED', 'fit_converged': fit['converged'],
             'covariance_symmetry_max_abs': symmetry, 'minimum_static_covariance_eigenvalue': minimum_covariance_eigenvalue,
             'all_response_entries_finite': bool(np.isfinite(C).all()), 'top5_by_time': [model['canonical'][o[:5]].tolist() for o in orders]}
    (out/'numerical-checks.json').write_text(json.dumps(check, indent=2))
    receipt = {'backend': backend_name, 'numpy_version': np.__version__, 'platform': platform.platform(),
               'total_seconds_including_load_fit_covariance_readout_save': time.perf_counter()-begin,
               'fit_seconds': fit['seconds'], 'moment_timings': timings,
               'peak_process_rss_platform_units': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
               'rss_units': 'bytes on macOS; KiB on Linux', 'block_edges': block,
               'dense_covariance_parameter_hessian_allocated': False, 'source_sha256': digest(__file__)}
    if xp.__name__ == 'cupy':
        receipt.update(cupy_version=xp.__version__, cuda_runtime=xp.cuda.runtime.runtimeGetVersion(),
                       gpu_name=xp.cuda.runtime.getDeviceProperties(0)['name'].decode(),
                       device_memory_pool_total_bytes=xp.get_default_memory_pool().total_bytes())
    (out/'total-cost-receipt.json').write_text(json.dumps(receipt, indent=2))
    return check, receipt


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--model', required=True); parser.add_argument('--scales', required=True)
    parser.add_argument('--output', required=True); parser.add_argument('--label', default='protein')
    parser.add_argument('--backend', default='numpy'); parser.add_argument('--block', default=32, type=int)
    args = parser.parse_args()
    result = run_model(args.model, args.scales, args.output, args.backend, args.label, args.block)
    print(json.dumps(result, indent=2))
