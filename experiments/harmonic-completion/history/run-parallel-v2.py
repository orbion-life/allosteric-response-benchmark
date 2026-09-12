"""Locked additive response correction on a fully resolved nonlinear network.

The three-node geometry has three internal Cartesian Hessian modes. There is
no protein label, molecular trajectory, fitting or stochastic simulation here.
"""
from pathlib import Path
import argparse, datetime, hashlib, json, resource, sys, time
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
from scipy.linalg import eigh
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import expm_multiply
from scipy.special import expit, logsumexp

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'vendor'))
from prepare import polynomial, gaussian_sd
from run_pilot import harmonic_covariances

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def contract(a, b):
    return np.einsum('ij,jk->ik', a, b, optimize=False)

def static_model(protocol, kappa):
    r0 = np.array(protocol['geometry_A'], dtype=float)
    edges = np.array(protocol['edges_zero_based'], dtype=int)
    hess = np.zeros((r0.size, r0.size))
    for i, j in edges:
        v = r0[i] - r0[j]
        v /= np.linalg.norm(v)
        v = np.outer(v, v) * kappa
        for a, b, sign in [(i, i, 1), (j, j, 1), (i, j, -1), (j, i, -1)]:
            hess[3*a:3*a+3, 3*b:3*b+3] += sign*v
    vals, B = eigh(hess)
    keep = vals > 1e-8*kappa
    assert keep.sum() == protocol['all_internal_modes']
    lam, B = vals[keep], B[:, keep]
    for c in B.T:
        if c[np.argmax(np.abs(c))] < 0:
            c *= -1
    A, powers = polynomial(r0, edges, B, kappa=kappa)
    sd, _, _ = gaussian_sd(A, powers, lam)
    return r0, edges, lam, B, sd

def grid(r0, edges, lam, B, kappa, d, n, extent, kind, times, sd):
    started = time.monotonic()
    axis = np.linspace(-extent/np.sqrt(lam[0]), extent/np.sqrt(lam[0]), n)
    delta = axis[1]-axis[0]
    inds = np.indices((n,)*d).reshape(d, -1).T
    q = axis[inds]
    G = len(q)
    rr = r0[None] + contract(q, B[:, :d].T).reshape(G, len(r0), 3)
    E = np.zeros((G, len(r0)))
    U = np.zeros(G)
    strain = np.zeros(G)
    degree = np.bincount(edges.ravel(), minlength=len(r0))
    for i, j in edges:
        rest2 = np.sum((r0[i]-r0[j])**2)
        length2 = np.sum((rr[:, i]-rr[:, j])**2, axis=1)
        u = kappa*(length2-rest2)**2/(8*rest2)
        U += u
        E[:, i] += u/degree[i]
        E[:, j] += u/degree[j]
        strain = np.maximum(strain, np.abs(np.sqrt(length2/rest2)-1))
    if kind == 'harmonic':
        U = .5*np.sum(lam[:d]*q*q, axis=1)
    pi = np.exp(-U-logsumexp(-U))
    f = np.sqrt(pi)[:, None]*(E-np.einsum('i,ij->j', pi, E, optimize=False))
    row, col, data = [], [], []
    diag = np.zeros(G)
    flat = np.arange(G).reshape((n,)*d)
    for dim in range(d):
        lo, hi = [slice(None)]*d, [slice(None)]*d
        lo[dim], hi[dim] = slice(0, -1), slice(1, None)
        x, y = flat[tuple(lo)].ravel(), flat[tuple(hi)].ravel()
        xy = 2/delta**2*expit(-(U[y]-U[x]))
        yx = 2/delta**2*expit(U[y]-U[x])
        h = -np.sqrt(xy*yx)
        row.extend([x, y]); col.extend([y, x]); data.extend([h, h])
        np.add.at(diag, x, xy); np.add.at(diag, y, yx)
    row.append(np.arange(G)); col.append(np.arange(G)); data.append(diag)
    H = coo_matrix((np.concatenate(data), (np.concatenate(row), np.concatenate(col))), shape=(G,G)).tocsr()
    cov = contract(f.T, f)
    norm = np.outer(sd, sd)
    delayed = []
    for t in times:
        ft = expm_multiply(-t*H, f, traceA=-t*diag.sum())
        delayed.append(contract(f.T, ft)/norm)
    delayed = np.array(delayed)
    C = delayed-cov/norm
    assert np.isfinite(C).all()
    stationarity = float(np.max(np.abs(H@np.sqrt(pi))))
    assert stationarity < 1e-9
    meta = dict(d=d,n=n,extent=extent,kind=kind,kappa=kappa,states=G,
                seconds=time.monotonic()-started,stationarity_residual=stationarity,
                boundary_probability=float(pi[np.any((inds==0)|(inds==n-1),axis=1)].sum()),
                strain_probability_gt_0_2=float(pi[strain>.2].sum()),
                mean_max_contact_strain=float(np.sum(pi*strain)))
    return dict(C=C, equilibrium=cov/norm, delayed=delayed, harmonic_sd=sd,
                eigenvalues=lam, basis=B, times=times), meta

def worker(job):
    output, protocol, kappa, d, n, extent, kind = job
    r0,edges,lam,B,sd=static_model(protocol,kappa)
    times=np.array(protocol['times_over_tau'])/lam[0]
    arrays,meta=grid(r0,edges,lam,B,kappa,d,n,extent,kind,times,sd)
    name=f'{kind}-k{kappa}-d{d}-n{n}-e{extent}'
    np.savez_compressed(output/(name+'.npz'),**arrays)
    (output/(name+'.json')).write_text(json.dumps(meta,indent=2)+'\n')
    meta['worker_peak_bytes']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024)
    return name,meta

def run(output, workers=4, resume=False):
    output.mkdir(parents=True, exist_ok=resume)
    protocol = json.loads((ROOT/'preanalysis.json').read_text())
    execution = json.loads((ROOT/'execution-amendment.json').read_text())
    assert sha(ROOT/'preanalysis.json') == (ROOT/'preanalysis.sha256').read_text().strip()
    started = time.monotonic()
    receipt = dict(status='RUNNING',started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                   protocol_sha256=sha(ROOT/'preanalysis.json'),script_sha256=sha(__file__),cases=[])
    if resume:
        old=json.loads((output/'run.json').read_text())
        assert old['protocol_sha256']==receipt['protocol_sha256']
        if not (output/'previous-execution.json').exists():
            (output/'previous-execution.json').write_text(json.dumps(old,indent=2)+'\n')
        # Only a completed pair of files is a checkpoint. The numeric grid
        # implementation is unchanged; the execution amendment preserves its AST.
    receipt.update(workers=workers,execution_amendment_sha256=sha(ROOT/'execution-amendment.json'))
    (output/'run.json').write_text(json.dumps(receipt,indent=2)+'\n')
    jobs=[]
    for kappa in protocol['kappa']:
        r0,edges,lam,B,sd = static_model(protocol,kappa)
        times = np.array(protocol['times_over_tau'])/lam[0]
        for d in protocol['retained_dimensions']:
            cov, delayed = [], []
            for t in times:
                a,b=harmonic_covariances(r0,edges,B[:,:d],lam[:d],t)
                cov.append(a*kappa**2/np.outer(sd,sd));delayed.append(b*kappa**2/np.outer(sd,sd))
            np.savez_compressed(output/f'analytic-k{kappa}-d{d}.npz',C=np.array(delayed)-np.array(cov),equilibrium=cov[0],delayed=delayed,harmonic_sd=sd,eigenvalues=lam,basis=B,times=times)
            cases=[(n,protocol['domain_sigma_slowest']) for n in protocol['grid_n']]
            cases.append((protocol['domain_check']['grid_n'],protocol['domain_check']['extent']))
            for n,extent in cases:
                for kind in ['biquadratic','harmonic']:
                    name=f'{kind}-k{kappa}-d{d}-n{n}-e{extent}'
                    if resume and (output/(name+'.npz')).exists() and (output/(name+'.json')).exists():
                        meta=json.loads((output/(name+'.json')).read_text())
                        assert all(meta[k]==v for k,v in dict(kappa=kappa,d=d,n=n,extent=extent,kind=kind).items())
                        with np.load(output/(name+'.npz')) as saved:
                            assert saved['C'].shape==(3,3,3) and all(np.isfinite(saved[k]).all() for k in saved.files)
                        receipt['cases'].append(meta)
                    else:
                        jobs.append((output,protocol,kappa,d,n,extent,kind))
    with ProcessPoolExecutor(max_workers=workers) as pool:
        pending=[pool.submit(worker,job) for job in jobs]
        try:
            for future in as_completed(pending,timeout=execution['max_seconds_per_execution']):
                name,meta=future.result()
                assert meta['worker_peak_bytes']*workers<=execution['peak_total_budget_bytes']
                receipt['cases'].append(meta)
                (output/'run.json').write_text(json.dumps(receipt,indent=2)+'\n')
                print(name,round(meta['seconds'],3),flush=True)
        except BaseException:
            for future in pending:future.cancel()
            for process in pool._processes.values():process.terminate()
            raise
    receipt.update(status='COMPLETE',seconds=time.monotonic()-started,peak_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024),completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
    (output/'run.json').write_text(json.dumps(receipt,indent=2)+'\n')
    (output/'freeze.json').write_text(json.dumps({p.name:sha(p) for p in sorted(output.iterdir()) if p.is_file()},indent=2)+'\n')

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,default=ROOT/'results');parser.add_argument('--workers',type=int,choices=range(1,5),default=4);parser.add_argument('--resume',action='store_true')
    args=parser.parse_args()
    try:
        run(args.output,args.workers,args.resume)
    except Exception as exc:
        if args.output.exists():
            (args.output/'failure.json').write_text(json.dumps({'type':type(exc).__name__,'message':str(exc)},indent=2)+'\n')
        raise
