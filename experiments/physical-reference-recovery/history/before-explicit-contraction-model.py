"""Independent three-node model and polynomial algebra; no saved responses used."""
from collections import defaultdict
from functools import lru_cache
from itertools import product
import json
import os
from pathlib import Path
import numpy as np
from scipy.linalg import eigh

ROOT = Path(__file__).resolve().parent


def runtime_metadata():
    from importlib.metadata import version
    from threadpoolctl import threadpool_info
    pools = threadpool_info()
    if any(p['num_threads'] != 1 for p in pools):
        raise RuntimeError('Every numerical threadpool must have exactly one thread.')
    clean = [{k: p.get(k) for k in ('user_api', 'internal_api', 'num_threads', 'version', 'threading_layer', 'architecture')} for p in pools]
    return dict(versions={name: version(name) for name in ('numpy', 'scipy', 'psutil', 'threadpoolctl')},
                threadpools=clean,
                thread_environment={k: os.environ.get(k) for k in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS')})


def protocol():
    return json.loads((ROOT / 'protocol.json').read_text())


def multiply(a, b):
    out = defaultdict(float)
    for p, u in a.items():
        for q, v in b.items():
            out[tuple(x+y for x, y in zip(p, q))] += u*v
    return dict(out)


def build_model(geometry=None, kappa=1.):
    p = protocol()
    r0 = np.asarray(p['geometry_A'] if geometry is None else geometry, dtype=float)
    edges = np.asarray(p['edges_zero_based'], dtype=int)
    H = np.zeros((9, 9))
    for i, j in edges:
        a = r0[i]-r0[j]
        block = kappa*np.outer(a, a)/np.dot(a, a)
        for u, v, sign in [(i, i, 1), (j, j, 1), (i, j, -1), (j, i, -1)]:
            H[3*u:3*u+3, 3*v:3*v+3] += sign*block
    eigenvalues, basis = eigh(H)
    keep = eigenvalues > 1e-8*kappa
    lam, B = eigenvalues[keep], basis[:, keep]
    if len(lam) != 3:
        raise ValueError('Expected all three positive internal modes.')
    for column in B.T:
        if column[np.argmax(np.abs(column))] < 0:
            column *= -1
    beta, mu = p['beta'], p['mu']
    Bx = B/np.sqrt(beta*lam)[None, :]
    degree = np.bincount(edges.ravel(), minlength=3)
    polys = [defaultdict(float) for _ in range(3)]
    edge_polys = []
    for i, j in edges:
        a = r0[i]-r0[j]
        D = Bx[3*i:3*i+3]-Bx[3*j:3*j+3]
        strain = defaultdict(float)
        for axis in range(3):
            z = tuple(int(k == axis) for k in range(3))
            strain[z] += 2*np.dot(a, D[:, axis])
        for axis in range(3):
            for other in range(3):
                z = tuple(int(k == axis)+int(k == other) for k in range(3))
                strain[z] += np.dot(D[:, axis], D[:, other])
        edge = {z: val*kappa/(8*np.dot(a, a)) for z, val in multiply(strain, strain).items()}
        edge_polys.append(edge)
        for node in (i, j):
            for z, val in edge.items():
                polys[node][z] += val/degree[node]
    powers = np.array(sorted(set().union(*(x.keys() for x in polys))), dtype=int)
    coefficients = np.array([[poly.get(tuple(z), 0.) for z in powers] for poly in polys])
    return dict(r0=r0, edges=edges, eigenvalues=lam, basis=B, Bx=Bx,
                powers=powers, coefficients=coefficients, kappa=kappa, beta=beta,
                mu=mu, degree=degree, times=np.array(p['times_over_tau'])/(mu*lam[0]))


def energies(model, x):
    """Evaluate original contact law directly, independently of coefficients."""
    x = np.asarray(x)
    rr = model['r0'][None]+(x @ model['Bx'].T).reshape(-1, 3, 3)
    E = np.zeros((len(x), 3))
    U = np.zeros(len(x))
    for i, j in model['edges']:
        rest2 = np.sum((model['r0'][i]-model['r0'][j])**2)
        u = model['kappa']*(np.sum((rr[:, i]-rr[:, j])**2, axis=1)-rest2)**2/(8*rest2)
        U += u
        E[:, i] += u/model['degree'][i]
        E[:, j] += u/model['degree'][j]
    return U, E


def polynomial_values(model, x):
    monomials = np.prod(np.asarray(x)[:, None, :]**model['powers'][None, :, :], axis=2)
    return monomials @ model['coefficients'].T


@lru_cache(None)
def univariate_moment(n):
    if n == 0:
        return 1.
    if n % 2:
        return 0.
    return (n-1)*univariate_moment(n-2)


def bivariate_moments(rho):
    """Wick/Stein recurrence for E[X^a Y^b], a,b<=4, unit variances."""
    @lru_cache(None)
    def moment(a, b):
        if a == 0:
            return univariate_moment(b)
        total = (a-1)*moment(a-2, b) if a >= 2 else 0.
        if b:
            total += b*rho*moment(a-1, b-1)
        return total
    return np.array([[moment(a, b) for b in range(5)] for a in range(5)])


def contract_moments(model, moments, first):
    p, A = model['powers'], model['coefficients']
    kernel = np.ones((len(p), len(p)))
    means = np.ones(len(p))
    for k in range(3):
        kernel *= moments[k][p[:, k, None], p[None, :, k]]
        means *= first[k][p[:, k]]
    m = A @ means
    return A @ kernel @ A.T-np.outer(m, m), m


def analytic_gaussian(model):
    first = [np.array([univariate_moment(k) for k in range(5)])]*3
    cov, mean = contract_moments(model, [bivariate_moments(1.)]*3, first)
    sd = np.sqrt(np.diag(cov))
    if not np.all(np.isfinite(sd) & (sd > 0)):
        raise ValueError('Nonpositive or invalid harmonic scales.')
    delayed = []
    for t in model['times']:
        moments = [bivariate_moments(np.exp(-model['mu']*lam*t)) for lam in model['eigenvalues']]
        value, _ = contract_moments(model, moments, first)
        delayed.append(value)
    norm = np.outer(sd, sd)
    G0, K = cov/norm, np.array(delayed)/norm
    return dict(G0=G0, K=K, C=K-G0, sd=sd, mean=mean)


def grid_shape(faces, spacing):
    ratio = 2*np.asarray(faces)/spacing
    if not np.allclose(ratio, np.round(ratio), atol=1e-12, rtol=0):
        raise ValueError('Faces must be integer cell multiples.')
    return tuple(int(x) for x in np.round(ratio))


def known_minimum(model):
    displacement = -2*(model['r0']-model['r0'].mean(0))
    q = model['basis'].T @ displacement.ravel()
    x = np.sqrt(model['beta']*model['eigenvalues'])*q
    rebuilt = (model['basis'] @ q).reshape(3, 3)
    U, _ = energies(model, x[None])
    faces = np.ceil(np.maximum(4, np.abs(x)+2)).astype(int)
    return dict(q=q.tolist(), x=x.tolist(), faces=faces.tolist(),
                reconstruction_error=float(np.max(np.abs(rebuilt-displacement))),
                contact_energy=float(U[0]))
