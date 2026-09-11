"""Independent dense Frechet verification. Does not import experiment functions."""
from pathlib import Path
import json, itertools, time
import numpy as np
from numpy.polynomial.hermite import hermgauss
from scipy.linalg import expm_frechet

def mm(a,b): return np.einsum("ij,jk->ik",a,b)
def vm(a,b): return np.einsum("i,ij->j",a,b)
def mv(a,b): return np.einsum("ij,j->i",a,b)
def expm(a):
    steps=max(0,int(np.ceil(np.log2(max(1.,np.linalg.norm(a,np.inf))))))
    b=a/(2**steps);term=np.eye(len(a));total=term.copy()
    for degree in range(1,200):
        term=mm(term,b)/degree;total+=term
        if np.max(np.abs(term))<1e-18:break
    for _ in range(steps):total=mm(total,total)
    return total

np.seterr(all="raise")
from scipy.special import expit

import os
HERE=Path(os.environ.get('PULSAR_OUTPUT_DIR',Path(__file__).resolve().parent/'results'))
original=json.loads((HERE/'results.json').read_text())
r0=np.array(original['parameters']['positions']);edges=original['parameters']['edges_zero_based']
B=np.array(original['geometry']['B']);base_eigen=np.array(original['geometry']['retained_eigenvalues'])
base_k=original['parameters']['stiffness'];n=8; beta=1.;mu=1.
degrees=np.zeros(len(r0))
for i,j in edges: degrees[i]+=1;degrees[j]+=1

def energy(q,k):
    positions=r0+mv(B,q).reshape(-1,3)
    E=np.zeros(len(r0));u=0.
    for i,j in edges:
        v0=r0[i]-r0[j];v=positions[i]-positions[j]
        edge=k*(np.sum(v*v)-np.sum(v0*v0))**2/(8*np.sum(v0*v0))
        u+=edge;E[i]+=edge/degrees[i];E[j]+=edge/degrees[j]
    return u,E

all_results=[];began=time.time()
for k in [10.,100.,1000.]:
    lam=base_eigen*k/base_k;t=1/lam[0]
    half=4/np.sqrt(lam[0]);axis=np.linspace(-half,half,n);dx=axis[1]-axis[0]
    points=np.array(list(itertools.product(axis,axis)))
    evaluated=[energy(q,k) for q in points];fullU=np.array([p[0] for p in evaluated]);E=np.array([p[1] for p in evaluated])
    # Gaussian quadrature independently checks normalization of degree-four observables.
    gh,gw=hermgauss(10);G_E=[];G_W=[]
    for a,b in itertools.product(range(10),repeat=2):
        q=np.sqrt(2/lam)*[gh[a],gh[b]];G_E.append(energy(q,k)[1]);G_W.append(gw[a]*gw[b]/np.pi)
    G_E=np.array(G_E);G_W=np.array(G_W);gm=vm(G_W,G_E)
    ghcov=mm((G_E-gm).T,(G_W[:,None]*(G_E-gm)));norm=np.sqrt(np.outer(np.diag(ghcov),np.diag(ghcov)))
    for model in ['harmonic','biquadratic']:
        U=fullU if model=='biquadratic' else .5*np.sum(points**2*lam,axis=1)
        pi=np.exp(-beta*(U-U.min()));pi/=pi.sum()
        L=np.zeros((n*n,n*n));derivatives=np.zeros((5,n*n,n*n))
        for ix in range(n*n):
            a,b=divmod(ix,n)
            for aa,bb in [(a-1,b),(a+1,b),(a,b-1),(a,b+1)]:
                if aa<0 or aa>=n or bb<0 or bb>=n:continue
                j=aa*n+bb;f=expit(-beta*(U[j]-U[ix]));rate=2*mu/(beta*dx**2)*f
                L[ix,j]=rate;L[ix,ix]-=rate
                for sender in range(5):
                    value=-2*mu/dx**2*(E[j,sender]-E[ix,sender])*f*(1-f)
                    derivatives[sender,ix,j]=value;derivatives[sender,ix,ix]-=value
        centered=E-vm(pi,E);prop=expm(t*L)
        Rcov=-beta*mm(centered.T,(pi[:,None]*(centered-mm(prop,centered))))
        Rf=np.column_stack([vm(vm(pi,expm_frechet(t*L,t*d,compute_expm=False)),E) for d in derivatives])
        C=Rcov/(beta*norm)
        item={'stiffness':k,'model':model,'n':n,'extent':4,'horizon':t,'C_receiver4_sender3':float(C[4,3]),
              'frechet_covariance_max_error':float(np.max(np.abs(Rf-Rcov))),
              'frechet_C_max_error':float(np.max(np.abs((Rf-Rcov)/(beta*norm)))),
              'R':Rcov.tolist(),'C':C.tolist(),'normalization_covariance':ghcov.tolist()}
        if k==base_k:
            ref=next(c for c in original['cases'] if c['n']==n and c['extent_slowest_harmonic_sigma']==4 and c['model']==model)
            item['recorded_response_max_difference']=float(np.max(np.abs(Rcov-np.array(ref['R']))))
            item['recorded_gaussian_covariance_max_difference']=float(np.max(np.abs(ghcov-np.array(original['polynomial']['analytic_reduced_harmonic_covariance']))))
        assert item['frechet_C_max_error']<1e-10
        all_results.append(item)
out={'method':'Independent scalar energy evaluation, dense backward generator, Taylor scaling/squaring exponential and SciPy analytic Frechet derivative; Gauss-Hermite normalization; strict floating-point checking and einsum contractions. No import of original experiment code.',
     'cases':all_results,'seconds':time.time()-began,'status':'all six cases passed normalized error <1e-10'}
(HERE/'independent-response-check.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'status':out['status'],'max_C_error':max(x['frechet_C_max_error'] for x in all_results),'seconds':out['seconds']}))
