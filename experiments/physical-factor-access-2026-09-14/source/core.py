"""Explicit physical contact factors and normalized Hermite coefficients.

All dense covariance, coordinate transforms and snapshot whitening are charged
classical inputs. A selected coefficient is not a complete quantum access oracle.
"""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[key]='4'
import sys, math, itertools, time, json, hashlib
from pathlib import Path
import numpy as np
from scipy.linalg import eigh, solve, expm
from scipy.linalg.blas import dgemm
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import LinearOperator, expm_multiply

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'vendor'))
import gaussian_protein as gp
from fixed_adapter import OUKernel

def mm(a,b):
    a,b=np.asarray(a),np.asarray(b)
    if a.ndim==b.ndim==2:return dgemm(1.,a,b)
    return np.einsum({(1,2):'i,ij->j',(2,1):'ij,j->i',(1,1):'i,i->'}[(a.ndim,b.ndim)],a,b,optimize=False)

def dump(path,data):Path(path).write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def force_factor(model,Sigma,beta=1.,mu=1.,kappa=1.):
    tic=time.perf_counter();r=model['r0'];edges=model['edges'];B0=model['B'];n=len(r)
    cart=mm(mm(B0,Sigma),B0.T).reshape(n,3,n,3).transpose(0,2,1,3)
    i,j=edges.T;S=cart[i,i]-cart[i,j]-cart[j,i]+cart[j,j];S=(S+S.transpose(0,2,1))/2
    a=r[i]-r[j];l2=np.sum(a*a,axis=1);tr=np.trace(S,axis1=1,axis2=2)
    blocks=kappa*(np.einsum('ei,ej->eij',a,a)+S+tr[:,None,None]*np.eye(3)/2)/l2[:,None,None]
    vals,vecs=np.linalg.eigh(blocks)
    if vals.min() < -1e-12:raise ArithmeticError('Non-PSD local contact block')
    roots=np.einsum('eik,ek,ejk->eij',vecs,np.sqrt(np.maximum(vals,0)),vecs)
    rows=[];cols=[];data=[]
    for e,(u,v) in enumerate(edges):
        for aa in range(3):
            for bb in range(3):
                for site,sign in [(u,1),(v,-1)]:
                    rows.append(3*e+aa);cols.append(3*site+bb);data.append(sign*roots[e,aa,bb])
    L=coo_matrix((data,(rows,cols)),shape=(3*len(edges),3*n)).tocsr();L.eliminate_zeros()
    generation_seconds=time.perf_counter()-tic
    tic=time.perf_counter();Kcart=(L.T@L).tocsr();K=mm(B0.T,Kcart@B0);K=(K+K.T)/2
    assemble_seconds=time.perf_counter()-tic
    tic=time.perf_counter();F=np.sqrt(mu)*(L@B0);factor_projection_seconds=time.perf_counter()-tic
    tic=time.perf_counter();FtF=mm(F.T,F);factor_identity=float(np.max(abs(FtF-mu*K)));identity_seconds=time.perf_counter()-tic
    tic=time.perf_counter();Sinv=solve(Sigma,np.eye(len(Sigma)),assume_a='pos');inverse_seconds=time.perf_counter()-tic
    err=beta*K-Sinv;err=(err+err.T)/2
    drift_error=float(np.max(abs(eigh(mu*err/beta,eigvals_only=True))))
    sigma_physical=solve(beta*K,np.eye(len(K)),assume_a='pos')
    stats={'beta':beta,'mu':mu,'kappa':kappa,'contacts':len(edges),'coordinate_dimension':len(K),
      'cartesian_factor_shape':list(L.shape),'cartesian_factor_nnz':L.nnz,'cartesian_factor_bytes':L.data.nbytes+L.indices.nbytes+L.indptr.nbytes,
      'local_block_min_eigenvalue':float(vals.min()),'local_root_reconstruction_maxabs':float(np.max(abs(np.einsum('eji,ejk->eik',roots,roots)-blocks))),
      'factor_identity_maxabs':factor_identity,'factor_identity_relative_frobenius':float(np.linalg.norm(FtF-mu*K)/np.linalg.norm(mu*K)),
      'stationarity_maxabs':float(np.max(abs(err))),'stationarity_relative_frobenius':float(np.linalg.norm(err)/np.linalg.norm(Sinv)),
      'coordinate_drift_mismatch_spectral_norm':drift_error,'factor_consistent_covariance_maxabs':float(np.max(abs(sigma_physical-Sigma))),
      'dense_B0_bytes':B0.nbytes,'dense_projected_factor_shape':list(F.shape),'dense_projected_factor_bytes':F.nbytes,
      'force_factor_frobenius_squared':float(np.sum(F*F)),'cartesian_factor_frobenius_squared':float(np.sum(L.data*L.data)),
      'timing_seconds':{'local_covariance_blocks_and_sparse_factor':generation_seconds,'coordinate_drift_assembly':assemble_seconds,'explicit_dense_factor_projection':factor_projection_seconds,'independent_factor_identity':identity_seconds,'reference_covariance_inverse':inverse_seconds}}
    return dict(L=L,F=F,K=K,Gamma=mu*K,Sigma_physical=sigma_physical,blocks=blocks,S=S,stats=stats)

def occupations(d,max_degree=4):
    return [tuple(x) for k in range(1,max_degree+1) for x in itertools.combinations_with_replacement(range(d),k)]

class Hermite:
    """Coefficient oracle implemented classically, with explicit dense setup."""
    def __init__(self,model,Sigma,sd,mu=1.,beta=1.):
        tic=time.perf_counter();self.model=model;self.sd=np.asarray(sd);self.n=len(sd);self.edges=model['edges'];self.d=len(Sigma)
        self.values,self.V=eigh(Sigma);self.rates=mu/(beta*self.values)
        self.modes=mm(model['B'],self.V);i,j=self.edges.T
        self.A=(self.modes.reshape(self.n,3,self.d)[i]-self.modes.reshape(self.n,3,self.d)[j])*np.sqrt(self.values)
        self.a=model['r0'][i]-model['r0'][j];self.c=1/(8*np.sum(self.a*self.a,axis=1))
        self.S=np.einsum('ead,ebd->eab',self.A,self.A);self.tr=np.trace(self.S,axis1=1,axis2=2)
        self.degree=np.bincount(self.edges.ravel(),minlength=self.n)
        self.ell=np.einsum('ea,ead->ed',self.a,self.A)
        self.setup_seconds=time.perf_counter()-tic

    def coeff(self,indices):
        ix=tuple(indices);k=len(ix);fac=math.sqrt(math.prod(math.factorial(ix.count(j)) for j in set(ix)))
        aa=[self.A[:,:,j] for j in ix];ll=[self.ell[:,j] for j in ix]
        q=lambda a,b:np.einsum('ei,ei->e',aa[a],aa[b])
        if k==1:
            val=4*self.c*(self.tr*ll[0]+2*np.einsum('ei,eij,ej->e',self.a,self.S,aa[0]))
        elif k==2:
            val=(2/fac)*self.c*(4*ll[0]*ll[1]+2*self.tr*q(0,1)+4*np.einsum('ei,eij,ej->e',aa[0],self.S,aa[1]))
        elif k==3:
            val=(8/fac)*self.c*(ll[0]*q(1,2)+ll[1]*q(0,2)+ll[2]*q(0,1))
        elif k==4:
            val=(8/fac)*self.c*(q(0,1)*q(2,3)+q(0,2)*q(1,3)+q(0,3)*q(1,2))
        else:raise ValueError(k)
        v=np.bincount(self.edges[:,0],weights=val,minlength=self.n)+np.bincount(self.edges[:,1],weights=val,minlength=self.n)
        return v/(self.degree*self.sd)

    def row(self,indices,W,nodes):
        h=self.coeff(indices);omega=float(np.sum(self.rates[list(indices)]))
        features=np.exp(-np.asarray(nodes)*omega)[:,None]*h[None,:]
        psi=mm(features.ravel(),W)
        return psi,np.sqrt(omega)*psi,omega

def pair_kernel(kernel,t,i,j):
    """Exact Gaussian scalar covariance without evaluating all residue pairs."""
    edges=kernel.edges;ei=np.where((edges==i).any(axis=1))[0];ej=np.where((edges==j).any(axis=1))[0]
    modes=kernel.modes.reshape(kernel.n,3,kernel.d)
    U=modes[edges[ei,0]]-modes[edges[ei,1]];V=modes[edges[ej,0]]-modes[edges[ej,1]]
    C=np.einsum('ead,d,fbd->efab',U,kernel.values*np.exp(-t/kernel.values),V,optimize=True)
    raw=gp.contracted_covariance(kernel.a[ei,None],kernel.a[None,ej],kernel.S[ei,None],kernel.S[None,ej],C)
    return float(np.sum(raw*kernel.c[ei,None]*kernel.c[None,ej])/(kernel.degree[i]*kernel.degree[j]*kernel.sd[i]*kernel.sd[j]))

def response_from_cross(kernel,Q):
    """All normalized residue covariances from a supplied Cartesian lag covariance."""
    n=kernel.n;Q=Q.reshape(n,3,n,3).transpose(0,2,1,3);edges=kernel.edges;m=len(edges);KW=np.empty((m,n))
    u,v=edges[None,:,0],edges[None,:,1]
    for start in range(0,m,32):
        ix=np.arange(start,min(start+32,m));x,y=edges[ix,0,None],edges[ix,1,None]
        C=Q[x,u]-Q[x,v]-Q[y,u]+Q[y,v]
        K=gp.contracted_covariance(kernel.a[ix,None],kernel.a[None],kernel.S[ix,None],kernel.S[None],C)*kernel.c[ix,None]*kernel.c[None]
        KW[ix]=(kernel.W@K.T).T
    return (kernel.W@KW)/np.outer(kernel.sd,kernel.sd)

def drift_action(factor,model):
    B0=model['B'];L=factor['L'];mu=factor['stats']['mu'];d=B0.shape[1]
    def act(x):
        vec=np.asarray(x).ndim==1;x=np.asarray(x).reshape(d,-1)
        y=mu*mm(B0.T,L.T@(L@mm(B0,x)))
        return y.ravel() if vec else y
    return LinearOperator((d,d),matvec=act,rmatvec=act,matmat=act,dtype=float)

def physical_small_model(full):
    """First three real ABL residues, inherited contacts, first2 positive modes."""
    r=full['r0'][:3].copy();edges=full['edges'][(full['edges']<3).all(axis=1)].copy();K=np.zeros((9,9))
    for i,j in edges:
        a=r[i]-r[j];block=np.outer(a,a)/np.dot(a,a)
        for u,su in [(i,1),(j,-1)]:
            for v,sv in [(i,1),(j,-1)]:K[3*u:3*u+3,3*v:3*v+3]+=su*sv*block
    lam,B=eigh(K);keep=np.flatnonzero(lam>1e-8)[:2]
    model={'r0':r,'edges':edges,'B':B[:,keep],'eigenvalues':lam[keep],'canonical':full['canonical'][:3]}
    obj=gp.CovarianceObjective(r,edges,model['B'],model['eigenvalues']);Y,fit=gp.newton_cg(obj)
    if not fit['converged']:raise ArithmeticError('Physical small covariance fit failed')
    Sigma=Y/np.sqrt(np.outer(model['eigenvalues'],model['eigenvalues']))
    harmonic=OUKernel(model,np.diag(1/model['eigenvalues']),np.ones(3));sd=np.sqrt(np.diag(harmonic.kernel(0,False)))
    return model,Sigma,sd,fit
