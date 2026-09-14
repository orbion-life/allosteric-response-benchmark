"""Exact Gaussian-OU semigroup snapshots, distinct from nonlinear equilibrium."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[key]='1'
from pathlib import Path
import sys,json,time,resource,hashlib,itertools,argparse,warnings
import numpy as np
from scipy.linalg import eigh
from scipy.sparse import coo_matrix
from numpy.polynomial.hermite import hermgauss
import gaussian_kernel as gp
warnings.filterwarnings('error',category=RuntimeWarning)
ROOT=Path(__file__).resolve().parent

def mm(a,b):
    if a.ndim==2 and b.ndim==2 and not (np.iscomplexobj(a) or np.iscomplexobj(b)):
        from scipy.linalg.blas import dgemm
        return dgemm(1.,a,b)
    return np.einsum({(1,2):'i,ij->j',(2,1):'ij,j->i',(2,2):'ij,jk->ik'}[(a.ndim,b.ndim)],a,b,optimize=False)

def dump(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

class OUKernel:
    def __init__(self,model,Sigma,sd):
        self.model=model;self.Sigma=Sigma;self.sd=sd;self.n=len(model['r0']);self.d=len(Sigma)
        self.tau=1/model['eigenvalues'][0];self.beta=1.;self.mu=1.
        self.values,self.vectors=eigh(Sigma);self.modes=mm(model['B'],self.vectors)
        self.cart=mm(self.modes*self.values,self.modes.T)
        self.edges=model['edges'];self.r0=model['r0'];self.degree=np.bincount(self.edges.ravel(),minlength=self.n)
        self.a=self.r0[self.edges[:,0]]-self.r0[self.edges[:,1]];self.l2=np.sum(self.a*self.a,axis=1);self.c=1/(8*self.l2)
        self.W=coo_matrix((1/self.degree[self.edges.ravel()],(self.edges.ravel(),np.repeat(np.arange(len(self.edges)),2))),shape=(self.n,len(self.edges))).tocsr()
        V=self.cart.reshape(self.n,3,self.n,3).transpose(0,2,1,3);i,j=self.edges.T
        self.S=V[i,i]-V[i,j]-V[j,i]+V[j,j]
        tr=np.trace(self.S,axis1=1,axis2=2)
        self.means=self.W@(self.c*(4*np.einsum('ei,eij,ej->e',self.a,self.S,self.a)+2*np.sum(self.S*self.S,axis=(1,2))+tr*tr))
        self.cache={}

    def kernel(self,t,derivative=True):
        """Return normalized K(t) and dK/dt with fixed Sigma, original Ei."""
        decay=np.exp(-t/self.values)
        Q=mm(self.modes*(self.values*decay),self.modes.T)
        if derivative:
            Qprime=-mm(self.modes*decay,self.modes.T)
            Q=Q.astype(complex)+1j*1e-20*Qprime
        Q=Q.reshape(self.n,3,self.n,3).transpose(0,2,1,3)
        m=len(self.edges);KW=np.empty((m,self.n),dtype=Q.dtype)
        u,v=self.edges[None,:,0],self.edges[None,:,1]
        for start in range(0,m,32):
            ix=np.arange(start,min(start+32,m));x,y=self.edges[ix,0,None],self.edges[ix,1,None]
            C=Q[x,u]-Q[x,v]-Q[y,u]+Q[y,v]
            K=gp.contracted_covariance(self.a[ix,None],self.a[None],self.S[ix,None],self.S[None],C)*self.c[ix,None]*self.c[None]
            KW[ix]=(self.W@K.T).T
        result=(self.W@KW)/np.outer(self.sd,self.sd)
        if derivative:return result.real,result.imag/1e-20
        return result

    def conditional(self,q,t,gradient=False):
        """P_t f(q) as an explicit conditional Gaussian quartic expectation."""
        decay=np.exp(-t/self.values);Bt=mm(self.modes*decay,self.vectors.T)
        mean_displacement=mm(q,Bt.T).reshape(len(q),self.n,3)
        noise=mm(self.modes*(self.values*(1-decay*decay)),self.modes.T).reshape(self.n,3,self.n,3).transpose(0,2,1,3)
        E=np.zeros((len(q),self.n));dE=np.zeros((len(q),self.n,self.d)) if gradient else None
        for e,(i,j) in enumerate(self.edges):
            v=self.a[e]+mean_displacement[:,i]-mean_displacement[:,j]
            S=noise[i,i]-noise[i,j]-noise[j,i]+noise[j,j]
            tr=np.trace(S);d=np.sum(v*v,axis=1)-self.l2[e]
            u=self.c[e]*(d*d+2*d*tr+4*np.einsum('si,ij,sj->s',v,S,v)+tr*tr+2*np.sum(S*S))
            E[:,i]+=u/self.degree[i];E[:,j]+=u/self.degree[j]
            if gradient:
                T=Bt[3*i:3*i+3]-Bt[3*j:3*j+3]
                g=mm(4*self.c[e]*((d+tr)[:,None]*v+2*mm(v,S)),T)
                dE[:,i]+=g/self.degree[i];dE[:,j]+=g/self.degree[j]
        value=(E-self.means)/self.sd
        return (value,dE/self.sd[None,:,None]) if gradient else value

def load(name):
    if name=='triangle':
        sys.path.insert(0,str(ROOT.parent/'galerkin'))
        from run import triangle
        a,p=triangle();model={'r0':a.r0,'edges':a.edges,'B':a.B,'eigenvalues':a.lam};sd=a.sd
        obj=gp.CovarianceObjective(a.r0,a.edges,a.B,a.lam);Y,fit=gp.newton_cg(obj)
        Sigma=Y/np.sqrt(a.lam[:,None]*a.lam[None,:]);source={'origin':'same development triangle; all3 positive modes; fixed-centroid Gaussian fit','fit':fit}
    else:
        base=ROOT.parent/'galerkin/inputs'/name
        model=dict(np.load(base/'static-model.npz'));sd=np.load(base/'all-mode-harmonic-scales.npz')['harmonic_sd'];Sigma=np.load(base/'fitted-covariance.npz')['covariance']
        source={'model':str(base/'static-model.npz'),'model_sha256':sha(base/'static-model.npz'),'covariance_sha256':sha(base/'fitted-covariance.npz')}
    return OUKernel(model,Sigma,sd),source

def whiten(M,n):
    """Seed-first Gram whitening; snapshot rank uses only Gram eigenvalues."""
    G=M[:n,:n];g,U=eigh(G)
    if g[0]<=0:raise ValueError('Seed covariance is not positive definite')
    Q0=U/np.sqrt(g)
    cross=mm(Q0.T,M[:n,n:]);Schur=M[n:,n:]-mm(cross.T,cross);Schur=(Schur+Schur.T)/2
    values,V=eigh(Schur);threshold=max(1e-30,1e-10*max(0,float(values[-1])))
    keep=values>threshold;Z=V[:,keep]/np.sqrt(values[keep])
    W0=np.zeros((len(M),n));W0[:n]=Q0
    Wrest=np.vstack([-mm(Q0,mm(cross,Z)),Z]);W=np.concatenate([W0,Wrest],axis=1)
    mass=mm(W.T,mm(M,W));mass=(mass+mass.T)/2
    v,R=eigh(mass)
    if v[0]<=0:raise ValueError('Whitened mass lost positive definiteness')
    W=mm(W,mm(R/np.sqrt(v),R.T))
    return W,{'rank':W.shape[1],'seed_rank':n,'snapshot_rank':int(keep.sum()),'Gram_cutoff':threshold,'minimum_seed_eigenvalue':float(g[0]),'minimum_Schur_eigenvalue':float(values[0]),'mass_error':float(np.max(abs(mm(W.T,mm(M,W))-np.eye(W.shape[1]))))}

def verify_conditionals(kernel):
    q=np.random.default_rng(7281).normal(size=(3,kernel.d))*.1;t=.1*kernel.tau
    value,gradient=kernel.conditional(q,t,True);h=1e-5;finite=np.empty_like(gradient)
    for k in range(kernel.d):
        step=np.eye(kernel.d)[k]*h;finite[:,:,k]=(kernel.conditional(q+step,t)-kernel.conditional(q-step,t))/(2*h)
    # Exact-degree Gaussian quadrature independently checks conditional values.
    x,w=hermgauss(5);indices=np.array(list(itertools.product(range(5),repeat=kernel.d)));nodes=np.sqrt(2)*x[indices];weights=np.prod(w[indices],axis=1)/np.pi**(kernel.d/2)
    decay=np.exp(-t/kernel.values);P=mm(kernel.vectors*decay,kernel.vectors.T)
    root=kernel.vectors*np.sqrt(kernel.values*(1-decay*decay))
    independent=[]
    for point in q:
        samples=mm(nodes,root.T)+mm(point,P.T)
        independent.append(mm(weights,kernel.conditional(samples,0)))
    return {'conditional_gradient_max_abs':float(np.max(abs(finite-gradient))),'conditional_Gauss_Hermite_max_abs':float(np.max(abs(value-np.array(independent))))}

def run(name):
    start=time.perf_counter();out=ROOT/name;out.mkdir(exist_ok=False);kernel,source=load(name);N=kernel.n;tau=kernel.tau
    checks=verify_conditionals(kernel) if name=='triangle' else {}
    if checks and max(checks.values())>1e-8:raise ValueError('Conditional snapshot adapter failed')
    ratios=[0.,.1,.2,1.,1.1,2.];cache={};timings=[]
    for u in ratios:
        at=time.perf_counter();K,Kprime=kernel.kernel(u*tau);cache[u]=(K,Kprime);timings.append({'time_over_tau':u,'seconds':time.perf_counter()-at})
        print(name,'kernel',u,flush=True)
    nodes=[0.,.1,1.]
    M=np.block([[cache[round(a+b,8)][0] for b in nodes] for a in nodes]);A=np.block([[-cache[round(a+b,8)][1] for b in nodes] for a in nodes]);D=np.vstack([cache[a][0] for a in nodes]);M=(M+M.T)/2;A=(A+A.T)/2
    W,diagnostic=whiten(M,N);Hr=mm(W.T,mm(A,W));Hr=(Hr+Hr.T)/2;B=mm(W.T,D);values,Q=eigh(Hr)
    psd_allowance=1e-10*max(1,float(values[-1]));psd=values[0]>=-psd_allowance
    if not psd:raise ValueError('Projected Gaussian relaxation matrix is not PSD within allowance')
    static=mm(B.T,B);times=np.array([.1,1,10])*tau
    projected=np.array([mm(mm(B.T,mm(Q*np.exp(-t*values),Q.T)),B) for t in times]);ref10=kernel.kernel(times[-1],False);reference=np.array([cache[.1][0],cache[1.][0],ref10]);response=projected-static;Cref=reference-cache[0.][0]
    # Independent fourth-order time finite difference at the fixed .1tau check.
    h=1e-4*tau;t=.1*tau
    derivative=(-kernel.kernel(t+2*h,False)+8*kernel.kernel(t+h,False)-8*kernel.kernel(t-h,False)+kernel.kernel(t-2*h,False))/(12*h)
    checks['tau_scaled_derivative_max_abs']=float(tau*np.max(abs(derivative-cache[.1][1])))
    checks['static_containment']=float(np.max(abs(static-cache[0.][0])))
    checks['delayed_response']=float(np.max(abs(projected-reference)))
    checks['response']=float(np.max(abs(response-Cref)))
    checks['minimum_Hr_eigenvalue']=float(values[0]);checks['mass_orthogonality']=diagnostic['mass_error']
    np.savez_compressed(out/'result.npz',M=M,A=A,D=D,W=W,Hr=Hr,B=B,times=times,G0=static,K=projected,C=response,reference_G0=cache[0.][0],reference_K=reference,reference_C=Cref,kernel_times=np.array(ratios)*tau,kernels=np.array([cache[u][0] for u in ratios]),kernel_derivatives=np.array([cache[u][1] for u in ratios]),derivative_fd=derivative,sd=kernel.sd,Sigma=kernel.Sigma)
    numerical=checks['tau_scaled_derivative_max_abs']<=1e-8 and checks['static_containment']<=1e-10 and checks['mass_orthogonality']<=1e-8
    fidelity=max(checks['delayed_response'],checks['response'])<=.002
    result={'engineering_status':'COMPLETE','scientific_status':'PASS' if numerical and fidelity else 'FAIL','numerical_checks':'PASS' if numerical else 'FAIL','surrogate_response_fidelity':'PASS' if fidelity else 'FAIL','scope':'Exact fixed-centroid Gaussian-OU surrogate only; not original nonlinear density, physical accuracy, biological advantage or quantum advantage','source':source,'N':N,'physical_coordinates':kernel.d,'snapshot_nodes_over_tau':nodes,'basis':diagnostic,'checks':checks,'kernel_timings':timings,'wall_seconds_including_construction_reference_checks':time.perf_counter()-start,'peak_process_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'source_sha256':sha(__file__),'kernel_source_sha256':sha(ROOT/'gaussian_kernel.py'),'protocol_sha256':sha(ROOT/'protocol.json')}
    dump(out/'summary.json',result);print(json.dumps(result,indent=2),flush=True)
    return result

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('target',choices=['triangle','kras','abl']);args=parser.parse_args();run(args.target)
