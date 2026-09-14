"""Frozen six-node Gaussian basis construction and fresh query kernels."""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:
    os.environ[k]='3'
from pathlib import Path
import sys,time,json,hashlib,platform,resource
import numpy as np
from scipy.linalg import eigh,expm
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'vendor'))
import gaussian_protein as gp
from fixed_adapter import OUKernel

def sha(p):
    with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(p,x):Path(p).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')

def whiten(M,n):
    G=M[:n,:n];g,U=eigh(G)
    if g[0]<=0:raise ValueError('Seed covariance not positive definite')
    Q0=U/np.sqrt(g);cross=Q0.T@M[:n,n:]
    S=M[n:,n:]-cross.T@cross;S=(S+S.T)/2
    v,V=eigh(S);largest=float(eigh(M,eigvals_only=True,subset_by_index=[len(M)-1,len(M)-1])[0])
    threshold=max(1e-30,1e-10*largest);keep=v>threshold;Z=V[:,keep]/np.sqrt(v[keep])
    W0=np.zeros((len(M),n));W0[:n]=Q0
    W=np.concatenate([W0,np.vstack([-Q0@cross@Z,Z])],axis=1)
    mass=W.T@M@W;mass=(mass+mass.T)/2;g,R=eigh(mass)
    if g[0]<=0:raise ValueError('Whitened mass not positive')
    W=W@(R/np.sqrt(g))@R.T
    return W,{'rank':W.shape[1],'seed_rank':n,'snapshot_rank':int(keep.sum()),'cutoff':threshold,'full_Gram_largest_eigenvalue':largest,'minimum_Schur_eigenvalue':float(v[0]),'mass_error':float(np.max(abs(W.T@M@W-np.eye(W.shape[1]))))}

def run(target,out):
    import cupy as xp
    out=Path(out);out.mkdir(parents=True,exist_ok=False);start=time.perf_counter()
    p=json.loads((ROOT/'protocol.json').read_text());src=ROOT/'inputs'/target
    model=dict(np.load(src/'static-model.npz'));sigma=np.load(src/'fitted-covariance.npz')['covariance'];sd=np.load(src/'all-mode-harmonic-scales.npz')['harmonic_sd']
    tau=1/model['eigenvalues'][0];n=len(sd);scale=np.outer(sd,sd)
    values,vectors=xp.linalg.eigh(xp.asarray((sigma+sigma.T)/2));modes=xp.asarray(model['B'])@vectors
    cart=(modes*values[None,:])@modes.T
    def kernel(u,derivative=False):
        decay=xp.exp(-float(u*tau)/values);lag=(modes*(values*decay)[None,:])@modes.T
        if derivative:lag=lag.astype(xp.complex128)+1j*1e-20*(-(modes*decay[None,:])@modes.T)
        cov,_=gp.residue_covariance(model['r0'],model['edges'],cart,lag,xp,block=32);gp.sync(xp);z=gp.host(cov)/scale
        return (z.real,z.imag/1e-20) if derivative else z
    # Before construction, compare independent CPU OUKernel and derivative at 0.01tau.
    cpu=OUKernel(model,sigma,sd);at=time.perf_counter();kc,dc=cpu.kernel(.01*tau);kg,dg=kernel(.01,True)
    canary={'CPU_GPU_K_max_abs':float(np.max(abs(kc-kg))),'CPU_GPU_tau_derivative_max_abs':float(tau*np.max(abs(dc-dg))),'seconds':time.perf_counter()-at}
    refs=np.load(src/'reference-canary.npz');cache={}
    for i,u in enumerate([0.,.1,1.,10.]):
        z=kernel(u);cache[u]=z;ref=refs['G0'] if i==0 else refs['K'][i-1]
        canary[str(u)+'_saved_reference_max_abs']=float(np.max(abs(z-ref)))
    dump(out/'canary.json',canary)
    if max(v for k,v in canary.items() if k!='seconds')>1e-8:raise ValueError('Independent CPU/GPU parity failed')
    deriv={};cost=[]
    for u in p['basis_kernel_times_over_tau']:
        at=time.perf_counter();cache[u],deriv[u]=kernel(u,True);cost.append({'ratio':u,'seconds':time.perf_counter()-at});print(target,'basis',u,flush=True)
    np.savez_compressed(out/'basis-kernels.npz',times=np.array(p['basis_kernel_times_over_tau']),K=np.array([cache[u] for u in p['basis_kernel_times_over_tau']]),derivative=np.array([deriv[u] for u in p['basis_kernel_times_over_tau']]))
    at=time.perf_counter();nodes=p['basis_nodes_over_tau'];M=np.block([[cache[round(a+b,10)] for b in nodes] for a in nodes]);A=np.block([[-deriv[round(a+b,10)] for b in nodes] for a in nodes]);D=np.vstack([cache[a] for a in nodes]);M=(M+M.T)/2;A=(A+A.T)/2
    W,diag=whiten(M,n);H=W.T@A@W;H=(H+H.T)/2;B=W.T@D;lam,Q=eigh(H);static=B.T@B
    checks={'mass_error':diag['mass_error'],'static_error':float(np.max(abs(static-cache[0.]))),'minimum_H_eigenvalue':float(lam[0]),'maximum_H_eigenvalue':float(lam[-1]),'PSD_allowance':1e-10*max(1.,float(lam[-1]))}
    independent=[]
    for u in [.00031622776601683794,.03162277660168379,3.1622776601683795]:
        direct=B.T@expm(-u*tau*H)@B;spectral=(B.T@Q*np.exp(-u*tau*lam))@(Q.T@B)
        independent.append(float(np.max(abs(direct-spectral))))
    checks['independent_exponential_max_abs']=max(independent);reduction_seconds=time.perf_counter()-at
    np.savez_compressed(out/'operator.npz',M=M,A=A,D=D,W=W,Hr=H,B=B,lambda_=lam,Q=Q,G0=static,reference_G0=cache[0.],tau=tau)
    dump(out/'construction.json',{'basis':diag,'checks':checks,'reduction_seconds':reduction_seconds,'kernel_cost':cost})
    fresh=p['fresh_kernel_times_over_tau'];raw=np.lib.format.open_memmap(out/'fresh-kernels.npy',mode='w+',dtype='float64',shape=(len(fresh),n,n));cost=[]
    for i,u in enumerate(fresh):
        at=time.perf_counter();raw[i]=kernel(u);cost.append(time.perf_counter()-at)
        if i%12==0:print(target,'fresh',i+1,len(fresh),round(time.perf_counter()-start,1),flush=True)
    raw.flush();del raw
    receipt={'target':target,'N':n,'tau':float(tau),'wall_seconds':time.perf_counter()-start,'reduction_seconds':reduction_seconds,'kernel_seconds':cost,'rank':diag['rank'],'checks':checks,'canary':canary,'protocol_sha256':sha(ROOT/'protocol.json'),'source_sha256':sha(__file__),'input_files':{x.name:sha(x) for x in src.iterdir()},'fresh_kernel_sha256':sha(out/'fresh-kernels.npy'),'operator_sha256':sha(out/'operator.npz'),'peak_RSS_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,'GPU':xp.cuda.runtime.getDeviceProperties(0)['name'].decode(),'GPU_memory_pool_bytes':xp.get_default_memory_pool().total_bytes(),'numpy':np.__version__,'cupy':xp.__version__,'scipy':__import__('scipy').__version__,'python':platform.python_version(),'provider_task_id':os.environ.get('MODAL_TASK_ID')}
    dump(out/'receipt.json',receipt);return receipt
