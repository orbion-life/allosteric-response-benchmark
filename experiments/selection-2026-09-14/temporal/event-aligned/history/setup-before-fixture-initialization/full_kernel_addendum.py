"""New event-aligned full kernels; reuses fixed inputs without any fitting."""
from pathlib import Path
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[k]='3'
import sys,json,time,hashlib,platform,resource,numpy as np
ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent;sys.path.insert(0,str(BASE/'vendor'))
import gaussian_protein as gp

def sha(p):
 with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(p,x):Path(p).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def run(target,out,backend='cupy'):
 start=time.perf_counter();out=Path(out);out.mkdir(parents=True,exist_ok=False)
 protocol=json.loads((ROOT/'protocol.json').read_text());seal=json.loads((ROOT/'source-seal.json').read_text())
 for r in seal['files']:assert sha(BASE/r['file'])==r['sha256']
 manifest=json.loads((BASE/'input-manifest.json').read_text())
 for rec in manifest:
  if rec['target']==target:
   for f in rec['files']:assert sha(BASE/f['file'])==f['sha256']
 xp=gp.backend(backend);src=BASE/'inputs'/target;model=dict(np.load(src/'static-model.npz'));S=np.load(src/'fitted-covariance.npz')['covariance'];sd=np.load(src/'all-mode-harmonic-scales.npz')['harmonic_sd'];fixed=np.load(src/'fixed-operator.npz');tau=1/model['eigenvalues'][0]
 vals,V=xp.linalg.eigh(xp.asarray((S+S.T)/2));modes=gp.mm(xp.asarray(model['B']),V,xp);cart=gp.mm(modes*vals[None,:],modes.T,xp)
 def kernel(t):
  lag=cart if t==0 else gp.mm(modes*(vals*xp.exp(-float(t*tau)/vals))[None,:],modes.T,xp)
  K,_=gp.residue_covariance(model['r0'],model['edges'],cart,lag,xp,block=32);gp.sync(xp);return gp.host(K)/np.outer(sd,sd)
 parity=[]
 for i,t in enumerate([0.,.1,1.,10.]):
  K=kernel(t);ref=fixed['reference_G0'] if i==0 else fixed['reference_K'][i-1];parity.append(float(np.max(abs(K-ref))))
 dump(out/'canary.json',{'saved_CPU_max_abs':parity})
 if max(parity)>1e-8:raise ValueError('Known-query CPU parity failed')
 times=protocol['new_exact_kernel_times_over_tau'];n=len(sd);raw=np.lib.format.open_memmap(out/'new-full-kernels.npy',mode='w+',dtype='float64',shape=(len(times),n,n));cost=[]
 for i,t in enumerate(times):
  at=time.perf_counter();raw[i]=kernel(t);cost.append(time.perf_counter()-at)
  if i%20==0:raw.flush();print(target,i+1,len(times),round(time.perf_counter()-start,2),flush=True)
 raw.flush();del raw
 receipt={'status':'COMPLETE','target':target,'scope':protocol['scope'],'new_queries':len(times),'seconds':time.perf_counter()-start,'kernel_seconds':cost,'full_kernel_sha256':sha(out/'new-full-kernels.npy'),'bytes':(out/'new-full-kernels.npy').stat().st_size,'saved_CPU_parity_max_abs':max(parity),'source_sha256':sha(__file__),'protocol_sha256':sha(ROOT/'protocol.json'),'numpy':np.__version__,'python':platform.python_version(),'platform':platform.platform(),'peak_RSS_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024)}
 if backend=='cupy':receipt.update(GPU=xp.cuda.runtime.getDeviceProperties(0)['name'].decode(),cupy=xp.__version__,GPU_memory_pool_bytes=xp.get_default_memory_pool().total_bytes())
 dump(out/'receipt.json',receipt);return receipt
if __name__=='__main__':
 import argparse
 p=argparse.ArgumentParser();p.add_argument('--target',required=True);p.add_argument('--output',required=True);p.add_argument('--backend',default='numpy');a=p.parse_args();print(json.dumps(run(a.target,a.output,a.backend),indent=2))
