"""Portable fixed-covariance full Gaussian time characterization; never fits."""
from pathlib import Path
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[key]='3'
import sys,json,time,hashlib,platform,resource,gc
import numpy as np
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'vendor'))
import gaussian_protein as gp

def sha(p):
 with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(p,x):Path(p).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def run(target,out,backend='cupy'):
 out=Path(out);out.mkdir(parents=True,exist_ok=False);start=time.perf_counter()
 protocol=json.loads((ROOT/'protocol.json').read_text());manifest=json.loads((ROOT/'input-manifest.json').read_text())
 for rec in manifest:
  if rec['target']==target:
   for f in rec['files']:assert sha(ROOT/f['file'])==f['sha256']
 seal=json.loads((ROOT/'source-seal.json').read_text())
 for row in seal['files']:assert sha(ROOT/row['file'])==row['sha256']
 xp=gp.backend(backend)
 from check_independent import checks
 canary,fixtures=checks(xp)
 np.savez_compressed(out/'fixture-arrays.npz',**fixtures);dump(out/'fixture-check.json',canary)
 assert canary['passed'] and canary.get('cpu_gpu_passed',True)
 src=ROOT/'inputs'/target;model=dict(np.load(src/'static-model.npz'));sigma=np.load(src/'fitted-covariance.npz')['covariance']
 sd=np.load(src/'all-mode-harmonic-scales.npz')['harmonic_sd'];old=np.load(src/'fixed-operator.npz');tau=1/model['eigenvalues'][0]
 values,vectors=xp.linalg.eigh(xp.asarray((sigma+sigma.T)/2));modes=gp.mm(xp.asarray(model['B']),vectors,xp)
 cart=gp.mm(modes*values[None,:],modes.T,xp);scale=np.outer(sd,sd)
 def kernel(ratio):
  t=float(ratio*tau);lag=cart if ratio==0 else gp.mm(modes*(values*xp.exp(-t/values))[None,:],modes.T,xp)
  cov,_=gp.residue_covariance(model['r0'],model['edges'],cart,lag,xp,block=32)
  gp.sync(xp);return gp.host(cov)/scale
 cached={};parity=[]
 for i,u in enumerate([0.,.1,1.,10.]):
  begin=time.perf_counter();K=kernel(u);cached[u]=K
  ref=old['reference_G0'] if i==0 else old['reference_K'][i-1]
  parity.append({'time_over_tau':u,'saved_CPU_max_abs':float(np.max(abs(K-ref))),'seconds':time.perf_counter()-begin})
 dump(out/'saved-reference-canary.json',parity)
 if max(r['saved_CPU_max_abs'] for r in parity)>1e-8:raise ValueError('CPU reference parity failed before new queries')
 ratios=protocol['all_exact_kernel_times_over_tau'];n=len(sd)
 raw=np.lib.format.open_memmap(out/'full-kernels.npy',mode='w+',dtype='float64',shape=(len(ratios),n,n));cost=[]
 for i,u in enumerate(ratios):
  begin=time.perf_counter();raw[i]=cached[u] if u in cached else kernel(u);cost.append(time.perf_counter()-begin)
  if not np.isfinite(raw[i]).all():raise ValueError('Nonfinite kernel')
  if i%10==0:
   raw.flush();dump(out/'progress.json',{'target':target,'queries_complete':i+1,'total':len(ratios),'seconds':time.perf_counter()-start});print(target,i+1,len(ratios),round(time.perf_counter()-start,2),flush=True)
 raw.flush();del raw
 receipt={'status':'COMPLETE','target':target,'scope':protocol['scope'],'queries':len(ratios),'N':n,'tau':float(tau),'seconds':time.perf_counter()-start,'kernel_seconds':cost,'source_sha256':sha(__file__),'protocol_sha256':sha(ROOT/'protocol.json'),'numpy':np.__version__,'python':platform.python_version(),'platform':platform.platform(),'backend':backend,'peak_RSS_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024),'full_kernel_sha256':sha(out/'full-kernels.npy'),'bytes':(out/'full-kernels.npy').stat().st_size,'saved_CPU_parity_max_abs':max(r['saved_CPU_max_abs'] for r in parity),'provider_task_id':os.environ.get('MODAL_TASK_ID')}
 if backend=='cupy':receipt.update(cupy=xp.__version__,GPU=xp.cuda.runtime.getDeviceProperties(0)['name'].decode(),GPU_memory_pool_bytes=xp.get_default_memory_pool().total_bytes())
 dump(out/'receipt.json',receipt);return receipt
if __name__=='__main__':
 import argparse
 p=argparse.ArgumentParser();p.add_argument('--target',required=True);p.add_argument('--output',required=True);p.add_argument('--backend',default='numpy');a=p.parse_args();print(json.dumps(run(a.target,a.output,a.backend),indent=2))
