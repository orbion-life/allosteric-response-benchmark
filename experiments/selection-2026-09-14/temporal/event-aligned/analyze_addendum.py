"""Separate frozen event-aligned pulse metrics, preserving original results."""
from pathlib import Path
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[k]='3'
import sys,json,time,hashlib,numpy as np
from scipy.linalg import eigh,expm
ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent;sys.path.insert(0,str(BASE))
from analyze import diagnostics,mm

def sha(p):
 with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(p,x):Path(p).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def run(target):
 start=time.perf_counter();p=json.loads((ROOT/'protocol.json').read_text());oldp=json.loads((BASE/'protocol.json').read_text());out=ROOT/'results'/target;src=BASE/'inputs'/target
 oldpath=BASE/'results'/target/'full-kernels.npy';newpath=out/'new-full-kernels.npy';r=json.loads((out/'receipt.json').read_text())
 assert sha(oldpath)==p['prior_full_kernel_hashes'][target];assert sha(newpath)==r['full_kernel_sha256']
 old=np.load(oldpath,mmap_mode='r');new=np.load(newpath,mmap_mode='r');oldidx={float(t):i for i,t in enumerate(oldp['all_exact_kernel_times_over_tau'])};newidx={float(t):i for i,t in enumerate(p['new_exact_kernel_times_over_tau'])}
 def kernel(t):return old[oldidx[t]] if t in oldidx else new[newidx[t]]
 z=np.load(src/'fixed-operator.npz');meta=np.load(src/'ranking-metadata.npz');model=np.load(src/'static-model.npz');tau=1/model['eigenvalues'][0];H=z['Hr'];B=z['B'];v,Q=eigh(H);QB=mm(Q.T,B)
 def reducedK(t):return mm(QB.T*np.exp(-t*tau*v),QB)
 rows=[];verify=[]
 for group in p['groups']:
  d=group['duration_over_tau'];times=np.array(group['observation_times_over_tau']);lags=np.array(group['post_removal_lags_over_tau'])
  full=np.array([kernel(float(t))-kernel(float(lag)) for t,lag in zip(times,lags)])
  red=np.array([mm(QB.T*(np.exp(-t*tau*v)-np.exp(-lag*tau*v)),QB) for t,lag in zip(times,lags)])
  metric,arrays=diagnostics(full,red,times,meta,.002);metric['duration_over_tau']=d;metric['post_removal_lags_over_tau']=lags.tolist();metric['worst_post_removal_lag_over_tau']=float(lags[np.argmax(metric['error_by_time'])]);metric['roundoff_time_closure_max_abs']=float(np.max(abs(times-d-lags)));rows.append(metric)
  np.savez_compressed(out/f'pulse-{d:g}.npz',full=full,reduced=red,times_over_tau=times,post_removal_lags_over_tau=lags,duration_over_tau=d,**arrays)
  # Independent nonspectral exponential at the first positive lag, predeclared.
  t,lag=float(times[1]),float(lags[1]);ind=mm(mm(B.T,expm(-t*tau*H)-expm(-lag*tau*H)),B)
  verify.append({'duration_over_tau':d,'lag_over_tau':lag,'independent_exponential_max_abs':float(np.max(abs(ind-red[1])))})
 result={'target':target,'status':'COMPLETE','scope':p['scope'],'pulses':rows,'independent_verification':verify,'seconds':time.perf_counter()-start,'source_sha256':sha(__file__),'protocol_sha256':sha(ROOT/'protocol.json')};dump(out/'analysis.json',result);return result
if __name__=='__main__':print(json.dumps(run(sys.argv[1]),indent=2))
