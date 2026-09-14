"""Fixed-grid signed step/pulse comparison and ranking from saved full kernels."""
from pathlib import Path
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[k]='3'
import numpy as np,json,time,hashlib,sys,resource
from scipy.linalg import eigh,expm
from scipy.linalg.blas import dgemm
ROOT=Path(__file__).resolve().parent

def sha(p):
 with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(p,x):Path(p).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def mm(a,b):return dgemm(1.,a,b)
def rank(C,meta):
 score=np.sqrt(np.mean(C[:,:,meta['receiver']]**2,axis=-1));cand=np.flatnonzero(meta['candidate']);ids=meta['canonical']
 order=np.array([cand[np.lexsort((ids[cand],-row[cand]))] for row in score]);return score,order

def diagnostics(ref,reduced,times,meta,threshold):
 error=np.max(abs(reduced-ref),axis=(1,2));eligible=abs(ref)>threshold
 sign_count=int(np.count_nonzero(eligible));wrong=int(np.count_nonzero(eligible & (np.sign(ref)!=np.sign(reduced))))
 scores,orders=rank(reduced,meta);refs,reforders=rank(ref,meta)
 gaps=np.array([r[o[4]]-r[o[5]] for r,o in zip(refs,reforders)])
 peaks=[];mask=meta['candidate'][:,None]&meta['receiver'][None,:]
 for i,j in zip(*np.where(mask)):
  a=abs(ref[:,i,j]);b=abs(reduced[:,i,j]);p=int(np.argmax(a));q=int(np.argmax(b));height=float(a[p]);ties=int(np.count_nonzero(abs(a-height)<=1e-12))
  resolved=height>threshold and p not in [0,len(times)-1] and ties==1
  peaks.append((i,j,p,q,height,int(resolved),float(abs(times[p]-times[q]))))
 peaks=np.array(peaks);good=peaks[:,5].astype(bool) if len(peaks) else np.zeros(0,bool)
 peakmetric={'eligible_source_receiver_pairs':len(peaks),'resolved_reference_interior_peaks':int(good.sum()),'low_flat_or_boundary_reference_peaks':int(len(peaks)-good.sum()),'maximum_sampled_peak_time_error_over_tau':float(peaks[good,6].max()) if good.any() else None,'fraction_same_sampled_peak_index':float(np.mean(peaks[good,2]==peaks[good,3])) if good.any() else None,'scope':'Discrete sampled absolute signed-entry peak only, no continuous-time timing guarantee'}
 rankrows=[{'time_over_tau':float(t),'reference_top5':meta['canonical'][ro[:5]].tolist(),'reduced_top5':meta['canonical'][o[:5]].tolist(),'top5_membership':len(set(o[:5])&set(ro[:5]))/5,'top5_order_identical':bool(np.array_equal(o[:5],ro[:5])),'reference_gap_5_6':float(gap),'gap_above_0_004':bool(gap>2*threshold),'reference_max_score':float(rs.max()),'low_signal':bool(rs.max()<=threshold)} for t,o,ro,gap,rs in zip(times,orders,reforders,gaps,refs)]
 metrics={'maximum_absolute_response_error':float(error.max()),'worst_time_over_tau':float(times[int(np.argmax(error))]),'passing_times':int(np.count_nonzero(error<=threshold)),'total_times':len(times),'error_by_time':error.tolist(),'all_times_pass':bool(np.all(error<=threshold)),'sign_eligible_entries':sign_count,'sign_disagreements':wrong,'sign_disagreement_fraction':wrong/sign_count if sign_count else None,'minimum_top5_membership':min(x['top5_membership'] for x in rankrows),'identical_top5_order_times':sum(x['top5_order_identical'] for x in rankrows),'peak':peakmetric,'ranking':rankrows}
 return metrics,{'score':scores,'reference_score':refs,'order':orders,'reference_order':reforders,'peaks':peaks}

def run(target):
 start=time.perf_counter();protocol=json.loads((ROOT/'protocol.json').read_text());src=ROOT/'inputs'/target;out=ROOT/'results'/target
 receipt=json.loads((out/'receipt.json').read_text());assert sha(out/'full-kernels.npy')==receipt['full_kernel_sha256']
 z=np.load(src/'fixed-operator.npz');meta=np.load(src/'ranking-metadata.npz');model=np.load(src/'static-model.npz');tau=1/model['eigenvalues'][0]
 ratios=np.array(protocol['all_exact_kernel_times_over_tau']);times=np.array(protocol['observation_times_over_tau']);full=np.load(out/'full-kernels.npy',mmap_mode='r');index={float(t):i for i,t in enumerate(ratios)}
 H=z['Hr'];B=z['B'];v,Q=eigh(H);QB=mm(Q.T,B)
 def reduced(t):return mm(QB.T*np.expm1(-float(t*tau)*v),QB) if t>0 else np.zeros_like(z['G0'])
 R=np.stack([reduced(t) for t in ratios]);F=np.stack([full[index[float(t)]]-full[index[0.]] for t in times]);G=R[[index[float(t)] for t in times]]
 metric,arrays=diagnostics(F,G,times,meta,.002)
 np.savez_compressed(out/'step.npz',full=F,reduced=G,times_over_tau=times,**arrays)
 pulse=[]
 for d in protocol['pulse_durations_over_tau']:
  Pf=F.copy();Pr=G.copy()
  for i,t in enumerate(times):
   if t>=d:
    Pf[i]-=full[index[float(t-d)]]-full[index[0.]];Pr[i]-=R[index[float(t-d)]]
  m,a=diagnostics(Pf,Pr,times,meta,.002);m['duration_over_tau']=d;pulse.append(m)
  np.savez_compressed(out/f'pulse-{d:g}.npz',full=Pf,reduced=Pr,times_over_tau=times,duration_over_tau=d,**a)
 verification=[]
 for t in [.015,.03,.1,1.,10.]:
  direct=mm(mm(B.T,expm(-t*tau*H)),B)-mm(B.T,B)
  verification.append({'time_over_tau':t,'independent_exponential_max_abs':float(np.max(abs(direct-R[index[t]])))})
 original=[]
 for k,t in enumerate([.1,1.,10.]):original.append({'time_over_tau':t,'full_saved_max_abs':float(np.max(abs(F[list(times).index(t)]-z['reference_C'][k]))),'reduced_saved_max_abs':float(np.max(abs(G[list(times).index(t)]-z['C'][k])))})
 result={'target':target,'status':'COMPLETE','scope':protocol['scope'],'step':metric,'pulses':pulse,'numerical_verification':verification,'saved_reference':original,'rank':len(H),'seconds':time.perf_counter()-start,'peak_RSS_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024),'source_sha256':sha(__file__),'protocol_sha256':sha(ROOT/'protocol.json')}
 dump(out/'analysis.json',result);return result
if __name__=='__main__':
 print(json.dumps(run(sys.argv[1]),indent=2))
