"""Offline event-aligned source, protocol and complete saved-array verification."""
from pathlib import Path
import numpy as np,json,hashlib,argparse,sys,datetime
ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent

def sha(p):
 with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def protocol_checks():
 p=json.loads((ROOT/'protocol.json').read_text());s=json.loads((ROOT/'source-seal.json').read_text());old=json.loads((BASE/'protocol.json').read_text());checks={'sealed_sources':all(sha(BASE/x['file'])==x['sha256'] for x in s['files']),'prior_protocol':sha(BASE/'protocol.json')==p['previous_protocol_sha256'],'new_and_old_disjoint':set(p['new_exact_kernel_times_over_tau']).isdisjoint(old['all_exact_kernel_times_over_tau']),'same_error_allocation':p['threshold']==.002,'no_tuning':p['no_model_basis_or_outcome_tuning']}
 lags=np.geomspace(.001,30,61)
 checks['fixed_lag_rule']=all(g['post_removal_lags_over_tau']==[0.]+[float(s) for s in lags if g['duration_over_tau']+s<=30.] for g in p['groups'])
 checks['exact_time_construction']=all(np.array_equal(np.array(g['observation_times_over_tau']),g['duration_over_tau']+np.array(g['post_removal_lags_over_tau'])) for g in p['groups'])
 return checks

def verify_target(target):
 p=json.loads((ROOT/'protocol.json').read_text());o=json.loads((BASE/'protocol.json').read_text());out=ROOT/'results'/target;r=json.loads((out/'receipt.json').read_text());x=json.loads((out/'analysis.json').read_text())
 oldpath=BASE/'results'/target/'full-kernels.npy';newpath=out/'new-full-kernels.npy'
 hashes=sha(oldpath)==p['prior_full_kernel_hashes'][target] and sha(newpath)==r['full_kernel_sha256'];old=np.load(oldpath,mmap_mode='r');new=np.load(newpath,mmap_mode='r');oidx={float(t):i for i,t in enumerate(o['all_exact_kernel_times_over_tau'])};nidx={float(t):i for i,t in enumerate(p['new_exact_kernel_times_over_tau'])}
 def K(t):return old[oidx[t]] if t in oidx else new[nidx[t]]
 rows=[];meta=np.load(BASE/'inputs'/target/'ranking-metadata.npz');candidates=np.flatnonzero(meta['candidate']);ids=meta['canonical']
 for g,m in zip(p['groups'],x['pulses']):
  d=g['duration_over_tau'];z=np.load(out/f'pulse-{d:g}.npz');full=z['full'];red=z['reduced'];check=0.
  for i,(t,s) in enumerate(zip(g['observation_times_over_tau'],g['post_removal_lags_over_tau'])):check=max(check,float(np.max(abs(full[i]-(K(t)-K(s))))))
  e=np.max(abs(red-full),axis=(1,2));mask=abs(full)>.002;wrong=int(np.count_nonzero(mask & (np.sign(full)!=np.sign(red))))
  ranks=[]
  for arr in [red,full]:
   score=np.sqrt(np.mean(arr[:,:,meta['receiver']]**2,axis=-1));ranks.append(np.array([candidates[np.lexsort((ids[candidates],-r[candidates]))] for r in score]))
  rows.append({'duration_over_tau':d,'raw_full_identity_max_abs':check,'errors_agree':bool(np.allclose(e,m['error_by_time'],rtol=0,atol=1e-12)),'sign_count_agrees':wrong==m['sign_disagreements'],'rank_arrays_agree':bool(np.array_equal(ranks[0],z['order']) and np.array_equal(ranks[1],z['reference_order'])),'scientific_all_times_pass':bool(np.all(e<=.002))})
 ok=hashes and r['saved_CPU_parity_max_abs']<=1e-8 and max(z['independent_exponential_max_abs'] for z in x['independent_verification'])<=1e-8 and all(v['raw_full_identity_max_abs']<=1e-12 and v['errors_agree'] and v['sign_count_agrees'] and v['rank_arrays_agree'] for v in rows)
 return {'target':target,'implementation_checks_pass':ok,'input_hashes_pass':hashes,'pulses':rows}
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--mode',choices=['protocol','saved'],default='saved');p.add_argument('--output',required=True);a=p.parse_args();r={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'scope':'Implementation verification does not turn scientific temporal failures into passes','protocol':protocol_checks()}
 if a.mode=='saved':r['targets']=[verify_target(t) for t in ['kras','abl','myc','myh7']]
 r['passed']=all(r['protocol'].values()) and all(t['implementation_checks_pass'] for t in r.get('targets',[]));Path(a.output).write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2));sys.exit(0 if r['passed'] else 1)
