"""Response and candidate ranking comparison on frozen old and fresh queries."""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[k]='3'
from pathlib import Path
import sys,json,time,hashlib,resource
import numpy as np
from scipy.linalg.blas import dgemm
ROOT=Path(__file__).resolve().parent
P=json.loads((ROOT/'protocol.json').read_text());OLD=Path(P['original_root'])
def sha(p):
    with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(p,x):Path(p).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')

def run(target):
    start=time.perf_counter();out=ROOT/'results'/target;z=np.load(out/'operator.npz');receipt=json.loads((out/'receipt.json').read_text())
    assert sha(out/'operator.npz')==receipt['operator_sha256'];assert sha(out/'fresh-kernels.npy')==receipt['fresh_kernel_sha256']
    op=json.loads((OLD/'protocol.json').read_text());ep=json.loads((OLD/'event-aligned/protocol.json').read_text())
    oldpath=OLD/'results'/target/'full-kernels.npy';eventpath=OLD/'event-aligned/results'/target/'new-full-kernels.npy'
    er=json.loads((eventpath.parent/'receipt.json').read_text());assert sha(oldpath)==ep['prior_full_kernel_hashes'][target];assert sha(eventpath)==er['full_kernel_sha256']
    old=np.load(oldpath,mmap_mode='r');event=np.load(eventpath,mmap_mode='r');fresh=np.load(out/'fresh-kernels.npy',mmap_mode='r')
    oi={float(t):i for i,t in enumerate(op['all_exact_kernel_times_over_tau'])};ei={float(t):i for i,t in enumerate(ep['new_exact_kernel_times_over_tau'])};fi={float(t):i for i,t in enumerate(P['fresh_kernel_times_over_tau'])}
    meta=np.load(ROOT/'inputs'/target/'ranking-metadata.npz');cand=np.flatnonzero(meta['candidate']);rec=np.flatnonzero(meta['receiver']);ids=meta['canonical'];tau=float(z['tau']);lam=z['lambda_'];QB=dgemm(1.,z['Q'].T,z['B']);G0=old[oi[0.]];cache={}
    def K(t):
        if t in fi:return fresh[fi[t]]
        if t in oi:return old[oi[t]]
        return event[ei[t]]
    def C(t):
        if t not in cache:cache[t]=dgemm(1.,QB.T*np.expm1(-t*tau*lam),QB) if t else np.zeros_like(G0)
        return cache[t]
    def metric(times,lags=None):
        rows=[]
        for it,t in enumerate(times):
            t=float(t);f=K(t)-G0;r=C(t)
            if lags is not None:
                lag=float(lags[it]);f=f-(K(lag)-G0);r=r-C(lag)
            e=abs(r-f);fs=np.sqrt(np.mean(f[:,rec]**2,axis=1));rs=np.sqrt(np.mean(r[:,rec]**2,axis=1));fo=cand[np.lexsort((ids[cand],-fs[cand]))];ro=cand[np.lexsort((ids[cand],-rs[cand]))];membership=len(set(fo[:5])&set(ro[:5]))/5
            row={'time_over_tau':t,'max_entry_error':float(e.max()),'candidate_receiver_max_entry_error':float(e[np.ix_(cand,rec)].max()),'candidate_score_max_abs':float(abs(fs[cand]-rs[cand]).max()),'top5_set_agreement':membership,'top5_order_identical':bool(np.array_equal(fo[:5],ro[:5])),'reference_top5':ids[fo[:5]].tolist(),'reduced_top5':ids[ro[:5]].tolist(),'reference_gap_5_6':float(fs[fo[4]]-fs[fo[5]]),'reference_max_candidate_score':float(fs[cand].max()),'low_signal':bool(fs[cand].max()<=.002)}
            if lags is not None:row['post_removal_lag_over_tau']=float(lags[it])
            rows.append(row)
        passing=sum(r['max_entry_error']<=.002 for r in rows);eligible=[r for r in rows if not r['low_signal']]
        separated=[r for r in eligible if r['reference_gap_5_6']>.004]
        return {'query_count':len(rows),'all_entry_pass_count':passing,'all_entries_pass':passing==len(rows),'max_entry_error':max(r['max_entry_error'] for r in rows),'candidate_receiver_max_entry_error':max(r['candidate_receiver_max_entry_error'] for r in rows),'candidate_score_max_abs':max(r['candidate_score_max_abs'] for r in rows),'minimum_top5_set_agreement':min(r['top5_set_agreement'] for r in rows),'identical_top5_set_count':sum(r['top5_set_agreement']==1 for r in rows),'identical_top5_order_count':sum(r['top5_order_identical'] for r in rows),'non_low_signal_query_count':len(eligible),'non_low_signal_minimum_top5_set_agreement':min((r['top5_set_agreement'] for r in eligible),default=None),'non_low_signal_identical_top5_set_count':sum(r['top5_set_agreement']==1 for r in eligible),'non_low_signal_identical_top5_order_count':sum(r['top5_order_identical'] for r in eligible),'non_low_signal_minimum_reference_gap_5_6':min((r['reference_gap_5_6'] for r in eligible),default=None),'gap_above_0_004_query_count':len(separated),'gap_above_0_004_identical_top5_set_count':sum(r['top5_set_agreement']==1 for r in separated),'rows':rows}
    results={}
    results['fresh_step']=metric(P['fresh_step_times_over_tau']);results['fresh_removal']=[]
    for d in P['fresh_removal_durations_over_tau']:
        q=metric([d+u for u in P['fresh_post_removal_lags_over_tau']],P['fresh_post_removal_lags_over_tau']);q['duration_over_tau']=d;results['fresh_removal'].append(q)
    print(target,'fresh analyzed',flush=True)
    results['selected_times']=metric([.1,1.,10.]);results['historical_step']=metric(op['observation_times_over_tau']);results['historical_pulse']=[]
    for d in op['pulse_durations_over_tau']:
        q=metric(op['observation_times_over_tau'],[max(0.,float(t)-d) for t in op['observation_times_over_tau']]);q['duration_over_tau']=d;results['historical_pulse'].append(q)
    results['historical_event_aligned']=[]
    for g in ep['groups']:
        q=metric(g['observation_times_over_tau'],g['post_removal_lags_over_tau']);q['duration_over_tau']=g['duration_over_tau'];results['historical_event_aligned'].append(q)
    checks=receipt['checks'];num=checks['mass_error']<=1e-8 and checks['static_error']<=1e-10 and checks['minimum_H_eigenvalue']>=-checks['PSD_allowance'] and checks['independent_exponential_max_abs']<=1e-8
    allgroups=[results['fresh_step']]+results['fresh_removal']+[results['historical_step']]+results['historical_pulse']+results['historical_event_aligned']
    result={'target':target,'status':'COMPLETE','numerical_checks_pass':num,'all_sampled_fidelity_groups_pass':all(r['all_entries_pass'] for r in allgroups),'rank':len(lam),'N':len(G0),'checks':checks,'results':results,'seconds':time.perf_counter()-start,'peak_RSS_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024),'protocol_sha256':sha(ROOT/'protocol.json'),'source_sha256':sha(__file__),'inputs':{'old_full_kernel_sha256':sha(oldpath),'event_full_kernel_sha256':sha(eventpath),'fresh_full_kernel_sha256':sha(out/'fresh-kernels.npy')}}
    dump(out/'analysis.json',result);print(target,'COMPLETE',result['all_sampled_fidelity_groups_pass'],round(result['seconds'],2),flush=True);return result
if __name__=='__main__':run(sys.argv[1])
