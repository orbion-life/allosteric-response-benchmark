from pathlib import Path
import numpy as np,json,hashlib,datetime
ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent
p=json.loads((BASE/'protocol.json').read_text());old=set(p['all_exact_kernel_times_over_tau']);lags=np.geomspace(.001,30,61).tolist();groups=[]
for d in p['pulse_durations_over_tau']:
 times=[float(d)]+[float(d+s) for s in lags if d+s<=30.]
 groups.append({'duration_over_tau':d,'observation_times_over_tau':times,'post_removal_lags_over_tau':[0.]+[float(s) for s in lags if d+s<=30.]})
needed=set(x for g in groups for x in g['observation_times_over_tau']) | set(x for g in groups for x in g['post_removal_lags_over_tau']) | {0.}
new=sorted(needed-old)
protocol={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'FROZEN_BEFORE_EVENT_ALIGNED_NEW_KERNELS_OR_SCORES','scope':'Separate developmental event-aligned characterization after the absolute-grid sampling gap was identified; not fresh holdout and no model/basis tuning','previous_protocol_sha256':hashlib.sha256((BASE/'protocol.json').read_bytes()).hexdigest(),'reason':'Absolute observation times did not densely sample the early lag after each pulse removal. Preserve that original grid and all original failures, then characterize fixed lag locations after each removal.','groups':groups,'new_exact_kernel_times_over_tau':new,'reuse_kernel_times_over_tau':sorted(needed&old),'time_upper_limit_over_tau':30.,'threshold':.002,'fixed_diagnostics':'same signed-entry, ranking, low-signal and fixed tie-rule definitions as the original protocol; no continuous-time timing claim','no_model_basis_or_outcome_tuning':True,'resource_cap':{'GPU':'one A100-40GB','timeout_seconds':1200,'CPU_cores':3,'host_memory_GiB':12,'automatic_retries':0},'prior_full_kernel_hashes':{t:json.loads((BASE/'results'/t/'receipt.json').read_text())['full_kernel_sha256'] for t in p['targets']}}
(ROOT/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
print({'groups':[len(x['observation_times_over_tau']) for x in groups],'new_kernels':len(new),'reused_kernels':len(needed&old)})
