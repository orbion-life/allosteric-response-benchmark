"""Compact reporting of frozen metrics; supplementary descriptions do not change gates."""
from pathlib import Path
import json,numpy as np,hashlib,datetime
ROOT=Path(__file__).resolve().parent

def main():
 p=json.loads((ROOT/'protocol.json').read_text());times=np.array(p['observation_times_over_tau']);rows=[]
 for target in p['targets']:
  out=ROOT/'results'/target;x=json.loads((out/'analysis.json').read_text());gpu=json.loads((out/'receipt.json').read_text())
  def compact(m,d=None):
   errs=np.array(m['error_by_time']);rank=m['ranking'];signal=[r for r in rank if not r['low_signal']]
   after=times>d if d is not None else np.ones(len(times),bool)
   return {'duration_over_tau':d,'max_absolute_error':m['maximum_absolute_response_error'],'worst_time_over_tau':m['worst_time_over_tau'],'passing_times':m['passing_times'],'total_times':len(times),'max_error_strictly_after_removal':float(errs[after].max()) if d is not None else None,'sign_disagreements':m['sign_disagreements'],'sign_eligible_entries':m['sign_eligible_entries'],'all_top5_order_identical_times':m['identical_top5_order_times'],'high_signal_times':len(signal),'high_signal_min_top5_membership':min([r['top5_membership'] for r in signal],default=None),'high_signal_identical_top5_order_times':sum(r['top5_order_identical'] for r in signal),'high_signal_gap_certifiable_times':sum(r['gap_above_0_004'] for r in signal),'sampled_peak':m['peak']}
  rows.append({'target':target,'rank':x['rank'],'step':compact(x['step']),'pulses':[compact(m,m['duration_over_tau']) for m in x['pulses']],'full_GPU_seconds':gpu['seconds'],'analysis_seconds':x['seconds'],'GPU_saved_CPU_parity':gpu['saved_CPU_parity_max_abs'],'independent_reduced_exponential_max_abs':max(z['independent_exponential_max_abs'] for z in x['numerical_verification'])})
 r={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'COMPLETE','scope':p['scope'],'observation_times':len(times),'full_kernel_times':len(p['all_exact_kernel_times_over_tau']),'rows':rows,'tiny_field':json.loads((ROOT/'tiny-field/receipt.json').read_text()),'rating_implication':'Measured characterization reveals the current temporal limitations; it does not establish improved biology or quantum advantage.'}
 (ROOT/'summary.json').write_text(json.dumps(r,indent=2)+'\n')
 print(json.dumps(r,indent=2))
if __name__=='__main__':main()
