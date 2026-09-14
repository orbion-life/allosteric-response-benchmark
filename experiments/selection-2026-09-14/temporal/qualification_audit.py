"""Explicit post-execution interpretation of candidate-specific signal floors.

Preserves frozen global-score flags and all whole-matrix pass/fail decisions.
"""
from pathlib import Path
import numpy as np,json,datetime
ROOT=Path(__file__).resolve().parent

def main():
 rows=[]
 for label,base in [('original-grid',ROOT),('event-aligned',ROOT/'event-aligned')]:
  for target in ['kras','abl','myc','myh7']:
   p=base/'results'/target/'analysis.json'
   if not p.exists():continue
   x=json.loads(p.read_text());meta=np.load(ROOT/'inputs'/target/'ranking-metadata.npz');cand=meta['candidate'];receiver=meta['receiver']
   endpoints=([('step',x['step'])] if 'step' in x else [])+[(f"pulse-{m['duration_over_tau']:g}",m) for m in x['pulses']]
   for name,m in endpoints:
    z=np.load(base/'results'/target/f'{name}.npz');refs=z['reference_score'];orders=z['reference_order'];eligible=refs[:,cand].max(axis=1)>.002;selected=[v for v,k in zip(m['ranking'],eligible) if k]
    full=z['full'][:,cand][:,:,receiver];red=z['reduced'][:,cand][:,:,receiver];errors=np.max(abs(full-red),axis=(1,2));mask=abs(full)>.002
    peak=float(np.max(abs(full)));norm=float(np.linalg.norm(full));errornorm=float(np.linalg.norm(red-full));relative_status='DESCRIPTIVE' if peak>.002 else 'LOW_SIGNAL_INDETERMINATE'
    rows.append({'grid':label,'target':target,'endpoint':name,'qualification_scope':'Additional candidate-only descriptive audit; frozen global-score flags and whole-matrix decisions retained','original_global_high_signal_times':sum(not v['low_signal'] for v in m['ranking']),'candidate_high_signal_times':int(eligible.sum()),'candidate_high_signal_min_top5_membership':min([v['top5_membership'] for v in selected],default=None),'candidate_high_signal_order_identical':sum(v['top5_order_identical'] for v in selected),'candidate_gap_above_0_004_times':sum(v['gap_above_0_004'] for v in selected),'max_candidate_to_receiver_entry_error':float(errors.max()),'candidate_receiver_reference_peak_magnitude':peak,'candidate_receiver_reference_Frobenius_norm':norm,'candidate_receiver_error_to_peak_magnitude':float(errors.max())/peak if peak>.002 else None,'candidate_receiver_relative_Frobenius_error':errornorm/norm if peak>.002 else None,'relative_error_status':relative_status,'candidate_to_receiver_error_by_time':errors.tolist(),'candidate_receiver_sign_eligible_entries':int(mask.sum()),'candidate_receiver_sign_disagreements':int(np.count_nonzero(mask & (np.sign(full)!=np.sign(red))))})
 result={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'interpretive_correction':'The frozen low_signal flag uses the maximum score across all residues. This audit separately quantifies eligible candidate scores; it does not revise flags, outputs, thresholds or original scientific decisions.','relative_Frobenius_scope':'Unweighted Frobenius norm over all listed sample times and fixed candidate-by-receiver entries, not a continuous-time integral. Error-to-peak uses the maximum full reference magnitude on that same submatrix/time grid. Near-zero global reference peak (<=0.002) is indeterminate.','rows':rows}
 (ROOT/'candidate-qualification-audit.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
