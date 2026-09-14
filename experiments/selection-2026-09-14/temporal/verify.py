"""Offline protocol or complete saved-array audit; does not rerun GPU campaign."""
from pathlib import Path
import json,hashlib,argparse,datetime,sys
import numpy as np
ROOT=Path(__file__).resolve().parent

def sha(p):
 with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def verify_protocol():
 p=json.loads((ROOT/'protocol.json').read_text());seal=json.loads((ROOT/'source-seal.json').read_text())
 checks={'sealed_sources':all(sha(ROOT/x['file'])==x['sha256'] for x in seal['files']),'targets':p['targets']==['kras','abl','myc','myh7'],'log_grid_contained':all(float(x) in p['observation_times_over_tau'] for x in np.geomspace(.001,30,61)),'pulse_times_complete':all(float(t-d) in p['all_exact_kernel_times_over_tau'] for d in p['pulse_durations_over_tau'] for t in p['observation_times_over_tau'] if t>=d),'strictly_increasing':bool(np.all(np.diff(p['observation_times_over_tau'])>0) and np.all(np.diff(p['all_exact_kernel_times_over_tau'])>0)),'zero_included':p['observation_times_over_tau'][0]==0,'fixed_error':p['max_step_error']==.002 and p['max_pulse_error']==.002,'no_refit':p['no_refit_or_basis_changes']}
 dep=json.loads((ROOT/'tiny-dependency-receipt.json').read_text());checks['tiny_model_dependency']=sha(ROOT/dep['file'])==dep['sha256']
 return checks

def verify_results(target):
 p=json.loads((ROOT/'protocol.json').read_text());src=ROOT/'inputs'/target;out=ROOT/'results'/target
 manifest=next(x for x in json.loads((ROOT/'input-manifest.json').read_text()) if x['target']==target)
 checks={'input_hashes':all(sha(ROOT/f['file'])==f['sha256'] for f in manifest['files'])}
 receipt=json.loads((out/'receipt.json').read_text());checks['full_kernel_hash']=sha(out/'full-kernels.npy')==receipt['full_kernel_sha256']
 checks['canary_pass']=receipt['saved_CPU_parity_max_abs']<=1e-8
 R=np.load(out/'full-kernels.npy',mmap_mode='r');nodes=np.array(p['all_exact_kernel_times_over_tau']);idx={float(t):i for i,t in enumerate(nodes)}
 step=np.load(out/'step.npz');summary=json.loads((out/'analysis.json').read_text());times=step['times_over_tau'];reference=step['full'];reduced=step['reduced'];worst=0.
 for j,t in enumerate(times):worst=max(worst,float(np.max(abs(reference[j]-(R[idx[float(t)]]-R[idx[0.]])))))
 checks['step_full_identity']=worst<=1e-12
 error=np.max(abs(reference-reduced),axis=(1,2));checks['step_metrics']=bool(np.allclose(error,summary['step']['error_by_time'],rtol=0,atol=1e-12))
 checks['independent_exponential']=max(r['independent_exponential_max_abs'] for r in summary['numerical_verification'])<=1e-8
 # This audit recomputes the pulse reference and declared maxima from saved raw arrays.
 # Reduced propagation is checked by the independently calculated exponential receipts.
 per_pulse=[]
 for d,m in zip(p['pulse_durations_over_tau'],summary['pulses']):
  z=np.load(out/f'pulse-{d:g}.npz');full=z['full'];red=z['reduced'];check=0.
  for j,t in enumerate(times):
   ref=R[idx[float(t)]]-R[idx[0.]]
   if t>=d:ref=R[idx[float(t)]]-R[idx[float(t-d)]]
   check=max(check,float(np.max(abs(full[j]-ref))))
  errors=np.max(abs(full-red),axis=(1,2))
  per_pulse.append({'duration_over_tau':d,'full_identity_max_abs':check,'metrics_agree':bool(np.allclose(errors,m['error_by_time'],rtol=0,atol=1e-12)),'scientific_all_times_pass':bool(np.all(errors<=.002))})
 checks['pulse_identities_and_metrics']=all(x['full_identity_max_abs']<=1e-12 and x['metrics_agree'] for x in per_pulse)
 return {'target':target,'checks':checks,'all_implementation_checks_pass':all(checks.values()),'step_scientific_all_times_pass':bool(np.all(error<=.002)),'pulse_results':per_pulse}

def verify_tiny():
 out=ROOT/'tiny-field';z=np.load(out/'raw.npz');receipt=json.loads((out/'receipt.json').read_text());eta=z['eta'];raw=z['raw_on_off_means'];calc=(raw[:,:,0]-raw[:,:,1])/(2*eta[None,:,None,None]);ref=z['pulse_reference'].transpose(1,0,2)[:,None];errors=np.max(abs(calc-ref),axis=(0,2,3));ind=(z['independent_forward_means'][:,:,0]-z['independent_forward_means'][:,:,1])/(2*eta[None,:,None,None])
 checks={'central_response':bool(np.max(abs(calc-z['central_response']))<=1e-12),'receipt_metrics':bool(np.allclose(errors,receipt['max_identity_error_by_amplitude'],rtol=0,atol=1e-12)),'independent_forward':bool(np.max(abs(calc-ind))<=1e-8),'pulse_gate':bool(errors[-1]<=.0005)}
 return checks
if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--mode',choices=['protocol','saved'],default='saved');a.add_argument('--output',required=True);args=a.parse_args()
 result={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'scope':'Offline implementation and provenance audit; a passing audit does not change scientific failures','protocol':verify_protocol(),'tiny_field':verify_tiny()}
 if args.mode=='saved':result['targets']=[verify_results(t) for t in ['kras','abl','myc','myh7']]
 result['passed']=all(result['protocol'].values()) and all(result['tiny_field'].values()) and all(x['all_implementation_checks_pass'] for x in result.get('targets',[]))
 Path(args.output).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2));sys.exit(0 if result['passed'] else 1)
