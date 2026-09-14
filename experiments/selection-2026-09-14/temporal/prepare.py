from pathlib import Path
import numpy as np,json,hashlib,shutil,datetime
ROOT=Path(__file__).resolve().parent
OLD=ROOT.parents[1]/'pulsar-phase1-execution-2026-09-13'
PRIOR=ROOT.parents[1]/'pulsar-feasibility-validation-2026-09-14'
NEW=ROOT.parents[1]/'pulsar-field-switch-execution-2026-09-13'
(ROOT/'vendor').mkdir(exist_ok=True);(ROOT/'inputs').mkdir(exist_ok=True)
for name in ['gaussian_protein.py','check_independent.py']:
 shutil.copy2(OLD/'gaussian'/name,ROOT/'vendor'/name)
for name in ['model.py','chebyshev.py']:
 shutil.copy2(NEW/'vendor'/name,ROOT/'vendor'/name)
shutil.copy2(NEW/'dynamics/core.py',ROOT/'vendor/finite_core.py')
records=[]
for target in ['kras','abl','myc','myh7']:
 source=OLD/'galerkin/inputs'/target if target in ['kras','abl'] else PRIOR/'representation/inputs'/target
 out=ROOT/'inputs'/target;out.mkdir(exist_ok=True)
 for name in ['static-model.npz','fitted-covariance.npz','all-mode-harmonic-scales.npz']:
  shutil.copy2(source/name,out/name)
 operator=PRIOR/'representation/global-gram-followup'/target/'result.npz'
 z=np.load(operator)
 np.savez_compressed(out/'fixed-operator.npz',**{k:z[k] for k in ['Hr','B','G0','reference_G0','reference_K','reference_C','C','times']})
 prior=np.load(OLD/'gaussian'/f'{target}-cpu/response-matrices.npz')
 np.savez_compressed(out/'ranking-metadata.npz',**{k:prior[k] for k in ['canonical','receiver','candidate']})
 records.append({'target':target,'original_operator_sha256':hashlib.sha256(operator.read_bytes()).hexdigest(),'original_operator':str(operator),'files':[{'file':str(p.relative_to(ROOT)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size} for p in sorted(out.glob('*.npz'))]})
base=sorted(set(np.geomspace(.001,30,61).tolist()+[0,.015,.03,.04,.08,.1,.3,.4,.8,1.,3.,4.,10.,30.]))
durations=[.03,.3,3.]
unions=sorted(set(base+[float(t-d) for d in durations for t in base if t>=d]))
protocol={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'FROZEN_BEFORE_NEW_TEMPORAL_SCORES','scope':'Dense fixed-operator Gaussian numerical characterization, not fresh temporal holdout, nonlinear fidelity or biological validation','targets':['kras','abl','myc','myh7'],'base_log_grid':{'count':61,'min':.001,'max':30.,'spacing':'numpy.geomspace'},'observation_times_over_tau':base,'all_exact_kernel_times_over_tau':unions,'pulse_durations_over_tau':durations,'step':'K(t)-K(0)','pulse':'C_step(t)-C_step(t-d), with second term exactly zero when t<d','time_scale':'inverse smallest positive native harmonic eigenvalue; no calibration to seconds','max_step_error':.002,'max_pulse_error':.002,'low_signal_floor':.002,'sign_eligibility':'absolute full reference entry >0.002; report disagreements and full error even if ineligible','top5':'RMS over fixed receivers; canonical tie breaking; report membership and order; reference fifth/sixth gap <=0.004 considered not certifiable at error allocation','sampled_peak':'absolute signed entry peak only for eligible source-receiver pairs with reference peak>0.002; interior unique peak required, ties within 1e-12 and boundary peaks unresolved; report discrete-grid errors only','resource_policy':{'local_threads':3,'memory_GiB':12,'remote_gpu':'one A100-40GB','remote_timeout_seconds':2700,'automatic_retries':0},'GPU_canary':'for each target compare K(0), K(0.1tau), K(tau), K(10tau) to saved CPU exact reference, <=1e-8 before new time outputs; independent Wick fixtures <=1e-8','no_refit_or_basis_changes':True,'failure_policy':'retain all fixed-grid scores and failure times; no responsive tuning','tiny_pulse':{'faces':[2.,2.,2.],'spacing':1.,'source_residues':[0,1,2],'eta':[.01,.005,.0025],'duration_over_tau':.3,'times_over_tau':[.03,.1,.3,.4,1.,3.],'acceptance_max_abs':.0005,'scope':'64-state original quartic finite box; independent dense eigensystem and forward-generator piecewise propagation, no physical convergence'}}
(ROOT/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
(ROOT/'input-manifest.json').write_text(json.dumps(records,indent=2)+'\n')
print({'observation_times':len(base),'exact_kernel_times':len(unions),'input_bytes':sum(p.stat().st_size for p in (ROOT/'inputs').rglob('*.npz'))})
