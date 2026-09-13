# Run this saved-array helper from the complete extracted raw stability/ tree.
"""Audit saved traces only; no basis, operator or reference construction."""
from pathlib import Path
import hashlib,json,time
import numpy as np
import diagnostic as d
ROOT=Path(__file__).resolve().parent
start=time.monotonic()
def read(p):return json.loads(p.read_text())
def maxdiff(a,b):return float(np.max(np.abs(a-b)))
cases=[];coverage=[]
for index in range(2):
 name=d.case_info(index)['name'];cs={'case':name,'input_comparison':read(ROOT/'results/comparison'/(name+'-input-comparison.json')),'comparisons':{},'traces':{}}
 for suffix in ['native','shared','system-input-effect']:
  comp=read(ROOT/'results/comparison'/(name+'-'+suffix+'.json'));full=[r for r in comp['rows'] if r['stage']=='full_V'];bad=[{'step':r['step'],'stage':r['stage'],'effective_ranks':r['gap']['effective_ranks']} for r in comp['rows'] if r.get('gap',{}).get('rank_deficient_input')]
  cs['comparisons'][suffix]={'threshold_crossings_all_stages':comp['thresholds'],'rank_deficient_stages':bad,'final_full_V':full[-1],'first_full_V_above':{str(t):next(({'step':r['step'],'rank':r['rank_a'],'maximum_gap':r['gap']['maximum_gap']} for r in full if r['gap']['maximum_gap']>t),None) for t in [1e-12,1e-10,1e-8,1e-6]}}
 for label,kind in [('pinned','native'),('system','native'),('system','shared')]:
  base=ROOT/('results' if label=='pinned' else 'trace-repair')/label/name/kind
  d.assert_trace(base);c=read(base/'trace-completeness.json');coverage.append({'path':str(base.relative_to(ROOT)),**c})
  steps=sorted(base.glob('step-*-svd.npz'));mins=[];metadata=[]
  for p in steps:
   z=np.load(p);mins.append(float(z['s'].min()/z['threshold']));metadata.append(read(base/(p.name.replace('-svd.npz','.json'))))
  r={'path':str(base.relative_to(ROOT)),'steps':len(steps),'minimum_singular_value_to_admission_threshold_ratio':min(mins),'prefix_cut_count':sum(x['prefix_cut'] for x in metadata),'eligible_widths':sorted(set(x['eligible_width'] for x in metadata)),'admitted_widths':sorted(set(x['admitted_width'] for x in metadata)),'evaluation':read(base/'evaluation.json')}
  if label=='system':
   parity=read(base/'prior-trace-bitwise-parity.json');assert all(parity.values());r['prior_files_bitwise_equal']=len(parity)
   a=np.load(base/'reference.npz');b=np.load(ROOT/'results'/label/name/kind/'reference.npz');assert set(a.files)==set(b.files);r['reference_bitwise_equal']=all(np.array_equal(a[k],b[k]) for k in a.files);r['reference_max_differences']={k:maxdiff(a[k],b[k]) for k in a.files}
  cs['traces'][label+'_'+kind]=r
 rank=cs['traces']['pinned_native']['evaluation']['rank'];a=np.load(ROOT/'results/pinned'/name/'native'/('checkpoint-%03d'%rank)/'outputs.npz')
 cs['final_differences']={}
 for kind in ['native','shared']:
  b=np.load(ROOT/'trace-repair/system'/name/kind/('checkpoint-%03d'%rank)/'outputs.npz');cs['final_differences'][kind]={k:maxdiff(a[k],b[k]) for k in ['C','bound']}
 cases.append(cs)
summary=read(ROOT/'results/summary.json')
for case,s in zip(cases,summary['cases']):
 case['fixed_V_max_C_difference']=max(x['C_max_difference'] for x in s['fixed_V_driver_comparisons']);case['fixed_V_max_bound_difference']=max(x['bound_max_difference'] for x in s['fixed_V_driver_comparisons']);case['same_W_max_admitted_gap']=max(x['admitted_gap']['maximum_gap'] for x in s['isolated_svd_comparisons']);case['same_W_keep_masks_equal']=all(x['keep_equal'] for x in s['isolated_svd_comparisons']);case['same_W_step']=s['isolated_svd_comparisons'][0]['step']
files=[p for p in ROOT.rglob('*') if p.is_file()];sizes={'total_bytes_before_inspection_receipt':sum(p.stat().st_size for p in files),'file_count':len(files),'fixed_V_driver_bytes':sum(p.stat().st_size for p in files if 'fixed-V-drivers' in p.parts),'results_bytes':sum(p.stat().st_size for p in files if p.is_relative_to(ROOT/'results')),'trace_repair_bytes':sum(p.stat().st_size for p in files if p.is_relative_to(ROOT/'trace-repair'))}
receipt={'status':'ALL SIX COMPLETE TRACES VERIFIED; ORIGINAL INCOMPLETE TRACES RETAINED','cases':cases,'trace_coverage':coverage,'sizes':sizes,'seconds':time.monotonic()-start,'method':'Only saved arrays and comparison receipts are inspected. No operator, basis or reference is constructed. No scientific settings change.'}
(ROOT/'inspection.json').write_text(json.dumps(receipt,indent=2,allow_nan=False)+'\n')
print(json.dumps({'status':receipt['status'],'seconds':receipt['seconds'],'sizes':sizes}))
