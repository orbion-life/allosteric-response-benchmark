"""Export exact completed receipts into portable CSV inputs; no new calculation."""
import json,csv,hashlib
from pathlib import Path
F=Path(__file__).resolve().parent;ROOT=F.parent;PHYSICAL=ROOT.parent/'physical'
sources={'harmonic-receipt.json':PHYSICAL/'harmonic-primary/calibration/receipt.json','nonlinear-h0.5.json':PHYSICAL/'nonlinear-primary/nonlinear-h0.5/receipt.json','nonlinear-h0.25.json':PHYSICAL/'nonlinear-primary/nonlinear-h0.25/receipt.json','nonlinear-watchdog.json':PHYSICAL/'nonlinear-primary/watchdog.json','physical-protocol.json':PHYSICAL/'protocol.json','stability-inspection.json':ROOT/'inspection.json','stability-summary.json':ROOT/'results/summary.json'}
records={}
for name,p in sources.items():
 data=p.read_bytes();(F/'sources'/name).write_bytes(data);records[name]={'source_path_from_v12_preparation':str(p.relative_to(ROOT.parent)),'sha256':hashlib.sha256(data).hexdigest()}
(F/'source-receipts.json').write_text(json.dumps(records,indent=2)+'\n')
def read(name):return json.loads((F/'sources'/name).read_text())
def write(name,rows):
 with (F/name).open('w',newline='') as stream:
  w=csv.DictWriter(stream,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
h=read('harmonic-receipt.json');write('harmonic-all-cases.csv',[{'case':c['name'],'geometry':'development','kappa':1,'A':c['faces'][0],'h':c['spacing'],**c['errors_against_analytic'],'max_component_error':max(c['errors_against_analytic'].values()),'screen_pass':c['all_components_within_0_001'],'equivalent_full_states':c['equivalent_full_states'],'full_grid_allocated':c['full_grid_allocated']} for c in h['cases']]);write('harmonic-all-successive-comparisons.csv',[{'first':c['first'],'second':c['second'],'kind':c['kind'],**c['errors'],'screen_pass':c['all_components_within_0_001']} for c in h['comparisons']])
s=read('stability-inspection.json');rows=[]
for c in s['cases']:
 for k,v in c['final_differences'].items():rows.append({'case':c['case'],'comparison':k,'scope':'final selected rank','C_difference':v['C'],'bound_difference':v['bound'],'basis_gap':c['comparisons'][k]['final_full_V']['gap']['maximum_gap']})
 rows.append({'case':c['case'],'comparison':'fixed_V','scope':'max over all checkpoints and driver comparisons','C_difference':c['fixed_V_max_C_difference'],'bound_difference':c['fixed_V_max_bound_difference'],'basis_gap':0})
write('stability-all-cases.csv',rows);write('stability-all-fixed-V-comparisons.csv',[{'case':c['case'],**x} for c in read('stability-summary.json')['cases'] for x in c['fixed_V_driver_comparisons']]);w=read('nonlinear-watchdog.json');rows=[];queries=[]
for spacing in [.5,.25]:
 n=read('nonlinear-h'+str(spacing)+'.json');outer=next(x['outer_worker_wall_seconds'] for x in w['jobs'] if x['name']=='nonlinear-h'+str(spacing));prep=sum(n['assembly_and_preprocessing'].values());taylor=sum(x['seconds'] for x in n['propagation']);cheby=sum(x['seconds'] for x in n['independent_chebyshev']);other=outer-prep-taylor-cheby;assert other>0
 rows.append({'h':spacing,'geometry':'development','kappa':1,'faces':'4,4,17','states':n['states'],'preparation_and_array_serialization_seconds':prep,'Taylor_three_queries_seconds':taylor,'Chebyshev_three_queries_seconds':cheby,'other_validation_startup_shutdown_seconds':other,'full_fresh_worker_seconds':outer,'worker_peak_rss_bytes':n['worker_peak_rss_bytes'],'independent_K_C_max_difference':n['independent_max_normalized_delayed_response_difference'],'boundary_mass':n['boundary_mass']})
 for method,key in [('Taylor','propagation'),('Chebyshev','independent_chebyshev')]:
  for i,q in enumerate(n[key]):queries.append({'h':spacing,'method':method,'time_over_tau':[.1,1,10][i],'seconds':q['seconds'],'vector_equivalent_products':q.get('vector_equivalent_products',q.get('counts',{}).get('vector_equivalents'))})
write('nonlinear-all-cases.csv',rows);write('nonlinear-all-queries.csv',queries)
print('Exact all-case CSV export complete.')
