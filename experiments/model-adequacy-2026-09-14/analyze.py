"""All fixed-time amplitudes and two-other-site rankings; no parameter selection."""
from pathlib import Path
import json,csv,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parent

def calculate():
    p=json.loads((ROOT/'protocol.json').read_text());times=np.array(p['all_times_sorted']);base=ROOT/'remote/all-five-grids-ten-times/data'
    z=[dict(np.load(base/f'case-{i}/response.npz')) for i in range(5)]
    receipts=[json.loads((base/f'case-{i}/receipt.json').read_text()) for i in range(5)]
    controls=dict(np.load(ROOT/'results/controls.npz'));cr=json.loads((ROOT/'results/controls-verification.json').read_text())
    for i,a in enumerate(z):
        assert np.array_equal(a['times_over_tau'],times)
        assert np.max(abs(a['sd']-controls['sd']))<=1e-10
        assert receipts[i]['faces']==p['cases'][i]['faces'] and receipts[i]['spacing']==p['cases'][i]['spacing']
    pairs=[('domain',0,1),('domain',1,2),('spacing',2,3),('spacing',3,4)]
    refine=[];empirical=np.zeros(10)
    for kind,i,j in pairs:
        gerr=float(np.max(abs(z[i]['G0']-z[j]['G0'])))
        for k,t in enumerate(times):
            ke=float(np.max(abs(z[i]['K_chebyshev'][k]-z[j]['K_chebyshev'][k])))
            ce=float(np.max(abs(z[i]['C_chebyshev'][k]-z[j]['C_chebyshev'][k])))
            refine.append(dict(kind=kind,first_grid=i,second_grid=j,time_over_tau=float(t),G0_error=gerr,K_error=ke,C_error=ce,pass_0_001=max(gerr,ke,ce)<=.001));empirical[k]+=ce
    propagation=[]
    for k,t in enumerate(times):
        dif=max(float(np.max(abs(a['K_chebyshev'][k]-a['K_uniformization'][k]))) for a in z)
        tail=max(sum(x['tail_bound_times_norms'] for x in r['propagation'] if abs(x['time_over_tau']-t)<1e-12) for r in receipts)
        empirical[k]+=dif+tail
        propagation.append(dict(time_over_tau=float(t),maximum_independent_difference=dif,maximum_sum_of_two_tail_allowances=tail,pass_1e_8=dif<=1e-8))
    static=[]
    for i,a in enumerate(z):
        old=json.loads((ROOT/f'inputs/historical/static-cell-{i}.json').read_text())
        static.append(float(np.max(abs(a['G0']-np.array(old['G0'])))))
    old=np.load(ROOT/'inputs/historical/controls.npz');ix=[list(times).index(t) for t in p['original_times']]
    historical=float(np.max(abs(z[-1]['C_chebyshev'][ix]-old['C_original'])))
    local=dict(np.load(ROOT/'local/pilot/response.npz'));pilot_difference=float(np.max(abs(local['C_chebyshev'][0]-z[0]['C_chebyshev'][0])))
    reference_pass=all(r['pass_0_001'] for r in refine) and all(r['pass_1e_8'] for r in propagation) and max(static)<=1e-8 and historical<=1e-8 and pilot_difference<=1e-8 and all(r['stationarity']<=1e-9 and r['independent_status']=='PASS' and all(v['tail_bound_times_norms']<=1e-9 for v in r['propagation']) for r in receipts)
    entries=[];orders=[];summaries={};ref=z[-1]['C_chebyshev']
    for name in ['Gaussian','harmonic']:
        C=controls['C_'+name];delta=abs(C-ref)
        numerical=2*cr['checks'][('gaussian' if name=='Gaussian' else 'harmonic')+'_Wick_maxabs']
        amplitude_counts=[]
        for k,t in enumerate(times):
            for sender in range(3):
                for receiver in range(3):
                    observed=float(delta[k,sender,receiver]);conditional=observed+float(empirical[k])+numerical
                    entries.append(dict(model=name,time_over_tau=float(t),sender=sender,receiver=receiver,reference_C=float(ref[k,sender,receiver]),model_C=float(C[k,sender,receiver]),signed_difference=float(C[k,sender,receiver]-ref[k,sender,receiver]),absolute_difference=observed,empirical_reference_envelope=float(empirical[k]),control_numerical_allowance=numerical,conditional_total_error=conditional,passes_0_002_observed=observed<=.002,passes_0_002_with_empirical_envelope=conditional<=.002))
            amplitude_counts.append(int(np.sum(delta[k]+empirical[k]+numerical<=.002)))
            for receiver in range(3):
                other=[i for i in range(3) if i!=receiver];u,v=other
                a=ref[k,other,receiver];b=C[k,other,receiver]
                ad=float(abs(a[0])-abs(a[1]));bd=float(abs(b[0])-abs(b[1]));combined=float(empirical[k])+numerical
                raw='tie' if ad==0 or bd==0 else 'agree' if ad*bd>0 else 'reverse'
                status='unresolved' if min(abs(ad),abs(bd))<=2*combined else 'agree' if ad*bd>0 else 'reverse'
                orders.append(dict(model=name,time_over_tau=float(t),receiver=receiver,other_site_1=u,other_site_2=v,reference_C_1=float(a[0]),reference_C_2=float(a[1]),model_C_1=float(b[0]),model_C_2=float(b[1]),reference_magnitude_gap=abs(ad),model_magnitude_gap=abs(bd),reference_preferred_site=None if ad==0 else u if ad>0 else v,model_preferred_site=None if bd==0 else u if bd>0 else v,raw_order_result=raw,conditional_result=status,combined_per_entry_empirical_allowance=combined))
        own=[x for x in orders if x['model']==name]
        summaries[name]={'maximum_C_discrepancy':float(delta.max()),'maximum_G0_discrepancy':float(np.max(abs(controls['G0_'+name]-z[-1]['G0']))),'maximum_K_discrepancy':float(np.max(abs(controls['K_'+name]-z[-1]['K_chebyshev']))),'maximum_C_discrepancy_location':{'time_over_tau':float(times[np.unravel_index(delta.argmax(),delta.shape)[0]]),'sender':int(np.unravel_index(delta.argmax(),delta.shape)[1]),'receiver':int(np.unravel_index(delta.argmax(),delta.shape)[2])},'total_entries':90,'observed_entries_within_0_002':int((delta<=.002).sum()),'entries_within_0_002_with_empirical_envelope':sum(amplitude_counts),'conditional_per_time_pass_counts':amplitude_counts,'raw_order_counts':{s:sum(x['raw_order_result']==s for x in own) for s in ['agree','reverse','tie']},'conditional_order_counts':{s:sum(x['conditional_result']==s for x in own) for s in ['agree','reverse','unresolved']},'original_time_raw_order_counts':{s:sum(x['raw_order_result']==s and x['time_over_tau'] in p['original_times'] for x in own) for s in ['agree','reverse','tie']},'amplitude_adequacy':'FAIL' if delta.max()>.002 else 'NOT_CERTIFIED','whole_space_fidelity_established':False}
    summary={'status':'COMPLETED_FIXED_MODEL_ADEQUACY_COMPARISON','protocol_sha256':hashlib.sha256((ROOT/'protocol.json').read_bytes()).hexdigest(),'all_times_over_tau':times.tolist(),'reference_numerical_screens':'PASS' if reference_pass else 'FAIL','maximum_domain_spacing_difference':max(max(r['G0_error'],r['K_error'],r['C_error']) for r in refine),'maximum_independent_propagation_difference':max(r['maximum_independent_difference'] for r in propagation),'static_independent_maximum':max(static),'historical_original_response_maximum_difference':historical,'local_CPU_A100_pilot_response_difference':pilot_difference,'empirical_reference_envelope_by_time':empirical.tolist(),'reference_error_scope':'Sum of observed domain/spacing differences plus independent propagation and tail allowances; empirical sensitivity envelope, not a certified whole-space error bound or proof of convergence. Ranking labels are conditional on this envelope.','ranking_rule':'Each receiver compares only two other sites; both reference and control magnitude gaps must exceed twice the combined empirical per-entry allowance; exact ties and overlapping intervals remain unresolved.','models':summaries,'new_fit':False,'new_basis':False,'biological_validation':False,'whole_space_dynamic_certificate':False,'scope':'Original nonlinear three-site finite reflecting model; no protein or experimental validation; fixed controls remain inadequate in normalized response amplitude.'}
    return summary,entries,orders,refine,propagation

def csv_write(path,rows):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

if __name__=='__main__':
    s,e,o,r,p=calculate();(ROOT/'results/summary.json').write_text(json.dumps(s,indent=2)+'\n');csv_write(ROOT/'results/all-response-entries.csv',e);csv_write(ROOT/'results/all-two-other-site-orderings.csv',o);csv_write(ROOT/'results/reference-refinement.csv',r);csv_write(ROOT/'results/independent-propagation.csv',p);print(json.dumps(s,indent=2))
