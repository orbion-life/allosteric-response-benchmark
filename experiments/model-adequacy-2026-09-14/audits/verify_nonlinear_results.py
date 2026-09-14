"""Independent saved-array audit; no propagation jobs or author-code imports.
Writes only to the required output directory. The refinement envelope is an empirical
sensitivity sum, not a certified error bound for unbounded dynamics.
"""
from pathlib import Path
import csv, hashlib, json, sys, time
import numpy as np
HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
OUT=HERE
INDEPENDENT=HERE/'independent-control-moments.npz'
DATA=ROOT/'remote/all-five-grids-ten-times/data'

def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read_npz(p):
    with np.load(p) as f:return {k:f[k].copy() for k in f.files}

def main():
    start=time.perf_counter();p=json.loads((ROOT/'protocol.json').read_text());frozen=json.loads((ROOT/'protocol-freeze.json').read_text())
    assert digest(ROOT/'protocol.json')==frozen['protocol_sha256']
    provider=json.loads((ROOT/'remote/all-five-grids-ten-times/provider-receipt.json').read_text())
    assert provider['status']=='COMPLETE'
    remote_mismatch=[]
    for r in provider['raw_files']:
        file=DATA/r['path']
        if not file.exists() or file.stat().st_size!=r['bytes'] or digest(file)!=r['sha256']:remote_mismatch.append(r['path'])
    assert not remote_mismatch,remote_mismatch
    z=[read_npz(DATA/f'case-{i}/response.npz') for i in range(5)]
    receipts=[json.loads((DATA/f'case-{i}/receipt.json').read_text()) for i in range(5)]
    control=read_npz(ROOT/'results/controls.npz');control_check=json.loads((ROOT/'results/controls-verification.json').read_text())
    times=np.array(p['all_times_sorted']);n=len(times)
    assert n==10 and all(t in times for t in p['original_times'])
    for i,case in enumerate(p['cases']):
        assert receipts[i]['faces']==case['faces'] and receipts[i]['spacing']==case['spacing']
        assert np.array_equal(z[i]['times_over_tau'],times)
        assert receipts[i]['source_sha256']==frozen['input_hashes']['source/dynamics/core.py']
        assert np.max(abs(z[i]['sd']-control['sd']))<1e-10
        assert receipts[i]['stationarity']<=p['numerical_gates']['stationarity']
        assert receipts[i]['independent_status']=='PASS'
    refinement=[];envelope=np.zeros(n)
    for kind,left,right in [('domain',0,1),('domain',1,2),('spacing',2,3),('spacing',3,4)]:
        static=float(abs(z[left]['G0']-z[right]['G0']).max())
        lag=np.max(abs(z[left]['K_chebyshev']-z[right]['K_chebyshev']),axis=(1,2))
        response=np.max(abs(z[left]['C_chebyshev']-z[right]['C_chebyshev']),axis=(1,2))
        assert max(static,lag.max(),response.max())<=p['numerical_gates']['successive_G0_K_C']
        envelope+=response
        refinement.append({'kind':kind,'left':left,'right':right,'static_max':static,'lag_max':float(lag.max()),'response_max':float(response.max())})
    propagation=np.max(np.stack([np.max(abs(a['K_chebyshev']-a['K_uniformization']),axis=(1,2)) for a in z]),axis=0)
    assert propagation.max()<=p['numerical_gates']['independent_propagation']
    tails=[]
    for t in times:
        sums=[]
        for r in receipts:
            terms=[a for a in r['propagation'] if abs(a['time_over_tau']-t)<1e-12]
            assert len(terms)==2
            assert {a['method'] for a in terms}=={'chebyshev','uniformization'}
            assert all(a['tail_bound_times_norms']<=1e-9 for a in terms)
            sums.append(sum(a['tail_bound_times_norms'] for a in terms))
        tails.append(max(sums))
    envelope+=propagation+tails
    ref=z[-1]['C_chebyshev'];offdiag=~np.eye(3,dtype=bool)
    model_rows=[];ordering_rows=[]
    with (ROOT/'results/all-two-other-site-orderings.csv').open() as file:author=list(csv.DictReader(file))
    assert len(author)==60
    for name in ['Gaussian','harmonic']:
        proposed=control['C_'+name];delta=abs(proposed-ref)
        numerical=2*control_check['checks'][('gaussian' if name=='Gaussian' else 'harmonic')+'_Wick_maxabs']
        counts=dict(agree=0,reverse=0,unresolved=0);raw=dict(agree=0,reverse=0,tie=0);unresolved=[];oldraw=dict(agree=0,reverse=0,tie=0)
        for k,t in enumerate(times):
            for receiver in range(3):
                others=np.array([i for i in range(3) if i!=receiver]);x=abs(ref[k,others,receiver]);y=abs(proposed[k,others,receiver])
                gx=float(x[0]-x[1]);gy=float(y[0]-y[1]);allow=float(envelope[k]+numerical)
                raw_result='tie' if gx==0 or gy==0 else ('agree' if gx*gy>0 else 'reverse')
                result='unresolved' if min(abs(gx),abs(gy))<=2*allow else raw_result
                raw[raw_result]+=1;counts[result]+=1
                if t in p['original_times']:oldraw[raw_result]+=1
                if result=='unresolved':unresolved.append({'time_over_tau':float(t),'receiver':receiver,'reference_gap':abs(gx),'model_gap':abs(gy),'threshold':2*allow})
                matched=[r for r in author if r['model']==name and float(r['time_over_tau'])==t and int(r['receiver'])==receiver]
                assert len(matched)==1 and matched[0]['conditional_result']==result and matched[0]['raw_order_result']==raw_result
                ordering_rows.append({'model':name,'time_over_tau':float(t),'receiver':receiver,'raw':raw_result,'conditional':result})
        model_rows.append({'model':name,'all_entry_count':delta.size,'offdiagonal_entry_count':delta[:,offdiag].size,
                           'minimum_C_discrepancy':float(delta.min()),'maximum_C_discrepancy':float(delta.max()),
                           'maximum_offdiagonal_C_discrepancy':float(delta[:,offdiag].max()),
                           'entries_within_0p002_observed':int(np.sum(delta<=.002)),
                           'entries_within_0p002_with_envelope':int(np.sum(delta+envelope[:,None,None]+numerical<=.002)),
                           'raw_order_counts':raw,'conditional_order_counts':counts,'original_times_raw':oldraw,'unresolved':unresolved})
    declared=json.loads((ROOT/'results/summary.json').read_text())
    assert np.max(abs(envelope-np.array(declared['empirical_reference_envelope_by_time'])))<1e-18
    for row in model_rows:
        original=declared['models'][row['model']]
        assert row['conditional_order_counts']==original['conditional_order_counts']
        assert row['raw_order_counts']==original['raw_order_counts']
        assert abs(row['maximum_C_discrepancy']-original['maximum_C_discrepancy'])<1e-15
    historical=read_npz(ROOT/'inputs/historical/controls.npz');oldidx=[list(times).index(t) for t in p['original_times']]
    old_error=float(abs(ref[oldidx]-historical['C_original']).max());assert old_error<1e-8
    independent=read_npz(INDEPENDENT)
    GH_error=max(float(abs(independent[name+'_'+kind]-control[kind+'_'+name]).max()) for name in ['Gaussian','harmonic'] for kind in ['G0','K','C'])
    receipt={'status':'PASS','scope':'Independent saved-array and download/provenance verification; finite reflecting triangle only. No new propagation, unbounded-dynamics certificate, or protein inference.',
             'time_count':n,'per_model_response_entry_count':90,'per_model_two_other_site_ordering_count':30,
             'same_fine_grid_operator_shared_by_two_propagators':True,'propagators':'Chebyshev/Bessel and Poisson uniformization; small CSR/expm tests supply independent implementation check',
             'refinement':refinement,'maximum_independent_difference':float(propagation.max()),
             'empirical_envelope':envelope.tolist(),'envelope_scope':'Sum of observed four successive response differences plus propagation discrepancy and both tail allowances; empirical only',
             'historical_original_time_response_error':old_error,'remote_manifest_mismatches':remote_mismatch,
             'model_results':model_rows,'independent_direct_Gauss_Hermite_control_max_error':GH_error,
             'source_hashes':{'audit_script':digest(__file__),'protocol':digest(ROOT/'protocol.json'),'provider_receipt':digest(ROOT/'remote/all-five-grids-ten-times/provider-receipt.json'),
                             'analyzed_summary':digest(ROOT/'results/summary.json'),'controls':digest(ROOT/'results/controls.npz'),
                             **{f'nonlinear_grid_{i}':digest(DATA/f'case-{i}/response.npz') for i in range(5)}},
             'warnings':['Ordering comparisons are conditional on an empirical allowance, not certified whole-space rankings.','Ninety entries and thirty orderings are repeated outputs from one three-site fixture, not independent biological samples.'],
             'seconds':time.perf_counter()-start,'python':sys.version,'numpy':np.__version__}
    (OUT/'nonlinear-results-check.json').write_text(json.dumps(receipt,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:receipt[k] for k in ['status','maximum_independent_difference','model_results','seconds']},indent=2))

if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True);parser.add_argument('--independent-controls',type=Path,default=INDEPENDENT);args=parser.parse_args()
    OUT=args.output;OUT.mkdir(parents=True,exist_ok=True);INDEPENDENT=args.independent_controls;main()
