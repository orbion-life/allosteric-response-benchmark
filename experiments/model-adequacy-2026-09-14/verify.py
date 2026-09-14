"""Read-only saved evidence verification; no GPU or provider call."""
from pathlib import Path
import json,hashlib,csv
import numpy as np
ROOT=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def verify():
    p=json.loads((ROOT/'protocol.json').read_text());freeze=json.loads((ROOT/'protocol-freeze.json').read_text())
    assert sha(ROOT/'protocol.json')==freeze['protocol_sha256']
    adaptations={x['path']:x for x in json.loads((ROOT/'PUBLIC-PACKAGING.json').read_text())['metadata_changes']}
    for name,h in freeze['input_hashes'].items():
        expected=h
        if name in adaptations:
            entry=adaptations[name];assert entry['original_sha256']==h
            expected=entry['public_sha256']
        assert sha(ROOT/name)==expected,name
    sourcehashes=json.loads((ROOT/'remote/all-five-grids-ten-times/provider-receipt.json').read_text())['source_sha256']
    for name,h in sourcehashes.items():
        local=ROOT/'protocol.json' if name=='study-protocol.json' else ROOT/'source'/name
        assert sha(local)==h,name
    fresh,control_receipt=__import__('controls').evaluate()
    saved=dict(np.load(ROOT/'results/controls.npz'))
    for name,x in fresh.items():assert np.allclose(saved[name],x,atol=1e-10,rtol=1e-12),name
    base=ROOT/'remote/all-five-grids-ten-times/data';grids=[dict(np.load(base/f'case-{i}/response.npz')) for i in range(5)]
    recs=[json.loads((base/f'case-{i}/receipt.json').read_text()) for i in range(5)]
    for a,r in zip(grids,recs):
        assert np.array_equal(a['times_over_tau'],p['all_times_sorted'])
        for method in ['chebyshev','uniformization']:
            assert np.allclose(a['C_'+method],a['K_'+method]-a['G0'],atol=1e-14,rtol=0)
        observed=np.max(abs(a['K_chebyshev']-a['K_uniformization']))
        assert abs(observed-r['independent_error'])<1e-14
        assert observed<=1e-8 and r['independent_status']=='PASS' and r['stationarity']<=1e-9
        assert all(v['tail_bound_times_norms']<=1e-9 for v in r['propagation'])
        assert len(r['propagation'])==20
    assert json.loads((ROOT/'local/small/receipt.json').read_text())['status']=='PASS'
    assert json.loads((base/'gpu-small-check/receipt.json').read_text())['status']=='PASS'
    report=json.loads((ROOT/'results/summary.json').read_text());ref=grids[-1]['C_chebyshev'];counts={}
    for model in ['Gaussian','harmonic']:
        c=saved['C_'+model];d=abs(c-ref);s=report['models'][model]
        assert abs(float(d.max())-s['maximum_C_discrepancy'])<1e-14
        assert s['total_entries']==90 and s['observed_entries_within_0_002']==int((d<=.002).sum())
        raw={'agree':0,'reverse':0,'tie':0}
        for k in range(10):
            for receiver in range(3):
                a,b=[i for i in range(3) if i!=receiver]
                x=abs(ref[k,a,receiver])-abs(ref[k,b,receiver]);y=abs(c[k,a,receiver])-abs(c[k,b,receiver])
                raw['tie' if x==0 or y==0 else 'agree' if np.sign(x)==np.sign(y) else 'reverse']+=1
        assert raw==s['raw_order_counts'];counts[model]=raw
    entries=list(csv.DictReader((ROOT/'results/all-response-entries.csv').open()));orders=list(csv.DictReader((ROOT/'results/all-two-other-site-orderings.csv').open()))
    assert len(entries)==180 and len(orders)==60
    assert len({(e['model'],e['time_over_tau'],e['sender'],e['receiver']) for e in entries})==180
    assert len({(e['model'],e['time_over_tau'],e['receiver']) for e in orders})==60
    from analyze import calculate
    replay=calculate()[0];assert replay==report
    manifest=ROOT/'PUBLIC-MANIFEST.json'
    assert manifest.exists(),'PUBLIC-MANIFEST.json is required'
    for row in json.loads(manifest.read_text())['files']:
        f=ROOT/row['path'];assert f.stat().st_size==row['bytes'] and sha(f)==row['sha256'],row['path']
    return {'status':'PASS','original_protocol_freeze_and_mapped_public_inputs':'PASS','metadata_hash_adaptations':len(adaptations),'original_solver_bytes':'PASS','fresh_control_Wick_check':control_receipt['checks'],'reference_grids':5,'propagations':100,'unique_amplitude_rows':len(entries),'unique_two_site_rows':len(orders),'independent_raw_order_counts':counts,'summary_replay_identical':True,'no_GPU_execution':True}
if __name__=='__main__':print(json.dumps(verify(),indent=2))
