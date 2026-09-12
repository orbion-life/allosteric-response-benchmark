#!/usr/bin/env python3
"""Check scientific invariants and compare a fresh calculation, not only file existence."""
from pathlib import Path
import json,sys
import numpy as np
from run import ROOT,sha,dump
def main():
    freeze=json.loads((ROOT/'prediction-freeze.json').read_text())
    assert all(sha(ROOT/p)==h for p,h in freeze['files'].items())
    repeat=ROOT/'checks/fresh-replay'; comparison=[]
    for path in sorted(ROOT.glob('model/*.npz'))+sorted(ROOT.glob('results/*.npz')):
        other=repeat/path.relative_to(ROOT);a=np.load(path);b=np.load(other);assert set(a.files)==set(b.files)
        differences=[]
        for k in a.files:
            assert a[k].shape==b[k].shape and a[k].dtype==b[k].dtype
            equal=np.array_equal(a[k],b[k]);err=float(np.max(np.abs(a[k].astype(float)-b[k].astype(float)))) if a[k].size else 0
            assert err<=1e-10,(str(path),k,err)
            differences.append({'array':k,'bitwise_equal_values':equal,'max_absolute_difference':err})
        comparison.append({'path':str(path.relative_to(ROOT)),'arrays':differences})
    checks=[]
    for p in sorted(ROOT.glob('results/*.npz')):
        z=np.load(p)
        if 'C' not in z.files:continue
        C=z['C'];assert C.shape==(252,252) and np.isfinite(C).all()
        asym=float(np.max(np.abs(C-C.T)));assert asym<1e-10
        mineig=float(np.linalg.eigvalsh(-z['R']).min());assert mineig>=-1e-10
        checks.append({'path':str(p.relative_to(ROOT)),'matrix_shape':list(C.shape),'max_asymmetry':asym,'smallest_eigenvalue_of_minus_response':mineig})
    lab=json.loads((ROOT/'evaluation/reference-labels.json').read_text());e=json.loads((ROOT/'evaluation/evaluation.json').read_text());unknown={r['canonical'] for r in lab['labels'] if r['state']=='unknown'}
    configs=np.load(ROOT/'evaluation/null-configurations.npz')['canonical_sets'];assert len(configs)==10000 and len(set(map(tuple,configs)))==10000
    assert not set(configs.ravel())&unknown
    assert not set(e['eligible_known_contacts'])&unknown
    for g in e['matched_null']['groups']:
        assert np.all(np.isin(configs,g['candidates']).sum(axis=1)==g['m'])
    model=np.load(ROOT/'model/static-model.npz');ids=model['canonical'];ii=np.flatnonzero(model['candidate'])
    for name,file in [('biquadratic_d2_n33','biquadratic-d2-n33-e4'),('biquadratic_d2_n65','biquadratic-d2-n65-e4')]:
        z=np.load(ROOT/f'results/{file}.npz');score=np.sqrt(np.mean(z['C'][:,model['receiver']]**2,axis=1));order=ii[np.lexsort((ids[ii],-score[ii]))]
        assert np.array_equal(score,z['score']) and np.array_equal(order,z['order'])
        assert ids[order[:5]].tolist()==[r['canonical'] for r in e['methods'][name]['top5']]
    dump(ROOT/'checks/verification.json',{'status':'passed','all_numeric_values_bitwise_equal_on_fresh_replay':all(a['bitwise_equal_values'] for f in comparison for a in f['arrays']),'archive_count':len(comparison),'comparison':comparison,'matrix_invariants':checks,'unknown_labels_excluded_from_null':True,'unique_matched_configurations':len(configs),'stratum_counts_preserved':True,'prediction_protocol_sha256':freeze['protocol_sha256'],'scope':'Repeatability in the recorded environment; this does not establish physical accuracy or biological validity.'})
    print(json.dumps({'status':'passed','archives':len(comparison),'bitwise_equal':all(a['bitwise_equal_values'] for f in comparison for a in f['arrays'])}))
if __name__=='__main__':main()
