#!/usr/bin/env python3
"""Check scientific invariants and compare a fresh calculation, not only file existence."""
from pathlib import Path
import json,sys
import numpy as np
from run import ROOT,sha,dump

ATOL = 1e-10
RAW_RTOL = 1e-10
# These quantities retain their unnormalised physical or polynomial scale.
# All other floating arrays, including C, score and equilibrium, use rtol=0.
RAW_SCALE_FIELDS = frozenset({
    'monomial_mean', 'monomial_centered', 'monomial_covariance',
    'monomial_time_covariance', 'covariance', 'time_covariance',
    'R', 'sd', 'harmonic_sd',
})
EXACT_FIELDS = frozenset({'order', 'canonical', 'receiver', 'candidate',
                          'edges', 'degree', 'powers', 'dimensions',
                          'H_indices', 'H_indptr'})


def compare_array(key, reference, actual):
    """Compare actual against the frozen reference using a field-specific policy.

    The relative tolerance applies only to explicitly named unnormalised fields.
    Shapes, dtypes and discrete values must match exactly; nonfinite values fail.
    """
    reference, actual = np.asarray(reference), np.asarray(actual)
    assert reference.shape == actual.shape, (key, 'shape', reference.shape, actual.shape)
    assert reference.dtype == actual.dtype, (key, 'dtype', str(reference.dtype), str(actual.dtype))
    if reference.dtype.kind in 'fc':
        assert np.isfinite(reference).all() and np.isfinite(actual).all(), (key, 'nonfinite values')
    exact = key in EXACT_FIELDS or reference.dtype.kind not in 'fc'
    rtol = RAW_RTOL if key in RAW_SCALE_FIELDS and not exact else 0.0
    equal = bool(np.array_equal(reference, actual))
    err = float(np.max(np.abs(actual.astype(float) - reference.astype(float)))) if reference.size and reference.dtype.kind in 'biuf' else 0.0
    if exact:
        passed, policy = equal, 'exact'
    else:
        # Put the frozen reference second: NumPy scales rtol by this argument.
        passed = bool(np.allclose(actual, reference, atol=ATOL, rtol=rtol, equal_nan=False))
        policy = 'raw_scale_aware' if rtol else 'strict_absolute'
    diagnostic = {'array': key, 'policy': policy, 'atol': 0.0 if exact else ATOL,
                  'rtol': rtol, 'bitwise_equal_values': equal,
                  'max_absolute_difference': err}
    assert passed, diagnostic
    return diagnostic


def main():
    freeze=json.loads((ROOT/'prediction-freeze.json').read_text())
    assert all(sha(ROOT/p)==h for p,h in freeze['files'].items())
    repeat=ROOT/'checks/fresh-replay'; comparison=[]
    for path in sorted(ROOT.glob('model/*.npz'))+sorted(ROOT.glob('results/*.npz')):
        other=repeat/path.relative_to(ROOT);a=np.load(path);b=np.load(other);assert set(a.files)==set(b.files)
        differences=[]
        for k in a.files:
            try:
                differences.append(compare_array(k, a[k], b[k]))
            except AssertionError as exc:
                raise AssertionError((str(path.relative_to(ROOT)), str(exc))) from exc
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
        compare_array('score',z['score'],score)
        compare_array('order',z['order'],order)
        assert ids[order[:5]].tolist()==[r['canonical'] for r in e['methods'][name]['top5']]
    dump(ROOT/'checks/verification-current.json',{'status':'passed','verifier_sha256':sha(ROOT/'verify.py'),'comparison_policy':{'absolute_tolerance':ATOL,'raw_relative_tolerance':RAW_RTOL,'raw_scale_fields':sorted(RAW_SCALE_FIELDS),'exact_fields':sorted(EXACT_FIELDS),'default_floating_policy':'strict absolute','nonfinite_values':'rejected'},'historical_receipt':'checks/verification.json','all_numeric_values_bitwise_equal_on_fresh_replay':all(a['bitwise_equal_values'] for f in comparison for a in f['arrays']),'archive_count':len(comparison),'comparison':comparison,'matrix_invariants':checks,'unknown_labels_excluded_from_null':True,'unique_matched_configurations':len(configs),'stratum_counts_preserved':True,'prediction_protocol_sha256':freeze['protocol_sha256'],'scope':'Numerical replay agreement under the recorded field-specific tolerances; this does not establish physical accuracy or biological validity.'})
    print(json.dumps({'status':'passed','archives':len(comparison),'bitwise_equal':all(a['bitwise_equal_values'] for f in comparison for a in f['arrays'])}))
if __name__=='__main__':main()
