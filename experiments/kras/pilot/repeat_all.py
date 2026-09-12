"""Repeat the frozen pilot and diagnose numerical cross-platform agreement.

Scientific scripts are unchanged. Float calculations, discrete identities and
score-based ranks have separate rules. A fresh receipt precedes execution, and
complete diagnostics precede failure, preventing reuse of an old success record.
"""
from pathlib import Path
import argparse, datetime, json, platform, shutil, subprocess, sys, time
import numpy as np
from prepare import ROOT, sha

# Engineering comparison tolerances, not biological/discretization estimates.
# The moment normalization below was revised AFTER Linux run 34689513534;
# the existing 1e-10 dimensionless tolerance itself was not increased.
# The C tolerance is seven orders below the 1e-3 grid gate. Receiver RMS scores
# are 1-Lipschitz in the maximum C-entry norm. Other floats require unit scaling.
RESPONSE_ATOL = 1e-10
INTERMEDIATE_ATOL = 1e-12
INTERMEDIATE_RTOL = 1e-10
RESPONSE_FIELDS = {'C', 'equilibrium', 'score'}
MOMENT_FIELDS = {'monomial_covariance', 'monomial_time_covariance', 'monomial_centered'}
POLICY = {
    'normalized_response_and_score': {'absolute_tolerance': RESPONSE_ATOL, 'relative_tolerance': 0.0},
    'other_floating_fields': {'absolute_tolerance': INTERMEDIATE_ATOL, 'relative_tolerance': INTERMEDIATE_RTOL,
                              'scale': 'maximum absolute value of the two compared entries'},
    'monomial_moments': {'dimensionless_absolute_tolerance': RESPONSE_ATOL,
        'covariance_and_delayed_scale': 'common SD outer product; SD=sqrt(maximum of the two equilibrium covariance diagonals), reused for both matrices',
        'centered_vector_scale': 'common maximum absolute entry per column across the two centered arrays',
        'invalid_scale': 'missing, nonfinite or nonpositive equilibrium variances/column scales fail; no field is skipped'},
    'discrete_fields': 'exact dtype, shape and values; order uses the separate score-aware rule',
    'rank_rule': 'same unique eligible node IDs; each order exactly sorted by its own finite scores/canonical tie-break; every cross-run inverted pair must have a score gap <=2e-10 in both runs',
    'rank_pair_tolerance': 2*RESPONSE_ATOL,
    'rationale': 'The absolute C/score tolerance is seven orders below the frozen 0.001 grid/domain gate. RMS-score error is bounded by maximum C-entry error. A score gap exceeding twice the numerical tolerance must not reverse. Covariance entries near zero can result from cancellation, so their own entry magnitude is not a suitable scale: equilibrium SD products give dimensionless moment errors, including delayed covariance. Centered arrays use their column amplitude. Other dimensional intermediates retain relative scaling.',
    'revision_history': ['Initial diagnostic policy: absolute C/score tolerance 1e-10 and entry-relative intermediate tolerance 1e-10.',
        'After inspecting run34689513534: use common equilibrium fluctuation scales for monomial covariance/time covariance and common column amplitudes for centered arrays; retain the existing 1e-10 dimensionless criterion. The previous failed receipt and Linux arrays remain part of the audit history.'],
    'interpretation': 'Passing classifies numerical agreement under declared scales. It does not prove the origin of every residual or establish physical/biological model accuracy. Actual C, equilibrium and score comparisons are unchanged.',
}


def bitwise_equal(a, b):
    return a.dtype == b.dtype and a.shape == b.shape and a.tobytes() == b.tobytes()


def numeric_field(a, b, name, reference=None, repeated=None):
    item = {'shape': list(a.shape), 'reference_dtype': str(a.dtype), 'repeated_dtype': str(b.dtype),
            'bitwise_identical': bitwise_equal(a, b)}
    if a.shape != b.shape or a.dtype != b.dtype:
        return dict(item, passed=False, reason='shape or dtype differs', repeated_shape=list(b.shape))
    if a.dtype.kind not in 'fc':
        differing = np.argwhere(a != b)
        return dict(item, passed=bool(np.array_equal(a,b)), rule='exact discrete identity',
                    differing_entries=int(len(differing)),
                    first_differing_index=differing[0].tolist() if len(differing) else None)
    if not (np.isfinite(a).all() and np.isfinite(b).all()):
        return dict(item, passed=False, reason='nonfinite floating value',
                    reference_nonfinite=int(np.count_nonzero(~np.isfinite(a))),
                    repeated_nonfinite=int(np.count_nonzero(~np.isfinite(b))))
    difference = np.abs(a-b)
    scaling = {}
    if name in MOMENT_FIELDS:
        if reference is None or repeated is None:
            return dict(item,passed=False,reason='moment comparison requires both complete array records')
        if name == 'monomial_centered':
            if a.ndim != 2:
                return dict(item,passed=False,reason='centered monomials must be a two-dimensional array')
            scale = np.maximum(np.max(abs(a),axis=0),np.max(abs(b),axis=0))[None,:]
            label = 'common maximum absolute centered-vector entry per column'
        else:
            covariances = [record.get('monomial_covariance') for record in (reference,repeated)]
            if any(cov is None or cov.ndim!=2 or cov.shape!=a.shape or cov.shape[0]!=cov.shape[1]
                   or not np.isfinite(cov).all() for cov in covariances):
                return dict(item,passed=False,reason='missing, incompatible or nonfinite equilibrium monomial covariance')
            diagonals = [np.diag(cov) for cov in covariances]
            if any(np.any(diagonal<=0) for diagonal in diagonals):
                return dict(item,passed=False,reason='nonpositive equilibrium monomial variance')
            sd = np.sqrt(np.maximum(*diagonals))
            scale = np.outer(sd,sd)
            label = 'common equilibrium monomial SD outer product for both equal-time and delayed covariance'
        if not np.isfinite(scale).all() or np.any(scale<=0):
            return dict(item,passed=False,reason='nonfinite or nonpositive monomial comparison scale')
        atol, rtol = 0.0, 0.0
        allowed = RESPONSE_ATOL*scale
        scaling = {'normalization':label,'dimensionless_absolute_tolerance':RESPONSE_ATOL,
                   'normalization_scale_min':float(np.min(scale)),'normalization_scale_max':float(np.max(scale)),
                   'max_dimensionless_difference':float(np.max(difference/scale)) if a.size else 0.0}
    else:
        atol, rtol = (RESPONSE_ATOL,0.0) if name in RESPONSE_FIELDS else (INTERMEDIATE_ATOL,INTERMEDIATE_RTOL)
        allowed = atol + rtol*np.maximum(np.abs(a),np.abs(b))
    ratio = difference/allowed
    index = np.unravel_index(int(np.argmax(difference)),a.shape) if a.size else ()
    normalized_index = np.unravel_index(int(np.argmax(ratio)),a.shape) if a.size else ()
    def value(array, position):
        x = array[position].item()
        return {'real':float(x.real),'imag':float(x.imag)} if isinstance(x,complex) else float(x)
    return dict(item, passed=bool(np.all(difference<=allowed)), absolute_tolerance=atol,
                relative_tolerance=rtol, max_absolute_difference=float(np.max(difference)) if a.size else 0.0,
                max_tolerance_ratio=float(np.max(ratio)) if a.size else 0.0,
                violating_entries=int(np.count_nonzero(difference>allowed)),
                worst_absolute_index=[int(i) for i in index], worst_ratio_index=[int(i) for i in normalized_index],
                reference_at_worst_absolute=value(a,index) if a.size else None,
                repeated_at_worst_absolute=value(b,index) if a.size else None,**scaling)


def compare_order(reference, repeated, reference_model, repeated_model):
    """Allow only numerically unresolved reversals, never changed node IDs."""
    result = {'passed':False,'rule':POLICY['rank_rule'],'inverted_pairs':[]}
    if 'score' not in reference or 'score' not in repeated:
        return dict(result,reason='order requires its corresponding scores')
    orders = [reference['order'],repeated['order']]
    scores = [reference['score'],repeated['score']]
    for label,data,model in zip(('reference','repeated'),(reference,repeated),(reference_model,repeated_model)):
        order,score = data['order'],data['score']
        eligible = np.flatnonzero(model['candidate'])
        canonical = model['canonical']
        if (order.dtype.kind not in 'iu' or order.ndim!=1 or score.shape!=canonical.shape
                or not np.isfinite(score).all() or not np.array_equal(np.sort(order),eligible)):
            return dict(result,reason=f'{label} ranking is not a unique permutation of all eligible node IDs')
        expected = eligible[np.lexsort((canonical[eligible],-score[eligible]))]
        if not np.array_equal(order,expected):
            return dict(result,reason=f'{label} ranking disagrees with its own scores/canonical tie-break')
        if 'C' in data:
            recalculated = np.sqrt(np.mean(data['C'][:,model['receiver']]**2,axis=1))
            error = float(np.max(abs(recalculated-score)))
            result[label+'_score_reconstruction_error'] = error
            if error>RESPONSE_ATOL:
                return dict(result,reason=f'{label} scores do not reconstruct from C')
    if orders[0].dtype!=orders[1].dtype or not np.array_equal(np.sort(orders[0]),np.sort(orders[1])):
        return dict(result,reason='eligible identities or ranking dtype differ across runs')
    if not np.array_equal(reference_model['canonical'],repeated_model['canonical']):
        return dict(result,reason='canonical identities differ across runs')
    canonical = reference_model['canonical']
    positions = {int(node):rank for rank,node in enumerate(orders[1])}
    for first_position,first in enumerate(orders[0]):
        for second in orders[0][first_position+1:]:
            if positions[int(first)]>positions[int(second)]:
                gaps = [float(abs(score[first]-score[second])) for score in scores]
                result['inverted_pairs'].append({'node_indices':[int(first),int(second)],
                    'canonical_residues':[int(canonical[first]),int(canonical[second])],
                    'reference_score_gap':gaps[0],'repeated_score_gap':gaps[1],
                    'allowed_gap':2*RESPONSE_ATOL,'numerically_unresolved':all(gap<=2*RESPONSE_ATOL for gap in gaps)})
    result.update(passed=all(row['numerically_unresolved'] for row in result['inverted_pairs']),
                  bitwise_identical=bitwise_equal(*orders),
                  reference_top5_canonical=canonical[orders[0][:5]].tolist(),
                  repeated_top5_canonical=canonical[orders[1][:5]].tolist(),
                  top5_set_identical=set(orders[0][:5])==set(orders[1][:5]))
    return result


def load_npz(path):
    with np.load(path,allow_pickle=False) as data:
        return {name:data[name] for name in data.files}


def compare_directories(reference, repeated):
    """Collect every field diagnostic before deciding whether to raise."""
    files, failures, worst = [], [], {}
    try:
        models = [load_npz(folder/'model/static-model.npz') for folder in (reference,repeated)]
    except Exception as error:
        return {'comparison_passed':False,'failures':[{'scope':'static model load','error_type':type(error).__name__}],
                'result_files':[],'worst_by_floating_field':{},'all_arrays_bitwise_identical':False}
    expected_paths = sorted((reference/'model').glob('*.npz'))+sorted((reference/'results').glob('*.npz'))
    for original in expected_paths:
        relative = original.relative_to(reference)
        item = {'file':relative.as_posix(),'fields':{}}
        try:
            a,b = load_npz(original),load_npz(repeated/relative)
            item['missing_fields'] = sorted(set(a)-set(b))
            item['unexpected_fields'] = sorted(set(b)-set(a))
            for name in sorted(set(a)&set(b)):
                detail = compare_order(a,b,*models) if name=='order' else numeric_field(a[name],b[name],name,a,b)
                item['fields'][name] = detail
                if 'max_tolerance_ratio' in detail:
                    entry = {'file':relative.as_posix(),**detail}
                    if name not in worst or entry['max_tolerance_ratio']>worst[name]['max_tolerance_ratio']:
                        worst[name] = entry
            item['passed'] = not item['missing_fields'] and not item['unexpected_fields'] and all(x['passed'] for x in item['fields'].values())
            item['all_arrays_bitwise_identical'] = set(a)==set(b) and all(bitwise_equal(a[k],b[k]) for k in a)
        except Exception as error:
            item.update(passed=False,error_type=type(error).__name__,all_arrays_bitwise_identical=False)
        if not item['passed']:
            failures.append({'file':relative.as_posix(),'fields':[k for k,v in item['fields'].items() if not v['passed']],
                             'error_type':item.get('error_type'),'missing_fields':item.get('missing_fields',[]),
                             'unexpected_fields':item.get('unexpected_fields',[])})
        files.append(item)
    return {'comparison_passed':bool(files) and not failures,'result_files':files,'failures':failures,
            'worst_by_floating_field':worst,
            'all_arrays_bitwise_identical':bool(files) and all(x['all_arrays_bitwise_identical'] for x in files)}


def write_receipt(report):
    destination = ROOT/'full-repeat-verification.json'
    temporary = destination.with_suffix('.json.tmp')
    temporary.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    temporary.replace(destination)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--compare-only',action='store_true',help='Compare an existing repeat without rerunning scientific work')
    parser.add_argument('--candidate-dir',type=Path,default=ROOT/'repeat-all')
    args = parser.parse_args()
    work = args.candidate_dir.resolve()
    if work==ROOT.resolve():parser.error('The repeat directory must differ from the frozen pilot directory')
    begin = time.monotonic()
    report = {'status':'running','repeated_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
              'scope':'fresh preparation, all24 finite-grid cases, eight analytic cases and graph/geometric controls; or explicit comparison-only replay',
              'comparison_only':args.compare_only,'comparison_policy':POLICY,'runs':[],
              'environment':{'python':platform.python_version(),'numpy':np.__version__,
                             'system':platform.system(),'architecture':platform.machine()},
              'cross_platform_guarantee':False,'comparator_sha256':sha(__file__)}
    write_receipt(report)
    try:
        if not args.compare_only:
            work.mkdir(parents=True,exist_ok=True)
            (work/'model').mkdir(exist_ok=True)
            shutil.copytree(ROOT/'raw',work/'raw',dirs_exist_ok=True)
            for name in ['prepare.py','run_pilot.py','preanalysis-protocol.json','preanalysis-protocol.sha256']:
                shutil.copyfile(ROOT/name,work/name)
            for script in ['prepare.py','run_pilot.py']:
                process = subprocess.run([sys.executable,'-W','error',script],cwd=work,capture_output=True,text=True)
                (work/(script+'.log')).write_text(process.stdout+process.stderr)
                report['runs'].append({'script':script,'code_sha256':sha(ROOT/script),'returncode':process.returncode})
                write_receipt(report)
                if process.returncode:raise RuntimeError(f'{script} failed; see the phase log')
        report.update(compare_directories(ROOT,work))
        report['status'] = 'passed' if report['comparison_passed'] else 'failed'
    except Exception as error:
        report.update(status='failed',execution_error_type=type(error).__name__)
    report['seconds'] = time.monotonic()-begin
    write_receipt(report)
    print(json.dumps({'status':report['status'],'files':len(report.get('result_files',[])),
                      'all_identical':report.get('all_arrays_bitwise_identical',False),
                      'failed_files':report.get('failures',[]),'seconds':report['seconds']}))
    if report['status']!='passed':
        raise RuntimeError('Pilot replay failed; full-repeat-verification.json contains complete available diagnostics')


if __name__=='__main__':main()
