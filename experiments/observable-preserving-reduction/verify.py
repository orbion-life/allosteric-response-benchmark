"""Portable semantic verification; scientific rejection is not software failure.

Cross-platform comparisons use the declared 1e-9 absolute allowance on physical
responses, covariances and bounds. Raw reduced-coordinate arrays have basis
gauge freedom; check their shapes, positivity and physical reconstruction instead.
The archived local audit separately compared every array on the original host.
"""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import math
import platform
from pathlib import Path

import numpy as np
import scipy
from scipy.linalg import eigh

ROOT = Path(__file__).resolve().parent
TOL = 1e-9
GAUGE_FIELDS = {'Hr', 'coefficients', 'residual_Gram'}
PHYSICAL_FIELDS = {'C', 'delayed', 'equilibrium', 'bound', 'one_sided_bound',
                   'two_sided_bound', 'uniform_bound'}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_json(path):
    def reject(value):
        raise ValueError('Nonfinite JSON value: ' + value)
    return json.loads(Path(path).read_text(), parse_constant=reject)


def read_arrays(path):
    with np.load(path, allow_pickle=False) as archive:
        return {key: archive[key].copy() for key in archive.files}


def product(a, b):
    return np.einsum('ij,jk->ik', a, b, optimize=False)


def physical_comparison(expected, actual, tolerance=TOL):
    """Return full diagnostics; use absolute units, not relative entry scaling."""
    result = {'shape_reference': list(expected.shape), 'shape_replay': list(actual.shape),
              'dtype_reference': str(expected.dtype), 'dtype_replay': str(actual.dtype)}
    if expected.shape != actual.shape or expected.dtype != actual.dtype:
        return dict(result, passed=False, reason='shape or dtype mismatch')
    if expected.dtype.kind not in 'biufc':
        return dict(result, passed=False, reason='non-numeric array')
    if not np.isfinite(expected).all() or not np.isfinite(actual).all():
        return dict(result, passed=False, reason='nonfinite array')
    if expected.dtype.kind in 'biu':
        return dict(result, passed=bool(np.array_equal(expected, actual)), policy='exact discrete identity')
    difference = np.abs(expected - actual)
    maximum = float(difference.max()) if difference.size else 0.
    index = list(np.unravel_index(int(difference.argmax()), difference.shape)) if difference.size else []
    return dict(result, passed=maximum <= tolerance, policy='absolute physical-array tolerance',
                tolerance=tolerance, max_absolute_difference=maximum,
                worst_index=[int(x) for x in index])


def gauge_consistency(arrays, prefix, times):
    """Check physical reconstruction without requiring a particular basis gauge."""
    h = arrays[prefix + 'Hr']; b = arrays[prefix + 'coefficients']
    residual = arrays[prefix + 'residual_Gram']; m = b.shape[0]
    if m != int(prefix[4:-1]) or h.shape != (m, m) or residual.shape != (m, m) or b.shape != (m, 3):
        return {'passed': False, 'reason': 'reduced-coordinate shape mismatch'}
    if not all(np.isfinite(a).all() for a in (h, b, residual)):
        return {'passed': False, 'reason': 'nonfinite reduced-coordinate array'}
    hs = float(np.max(np.abs(h-h.T))); rs = float(np.max(np.abs(residual-residual.T)))
    eigenvalues, u = eigh((h+h.T)/2)
    residual_eigenvalues = eigh((residual+residual.T)/2, eigvals_only=True)
    hscale = max(1., float(np.max(np.abs(eigenvalues))))
    rscale = max(1., float(np.max(np.abs(residual_eigenvalues))))
    static = product(b.T, b); rotated = product(u.T, b)
    delayed = np.array([product(rotated.T*np.exp(-t*np.maximum(eigenvalues, 0)), rotated) for t in times])
    checks = {name: physical_comparison(arrays[prefix+name], value) for name, value in
              [('equilibrium', static), ('delayed', delayed), ('C', delayed-static)]}
    return {'passed': bool(hs <= TOL and rs <= TOL
                           and eigenvalues[0] >= -1e-10*hscale
                           and residual_eigenvalues[0] >= -TOL*rscale
                           and all(c['passed'] for c in checks.values())),
            'policy': 'basis-gauge invariant reconstruction and PSD checks; no raw basis comparison',
            'Hr_symmetry_error': hs, 'residual_Gram_symmetry_error': rs,
            'minimum_Hr_eigenvalue': float(eigenvalues[0]),
            'minimum_residual_Gram_eigenvalue': float(residual_eigenvalues[0]),
            'reconstruction': checks}


class Review:
    def __init__(self):
        self.failures = []

    def require(self, condition, context, **detail):
        if not condition:
            self.failures.append({'context': context, **detail})
        return bool(condition)


def expected_cases(protocol):
    return {(g['id'], k, n, protocol['domain_sigma_slowest']): False
            for g in protocol['geometries'] for k in protocol['kappa'] for n in protocol['grid_n']} | {
                (g['id'], k, protocol['domain_check']['n'], protocol['domain_check']['extent']): True
                for g in protocol['geometries'] for k in protocol['kappa']}


def campaign(folder, protocol, review, label):
    run = read_json(folder/'run.json'); watchdog = read_json(folder/'watchdog.json')
    review.require(run['status'] == watchdog['status'] == 'COMPLETE', label+': completion')
    review.require(watchdog['returncode'] == 0 and watchdog['limit_reason'] is None, label+': watchdog')
    review.require(watchdog['limits'] == protocol['resources'], label+': resource protocol')
    review.require(0 <= watchdog['wall_seconds'] <= protocol['resources']['max_total_wall_seconds'], label+': wall cap')
    review.require(0 <= watchdog['sampled_peak_aggregate_bytes'] <= protocol['resources']['max_aggregate_bytes'], label+': sampled memory cap')
    review.require(run['protocol_sha256'] == sha(ROOT/'preanalysis.json'), label+': protocol hash')
    review.require(run['script_sha256'] == sha(ROOT/'operator_study.py'), label+': scientific source hash')
    declarations=[protocol['created_utc']]
    declarations.extend(protocol[k]['date'] for k in ['pre_execution_metadata_correction','pre_execution_numeric_clarification'] if k in protocol)
    review.require(max(datetime.datetime.fromisoformat(x) for x in declarations)
                   < datetime.datetime.fromisoformat(run['started_utc']), label+': protocol recorded before run')
    expected = expected_cases(protocol)
    expected_names={f'{g}-k{k}-n{n}-e{e}' for g,k,n,e in expected}
    review.require({p.stem for p in folder.glob('*.npz')} == expected_names,label+': exact archive coverage')
    identities = [(x['geometry'], x['kappa'], x['n'], x['extent']) for x in run['cases']]
    review.require(len(identities) == len(set(identities)) == len(expected) and set(identities) == set(expected), label+': exact case coverage')
    outputs = {}; summary = []; gauges = []; checkpoint_count = 0
    for identity in expected:
        geometry, kappa, n, extent = identity
        name = f'{geometry}-k{kappa}-n{n}-e{extent}'
        try:
            c = read_json(folder/(name+'.json')); z = read_arrays(folder/(name+'.npz')); outputs[name] = z
            matches = [item for item in run['cases'] if (item['geometry'], item['kappa'], item['n'], item['extent']) == identity]
            review.require(len(matches) == 1 and c == matches[0], label+': '+name+': run/case receipt equality')
            review.require(c['reference_only'] == expected[identity] and c['states'] == n**3, label+': '+name+': model identity')
            review.require(np.array_equal(z['times'], c['times']) and z['times'].shape == (3,), label+': '+name+': query times')
            review.require(z['harmonic_sd'].shape == (3,) and np.all(z['harmonic_sd'] > 0), label+': '+name+': positive common scales')
            for key, value in z.items():
                review.require(value.dtype.kind in 'biufc' and np.isfinite(value).all(), label+': '+name+': finite '+key)
            for key, shape in [('reference_C', (3,3,3)), ('reference_delayed', (3,3,3)), ('reference_equilibrium', (3,3))]:
                review.require(z[key].shape == shape, label+': '+name+': '+key+' shape')
            ref_reconstruction = physical_comparison(z['reference_C'], z['reference_delayed']-z['reference_equilibrium'])
            review.require(ref_reconstruction['passed'], label+': '+name+': reference reconstruction', diagnostic=ref_reconstruction)
            if c['reference_only']:
                review.require(set(z) == {'reference_C','reference_delayed','reference_equilibrium','times','harmonic_sd'}, label+': '+name+': reference-only fields')
                continue
            rows = []; checkpoints = c['checkpoints']
            ranks = [row['rank'] for row in checkpoints]
            review.require(ranks == sorted(set(ranks)) and len(ranks)>0 and ranks[-1] <= protocol['basis']['maximum_rank'], label+': '+name+': rank progression')
            for row in checkpoints:
                prefix = f"rank{row['rank']}_"; checkpoint_count += 1
                static = float(np.max(np.abs(z[prefix+'equilibrium']-z['reference_equilibrium'])))
                delayed = float(np.max(np.abs(z[prefix+'delayed']-z['reference_delayed'])))
                error = np.abs(z[prefix+'C']-z['reference_C']); bound = z[prefix+'bound']
                metrics = {'static_max_error': static, 'delayed_max_error': delayed,
                           'response_max_error': float(error.max()), 'max_response_bound': float(bound.max()),
                           'minimum_bound_slack': float((bound-error).min())}
                for key, value in metrics.items():
                    review.require(math.isfinite(row[key]) and abs(row[key]-value) <= TOL, label+': '+name+': '+prefix+key, saved=row[key], recomputed=value)
                review.require(np.all(bound >= 0), label+': '+name+': nonnegative response bound')
                coverage = bool(np.all(error <= bound+protocol['gates']['bound_comparison_arithmetic_allowance']))
                stop = bool(bound.max() <= protocol['gates']['calculated_response_error_bound']
                            and static <= protocol['gates']['normalized_static_error']
                            and row['orthogonality_error'] <= protocol['gates']['basis_orthogonality'])
                review.require(coverage == row['actual_error_within_calculated_bound'], label+': '+name+': bound coverage flag')
                review.require(stop == row['bound_stop_pass'], label+': '+name+': pre-reference selection flag')
                review.require(0 <= row['operator_seconds_to_checkpoint'] <= c['operator_seconds'], label+': '+name+': cumulative checkpoint cost')
                review.require(abs(row['dimension_reduction']-n**3/row['rank']) <= TOL, label+': '+name+': dimension accounting')
                g = gauge_consistency(z, prefix, z['times']); gauges.append({'case': name, 'rank': row['rank'], **g})
                review.require(g['passed'], label+': '+name+': '+prefix+'gauge consistency', diagnostic=g)
                rows.append({'rank': row['rank'], 'stop': stop, 'coverage': coverage, **metrics})
            last = rows[-1]
            review.require(not any(row['stop'] for row in rows[:-1]), label+': '+name+': first passing checkpoint')
            review.require(c['selected_by_bound'] == last['stop'] and c['selected_rank'] == last['rank'], label+': '+name+': selected rank')
            finite = bool(last['stop'] and last['response_max_error'] <= protocol['gates']['normalized_response_error']
                          and last['delayed_max_error'] <= protocol['gates']['normalized_delayed_error'] and last['coverage'])
            for field in ['common_setup_seconds','operator_seconds','reference_seconds','total_operator_seconds','total_reference_seconds']:
                review.require(math.isfinite(c[field]) and c[field] >= 0, label+': '+name+': finite nonnegative cost '+field)
            for total, route in [('total_operator_seconds','operator_seconds'),('total_reference_seconds','reference_seconds')]:
                review.require(abs(c[total]-c['common_setup_seconds']-c[route]) <= 1e-12*max(1.,c[total]), label+': '+name+': setup charged once '+total)
            ratio = c['total_operator_seconds']/c['total_reference_seconds']
            review.require(abs(ratio-c['total_time_ratio']) <= 1e-12*max(1.,ratio), label+': '+name+': time ratio')
            cost = ratio <= protocol['gates']['total_wall_time_ratio_at_same_three_time_queries']
            practical = bool(finite and n**3/last['rank'] >= protocol['gates']['operator_dimension_reduction_factor'] and cost)
            review.require(c['finite_fidelity_pass'] == finite and c['practical_pass'] == practical, label+': '+name+': scientific decisions')
            summary.append({'case':name,'ranks':ranks,'selected_rank':last['rank'],'selected_by_bound':last['stop'],
                            'finite_fidelity_pass':finite,'cost_only_pass':bool(cost),'practical_pass':practical,
                            'time_ratio':ratio,'response_error':last['response_max_error'],'bound':last['max_response_bound']})
        except Exception as exc:
            review.require(False, label+': '+name+': readable complete case', exception=type(exc).__name__+': '+str(exc))
    discretization = []
    for g in protocol['geometries']:
        for k in protocol['kappa']:
            try:
                stem=f"{g['id']}-k{k}"
                a=outputs[stem+'-n17-e4']['reference_C']; b=outputs[stem+'-n33-e4']['reference_C']; c=outputs[stem+'-n41-e5']['reference_C']
                grid=float(np.max(np.abs(a-b))); domain=float(np.max(np.abs(b-c)))
                discretization.append({'case':stem,'grid_max':grid,'domain_max':domain,
                                       'grid_pass':grid <= protocol['gates']['grid_response_difference'],
                                       'domain_pass':domain <= protocol['gates']['domain_response_difference']})
            except Exception as exc:
                review.require(False,label+': discretization '+g['id'],exception=str(exc))
    receipt={'operator_cases':summary,'checkpoint_count':checkpoint_count,'discretization':discretization,
             'finite_fidelity_pass_count':sum(x['finite_fidelity_pass'] for x in summary),
             'cost_only_pass_count':sum(x['cost_only_pass'] for x in summary),
             'practical_pass_count':sum(x['practical_pass'] for x in summary),
             'grid_pass_count':sum(x['grid_pass'] for x in discretization),
             'domain_pass_count':sum(x['domain_pass'] for x in discretization),
             'run_wall_seconds':run['wall_seconds'],'watchdog':watchdog,'gauge_checks':gauges}
    return receipt,outputs


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--replay',type=Path,default=ROOT/'replay-independent')
    parser.add_argument('--output',type=Path,default=ROOT/'checks/verification-current.json')
    args=parser.parse_args(); review=Review()
    # Refuse to overwrite any frozen evidence even when verification fails.
    frozen=read_json(ROOT/'frozen-evidence-sha256.json')
    if args.output.resolve() in {(ROOT/name).resolve() for name in frozen}:
        parser.error('--output would overwrite frozen evidence')
    result={'status':'RUNNING','absolute_physical_array_tolerance':TOL,
            'source_sha256':sha(ROOT/'operator_study.py'),'protocol_sha256':sha(ROOT/'preanalysis.json'),
            'verifier_sha256':sha(__file__),'software':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__},
            'performance_policy':'Recompute cost decisions; do not require a timing or practical pass or cross-host timing agreement.',
            'basis_policy':'Raw Hr, coefficients and residual_Gram are gauge-dependent. Validate physical reconstruction/PSD; same-host full-array parity is archived separately.'}
    try:
        review.require(args.replay.resolve() != (ROOT/'results').resolve(),'replay must differ from primary evidence directory')
        for name,expected_hash in frozen.items():
            review.require((ROOT/name).is_file() and sha(ROOT/name)==expected_hash,'frozen evidence hash',file=name)
        for name,expected_hash in read_json(ROOT/'source-receipts.json').items():
            review.require(sha(ROOT/name)==expected_hash,'scientific source receipt',file=name)
        protocol=read_json(ROOT/'preanalysis.json')
        review.require(sha(ROOT/'preanalysis.json')==(ROOT/'preanalysis.sha256').read_text().strip(),'protocol digest')
        review.require(TOL==protocol['gates']['bound_comparison_arithmetic_allowance'],'declared numerical allowance')
        primary,pa=campaign(ROOT/'results',protocol,review,'primary')
        replay,ra=campaign(args.replay,protocol,review,'replay')
        result.update(primary=primary,replay=replay,frozen_hash_count=len(frozen))
        comparisons=[]; physical_count=0; gauge_count=0; worst=0.
        review.require(set(pa)==set(ra),'archive coverage equality')
        for name in sorted(pa.keys() & ra.keys()):
            a,b=pa[name],ra[name]; review.require(set(a)==set(b),'field coverage equality',case=name)
            for key in sorted(a.keys() & b.keys()):
                suffix=key.split('_',1)[1] if key.startswith('rank') else None
                if suffix in GAUGE_FIELDS:
                    review.require(a[key].shape==b[key].shape and a[key].dtype==b[key].dtype,
                                   'gauge array shape/dtype equality',case=name,field=key)
                    gauge_count+=1
                    continue
                allowed=key in {'reference_C','reference_delayed','reference_equilibrium','times','harmonic_sd'} or suffix in PHYSICAL_FIELDS
                review.require(allowed,'explicit field policy',case=name,field=key)
                comparison=physical_comparison(a[key],b[key]); physical_count+=1
                comparisons.append({'case':name,'field':key,**comparison})
                worst=max(worst,comparison.get('max_absolute_difference',0.))
                review.require(comparison['passed'],'semantic physical array',case=name,field=key,diagnostic=comparison)
        for a,b in zip(primary['operator_cases'],replay['operator_cases']):
            for key in ['case','ranks','selected_rank','selected_by_bound','finite_fidelity_pass']:
                review.require(a[key]==b[key],'rank/selection/scientific decision equality',case=a['case'],field=key,primary=a[key],replay=b[key])
        for a,b in zip(primary['discretization'],replay['discretization']):
            for key in ['case','grid_pass','domain_pass']:
                review.require(a[key]==b[key],'discretization decision equality',case=a['case'],field=key,primary=a[key],replay=b[key])
        result.update(archives_compared=len(pa.keys() & ra.keys()),physical_arrays_compared=physical_count,
                      gauge_arrays_checked_semantically=gauge_count,maximum_physical_absolute_difference=worst,
                      array_comparisons=comparisons)
    except Exception as exc:
        review.require(False,'verification could not finish',exception=type(exc).__name__+': '+str(exc))
    result.update(status='FAIL' if review.failures else 'PASS',failures=review.failures,
                  meaning='PASS means reproducible software/evidence under this comparison policy; scientific failures remain in the result.')
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:result[k] for k in ['status','archives_compared','physical_arrays_compared','gauge_arrays_checked_semantically','maximum_physical_absolute_difference'] if k in result}))
    if review.failures:
        raise SystemExit(1)


if __name__=='__main__':
    main()
