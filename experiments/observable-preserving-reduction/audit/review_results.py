"""Read-only independent audit of the frozen operator campaign and its replay."""
from pathlib import Path
import hashlib, io, json, platform, subprocess
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
ALLOWANCE = 1e-9


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dump(name, value):
    (ROOT / 'audit' / name).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def audit_campaign(folder, protocol):
    run = json.loads((folder / 'run.json').read_text())
    watchdog = json.loads((folder / 'watchdog.json').read_text())
    assert run['status'] == watchdog['status'] == 'COMPLETE'
    assert run['protocol_sha256'] == sha(ROOT / 'preanalysis.json')
    assert run['script_sha256'] == sha(ROOT / 'operator_study.py')
    assert watchdog['returncode'] == 0 and watchdog['limit_reason'] is None
    assert watchdog['wall_seconds'] < protocol['resources']['max_total_wall_seconds']
    assert watchdog['sampled_peak_aggregate_bytes'] < protocol['resources']['max_aggregate_bytes']
    expected = {(g['id'], k, n, 4) for g in protocol['geometries'] for k in protocol['kappa'] for n in protocol['grid_n']}
    expected |= {(g['id'], k, 41, 5) for g in protocol['geometries'] for k in protocol['kappa']}
    assert {(c['geometry'], c['kappa'], c['n'], c['extent']) for c in run['cases']} == expected
    assert len(run['cases']) == len(expected) == 27
    records, checkpoints = [], 0
    for c in run['cases']:
        name = f"{c['geometry']}-k{c['kappa']}-n{c['n']}-e{c['extent']}"
        assert c == json.loads((folder / (name + '.json')).read_text())
        z = np.load(folder / (name + '.npz'))
        assert all(np.isfinite(z[k]).all() for k in z.files)
        assert c['states'] == c['n'] ** 3
        np.testing.assert_array_equal(z['times'], c['times'])
        if c['reference_only']:
            assert (c['n'], c['extent']) == (41, 5)
            continue
        checkpoint_passes = []
        for r in c['checkpoints']:
            prefix = f"rank{r['rank']}_"
            error = np.abs(z[prefix + 'C'] - z['reference_C'])
            delay_error = np.abs(z[prefix + 'delayed'] - z['reference_delayed'])
            static_error = np.abs(z[prefix + 'equilibrium'] - z['reference_equilibrium'])
            bound = z[prefix + 'bound']
            assert np.all(bound >= 0)
            metrics = {
                'response_max_error': float(error.max()),
                'delayed_max_error': float(delay_error.max()),
                'static_max_error': float(static_error.max()),
                'max_response_bound': float(bound.max()),
                'minimum_bound_slack': float((bound - error).min()),
            }
            for k, v in metrics.items():
                assert abs(r[k] - v) <= 1e-14, (name, r['rank'], k, r[k], v)
            covers = bool(np.all(error <= bound + ALLOWANCE))
            assert covers == r['actual_error_within_calculated_bound']
            stop = bool(metrics['max_response_bound'] <= protocol['gates']['calculated_response_error_bound']
                        and metrics['static_max_error'] <= protocol['gates']['normalized_static_error']
                        and r['orthogonality_error'] <= protocol['gates']['basis_orthogonality'])
            assert stop == r['bound_stop_pass']
            checkpoint_passes.append(stop)
            assert r['operator_seconds_to_checkpoint'] <= c['operator_seconds']
            assert r['rank'] <= protocol['basis']['maximum_rank']
            checkpoints += 1
        assert not any(checkpoint_passes[:-1]), name
        assert c['selected_by_bound'] == checkpoint_passes[-1]
        last = c['checkpoints'][-1]
        assert c['selected_rank'] == last['rank']
        if not c['selected_by_bound']:
            assert c['selected_rank'] == protocol['basis']['maximum_rank'], name
        for field, other in [('total_operator_seconds', 'operator_seconds'), ('total_reference_seconds', 'reference_seconds')]:
            assert abs(c[field] - c['common_setup_seconds'] - c[other]) < 1e-12
        ratio = c['total_operator_seconds'] / c['total_reference_seconds']
        assert abs(ratio - c['total_time_ratio']) < 1e-12
        finite = bool(c['selected_by_bound']
                      and last['response_max_error'] <= protocol['gates']['normalized_response_error']
                      and last['delayed_max_error'] <= protocol['gates']['normalized_delayed_error']
                      and last['actual_error_within_calculated_bound'])
        practical = bool(finite and c['states'] / last['rank'] >= protocol['gates']['operator_dimension_reduction_factor']
                         and ratio <= protocol['gates']['total_wall_time_ratio_at_same_three_time_queries'])
        assert finite == c['finite_fidelity_pass'] and practical == c['practical_pass']
        records.append({'case': name, 'rank': last['rank'], 'bound_selection_pass': c['selected_by_bound'],
                        'response_error': last['response_max_error'], 'bound': last['max_response_bound'],
                        'finite_fidelity_pass': finite, 'practical_pass': practical,
                        'total_time_ratio': ratio, 'operator_total_seconds': c['total_operator_seconds'],
                        'reference_total_seconds': c['total_reference_seconds'],
                        'static_error': last['static_max_error']})
    discretization = []
    for g in protocol['geometries']:
        for k in protocol['kappa']:
            label = f"{g['id']}-k{k}"
            fine = np.load(folder / f'{label}-n33-e4.npz')['reference_C']
            coarse = np.load(folder / f'{label}-n17-e4.npz')['reference_C']
            larger = np.load(folder / f'{label}-n41-e5.npz')['reference_C']
            grid = np.max(np.abs(fine - coarse), axis=(1, 2))
            domain = np.max(np.abs(larger - fine), axis=(1, 2))
            discretization.append({'case': label, 'grid_by_time': grid.tolist(), 'domain_by_time': domain.tolist(),
                                   'grid_max': float(grid.max()), 'domain_max': float(domain.max()),
                                   'grid_pass': bool(grid.max() <= protocol['gates']['grid_response_difference']),
                                   'domain_pass': bool(domain.max() <= protocol['gates']['domain_response_difference'])})
    return {'status': 'PASS', 'numerical_problem_count': 27, 'operator_problem_count': len(records),
            'checkpoint_count': checkpoints, 'selected_cases': records, 'discretization': discretization,
            'finite_fidelity_pass_count': sum(r['finite_fidelity_pass'] for r in records),
            'practical_pass_count': sum(r['practical_pass'] for r in records),
            'all_selected_observed_response_errors_pass': all(r['response_error'] <= .002 for r in records),
            'run': {k: v for k, v in run.items() if k != 'cases'}, 'watchdog': watchdog}


def main():
    protocol = json.loads((ROOT / 'preanalysis.json').read_text())
    assert protocol['gates']['bound_comparison_arithmetic_allowance'] == ALLOWANCE
    before = json.loads((ROOT / 'audit/independent-before-replay-hashes.json').read_text())
    for name, h in before.items():
        assert sha(ROOT / name) == h, name
    primary = audit_campaign(ROOT / 'results', protocol)
    replay = audit_campaign(ROOT / 'replay-independent', protocol)
    numeric, arrays, max_error = [], 0, 0.
    for p in sorted((ROOT / 'results').glob('*.npz')):
        a = np.load(p); b = np.load(ROOT / 'replay-independent' / p.name)
        assert set(a.files) == set(b.files)
        fields = []
        for key in a.files:
            assert a[key].dtype == b[key].dtype and a[key].shape == b[key].shape
            assert np.isfinite(a[key]).all() and np.isfinite(b[key]).all()
            difference = float(np.max(np.abs(a[key] - b[key]))) if a[key].size else 0.
            assert difference <= ALLOWANCE, (p.name, key, difference)
            if a[key].dtype.kind in 'biu':
                assert np.array_equal(a[key], b[key])
            fields.append({'array': key, 'max_absolute_difference': difference,
                           'values_exact': bool(np.array_equal(a[key], b[key]))})
            arrays += 1; max_error = max(max_error, difference)
        numeric.append({'archive': p.name, 'fields': fields})
    timing = []
    for a, b in zip(primary['selected_cases'], replay['selected_cases']):
        assert a['case'] == b['case'] and a['rank'] == b['rank']
        assert a['bound_selection_pass'] == b['bound_selection_pass']
        assert a['finite_fidelity_pass'] == b['finite_fidelity_pass']
        timing.append({'case': a['case'], 'primary_time_ratio': a['total_time_ratio'],
                       'replay_time_ratio': b['total_time_ratio'],
                       'primary_practical_pass': a['practical_pass'], 'replay_practical_pass': b['practical_pass']})
    # Read immutable release content via Git, not a possibly modified working copy.
    repo = ROOT.parents[3] / 'pulsar-response-open-source'
    tag = subprocess.check_output(['git', '-C', str(repo), 'rev-parse', 'v0.3.0'], text=True).strip()
    commit = subprocess.check_output(['git', '-C', str(repo), 'rev-parse', 'v0.3.0^{}'], text=True).strip()
    original = []
    for k in protocol['kappa']:
        for n in protocol['grid_n']:
            rel = f'experiments/harmonic-completion/results/biquadratic-k{k}-d3-n{n}-e4.npz'
            raw = subprocess.check_output(['git', '-C', str(repo), 'show', f'v0.3.0:{rel}'])
            old = np.load(io.BytesIO(raw)); new = np.load(ROOT / 'results' / f'development-k{k}-n{n}-e4.npz')
            diff = {}
            for field in ['C', 'equilibrium', 'delayed']:
                value = float(np.max(np.abs(old[field] - new['reference_' + field])))
                assert value <= ALLOWANCE, (rel, field, value)
                diff[field] = value
            assert np.array_equal(old['times'], new['times'])
            assert np.array_equal(old['harmonic_sd'], new['harmonic_sd'])
            original.append({'release_path': rel, 'sha256': hashlib.sha256(raw).hexdigest(), 'max_differences': diff})
    dump('independent-primary-review.json', primary)
    dump('independent-replay-review.json', replay)
    dump('independent-numeric-comparison.json', {'status': 'PASS', 'absolute_tolerance': ALLOWANCE,
         'archives': len(numeric), 'arrays': arrays, 'maximum_absolute_difference': max_error,
         'archive_comparisons': numeric, 'timing_compared_separately': timing,
         'frozen_source_and_primary_files_unchanged': len(before),
         'software': {'python': platform.python_version(), 'numpy': np.__version__}})
    dump('independent-v0.3.0-parity.json', {'status': 'PASS', 'tag_object': tag, 'commit': commit,
         'archives': original, 'absolute_tolerance': ALLOWANCE})
    print(json.dumps({'status': 'PASS', 'arrays': arrays, 'max_replay_error': max_error,
                     'finite_fidelity_passes': primary['finite_fidelity_pass_count'],
                     'primary_practical_passes': primary['practical_pass_count'],
                     'replay_practical_passes': replay['practical_pass_count'],
                     'grid_passes': sum(r['grid_pass'] for r in primary['discretization']),
                     'domain_passes': sum(r['domain_pass'] for r in primary['discretization'])}))


if __name__ == '__main__':
    main()
