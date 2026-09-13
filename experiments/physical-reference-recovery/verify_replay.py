"""Compare fresh bounded-run evidence with the frozen physical archive.

Absolute 1e-9 is a numerical replay allowance, not the physical 0.001 screen,
independent-solver 1e-8 gate, or harmonic implementation 1e-10 gate. All 82
stored arrays are checked, including the implementation's sign-canonical mode
bases. No timing equality or nonlinear convergence success is required.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent
REPLAY_ATOL = 1e-9
PHYSICAL_SCREEN = 0.001
COMPONENTS = ('G0', 'K', 'C')
FINAL_PAIRS = {
    ('harmonic-A6-h0.0625', 'harmonic-A6-h0.03125', 'fine_grid'),
    ('harmonic-A6-h0.03125', 'harmonic-A6-h0.015625', 'fine_grid'),
    ('harmonic-A6-h0.015625', 'harmonic-A7-h0.015625', 'fine_domain'),
    ('harmonic-A7-h0.015625', 'harmonic-A8-h0.015625', 'fine_domain'),
}


def compare_array(actual, expected):
    """Return diagnostics even when a malformed array cannot be subtracted."""
    a, b = np.asarray(actual), np.asarray(expected)
    out = dict(shape=list(a.shape), expected_shape=list(b.shape),
               dtype=str(a.dtype), expected_dtype=str(b.dtype), pass_check=False)
    if a.shape != b.shape or a.dtype.kind != b.dtype.kind:
        out['reason'] = 'shape or numeric-kind mismatch'
        return out
    if a.dtype.kind not in 'fiu' or a.size == 0 or not np.isfinite(a).all() or not np.isfinite(b).all():
        out['reason'] = 'nonfinite, empty or nonnumeric array'
        return out
    delta = float(np.max(np.abs(a.astype(float)-b.astype(float))))
    out.update(maximum_absolute_difference=delta,
               pass_check=bool(np.array_equal(a, b)) if a.dtype.kind in 'iu' else delta <= REPLAY_ATOL)
    return out


def run(harmonic, nonlinear, reference=ROOT):
    reference, harmonic, nonlinear = map(Path, (reference, harmonic, nonlinear))
    checks, array_checks, screens, resources = [], [], [], []
    data = {}

    def check(name, passed, **details):
        checks.append(dict(name=name, pass_check=bool(passed), **details))

    def guarded(name, function):
        try:
            function()
        except (OSError, ValueError, KeyError, TypeError, IndexError, AssertionError, StopIteration, RuntimeWarning, FloatingPointError) as error:
            check(name, False, error=f'{type(error).__name__}: {error}')

    def read_json(path):
        return json.loads(path.read_text())

    def scalar_check(name, value, limit):
        valid = (type(value) in (int, float) and np.isfinite(value) and 0 <= value <= limit)
        check(name, valid, value=value if type(value) in (int, float) and np.isfinite(value) else repr(value), maximum=limit)

    def differences(a, b):
        for key in COMPONENTS:
            if a[key].shape != b[key].shape or not np.isfinite(a[key]).all() or not np.isfinite(b[key]).all():
                raise ValueError(f'Invalid shape or nonfinite physical component: {key}')
        return {k: float(np.max(np.abs(a[k]-b[k]))) for k in COMPONENTS}

    def hashes():
        manifest = read_json(reference/'execution-source-hashes.json')['sha256']
        for filename, expected in manifest.items():
            check(f'source/{filename}', hashlib.sha256((reference/filename).read_bytes()).hexdigest() == expected)
    guarded('source-manifest', hashes)
    protocol_hash = hashlib.sha256((reference/'protocol.json').read_bytes()).hexdigest()

    def load_archive(label, expected_dir, actual_dir):
        expected_paths = sorted(expected_dir.glob('*.npz'))
        check(f'{label}/file-set', bool(expected_paths) and
              {p.name for p in expected_paths} == {p.name for p in actual_dir.glob('*.npz')})
        for expected_path in expected_paths:
            name = f'{label}/{expected_path.name}'
            def load_one():
                with np.load(expected_path, allow_pickle=False) as z:
                    expected = {k: z[k] for k in z.files}
                with np.load(actual_dir/expected_path.name, allow_pickle=False) as z:
                    actual = {k: z[k] for k in z.files}
                check(f'{name}/keys', set(actual) == set(expected))
                for key in sorted(expected):
                    result = compare_array(actual[key], expected[key])
                    array_checks.append(dict(file=name, field=key, **result))
                    check(f'{name}/{key}', result['pass_check'])
                data[name] = actual
            guarded(name, load_one)

    load_archive('harmonic', reference/'harmonic-primary/calibration', harmonic/'calibration')
    for name in ('nonlinear-h0.5', 'nonlinear-h0.25'):
        load_archive(name, reference/'nonlinear-primary'/name, nonlinear/name)
    check('all-82-archived-arrays-compared', len(array_checks) == 82)

    def watchdog(stage, actual_dir, expected_dir):
        a, b = read_json(actual_dir/'watchdog.json'), read_json(expected_dir/'watchdog.json')
        for key in ('status', 'stage', 'workers', 'BLAS_threads', 'cap_seconds', 'memory_cap_bytes', 'protocol_sha256'):
            check(f'{stage}/watchdog/{key}', a[key] == b[key])
        check(f'{stage}/watchdog/complete', a['status'] == 'COMPLETE')
        scalar_check(f'{stage}/wall-cap', a['wall_seconds'], b['cap_seconds'])
        scalar_check(f'{stage}/aggregate-memory-cap', a['sampled_aggregate_peak_bytes'], b['memory_cap_bytes'])
        check(f'{stage}/jobs', [j['name'] for j in a['jobs']] == [j['name'] for j in b['jobs']])
        for job in a['jobs']:
            check(f'{stage}/{job["name"]}/exit', job['return_code'] == 0 and job['watchdog_stop'] is None)
        resources.append(dict(stage=stage, wall_seconds=a['wall_seconds'],
                              sampled_aggregate_peak_bytes=a['sampled_aggregate_peak_bytes'],
                              cap_seconds=b['cap_seconds'], memory_cap_bytes=b['memory_cap_bytes']))

    guarded('harmonic/watchdog', lambda: watchdog('harmonic', harmonic, reference/'harmonic-primary'))
    guarded('nonlinear/watchdog', lambda: watchdog('nonlinear', nonlinear, reference/'nonlinear-primary'))

    def receipt_identity(label, a, b):
        check(f'{label}/complete', a['status'] == 'COMPLETE')
        check(f'{label}/protocol', a['protocol_sha256'] == b['protocol_sha256'] == protocol_hash)
        check(f'{label}/source-hashes', a['source_hashes'] == b['source_hashes'])
        check(f'{label}/random-seed', a['numpy_random_seed'] == b['numpy_random_seed'])

    def harmonic_screens():
        a = read_json(harmonic/'calibration/receipt.json')
        b = read_json(reference/'harmonic-primary/calibration/receipt.json')
        receipt_identity('harmonic', a, b)
        check('harmonic/implementation-milestone', a['harmonic_implementation_milestone_pass'] is True)
        scalar_check('harmonic/worker-wall-cap', a['wall_seconds'], 60)
        scalar_check('harmonic/worker-memory-cap', a['worker_peak_rss_bytes'], 512_000_000)
        for field in ('exact_method_difference', 'tensor_full_difference'):
            check(f'harmonic/{field}/keys', set(a[field]) == set(COMPONENTS))
            for key in COMPONENTS:
                scalar_check(f'harmonic/{field}/{key}', a[field][key], 1e-10)
        scalar_check('harmonic/relative-scale-difference', a['relative_scale_difference'], 1e-10)
        check('harmonic/case-order', [c['name'] for c in a['cases']] == [c['name'] for c in b['cases']])
        analytic = data['harmonic/analytic.npz']
        for actual, expected in zip(a['cases'], b['cases'], strict=True):
            name = actual['name']
            for field in ('faces', 'spacing', 'equivalent_full_states', 'full_grid_allocated'):
                check(f'{name}/{field}', actual[field] == expected[field])
            error = differences(data[f'harmonic/{name}.npz'], analytic)
            passed = all(v <= PHYSICAL_SCREEN for v in error.values())
            check(f'{name}/screen', actual['all_components_within_0_001'] is passed and
                  passed == expected['all_components_within_0_001'])
            check(f'{name}/reported-error-keys', set(actual['errors_against_analytic']) == set(COMPONENTS))
            for key in COMPONENTS:
                scalar_check(f'{name}/reported-error/{key}', abs(actual['errors_against_analytic'][key]-error[key]), REPLAY_ATOL)
            screens.append(dict(kind='harmonic-case', case=name, errors=error, screen_pass=passed,
                                archived_screen_pass=expected['all_components_within_0_001']))
        identities = lambda receipt: [(c['first'], c['second'], c['kind']) for c in receipt['comparisons']]
        check('harmonic/comparison-order', identities(a) == identities(b))
        final_seen = set()
        for actual, expected in zip(a['comparisons'], b['comparisons'], strict=True):
            pair = (actual['first'], actual['second'], actual['kind'])
            error = differences(data[f'harmonic/{pair[1]}.npz'], data[f'harmonic/{pair[0]}.npz'])
            passed = all(v <= PHYSICAL_SCREEN for v in error.values())
            check(f'harmonic/pair/{pair}', actual['all_components_within_0_001'] is passed and
                  passed == expected['all_components_within_0_001'])
            for key in COMPONENTS:
                scalar_check(f'harmonic/pair-reported-error/{pair}/{key}', abs(actual['errors'][key]-error[key]), REPLAY_ATOL)
            if pair in FINAL_PAIRS:
                final_seen.add(pair)
                check(f'harmonic/final-pair/{pair}', passed)
            screens.append(dict(kind=pair[2], first=pair[0], second=pair[1], errors=error,
                                screen_pass=passed, archived_screen_pass=expected['all_components_within_0_001']))
        check('harmonic/all-four-final-pairs', final_seen == FINAL_PAIRS)
        final = next(c for c in a['cases'] if c['name'] == 'harmonic-A8-h0.015625')
        check('harmonic/final-analytic-screen', final['all_components_within_0_001'] is True)

    guarded('harmonic/receipt-and-screens', harmonic_screens)

    def nonlinear_checks():
        for name in ('nonlinear-h0.5', 'nonlinear-h0.25'):
            a = read_json(nonlinear/name/'receipt.json')
            b = read_json(reference/'nonlinear-primary'/name/'receipt.json')
            receipt_identity(name, a, b)
            for key in ('faces', 'spacing', 'shape', 'states', 'full_positive_coordinates'):
                check(f'{name}/{key}', a[key] == b[key])
            scalar_check(f'{name}/stationarity', a['stationarity_residual'], 1e-9)
            scalar_check(f'{name}/matrix-action', a['csr_stencil_scaled_action_error'], 1e-12)
            z = data[f'{name}/response.npz']
            error = max(float(np.max(np.abs(z[k]-z[k+'_chebyshev']))) for k in ('K', 'C'))
            scalar_check(f'{name}/independent-method', error, 1e-8)
            scalar_check(f'{name}/response-identity', float(np.max(np.abs(z['C']-(z['K']-z['G0'])))), 1e-14)
            check(f'{name}/recorded-independent-pass', a['independent_comparison_pass'] is True)
            check(f'{name}/all-three-independent-times', len(a['independent_chebyshev']) == 3)
            for index, c in enumerate(a['independent_chebyshev']):
                scalar_check(f'{name}/tail/{index}', c['maximum_normalized_covariance_tail_bound'], 1e-9)
        error = differences(data['nonlinear-h0.5/response.npz'], data['nonlinear-h0.25/response.npz'])
        pattern = {k: error[k] <= PHYSICAL_SCREEN for k in COMPONENTS}
        check('nonlinear/old-box-component-pattern', pattern == dict(G0=True, K=False, C=False))
        screens.append(dict(kind='nonlinear-old-box-refinement', errors=error,
                            component_screen_pass=pattern, overall_screen_pass=all(pattern.values()),
                            interpretation='Expected G0 PASS, delayed K FAIL and response C FAIL. This failure is preserved, not waived.'))

    guarded('nonlinear/receipt-and-screens', nonlinear_checks)
    failed = [c for c in checks if not c['pass_check']]
    finite_deltas = [c['maximum_absolute_difference'] for c in array_checks if 'maximum_absolute_difference' in c]
    result = dict(status='FAIL' if failed else 'PASS', numerical_replay_absolute_tolerance=REPLAY_ATOL,
                physical_screen=PHYSICAL_SCREEN, independent_solver_gate=1e-8,
                arrays_compared=len(array_checks), maximum_array_absolute_difference=max(finite_deltas, default=None),
                checks=checks, failed_checks=failed, array_checks=array_checks, scientific_screens=screens,
                resource_checks=resources, scope='All archived arrays and named screen outcomes are compared. Mode bases use the original sign convention; no gauge transformation is fitted. Timings are checked only against declared caps, not against archived durations. Expected scientific failures remain failures. No new propagation is computed.')

    def json_safe(value):
        if isinstance(value, dict):
            return {k: json_safe(v) for k, v in value.items()}
        if isinstance(value, list):
            return [json_safe(v) for v in value]
        if isinstance(value, float) and not np.isfinite(value):
            return repr(value)
        return value
    return json_safe(result)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--harmonic', type=Path, required=True, help='Fresh harmonic watchdog directory containing calibration/')
    parser.add_argument('--nonlinear', type=Path, required=True, help='Fresh nonlinear watchdog directory containing both cases')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    receipt = run(args.harmonic, args.nonlinear)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as stream:
        json.dump(receipt, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps({k: receipt[k] for k in ('status', 'arrays_compared', 'maximum_array_absolute_difference', 'failed_checks')}, indent=2))
    if receipt['status'] != 'PASS':
        raise SystemExit(1)
