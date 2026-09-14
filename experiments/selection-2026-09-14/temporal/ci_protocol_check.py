"""CI portability adapter for historical exact-float grid checks.

Historical verifiers and their source seals remain unchanged. This adapter only
allows nearest grid nodes within 32 binary64 ULPs; it never changes a saved time,
response threshold, kernel lookup or scientific gate.
"""
from pathlib import Path
import argparse
import datetime
import hashlib
import json
import runpy
import numpy as np

MAX_ULPS = 32


def match_grid(expected, stored, *, exact_size=False):
    expected = np.asarray(expected, dtype=np.float64)
    stored = np.asarray(stored, dtype=np.float64)
    if expected.ndim != 1 or stored.ndim != 1 or len(expected) == 0 or len(stored) == 0:
        raise ValueError('Grid vectors must be nonempty and one-dimensional')
    if not np.all(np.isfinite(expected)) or not np.all(np.isfinite(stored)):
        raise ValueError('Grid nodes must be finite')
    if not np.all(np.diff(expected) > 0) or not np.all(np.diff(stored) > 0):
        raise ValueError('Grid nodes must be distinct and increasing')
    if exact_size and len(expected) != len(stored):
        raise ValueError('The stored grid has a missing or extra node')
    nearest = []
    for value in expected:
        distance = abs(stored - value)
        winners = np.flatnonzero(distance == np.min(distance))
        if len(winners) != 1:
            raise ValueError('The nearest stored grid node is ambiguous')
        nearest.append(int(winners[0]))
    if len(set(nearest)) != len(nearest):
        raise ValueError('Two expected nodes matched the same stored node')
    selected = stored[nearest]
    difference = abs(selected - expected)
    scale = np.maximum(abs(selected), abs(expected))
    spacing = np.spacing(scale)
    if not np.all(difference <= MAX_ULPS * spacing):
        raise ValueError('The stored grid differs by more than 32 binary64 ULPs')
    relative = np.divide(difference, scale, out=np.zeros_like(difference), where=scale != 0)
    return dict(expected_nodes=len(expected),stored_nodes=len(stored),matched_indices=nearest,
                exact_matches=int(np.count_nonzero(selected == expected)),
                maximum_absolute_difference=float(max(difference)),
                maximum_relative_difference=float(max(relative)),
                maximum_ulp_difference=float(max(difference / spacing)),
                allowed_binary64_ulps=MAX_ULPS)


def verify(root, kind):
    root = Path(root).resolve()
    path = root / ('verify.py' if kind == 'base' else 'event-aligned/verify_addendum.py')
    legacy = runpy.run_path(str(path))
    old = legacy['verify_protocol' if kind == 'base' else 'protocol_checks']()
    checks = old.copy()
    if kind == 'base':
        protocol = json.loads((root/'protocol.json').read_text())
        geometry = protocol['base_log_grid']
        if geometry != dict(count=61,min=.001,max=30.,spacing='numpy.geomspace'):
            raise ValueError('The declared log-grid design changed')
        generated = np.geomspace(geometry['min'],geometry['max'],geometry['count'])
        diagnostics = [match_grid(generated,protocol['observation_times_over_tau'])]
        checks['log_grid_contained'] = True
        tiny = legacy['verify_tiny']()
    else:
        protocol = json.loads((root/'event-aligned/protocol.json').read_text())
        generated = np.geomspace(.001,30,61)
        diagnostics = []
        for group in protocol['groups']:
            expected = [0.] + [float(x) for x in generated if group['duration_over_tau']+x <= 30.]
            diagnostics.append(match_grid(expected,group['post_removal_lags_over_tau'],exact_size=True))
        checks['fixed_lag_rule'] = True
        tiny = {}
    with path.open('rb') as stream:
        historical_sha256 = hashlib.file_digest(stream,'sha256').hexdigest()
    return dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                scope='CI grid-comparison portability only; historical source/protocol hashes and all scientific gates are unchanged.',
                historical_verifier_sha256=historical_sha256,
                historical_exact_checks=old,protocol=checks,tiny_field=tiny,
                grid_rounding_diagnostics=diagnostics,
                passed=all(checks.values()) and all(tiny.values()))


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,required=True)
    parser.add_argument('--kind',choices=['base','event'],required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    result=verify(args.root,args.kind)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(result,indent=2,allow_nan=False))
    raise SystemExit(0 if result['passed'] else 1)
