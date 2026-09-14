"""Independently verify receiver scoring and reduced responses from full saved arrays.

This verifier does not import the production selector or scoring implementation.
It diagonalizes the saved reduced generator again, selects receiver columns with
np.take and ranks with Python's tuple sort. It never fits a Gaussian model.
"""
import os
for key in ['OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS']:
    os.environ[key] = '3'
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np
from scipy.linalg import eigh
from scipy.linalg.blas import dgemm

ROOT = Path(__file__).resolve().parent
TARGETS = ['kras', 'abl', 'myc', 'myh7']

def digest(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()

def groups(results):
    for name, value in results.items():
        for group in value if isinstance(value, list) else [value]:
            yield name, group

def rank_scores(matrix, receivers, candidates, ids):
    scores = np.linalg.norm(np.take(matrix, receivers, axis=1), axis=1) / np.sqrt(len(receivers))
    order = sorted(candidates, key=lambda index: (-scores[index], ids[index]))
    return scores, order

def verify_target(target, data_root, reference_root, analysis_root):
    directory = data_root/'results'/target
    metadata_path = data_root/'inputs'/target/'ranking-metadata.npz'
    with np.load(metadata_path, allow_pickle=False) as metadata:
        raw = metadata['receiver']; candidate = metadata['candidate']; ids = metadata['canonical']
    # Check the actual stored types before any conversion or numerical work.
    assert raw.dtype.kind in 'iu' and raw.ndim == 1 and len(raw) > 0
    assert len(set(raw.tolist())) == len(raw) and min(raw) >= 0 and max(raw) < len(ids)
    assert candidate.dtype.kind == 'b' and candidate.shape == ids.shape
    receivers = raw.tolist()
    candidates = [i for i, included in enumerate(candidate) if included]
    report = json.loads((analysis_root/'results'/target/'analysis.json').read_text())
    assert report['receiver_semantics']['receiver_indices'] == receivers
    assert report['receiver_semantics']['metadata_sha256'] == digest(metadata_path)
    with np.load(directory/'operator.npz', allow_pickle=False) as operator:
        eigenvalues, eigenvectors = eigh(operator['Hr'], driver='evd')
        projections = dgemm(1., eigenvectors.T, operator['B'])
        tau = float(operator['tau'])
    protocol = json.loads((ROOT/'protocol.json').read_text())
    old_protocol = json.loads((reference_root/'protocol.json').read_text())
    event_protocol = json.loads((reference_root/'event-aligned/protocol.json').read_text())
    paths = [directory/'fresh-kernels.npy', reference_root/'results'/target/'full-kernels.npy', reference_root/'event-aligned/results'/target/'new-full-kernels.npy']
    arrays = [np.load(path, mmap_mode='r') for path in paths]
    clocks = [protocol['fresh_kernel_times_over_tau'], old_protocol['all_exact_kernel_times_over_tau'], event_protocol['new_exact_kernel_times_over_tau']]
    locations = [{float(time): index for index, time in enumerate(clock)} for clock in clocks]
    static = arrays[1][locations[1][0.0]]
    def full(time):
        for values, locations_one in zip(arrays, locations):
            if time in locations_one:
                return values[locations_one[time]] - static
        raise KeyError(time)
    cache = {}
    def reduced(time):
        if time not in cache:
            cache[time] = dgemm(1., projections.T * np.expm1(-time * tau * eigenvalues), projections)
        return cache[time]
    historical = json.loads((ROOT/'history/receiver-ranking-before-2026-09-15/results'/target/'analysis.json').read_text())
    total = fresh_total = informative = matching_order = matching_set = 0
    max_error_discrepancy = max_score_discrepancy = max_matrix_entry_error = 0.0
    mismatches = []
    fixture = None
    for (name, group), (old_name, old_group) in zip(groups(report['results']), groups(historical['results'])):
        assert name == old_name and group['query_count'] == old_group['query_count']
        for row, old_row in zip(group['rows'], old_group['rows']):
            time = row['time_over_tau']; lag = row.get('post_removal_lag_over_tau')
            f = full(time); r = reduced(time)
            if lag is not None:
                f = f - full(lag); r = r - reduced(lag)
            fs, fo = rank_scores(f, receivers, candidates, ids)
            rs, ro = rank_scores(r, receivers, candidates, ids)
            maxerr = float(np.max(abs(r-f)))
            score_error = float(np.max(abs(fs[candidates]-rs[candidates])))
            score_delta = max(abs(float(max(fs[candidates]))-row['reference_max_candidate_score']), abs(score_error-row['candidate_score_max_abs']))
            max_score_discrepancy = max(max_score_discrepancy, score_delta)
            max_error_discrepancy = max(max_error_discrepancy, abs(maxerr-row['max_entry_error']))
            assert abs(maxerr-row['max_entry_error']) < 1e-9
            assert score_delta < 1e-9
            assert abs(float(np.max(abs(r-f)[np.ix_(candidates,receivers)]))-row['candidate_receiver_max_entry_error']) < 1e-9
            # Unaffected full-matrix metrics retain their historical values.
            assert row['max_entry_error'] == old_row['max_entry_error']
            assert row['time_over_tau'] == old_row['time_over_tau']
            above = float(max(fs[candidates])) > .002
            assert above == (not row['low_signal'])
            # Exact order checks exclude low-signal numerical ties.
            if above:
                assert ids[fo[:5]].tolist() == row['reference_top5']
                assert ids[ro[:5]].tolist() == row['reduced_top5']
                assert (fo[:5] == ro[:5]) == row['top5_order_identical']
                assert len(set(fo[:5]) & set(ro[:5])) / 5 == row['top5_set_agreement']
            total += 1
            if name.startswith('fresh'):
                fresh_total += 1
                max_matrix_entry_error = max(max_matrix_entry_error, maxerr)
                if above:
                    informative += 1; matching_order += int(fo[:5] == ro[:5]); matching_set += int(set(fo[:5]) == set(ro[:5]))
                    if fo[:5] != ro[:5]:
                        mismatches.append(dict(group=name, duration_over_tau=group.get('duration_over_tau'), **row))
            # Preserve a compact actual-data receiver scoring fixture for CI.
            if name == 'selected_times' and time == 1.0:
                fixture = dict(reference=f, reduced=r, receiver=np.array(receivers), candidate=candidate, canonical=ids)
    expected = {'kras': (47, 47), 'abl': (47, 46), 'myc': (36, 36), 'myh7': (46, 46)}[target]
    assert (informative, matching_order) == expected
    assert fresh_total == 72 and max_matrix_entry_error <= .002
    hashes = {'metadata':digest(metadata_path), 'operator':digest(directory/'operator.npz'), 'fresh_kernels':digest(paths[0]), 'old_kernels':digest(paths[1]), 'event_kernels':digest(paths[2])}
    return dict(target=target,status='PASS',queries_verified=total,fresh_queries_verified=fresh_total,informative_queries=informative,identical_top5_order=matching_order,identical_top5_set=matching_set,maximum_fresh_matrix_error=max_matrix_entry_error,maximum_matrix_error_discrepancy=max_error_discrepancy,maximum_score_discrepancy=max_score_discrepancy,receiver_indices=receivers,receiver_canonical=ids[receivers].tolist(),mismatches=mismatches,source_hashes=hashes), fixture

if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-root',type=Path,default=ROOT)
    parser.add_argument('--reference-root',type=Path,default=ROOT/'reference')
    parser.add_argument('--analysis-root',type=Path,default=ROOT)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--write-fixtures',type=Path)
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=False)
    rows=[]
    for target in TARGETS:
        row, fixture=verify_target(target,args.data_root,args.reference_root,args.analysis_root)
        rows.append(row)
        (args.output/f'{target}.json').write_text(json.dumps(row,indent=2)+'\n')
        if args.write_fixtures:
            args.write_fixtures.mkdir(parents=True,exist_ok=True)
            np.savez_compressed(args.write_fixtures/f'{target}-tau.npz',**fixture)
        print(target,'verified',row['queries_verified'],'queries',flush=True)
    result=dict(status='PASS',scope='Independent saved-generator eigendecomposition and explicit receiver-column scores. All full-matrix error values checked; exact shortlist comparison applies above the signal floor.',targets=rows,queries_verified=sum(r['queries_verified'] for r in rows),fresh_queries_verified=sum(r['fresh_queries_verified'] for r in rows),informative_queries=sum(r['informative_queries'] for r in rows),identical_top5_order=sum(r['identical_top5_order'] for r in rows),identical_top5_set=sum(r['identical_top5_set'] for r in rows),verifier_sha256=digest(__file__))
    (args.output/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
