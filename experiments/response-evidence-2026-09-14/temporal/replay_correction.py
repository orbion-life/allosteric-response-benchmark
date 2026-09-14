"""Replay the corrected analysis from original saved arrays into a new directory."""
from pathlib import Path
import argparse
import json
from analyze import ROOT, OLD, P, run, dump


def summarize(output):
    historical = json.loads((ROOT/'history/receiver-ranking-before-2026-09-15/summary.json').read_text())
    for row in historical['rows']:
        analysis = json.loads((output/'results'/row['target']/'analysis.json').read_text())
        groups = [analysis['results']['fresh_step']] + analysis['results']['fresh_removal']
        row.update(
            fresh_candidate_receiver_error=max(g['candidate_receiver_max_entry_error'] for g in groups),
            fresh_informative_queries=sum(g['non_low_signal_query_count'] for g in groups),
            fresh_informative_identical_top5_set=sum(g['non_low_signal_identical_top5_set_count'] for g in groups),
            fresh_informative_identical_top5_order=sum(g['non_low_signal_identical_top5_order_count'] for g in groups),
            fresh_informative_gap_above_0_004=sum(g['gap_above_0_004_query_count'] for g in groups),
            fresh_gap_above_0_004_identical_set=sum(g['gap_above_0_004_identical_top5_set_count'] for g in groups),
        )
    historical['fresh_informative_top5_order'] = sum(r['fresh_informative_identical_top5_order'] for r in historical['rows'])
    historical['fresh_informative_top5_set'] = sum(r['fresh_informative_identical_top5_set'] for r in historical['rows'])
    historical['fresh_informative_queries'] = sum(r['fresh_informative_queries'] for r in historical['rows'])
    historical['analysis_correction_date'] = '2026-09-15'
    historical['status'] = 'COMPLETED_DEVELOPMENT_EXPERIMENT_WITH_RECEIVER_INDEX_CORRECTION'
    historical['correction'] = 'Receiver values are integer column indices. Candidate values are a boolean mask. All derived ranking fields were regenerated; original matrices, operators, construction timings and numerical qualification are unchanged.'
    dump(output/'summary.json', historical)
    return historical


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-root', type=Path, default=ROOT)
    parser.add_argument('--reference-root', type=Path, default=OLD)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Output already exists; choose a new destination')
    args.output.mkdir(parents=True)
    for target in P['targets']:
        run(target, args.data_root, args.reference_root, args.output)
    summary = summarize(args.output)
    print(f"Corrected informative ordered top five: {summary['fresh_informative_top5_order']}/{summary['fresh_informative_queries']}")
