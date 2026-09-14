"""Correct the ancillary receiver-index peak diagnostic without changing evidence.

This post-execution correction reads archived or regenerated waveforms. It never
imports the original diagnostic function, modifies inputs, or recomputes kernels.
"""
from pathlib import Path
import argparse
import datetime
import hashlib
import json
import os
import time

for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[key] = "3"
import numpy as np

HERE = Path(__file__).resolve().parent
DEFAULT_ROOT = HERE.parent.parent
THRESHOLD = 0.002
TIE_TOLERANCE = 1e-12


def digest(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def corrected_pairs(reference, reduced, times, candidate, receiver):
    assert candidate.dtype == np.bool_
    assert np.issubdtype(receiver.dtype, np.integer)
    assert receiver.ndim == 1 and len(set(receiver.tolist())) == len(receiver)
    assert np.all((0 <= receiver) & (receiver < reference.shape[2]))
    records = []
    for source_index in np.flatnonzero(candidate):
        for receiver_index in receiver:
            a = np.abs(reference[:, source_index, receiver_index])
            b = np.abs(reduced[:, source_index, receiver_index])
            reference_index = int(np.argmax(a))
            reduced_index = int(np.argmax(b))
            height = float(a[reference_index])
            ties = int(np.count_nonzero(np.abs(a - height) <= TIE_TOLERANCE))
            unique_interior = height > THRESHOLD and 0 < reference_index < len(times) - 1 and ties == 1
            records.append([int(source_index), int(receiver_index), reference_index,
                            reduced_index, height, int(unique_interior),
                            float(abs(times[reference_index] - times[reduced_index]))])
    return np.array(records, dtype=np.float64).reshape(-1, 7)


def summary(records):
    eligible = records[:, 5].astype(bool)
    return {
        "candidate_receiver_pair_count": len(records),
        "unique_interior_reference_sampled_maxima_under_fixed_tie_rule": int(eligible.sum()),
        "low_signal_tied_or_boundary_reference_sampled_maxima": int((~eligible).sum()),
        "maximum_sampled_peak_time_difference_over_tau": float(records[eligible, 6].max()) if eligible.any() else None,
        "fraction_same_sampled_peak_index": float(np.mean(records[eligible, 2] == records[eligible, 3])) if eligible.any() else None,
    }


def fixture_check():
    # Nonconsecutive receivers include odd and even indices. Each true pair has
    # its reference maximum at index 1 and reduced maximum at index 2.
    candidate = np.array([True, False, True, False, True, False, False])
    receiver = np.array([1, 6], dtype=np.int64)
    full = np.zeros((4, 7, 7)); reduced = np.zeros_like(full)
    for i in [0, 2, 4]:
        for j in [1, 6]:
            full[1, i, j] = 1 + i + j
            reduced[2, i, j] = 1 + i + j
    times = np.array([0., 0.1, 0.4, 1.])
    rows = corrected_pairs(full, reduced, times, candidate, receiver)
    expected = {(i, j) for i in [0, 2, 4] for j in [1, 6]}
    assert {(int(i), int(j)) for i, j in rows[:, :2]} == expected
    assert np.all(rows[:, 2] == 1) and np.all(rows[:, 3] == 2)
    assert np.allclose(rows[:, 6], 0.3, atol=1e-15, rtol=0)
    historical_mask = candidate[:, None] & receiver[None, :]
    assert set(zip(*np.where(historical_mask))) != expected
    return {"status": "PASS", "true_pair_count": 6,
            "nonconsecutive_odd_and_even_receiver_indices": receiver.tolist(),
            "historical_selector_fails_fixture": True}


def main(root, output):
    root = Path(root).resolve(); output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=False)
    started = time.perf_counter(); fixture = fixture_check()
    qualification = json.loads((root / "candidate-qualification-audit.json").read_text())
    qualification_rows = {(r["grid"], r["target"], r["endpoint"]): r for r in qualification["rows"]}
    rows = []; raw = {}; verification = []
    for grid, base in [("original-grid", root), ("event-aligned", root / "event-aligned")]:
        for target in ["kras", "abl", "myc", "myh7"]:
            metadata_path = root / "inputs" / target / "ranking-metadata.npz"
            with np.load(metadata_path) as metadata:
                candidate = metadata["candidate"].copy(); receiver = metadata["receiver"].copy()
                ids = metadata["canonical"].copy()
            candidates = np.flatnonzero(candidate)
            original = json.loads((base / "results" / target / "analysis.json").read_text())
            endpoints = ([("step", original["step"])] if "step" in original else [])
            endpoints += [(f"pulse-{r['duration_over_tau']:g}", r) for r in original["pulses"]]
            for endpoint, metrics in endpoints:
                waveform_path = base / "results" / target / f"{endpoint}.npz"
                with np.load(waveform_path) as z:
                    reference = z["full"]; reduced = z["reduced"]; times = z["times_over_tau"]
                    records = corrected_pairs(reference, reduced, times, candidate, receiver)
                    # Independent extraction: np.take receives the true axis indices.
                    distal_reference = np.take(np.take(reference, candidates, axis=1), receiver, axis=2)
                    distal_reduced = np.take(np.take(reduced, candidates, axis=1), receiver, axis=2)
                    vectorized_peak_indices = np.argmax(abs(distal_reference), axis=0).ravel()
                    assert np.array_equal(vectorized_peak_indices, records[:, 2])
                    assert len(records) == len(candidates) * len(receiver)
                    assert np.array_equal(records[:, 1].reshape(len(candidates), -1)[0], receiver)
                    primary_error = np.max(abs(reference - reduced), axis=(1, 2))
                    sign_mask = abs(reference) > THRESHOLD
                    sign_wrong = int(np.count_nonzero(sign_mask & (np.sign(reference) != np.sign(reduced))))
                    # Recompute RMS rankings from actual receiver columns.
                    rank_arrays = []
                    for matrix in [reference, reduced]:
                        scores = np.sqrt(np.mean(np.take(matrix, receiver, axis=2) ** 2, axis=2))
                        rank_arrays.append(np.array([candidates[np.lexsort((ids[candidates], -s[candidates]))] for s in scores]))
                    distal_error = np.max(abs(distal_reference - distal_reduced), axis=(1, 2))
                    peak = float(np.max(abs(distal_reference)))
                    relative_frobenius = float(np.linalg.norm(distal_reference - distal_reduced) / np.linalg.norm(distal_reference)) if peak > THRESHOLD else None
                    q = qualification_rows[(grid, target, endpoint)]
                    checks = {
                        "all_matrix_error_by_time_unchanged": bool(np.allclose(primary_error, metrics["error_by_time"], atol=1e-12, rtol=0)),
                        "all_matrix_gate_unchanged": bool(np.all(primary_error <= THRESHOLD)) == metrics["all_times_pass"],
                        "all_matrix_sign_counts_unchanged": int(sign_mask.sum()) == metrics["sign_eligible_entries"] and sign_wrong == metrics["sign_disagreements"],
                        "reference_RMS_ranking_unchanged": bool(np.array_equal(rank_arrays[0], z["reference_order"])),
                        "reduced_RMS_ranking_unchanged": bool(np.array_equal(rank_arrays[1], z["order"])),
                        "candidate_receiver_errors_unchanged": bool(np.allclose(distal_error, q["candidate_to_receiver_error_by_time"], atol=1e-12, rtol=0)),
                        "candidate_reference_peak_unchanged": abs(peak - q["candidate_receiver_reference_peak_magnitude"]) <= 1e-12,
                        "candidate_relative_Frobenius_unchanged": (relative_frobenius is None and q["candidate_receiver_relative_Frobenius_error"] is None) or (relative_frobenius is not None and abs(relative_frobenius - q["candidate_receiver_relative_Frobenius_error"]) <= 1e-12),
                    }
                    assert all(checks.values()), (grid, target, endpoint, checks)
                key = f"{grid}__{target}__{endpoint}"
                raw[key] = records
                rows.append({"grid": grid, "target": target, "endpoint": endpoint,
                             "corrected": summary(records), "withdrawn_historical_peak_metrics": metrics["peak"]})
                verification.append({"grid": grid, "target": target, "endpoint": endpoint,
                                     "checks": checks, "metadata_sha256": digest(metadata_path),
                                     "waveform_sha256": digest(waveform_path),
                                     "maximum_matrix_error": float(primary_error.max()),
                                     "maximum_candidate_receiver_error": float(distal_error.max()),
                                     "candidate_receiver_reference_peak": peak,
                                     "candidate_receiver_relative_Frobenius_error": relative_frobenius})
                print("checked", grid, target, endpoint, flush=True)
    np.savez_compressed(output / "corrected-pair-peaks.npz", **raw)
    scope = "Post-execution correction of an ancillary indexing defect; fixed data, basis, time grid, threshold and 1e-12 tie rule. Unique interior sampled maxima are not robust or continuous-time timing certificates."
    result = {"created_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(), "scope": scope,
              "source_sha256": digest(Path(__file__)), "error_threshold": THRESHOLD,
              "peak_tie_tolerance": TIE_TOLERANCE, "endpoint_count": len(rows),
              "pair_array_columns": ["source_zero_based_index", "actual_receiver_zero_based_index", "reference_sample_index", "reduced_sample_index", "reference_peak_magnitude", "unique_interior_reference_sampled_maximum_flag", "sampled_peak_time_difference_over_tau"],
              "rows": rows}
    (output / "corrected-peak-summary.json").write_text(json.dumps(result, indent=2) + "\n")
    audit = {"created_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(), "status": "PASS",
             "scope": "Independent corrected-pair fixture and vectorized extraction; direct recheck that every response gate, sign count, RMS ranking and candidate qualification error remains unchanged.",
             "fixture": fixture, "all_endpoint_checks_pass": True, "endpoints": verification,
             "seconds": time.perf_counter() - started}
    (output / "independent-verification.json").write_text(json.dumps(audit, indent=2) + "\n")
    print(json.dumps({k: v for k, v in audit.items() if k != "endpoints"}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT, help="Temporal tree with archived or regenerated waveforms")
    parser.add_argument("--output", type=Path, required=True, help="New correction output directory")
    args = parser.parse_args()
    main(args.root, args.output)
