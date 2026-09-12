#!/usr/bin/env python3
"""Amend reference observability without refitting or reranking predictions.

Only the output directory is written. The released pilot and protocol are read-only.
The legacy-universe analyses are retained solely to explain the historical change.
"""
from pathlib import Path
from itertools import combinations, product
from math import comb
import argparse
import csv
import hashlib
import importlib.util
import json
import platform
import sys

# Import the released parser without creating cache files beside read-only inputs.
sys.dont_write_bytecode = True

import numpy as np
import scipy
from scipy.stats import hypergeom, rankdata

HERE = Path(__file__).resolve().parent
SEED = 20260912
DRAWS = 10000
PINNED_COMMIT = "527712c281f1a366b9fb1524342a8daf31a06d7a"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def write_csv(path, rows):
    with Path(path).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def classify_distance(distance):
    """A missing distance is unknown; it is never an observed noncontact."""
    if not np.isfinite(distance):
        return "unknown"
    if distance < 0:
        raise ValueError("A finite heavy-atom distance cannot be negative.")
    return "contact" if distance <= 5.0 else "observed_noncontact"


def evaluation_rows(rows, labels):
    """Return a new table and preserve all original prediction eligibility fields."""
    if len({r["canonical"] for r in rows}) != len(rows):
        raise ValueError("Duplicate canonical residue identifiers.")
    if set(labels) != {r["canonical"] for r in rows}:
        raise ValueError("Every input residue requires an explicit label state.")
    allowed = {"contact", "observed_noncontact", "unknown"}
    if not set(labels.values()) <= allowed:
        raise ValueError("Invalid label state.")
    return [dict(r, label_state=labels[r["canonical"]],
                 reference_observed=labels[r["canonical"]] != "unknown",
                 evaluation_eligible=bool(r["eligible_6A"] and labels[r["canonical"]] != "unknown"))
            for r in rows]


def top_summary(order, ids, labels):
    top = [int(ids[i]) for i in order[:5]]
    states = [labels[i] for i in top]
    hits = states.count("contact")
    unknown = states.count("unknown")
    # Unknown predictions are not dropped, and no lower-ranked replacements are inserted.
    return {"top5": top, "top5_label_states": states, "hits_at_5A": hits,
            "top5_unknown_count": unknown,
            "precision_at_5": hits / 5 if len(top) == 5 and not unknown else None,
            "precision_at_5_lower_bound": hits / 5,
            "precision_at_5_upper_bound": (hits + unknown) / 5}


def matched_plan(rows, positive, universe_field, exposure, fixed_cuts=None):
    universe = [r for r in rows if r[universe_field]]
    ids = [r["canonical"] for r in universe]
    if len(ids) != len(set(ids)) or not universe:
        raise ValueError("The evaluation universe must contain distinct residues.")
    if not positive or not set(positive) <= set(ids):
        raise ValueError("At least one positive must lie inside the evaluation universe.")
    features = ["degree", "receiver_distance_A"]
    cuts = {key: [float(x) for x in np.quantile([r[key] for r in universe], [1/3, 2/3])]
            for key in features}
    # The published exposure-matched protocol removes repeated cut values.
    if exposure:
        cuts = {key: sorted(set(values)) for key, values in cuts.items()}
    if fixed_cuts is not None:
        cuts = fixed_cuts
    groups = {}
    for row in universe:
        key = tuple(int(np.searchsorted(cuts[k], row[k], side="right")) for k in features)
        if exposure:
            key += (int(row["rsa"] >= 0.20),)
        groups.setdefault(key, []).append(int(row["canonical"]))
    strata = []
    space = 1
    for key, candidates in sorted(groups.items()):
        m = sum(i in positive for i in candidates)
        space *= comb(len(candidates), m)
        strata.append({"key": list(key), "candidates": candidates, "n": len(candidates), "m": m})
    return {"candidate_count": len(universe), "positive_count": len(positive),
            "cut_policy": "historical126 cuts retained as sensitivity" if fixed_cuts is not None else "tertiles recomputed within this universe",
            "universe_field": universe_field, "matching_features": features + (["rsa_ge_0.20"] if exposure else []),
            "cuts": cuts, "groups": strata, "M": space, "smallest_exact_p": 1 / space,
            "first_Holm_threshold_in_planned_three_target_family": 0.05 / 3,
            "finite_space_supports_that_threshold": space >= 60}


def configurations(plan, positive, distinct):
    """Legacy sampling is preserved exactly; the newer protocol samples distinct sets."""
    groups = plan["groups"]
    if distinct and plan["M"] <= 100000:
        values = [tuple(sorted(i for block in blocks for i in block))
                  for blocks in product(*(combinations(g["candidates"], g["m"]) for g in groups))]
        return np.asarray(values, dtype=int), True
    rng = np.random.Generator(np.random.PCG64(SEED))
    seen = {tuple(sorted(positive))} if distinct else set()
    values = []
    while len(values) < DRAWS:
        chosen = []
        for group in groups:
            if group["m"]:
                chosen.extend(int(i) for i in rng.choice(group["candidates"], group["m"], replace=False))
        # Legacy summation order is retained for exact arithmetic and receipt parity.
        draw = tuple(sorted(chosen)) if distinct else tuple(chosen)
        if distinct and draw in seen:
            continue
        seen.add(draw)
        values.append(draw)
    return np.asarray(values, dtype=int), False


def statistic(scores, ids, positive, draws, exact):
    position = {int(i): k for k, i in enumerate(ids)}
    used = set(positive) | set(draws.ravel().tolist())
    if not used <= position.keys() or not all(np.isfinite(scores[position[i]]) for i in used):
        raise ValueError("Every evaluated residue requires a finite, fixed score.")
    observed = float(np.mean([scores[position[i]] for i in sorted(positive)]))
    indices = np.array([[position[int(i)] for i in row] for row in draws])
    null = np.mean(scores[indices], axis=1)
    extreme = int(np.sum(null >= observed))
    p = extreme / len(null) if exact else (1 + extreme) / (1 + len(null))
    mean = float(null.mean())
    return {"observed_mean": observed, "null_mean": mean,
            "observed_minus_null_mean": observed - mean,
            "observed_over_null_mean": observed / mean if mean else None,
            "null_q025": float(np.quantile(null, .025)), "null_median": float(np.median(null)),
            "null_q975": float(np.quantile(null, .975)), "extreme_count": extreme,
            "raw_p": float(p), "draw_count": len(null), "exact_enumeration": exact,
            "attainable_sampled_p_floor": 1 / len(null) if exact else 1 / (1 + len(null))}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--kras-root", type=Path, default=HERE.parent,
                        help="Directory containing the unchanged pilot/ and protocol/ folders.")
    parser.add_argument("--out", type=Path, default=HERE / "results")
    args = parser.parse_args(argv)
    kras = args.kras_root.resolve()
    pilot = kras / "pilot"
    protocol = kras / "protocol"
    out = args.out.resolve()
    if out == kras or kras in out.parents and (pilot == out or pilot in out.parents or protocol == out or protocol in out.parents):
        raise ValueError("Results cannot be written inside the frozen pilot or protocol.")
    out.mkdir(parents=True, exist_ok=True)
    expected = json.loads((HERE / "source-manifest.json").read_text())
    for rel, digest in expected["files"].items():
        if sha(kras / rel) != digest:
            raise ValueError("Released source checksum mismatch: " + rel)
    freeze = json.loads((pilot / "prediction-freeze.json").read_text())
    for rel, digest in freeze["files"].items():
        if sha(pilot / rel) != digest:
            raise ValueError("Frozen prediction checksum mismatch: " + rel)

    model = np.load(pilot / "model/static-model.npz")
    ids = model["canonical"]
    candidate = model["candidate"]
    idx = np.flatnonzero(candidate)
    rec = model["receiver"]
    input_rows = json.loads((protocol / "prepared/KRAS-input.json").read_text())["residues"]
    assert [r["canonical"] for r in input_rows] == ids.tolist()
    assert [r["canonical"] for r in input_rows if r["eligible_6A"]] == ids[idx].tolist()
    np.testing.assert_allclose([r["degree"] for r in input_rows], model["degree"], rtol=0, atol=0)
    # The protocol's coordinate parser and the original pilot retain slightly
    # different floating-point precision. Each historical null keeps its source.
    distance_parser_difference = float(np.max(np.abs(np.array([r["receiver_distance_A"] for r in input_rows]) - model["receiver_distance"])))
    np.testing.assert_allclose([r["receiver_distance_A"] for r in input_rows], model["receiver_distance"], rtol=0, atol=3e-6)

    mapping = json.loads((pilot / "raw/6OIM-mapping.json").read_text())["6oim"]["UniProt"]["P01116"]["mappings"]
    mapping = [m for m in mapping if m["chain_id"] == "A"]
    assert len(mapping) == 1
    assert mapping[0]["unp_start"] == mapping[0]["start"]["author_residue_number"] == 1
    assert mapping[0]["unp_end"] >= int(ids.max())
    prepare = load_module(pilot / "prepare.py", "released_pilot_prepare")
    atoms = prepare.atoms(pilot / "raw/6OIM.pdb")
    ligand = np.asarray([a["xyz"] for a in atoms if a["record"] == "HETATM" and a["resname"] == "MOV"])
    if not len(ligand):
        raise ValueError("The reference ligand has no heavy atoms.")
    labels, reference_rows = {}, []
    for residue, is_candidate in zip(ids, candidate):
        heavy = np.asarray([a["xyz"] for a in atoms if a["record"] == "ATOM" and a["auth_seq_id"] == residue])
        distance = float(np.sqrt(np.sum((heavy[:, None] - ligand[None]) ** 2, axis=-1)).min()) if len(heavy) else np.nan
        label = classify_distance(distance)
        labels[int(residue)] = label
        reference_rows.append({"canonical": int(residue), "prediction_candidate": bool(is_candidate),
                               "reference_heavy_atom_count": len(heavy), "reference_observed": bool(len(heavy)),
                               "MOV_distance_A": distance if np.isfinite(distance) else None,
                               "label_state": label, "contact_5A": None if label == "unknown" else label == "contact",
                               "evaluation_eligible": bool(is_candidate and len(heavy))})
    rows = evaluation_rows(input_rows, labels)
    unknown = sorted(i for i, state in labels.items() if state == "unknown")
    assert unknown == [105, 106, 107]
    # Independent direct inspection of deposit REMARK records corroborates no-atom labels.
    remarks = [line.rstrip() for line in (pilot / "raw/6OIM.pdb").read_text().splitlines()
               if line.startswith("REMARK 465") and line[19:20] == "A" and line[21:26].strip() in {str(i) for i in unknown}]
    assert len(remarks) == 3, "Missing-residue REMARK cross-check failed."
    positive = sorted(r["canonical"] for r in rows if r["evaluation_eligible"] and r["label_state"] == "contact")
    historical_labels = json.loads((protocol / "evaluation/KRAS-contacts.json").read_text())["eligible_canonical"]
    assert positive == historical_labels
    primary = np.load(pilot / "results/biquadratic-d2-n33-e4.npz")
    harmonic = np.load(pilot / "results/harmonic-d2-n33-e4.npz")
    hookean = np.load(pilot / "results/distance_hookean-d2-n33-e4.npz")
    full = np.load(pilot / "results/analytic-harmonic-d492.npz")
    structural = np.load(pilot / "results/structural-baselines.npz")
    rms = lambda matrix: np.sqrt(np.mean(matrix[:, rec] ** 2, axis=1))
    np.testing.assert_allclose(rms(primary["C"]), primary["score"], rtol=0, atol=0)
    methods = {"Biquadratic d2 grid": primary["score"], "Harmonic d2 same grid": harmonic["score"],
               "Distance-Hookean d2 same grid": hookean["score"], "Harmonic all 492 modes": full["score"],
               "Equilibrium d2 biquadratic": rms(primary["equilibrium"]),
               "Equilibrium all-mode harmonic": rms(full["equilibrium"]),
               "Graph heat kernel": structural["graph_heat_kernel"], "Contact degree": structural["degree"],
               "Receiver inverse distance": structural["receiver_inverse_distance"]}
    ohm = json.loads((pilot / "external/ohm-scores.json").read_text())
    methods["Ohm pinned upstream"] = np.array([ohm["score_by_canonical"][str(int(i))] for i in ids])
    historical = json.loads((pilot / "evaluation/evaluation.json").read_text())
    historical_methods = {r["method"]: r for r in historical["methods"]}
    rank_rows, top_results = [], []
    for method, score in methods.items():
        if not np.all(np.isfinite(score[idx])):
            raise ValueError("All prediction-candidate scores must be finite.")
        order = idx[np.lexsort((ids[idx], -score[idx]))]
        summary = {"method": method, **top_summary(order, ids, labels)}
        assert summary["top5"] == historical_methods[method]["top5"]
        assert summary["precision_at_5"] == historical_methods[method]["precision_at_5"]
        top_results.append(summary)
        for rank, i in enumerate(order, 1):
            rank_rows.append({"method": method, "prediction_rank": rank, "canonical": int(ids[i]),
                              "score": float(score[i]), "label_state": labels[int(ids[i])],
                              "evaluation_eligible": labels[int(ids[i])] != "unknown"})

    paired_scores = {}
    for method in ["Harmonic d2 same grid", "Distance-Hookean d2 same grid"]:
        delta = np.full(len(ids), np.nan)
        delta[idx] = (rankdata(primary["score"][idx], method="average") - rankdata(methods[method][idx], method="average")) / len(idx)
        paired_scores["Biquadratic minus " + method] = delta

    results, paired_results, plans, draws_to_save = [], [], {}, {}
    parity = {}
    for exposure in [False, True]:
        family = "degree_distance_exposure" if exposure else "legacy_degree_distance"
        for universe_name in ["historical126", "evaluable123", "evaluable123_fixed_historical_cuts"]:
            corrected = universe_name != "historical126"
            name = family + "__" + universe_name
            field = "evaluation_eligible" if corrected else "eligible_6A"
            family_rows = rows if exposure else [dict(r, receiver_distance_A=float(model["receiver_distance"][i])) for i, r in enumerate(rows)]
            fixed_cuts = plans[family + "__historical126"]["cuts"] if universe_name.endswith("fixed_historical_cuts") else None
            plan = matched_plan(family_rows, set(positive), field, exposure, fixed_cuts)
            draws, exact = configurations(plan, positive, distinct=exposure)
            if corrected:
                assert not set(unknown) & set(draws.ravel().tolist())
            plan.update({"seed": SEED, "draw_count": len(draws), "distinct_sampling": exposure,
                         "unique_configurations_evaluated": len({tuple(sorted(x)) for x in draws.tolist()}),
                         "observed_configuration_excluded_from_sample": exposure and not exact,
                         "exact_enumeration": exact, "attainable_sampled_p_floor": 1 / len(draws) if exact else 1 / (len(draws) + 1)})
            plans[name] = plan
            draws_to_save[name] = draws
            for method, score in methods.items():
                row = {"null_family": family, "universe": universe_name, "candidate_count": plan["candidate_count"],
                       "method": method, **statistic(score, ids, positive, draws, exact)}
                results.append(row)
                if not exposure and not corrected:
                    old = historical_methods[method]
                    np.testing.assert_allclose(row["null_mean"], old["matched_random_mean_score"], rtol=0, atol=0)
                    assert row["raw_p"] == old["matched_random_fraction_at_least_known_site"]
            for method, score in paired_scores.items():
                paired_statistic = statistic(score, ids, positive, draws, exact)
                paired_statistic.pop("observed_over_null_mean")
                paired_results.append({"null_family": family, "universe": universe_name, "candidate_count": plan["candidate_count"],
                                       "comparison": method, "percentile_rank_denominator": len(idx),
                                       **paired_statistic})
            if not exposure and not corrected:
                old_draws = np.load(pilot / "evaluation/reference-mapping.npz")["matched_random_nodes"]
                assert np.array_equal(ids[old_draws], draws)
                parity["legacy_null_draws_exactly_reproduced"] = True
            if exposure and fixed_cuts is None:
                # Compare our reusable implementation directly to the released protocol.
                published = load_module(protocol / "matched_null.py", "released_matched_null")
                pp_rows = [dict(r, eligible_6A=r[field]) for r in rows]
                pp_plan = published.plan(pp_rows, set(positive))
                pp_draws = np.array(list(published.configurations(pp_plan, positive, SEED)))
                assert plan["M"] == pp_plan["M"]
                assert np.array_equal(draws, pp_draws)
                parity[name + "_published_protocol_exact_draw_parity"] = True

    uniform = []
    fixed_cut_sensitivity = {}
    for family in ["legacy_degree_distance", "degree_distance_exposure"]:
        original_key = family + "__historical126"
        fixed_key = family + "__evaluable123_fixed_historical_cuts"
        fixed_cut_sensitivity[family] = {
            "historical_unknown_strata": [{"key": g["key"], "n": g["n"], "m": g["m"]}
                                           for g in plans[original_key]["groups"] if set(g["candidates"]) & set(unknown)],
            "unknown_residues_in_historical_sample": sorted(set(draws_to_save[original_key].ravel().tolist()) & set(unknown)),
            "fixed_cut_configs_equal_historical_configs": bool(np.array_equal(draws_to_save[original_key], draws_to_save[fixed_key])),
            "interpretation": "The unknown residues occupy a historical stratum with zero positives. Retaining historical cuts gives identical samples after excluding unknowns; amended re-tertiled statistics change through reassignment of observed residues to strata."}
    for count, name in [(len(idx), "historical126"), (sum(r["evaluation_eligible"] for r in rows), "evaluable123")]:
        p = len(positive)
        uniform.append({"universe": name, "candidate_count": count, "positive_count": p,
                        "draw_size": 5, "contact_fraction": p/count, "expected_hits_in_five": 5*p/count,
                        "probability_at_least_one_hit": float(hypergeom.sf(0, count, p, 5)),
                        "probability_at_least_four_hits": float(hypergeom.sf(3, count, p, 5)),
                        "precision_enrichment_for_fixed_four_of_five": .8/(p/count),
                        "interpretation": "Descriptive uniform draws without replacement; no degree, distance or exposure matching."})
    write_csv(out / "reference-labels.csv", reference_rows)
    write_json(out / "reference-labels.json", reference_rows)
    write_csv(out / "unchanged-rankings.csv", rank_rows)
    write_json(out / "top5.json", top_results)
    write_csv(out / "null-comparisons.csv", results)
    write_json(out / "null-comparisons.json", results)
    write_csv(out / "paired-comparisons.csv", paired_results)
    write_json(out / "paired-comparisons.json", paired_results)
    write_json(out / "null-plans.json", plans)
    write_json(out / "uniform-hit-baselines.json", uniform)
    np.savez_compressed(out / "matched-configurations.npz", **draws_to_save)
    write_json(out / "evaluation-input.json", {"residues": rows, "positive_canonical": positive})
    receipt = {"amendment": "KRAS reference observability correction", "date": "2026-09-12",
               "historical_release": "v0.2.0", "historical_commit": PINNED_COMMIT,
               "source_manifest_sha256": sha(HERE / "source-manifest.json"),
               "script_sha256": sha(Path(__file__)), "frozen_prediction_files_verified": len(freeze["files"]),
               "prediction_candidate_count": len(idx), "evaluable_candidate_count": sum(r["evaluation_eligible"] for r in rows),
               "unknown_reference_canonical": unknown, "reference_missing_residue_remarks": remarks,
               "positive_count": len(positive), "all_ten_top5_lists_unchanged": True,
               "original_prediction_ranks_retained_for_paired_controls": True,
               "maximum_input_distance_parser_difference_A": distance_parser_difference,
               "matching_distance_provenance": "Legacy null uses pilot/model/static-model.npz; exposure-matched null uses protocol/prepared/KRAS-input.json. Each original numerical convention is preserved.",
               "paired_percentile_rank_denominator": len(idx), "parity": parity,
               "fixed_cut_sensitivity": fixed_cut_sensitivity,
               "reference_contact_definition": "At least one deposited protein heavy atom lies within 5 angstrom of a deposited MOV heavy atom.",
               "observed_noncontact_definition": "At least one protein heavy atom is deposited, but none of the deposited atoms is within 5 angstrom of MOV; this is not a biological negative.",
               "limitations": ["The correction is retrospective and changes no fitted model.",
                               "Absence of deposited atoms is not evidence of no contact.",
                               "Partial side-chain disorder can still hide contacts at residues with some deposited atoms.",
                               "The comparison uses ligand-conditioned structures with different sequences and crystal states.",
                               "Conditional null sets are exchangeability assumptions, not experimental negative pockets or independent biological replicates.",
                               "The sampled p-value floor is resolution, not proof that the exact p-value reaches that floor.",
                               "Paired raw p-values do not complete the reserved six-test comparative family.",
                               "No nonlinear benefit, compression validity, noise guarantee, or prospective generalization follows from this correction."],
               "runtime": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__, "platform": platform.system()},
               "output_sha256": {p.name: sha(p) for p in sorted(out.iterdir()) if p.is_file() and p.name != "verification.json"}}
    write_json(out / "verification.json", receipt)
    print(json.dumps({"prediction_candidates": len(idx), "evaluable_candidates": receipt["evaluable_candidate_count"],
                      "unknown": unknown, "positive_count": len(positive), "null_plans": {k:v["M"] for k,v in plans.items()}}))


if __name__ == "__main__":
    main()
