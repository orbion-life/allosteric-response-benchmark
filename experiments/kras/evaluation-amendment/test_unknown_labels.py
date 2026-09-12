"""Regression tests for missing-data semantics and immutable prediction order."""
from pathlib import Path
from itertools import combinations
import copy
import json
import unittest

import numpy as np

import amend_evaluation as amendment

RESULTS = Path(__file__).resolve().parent / "results"


class MissingDataSemantics(unittest.TestCase):
    def test_missing_distance_is_explicitly_unknown(self):
        self.assertEqual(amendment.classify_distance(np.nan), "unknown")

    def test_geometric_contact_threshold_is_inclusive(self):
        self.assertEqual(amendment.classify_distance(5.0), "contact")
        self.assertEqual(amendment.classify_distance(5.000001), "observed_noncontact")

    def test_negative_distance_is_rejected(self):
        with self.assertRaises(ValueError):
            amendment.classify_distance(-0.1)

    def test_evaluation_mask_does_not_mutate_prediction_mask(self):
        rows = [{"canonical": 1, "eligible_6A": True}, {"canonical": 2, "eligible_6A": True}]
        saved = copy.deepcopy(rows)
        new = amendment.evaluation_rows(rows, {1: "contact", 2: "unknown"})
        self.assertEqual(rows, saved)
        self.assertTrue(new[1]["eligible_6A"])
        self.assertFalse(new[1]["evaluation_eligible"])
        self.assertEqual(new[1]["label_state"], "unknown")

    def test_incomplete_or_implicit_boolean_labels_are_rejected(self):
        rows = [{"canonical": 1, "eligible_6A": True}]
        for labels in [{}, {1: False}, {1: "negative"}]:
            with self.assertRaises(ValueError):
                amendment.evaluation_rows(rows, labels)

    def test_duplicate_input_residues_are_rejected(self):
        with self.assertRaises(ValueError):
            amendment.evaluation_rows([{"canonical": 1}, {"canonical": 1}], {1: "unknown"})

    def test_unknown_top_prediction_is_not_a_false_miss_or_replaced(self):
        ids = np.arange(1, 7)
        order = np.arange(6)
        labels = {1: "unknown", 2: "contact", 3: "contact", 4: "observed_noncontact", 5: "contact", 6: "contact"}
        result = amendment.top_summary(order, ids, labels)
        self.assertEqual(result["top5"], [1, 2, 3, 4, 5])
        self.assertEqual(result["hits_at_5A"], 3)
        self.assertEqual(result["top5_unknown_count"], 1)
        self.assertIsNone(result["precision_at_5"])
        self.assertEqual(result["precision_at_5_lower_bound"], .6)
        self.assertEqual(result["precision_at_5_upper_bound"], .8)

    def test_unknowns_cannot_enter_any_corrected_null_pool(self):
        plans = json.loads((RESULTS / "null-plans.json").read_text())
        draws = np.load(RESULTS / "matched-configurations.npz")
        for name, plan in plans.items():
            if "evaluable123" not in name:
                continue
            self.assertEqual(plan["candidate_count"], 123)
            pool = {i for g in plan["groups"] for i in g["candidates"]}
            self.assertFalse(pool & {105, 106, 107})
            self.assertFalse(set(draws[name].ravel()) & {105, 106, 107})

    def test_reference_missing_residues_have_null_contact_field(self):
        rows = json.loads((RESULTS / "reference-labels.json").read_text())
        missing = [r for r in rows if r["label_state"] == "unknown"]
        self.assertEqual([r["canonical"] for r in missing], [105, 106, 107])
        for row in missing:
            self.assertIsNone(row["contact_5A"])
            self.assertIsNone(row["MOV_distance_A"])
            self.assertEqual(row["reference_heavy_atom_count"], 0)
            self.assertTrue(row["prediction_candidate"])
            self.assertFalse(row["evaluation_eligible"])

    def test_ten_original_top_fives_and_precision_are_preserved(self):
        new = json.loads((RESULTS / "top5.json").read_text())
        old = {r["method"]: r for r in json.loads((RESULTS.parent / "historical/evaluation-v0.2.0.json").read_text())["methods"]}
        self.assertEqual(len(new), 10)
        for row in new:
            self.assertEqual(row["top5"], old[row["method"]]["top5"])
            self.assertEqual(row["precision_at_5"], old[row["method"]]["precision_at_5"])

    def test_corrected_paired_ranks_use_frozen_denominator(self):
        rows = json.loads((RESULTS / "paired-comparisons.json").read_text())
        self.assertTrue(all(r["percentile_rank_denominator"] == 126 for r in rows))
        for comparison in {r["comparison"] for r in rows}:
            self.assertEqual(len({r["observed_mean"] for r in rows if r["comparison"] == comparison}), 1)

    def test_exposure_null_configs_are_distinct_and_exclude_observed(self):
        plans = json.loads((RESULTS / "null-plans.json").read_text())
        positive = tuple(json.loads((RESULTS / "evaluation-input.json").read_text())["positive_canonical"])
        draws = np.load(RESULTS / "matched-configurations.npz")
        for name, plan in plans.items():
            if not name.startswith("degree_distance_exposure"):
                continue
            rows = {tuple(r) for r in draws[name]}
            self.assertEqual(len(rows), 10000)
            self.assertNotIn(positive, rows)
            for row in draws[name]:
                self.assertEqual(len(row), len(set(row)))
                self.assertEqual(len(row), 17)
                for group in plan["groups"]:
                    self.assertEqual(sum(i in group["candidates"] for i in row), group["m"])

    def test_single_configuration_exact_null_has_p_one(self):
        plan = {"M": 1, "groups": [{"candidates": [1, 2], "m": 2}]}
        draws, exact = amendment.configurations(plan, [1, 2], distinct=True)
        self.assertTrue(exact)
        result = amendment.statistic(np.array([2., 3.]), np.array([1, 2]), [1, 2], draws, exact)
        self.assertEqual(result["raw_p"], 1.)

    def test_small_exact_null_matches_manual_enumeration(self):
        plan = {"M": 6, "groups": [{"candidates": [1, 2, 3, 4], "m": 2}]}
        draws, exact = amendment.configurations(plan, [3, 4], distinct=True)
        self.assertTrue(exact)
        self.assertEqual(draws.tolist(), [list(x) for x in combinations([1, 2, 3, 4], 2)])
        result = amendment.statistic(np.arange(1., 5.), np.arange(1, 5), [3, 4], draws, exact)
        self.assertEqual(result["raw_p"], 1 / 6)

    def test_nonfinite_evaluated_score_is_rejected(self):
        with self.assertRaises(ValueError):
            amendment.statistic(np.array([np.nan, 1.]), np.array([1, 2]), [1], np.array([[2]]), False)

    def test_legacy_and_published_protocol_parity_were_verified(self):
        receipt = json.loads((RESULTS / "verification.json").read_text())
        self.assertEqual(receipt["frozen_prediction_files_verified"], 59)
        self.assertTrue(all(receipt["parity"].values()))

    def test_fixed_cuts_isolate_missing_residue_exclusion(self):
        receipt = json.loads((RESULTS / "verification.json").read_text())
        for value in receipt["fixed_cut_sensitivity"].values():
            self.assertEqual(value["unknown_residues_in_historical_sample"], [])
            self.assertTrue(value["fixed_cut_configs_equal_historical_configs"])
            self.assertTrue(all(g["m"] == 0 for g in value["historical_unknown_strata"]))


if __name__ == "__main__":
    unittest.main(verbosity=2)
