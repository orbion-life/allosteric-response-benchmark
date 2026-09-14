"""Adversarial checks for rounding-only grid correspondence."""
import unittest
import numpy as np
from ci_protocol_check import match_grid


class GridPortabilityTests(unittest.TestCase):
    def test_one_ulp_platform_rounding_preserves_unique_saved_nodes(self):
        expected = np.geomspace(.001,30,61)
        stored = np.nextafter(expected,np.inf)
        original = stored.copy()
        result = match_grid(expected,stored,exact_size=True)
        self.assertEqual(result['expected_nodes'],61)
        self.assertEqual(result['exact_matches'],0)
        self.assertEqual(result['maximum_ulp_difference'],1)
        np.testing.assert_array_equal(stored,original)

    def test_duplicate_missing_shifted_and_ambiguous_nodes_are_rejected(self):
        expected = np.array([1.,2.,3.])
        for stored in [np.array([1.,2.,2.,3.]),np.array([1.,3.]),np.array([1.,2.001,3.]),np.array([.5,1.5,2.,3.])]:
            with self.subTest(stored=stored), self.assertRaises(ValueError):
                match_grid(expected,stored)
        with self.assertRaises(ValueError):
            match_grid(expected,np.array([1.,2.,3.,4.]),exact_size=True)

    def test_beyond_rounding_budget_is_rejected(self):
        expected = np.array([1.,2.,3.]);stored=expected.copy()
        stored[1] += 33*np.spacing(stored[1])
        with self.assertRaises(ValueError):
            match_grid(expected,stored)


if __name__=='__main__':
    unittest.main()
