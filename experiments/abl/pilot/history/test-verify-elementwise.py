#!/usr/bin/env python3
"""Verify portability without allowing scientific response or rank corruption."""
import unittest
import numpy as np
from verify import compare_array, RAW_SCALE_FIELDS


class ComparisonPolicyTests(unittest.TestCase):
    def test_observed_linux_scale_roundoff_passes(self):
        # The Linux failure was 1.5916157281026244e-10 in monomial_mean.
        reference = np.array([4753.236120318028])
        actual = reference + 1.5916157281026244e-10
        self.assertGreater(float(abs(actual-reference)[0]), 1e-10)
        result = compare_array('monomial_mean', reference, actual)
        self.assertEqual(result['policy'], 'raw_scale_aware')
        self.assertFalse(result['bitwise_equal_values'])

    def test_raw_relative_corruption_fails_for_every_whitelisted_field(self):
        for key in RAW_SCALE_FIELDS:
            with self.subTest(key=key), self.assertRaises(AssertionError):
                compare_array(key, np.array([2.4e8]), np.array([2.4e8*(1+1e-7)]))

    def test_raw_near_zero_corruption_fails(self):
        with self.assertRaises(AssertionError):
            compare_array('monomial_mean', np.zeros(1), np.array([2e-10]))

    def test_normalized_fields_do_not_inherit_relative_tolerance(self):
        # A relative policy would accept this at 0.8, but the strict policy fails.
        for key in ('C', 'score', 'equilibrium', 'graph_heat_kernel'):
            with self.subTest(key=key), self.assertRaises(AssertionError):
                compare_array(key, np.array([0.8]), np.array([0.8+1.5e-10]))

    def test_default_floating_fields_are_strict(self):
        with self.assertRaises(AssertionError):
            compare_array('unlisted_future_field', np.array([2.4e8]), np.array([2.4e8+1e-3]))

    def test_small_normalized_roundoff_passes(self):
        result=compare_array('score', np.array([0.8]), np.array([0.8+5e-11]))
        self.assertEqual(result['policy'], 'strict_absolute')
        self.assertEqual(result['rtol'], 0)

    def test_rank_boolean_and_integer_changes_fail_exactly(self):
        fixtures=[('order', np.array([3, 2, 1]), np.array([2, 3, 1])),
                  ('candidate', np.array([True, False]), np.array([True, True])),
                  ('canonical', np.array([242, 243]), np.array([242, 244])),
                  ('degree', np.array([30.]), np.array([30.+1e-12]))]
        for key,reference,actual in fixtures:
            with self.subTest(key=key), self.assertRaises(AssertionError):
                compare_array(key,reference,actual)

    def test_nonfinite_values_fail_even_when_equal(self):
        for key in ('C','monomial_mean'):
            for value in (np.nan,np.inf,-np.inf):
                for reference,actual in ((np.array([value]),np.array([value])),
                                         (np.ones(1),np.array([value]))):
                    with self.subTest(key=key,value=value), self.assertRaises(AssertionError):
                        compare_array(key,reference,actual)

    def test_dtype_and_shape_changes_fail(self):
        for actual in (np.ones(1,dtype=np.float32),np.ones((1,1))):
            with self.assertRaises(AssertionError):
                compare_array('C',np.ones(1),actual)

    def test_real_primary_matrix_and_order_perturbations_fail(self):
        # Exercise the actual scientific output, rather than only toy values.
        from run import ROOT
        with np.load(ROOT/'results/biquadratic-d2-n33-e4.npz') as result:
            for key in ('C','score'):
                altered=result[key].copy()
                altered.flat[0]+=1e-7
                with self.subTest(key=key), self.assertRaises(AssertionError):
                    compare_array(key,result[key],altered)
            altered=result['order'].copy()
            altered[:2]=altered[1::-1]
            with self.assertRaises(AssertionError):
                compare_array('order',result['order'],altered)


if __name__=='__main__':
    unittest.main(verbosity=2)
