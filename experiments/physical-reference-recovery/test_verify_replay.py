"""Verifier tests use only compact archived arrays; no scientific worker runs."""
import json
from pathlib import Path
import shutil
import tempfile
import unittest
import numpy as np
from verify_replay import ROOT, compare_array, run


class ReplayVerifierTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.directory = Path(self.tmp.name)
        self.harmonic = self.directory/'harmonic'
        self.nonlinear = self.directory/'nonlinear'
        shutil.copytree(ROOT/'harmonic-primary', self.harmonic)
        shutil.copytree(ROOT/'nonlinear-primary', self.nonlinear)

    def verify(self):
        return run(self.harmonic, self.nonlinear)

    def mutate_array(self, path, key, operation):
        with np.load(path) as z:
            arrays = {name: z[name] for name in z.files}
        arrays[key] = operation(arrays[key])
        np.savez_compressed(path, **arrays)

    def mutate_receipt(self, path, operation):
        record = json.loads(path.read_text())
        operation(record)
        path.write_text(json.dumps(record))

    def test_archived_copy_passes_and_retains_scientific_failures(self):
        result = self.verify()
        self.assertEqual(result['status'], 'PASS', result['failed_checks'])
        self.assertEqual(result['arrays_compared'], 82)
        screens = result['scientific_screens']
        self.assertTrue(any(s.get('screen_pass') is False for s in screens))
        old = screens[-1]
        self.assertEqual(old['component_screen_pass'], dict(G0=True, K=False, C=False))
        self.assertFalse(old['overall_screen_pass'])

    def test_wrong_shape_fails_with_complete_diagnostic_receipt(self):
        self.mutate_array(self.nonlinear/'nonlinear-h0.5/response.npz', 'C', lambda x: x[:1])
        result = self.verify()
        self.assertEqual(result['status'], 'FAIL')
        self.assertEqual(result['arrays_compared'], 82)
        self.assertTrue(any('shape' in x.get('reason', '') for x in result['array_checks']))
        json.dumps(result, allow_nan=False)

    def test_changed_normalized_response_fails(self):
        self.mutate_array(self.harmonic/'calibration/harmonic-A8-h0.015625.npz', 'C', lambda x: x+1e-6)
        result = self.verify()
        self.assertEqual(result['status'], 'FAIL')
        self.assertGreater(result['maximum_array_absolute_difference'], 1e-9)

    def test_nan_array_is_rejected(self):
        self.mutate_array(self.harmonic/'calibration/harmonic-A8-h0.015625.npz', 'C', lambda x: x*np.nan)
        result = self.verify()
        self.assertEqual(result['status'], 'FAIL')
        json.dumps(result, allow_nan=False)

    def test_wrong_screen_cannot_relabel_a_historical_failure(self):
        def relabel(receipt):
            self.assertFalse(receipt['comparisons'][0]['all_components_within_0_001'])
            receipt['comparisons'][0]['all_components_within_0_001'] = True
        self.mutate_receipt(self.harmonic/'calibration/receipt.json', relabel)
        result = self.verify()
        self.assertEqual(result['status'], 'FAIL')
        self.assertTrue(any('harmonic/pair/' in x['name'] for x in result['failed_checks']))

    def test_failed_milestone_and_cap_violation_are_rejected(self):
        self.mutate_receipt(self.harmonic/'calibration/receipt.json',
                            lambda x: x.update(harmonic_implementation_milestone_pass=False))
        self.mutate_receipt(self.nonlinear/'watchdog.json', lambda x: x.update(wall_seconds=901))
        result = self.verify()
        self.assertEqual(result['status'], 'FAIL')
        names = {c['name'] for c in result['failed_checks']}
        self.assertIn('harmonic/implementation-milestone', names)
        self.assertIn('nonlinear/wall-cap', names)

    def test_missing_array_is_rejected_without_default_values(self):
        path = self.harmonic/'calibration/analytic.npz'
        with np.load(path) as z:
            arrays = {key: z[key] for key in z.files if key != 'times'}
        np.savez_compressed(path, **arrays)
        self.assertEqual(self.verify()['status'], 'FAIL')

    def test_below_allowance_array_difference_is_allowed(self):
        a = np.asarray([1.0, 2.0])
        self.assertTrue(compare_array(a+2e-10, a)['pass_check'])
        self.assertFalse(compare_array(a+2e-9, a)['pass_check'])

    def test_timing_variation_within_caps_is_not_a_physical_failure(self):
        self.mutate_receipt(self.nonlinear/'watchdog.json', lambda x: x.update(wall_seconds=100))
        self.assertEqual(self.verify()['status'], 'PASS')


if __name__ == '__main__':
    unittest.main()
