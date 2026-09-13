import unittest
import contextlib
import io
import json
from pathlib import Path
import tempfile
import numpy as np
from scipy.linalg import eigh
from model import analytic_gaussian, build_model, energies, known_minimum, polynomial_values
from calibration import independent_hermite, max_errors, small_full_grid, tensorized
from chebyshev import action, log_tail_bound, required_degree
from preflight import grid_record


class AlgebraChecks(unittest.TestCase):
    def test_grid_counts(self):
        self.assertEqual(grid_record([4, 4, 17], .125)['states'], 1114112)
        self.assertEqual(grid_record([6, 6, 19], .125)['states'], 2801664)

    def test_direct_quartic_observable(self):
        m = build_model()
        x = np.array([[0, 0, 0], [1., -2., .5], [-.5, .3, 2.], [4., 0., -14.]])
        _, direct = energies(m, x)
        np.testing.assert_allclose(polynomial_values(m, x), direct, rtol=1e-12, atol=1e-12)

    def test_known_minimum(self):
        result = known_minimum(build_model())
        self.assertEqual(result['faces'], [4, 4, 17])
        self.assertLess(result['contact_energy'], 1e-20)

    def test_exact_harmonic_methods(self):
        m = build_model()
        a, b = analytic_gaussian(m), independent_hermite(m)
        self.assertLess(max(max_errors(a, b).values()), 1e-10)
        np.testing.assert_allclose(a['sd'], b['sd'], rtol=1e-12, atol=1e-12)

    def test_tensor_full_grid(self):
        m = build_model()
        sd = analytic_gaussian(m)['sd']
        a, _ = tensorized(m, [3]*3, 1., sd)
        b, _ = small_full_grid(m, [3]*3, 1., sd)
        self.assertLess(max(max_errors(a, b).values()), 1e-10)

    def test_chebyshev_versus_dense_eigen(self):
        A = np.arange(64, dtype=float).reshape(8, 8)/64
        H = A.T @ A+np.diag(np.arange(8)/10)
        F = np.stack([np.ones(8), np.linspace(-1., 1., 8), np.arange(8)**2], axis=1)/10
        vals, V = eigh(H)
        upper = float(np.max(np.sum(np.abs(H), axis=1)))
        for t in (.1, 1., 10.):
            actual, meta = action(H, F, t, upper, covariance_tolerance=1e-11)
            expected = V @ (np.exp(-t*vals)[:, None]*(V.T @ F))
            error = np.max(np.abs(F.T @ (actual-expected)))
            self.assertLess(error, 1e-10)
            self.assertLessEqual(meta['maximum_normalized_covariance_tail_bound'], 1e-11)

    def test_degree_minimal_bound(self):
        for z in (.01, 1., 100., 10000.):
            m = required_degree(z, 1e-10)
            self.assertLessEqual(log_tail_bound(z, m), np.log(1e-10))
            if m:
                self.assertGreater(log_tail_bound(z, m-1), np.log(1e-10))

    def test_complete_tiny_receipt_serialization(self):
        from nonlinear_profile import run
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp)/'tiny'
            with contextlib.redirect_stdout(io.StringIO()):
                run(output, [2, 2, 2], 1.)
            receipt = json.loads((output/'receipt.json').read_text())
            self.assertEqual(receipt['shape'], [4, 4, 4])
            self.assertEqual(receipt['states'], 64)
            self.assertTrue(receipt['independent_comparison_pass'])
            self.assertEqual(len(receipt['propagation']), 3)
            self.assertEqual(len(receipt['independent_chebyshev']), 3)


if __name__ == '__main__':
    unittest.main()
