"""Numerical regression tests for the unchanged convex Gaussian objective."""
from decimal import Decimal, localcontext
import unittest
import numpy as np
import algebra_check as ac


class StableNewtonTests(unittest.TestCase):
    def test_sub_ulp_decrease_matches_high_precision_diagonal_reference(self):
        K = np.diag([1., 2., 4.])
        target = np.diag([1., .5, .25])
        C = target * (1 + 1e-8)
        D = target - C
        with localcontext() as context:
            context.prec = 70
            exact = Decimal(0)
            for k, c, d in zip(K.diagonal(), C.diagonal(), D.diagonal()):
                kd, cd, dd = map(Decimal.from_float, map(float, (k, c, d)))
                exact += (kd*dd - ((cd+dd)/cd).ln()) / 2
        actual = ac.objective_difference(C, D, K, [], gamma=0.)
        self.assertLess(actual, 0.)
        self.assertLess(abs(actual-float(exact)), 1e-22)
        # The decrease is below the spacing of the separately evaluated objective.
        self.assertLess(abs(actual), abs(np.spacing(ac.objective(C,K,[],gamma=0.))))

    def test_difference_agrees_away_from_roundoff_for_all_parameters(self):
        K, terms = ac.setup()
        rng = np.random.default_rng(20260912)
        for beta in (.5, 1., 2.):
            for gamma in (0., .2, 1.):
                A = rng.normal(size=(3,3))
                C = A@A.T + np.eye(3)
                D = .04 * rng.normal(size=(3,3)); D = (D+D.T)/2
                direct = ac.objective(C+D,K,terms,gamma,beta)-ac.objective(C,K,terms,gamma,beta)
                stable = ac.objective_difference(C,D,K,terms,gamma,beta)
                self.assertLess(abs(stable-direct), 2e-13)

    def test_indefinite_trial_is_rejected(self):
        K, terms = ac.setup(); C = np.linalg.inv(K)
        self.assertTrue(np.isinf(ac.objective_difference(C,-2*C,K,terms)))

    def test_rotated_coordinates_and_starts_reach_original_stationarity(self):
        K, terms = ac.setup(); harmonic = np.linalg.inv(K)
        rng = np.random.default_rng(20260912)
        for gamma in (0.,1.):
            reference,_ = ac.optimize(harmonic,K,terms,gamma)
            for _ in range(5):
                Q,_ = np.linalg.qr(rng.normal(size=(3,3)))
                Kr = Q.T@K@Q
                tr = [(a,T@Q,w) for a,T,w in terms]
                for scale in (.25,1.,4.):
                    C, history = ac.optimize(scale*(Q.T@harmonic@Q),Kr,tr,gamma)
                    self.assertLess(np.linalg.norm(ac.gradient(C,Kr,tr,gamma)),1e-11)
                    self.assertLess(np.max(np.abs(Q@C@Q.T-reference)),1e-10)
                    self.assertLess(len(history),100)
                    self.assertGreater(np.linalg.eigvalsh(C)[0],0.)
                    for step in history[:-1]:
                        self.assertLess(step['stable_objective_decrease'],0.)

    def test_harmonic_near_stationary_start_converges_without_tolerance_change(self):
        K=np.diag([1.,2.,4.]);target=np.linalg.inv(K)
        for scale in (1+1e-8,1-1e-8):
            C,history=ac.optimize(scale*target,K,[],0.)
            self.assertLess(np.linalg.norm(ac.gradient(C,K,[],0.)),1e-11)
            self.assertLess(np.max(np.abs(C-target)),1e-10)
            self.assertLessEqual(len(history),3)


if __name__=='__main__':
    unittest.main()
