# A proposed Gaussian closure that retains every internal coordinate

This directory contains a mathematical assessment, a proposed two-week experiment, and one **exploratory 3 × 3 algebra check**. No protein fit, nonlinear-response benchmark, ligand-label analysis or quantum experiment has been run here. The Gaussian closure is a candidate for testing, not a demonstrated repair of the failed representation tests.

The proposed model retains every internal displacement coordinate while approximating their joint distribution by one Gaussian. The covariance is selected by a strictly convex free-energy objective with a log determinant. An additional Ornstein–Uhlenbeck approximation supplies dynamics. The equilibrium variational principle does not validate those dynamics.

The full derivation and experiment plan are in `design.md`. The proposed parameters and gates are also recorded in `proposed-experiment.json`; their status remains proposed. `sources.json` records the directly read primary literature and its claim boundaries.

To repeat only the algebra check:

```sh
python3 -m pip install -r requirements.txt
python3 algebra_check.py
python3 coordinate_domain_check.py
```

The check uses one triangle with all three internal coordinates, κ = β = 1, and six deterministic optimizations: three positive-definite starts for the harmonic objective and three for the quartic objective. It checks the exact Gaussian energy against Gauss–Hermite quadrature, analytical derivatives, positive definiteness, stationarity, covariance contraction and the harmonic limit. `algebra-check.json` records the results and source hash. Runtime fields can vary between repeats.

This calculation verifies algebra and numerical implementation only. The fixed-centroid optimum has a nonzero gradient in the excluded centroid direction, which illustrates an approximation that must be tested. Its lower variational free energy does not establish improved covariance, response or biological prediction. No current or prospective success claim should be inferred from this directory.

The second script independently constructs a distant rigidly equivalent zero-energy configuration within the triangle's linear internal-coordinate chart. It demonstrates why local box refinement and removal of infinitesimal rigid modes do not certify a global nonlinear equilibrium or a physical single-basin distribution. It computes no response or state population.

The numerical line search was corrected without changing the objective or acceptance tolerances; [the correction record](solver-correction.md) preserves the original script and local results. Run `python -m unittest -v test_solver.py` before the algebra checks.
