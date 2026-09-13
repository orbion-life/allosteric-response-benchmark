# Executable harmonic calibration and nonlinear cost profiling

This experiment repairs the first milestone of the proposed v11 physical-reference plan. The old five-case nonlinear sequence requires 1,114,112, 1,843,200 and 2,801,664 states for its three largest development problems, exceeding its one-million-state cap. The old protocol and all old scientific outcomes remain unchanged.

**Completed:** the exact harmonic implementation and the predeclared finer harmonic screens pass. Two bounded nonlinear cost probes and one separately declared largest-grid construction-only capacity check are also complete. **Unresolved:** nonlinear grid/domain convergence and physical interpretation. See the [measured results and next allocation](results-and-recommendation.md). The later nine-hour campaign is proposed and has not run.

The new work separates a tractable harmonic calibration from a bounded nonlinear cost study. Both retain the original quartic contact-energy observables, all three positive internal coordinates, flat measure in the affine coordinate plane, constant mobility in those coordinates, and the prescribed physical times and harmonic normalization. The rescaled grid has the correctly transformed anisotropic mobility. No basin restriction is introduced.

## The first milestone has an affordable exact reference

The continuum harmonic process is Gaussian. `model.py` computes its static and delayed polynomial moments by a joint-Gaussian recurrence. `calibration.py` independently evaluates the same quantities through a Hermite expansion, obtaining coefficients from exact-degree Gaussian quadrature applied to the original geometric contact law. Agreement of these methods tests both the observable and its dynamics.

For a finite harmonic Cartesian box, the generator is a Kronecker sum of three one-dimensional generators. Its exponential factors exactly. The code contracts three 5-by-5 polynomial moment matrices for each time after solving the one-dimensional symmetric tridiagonal generators. It therefore returns the same finite-grid quartic response without constructing the full three-dimensional grid. A 216-state full-grid propagation checks this equivalence independently.

The complete prospective schedule includes all five v11 harmonic rows, three finer spacings at faces six, and domain faces seven and eight at the finest spacing. The largest one-dimensional problem has 1,024 points, with an 8-MB eigenvector matrix. The corresponding full Cartesian grid would have 1,073,741,824 configurations, which this implementation never allocates. This factorization is specific to the harmonic potential. The nonlinear generator is not separable.

The implementation milestone requires normalized analytic-method and finite-grid parity errors at most 10⁻¹⁰, a measured harmonic-stage time at most 60 seconds, and a worker peak at most 512 MB. All finite-grid errors against the analytic continuum reference are also reported. The distinct 0.001 grid/domain screens retain unsuccessful rows; implementation correctness does not automatically imply that a selected finite grid is accurate enough.

## The nonlinear profile includes independent propagation costs

The prescribed development box remains (4,4,17). Only its 17,408-state and 139,264-state cases are profiled. Each case starts in a fresh process and includes model preparation, equilibrium weights, observables, sparse operator assembly, three full propagation queries, an independent three-query calculation, verification and serialization. A one-worker watchdog caps the two-case nonlinear profile at 900 seconds and 4 GB aggregate resident memory. The harmonic stage has a separate 60-second allocation, giving 960 seconds in total. Dependency installation is excluded. Each worker seeds NumPy's global random state with 1729 before SciPy's possible randomized one-norm estimation; counts and timings are still reported rather than assumed portable.

The primary calculation uses the scaled Taylor matrix exponential action, with explicit counts of operator-vector and operator-block calls, including norm estimation. Its cost is determined by those actions as well as matrix size ([Al-Mohy and Higham, 2011](https://eprints.maths.manchester.ac.uk/1591/)).

The independent calculation uses the Chebyshev expansion of the imaginary-time exponential, with scaled Bessel coefficients. This construction is documented in the [QDYN Chebyshev propagator](https://ag-koch.gitpages.physik.fu-berlin.de/qdyn/v23.04/api/modules/cheby.html); coefficient evaluation uses [SciPy's scaled modified Bessel function](https://docs.scipy.org/doc/scipy-1.14.1/reference/generated/scipy.special.ive.html). The implementation and following truncation bound are independently specified here.

For a symmetric positive semidefinite generator with spectrum in [0,Λ], put z=tΛ/2 and X=I−2H/Λ. Then

$$
e^{-tH}=e^{-z}I_0(z)I+2\sum_{k\ge1}e^{-z}I_k(z)T_k(X).
$$

The Bessel generating function implies, for any s>0,

$$
2\sum_{k>m}e^{-z}I_k(z)
\le 2\exp\{z(\cosh s-1)-(m+1)s\}.
$$

The implementation chooses s=asinh((m+1)/z), then finds the smallest degree meeting the declared bound. It multiplies the operator bound by the actual largest pairwise product of normalized observable-vector norms. The resulting maximum omitted covariance contribution is at most 10⁻⁹. Floating-point coefficient, recurrence and spectral-bound errors remain distinct; this is not an interval-arithmetic certificate. All normalized components must agree with the primary calculation within 10⁻⁸. These are three fixed independent queries per profile case, with no adaptive method retries.

Both CSR and a slice-based matrix-free stencil are checked and timed for ten thin-block actions. Only CSR receives the complete propagation profile. The stencil action timing and storage formula inform a future matrix-free branch; they do not establish its full propagation performance. Because both operators are resident during this profile, its measured memory peak must not be presented as a matrix-free-only peak.

## Reproduce the bounded work

Use Python 3.12 and a new environment. Run these commands from this folder:

```sh
python3.12 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
.venv/bin/python -W error preflight.py --output replay-preflight.json
.venv/bin/python -W error -m unittest -v test_algebra
.venv/bin/python -W error run_bounded.py --stage calibration --output replay-harmonic
.venv/bin/python -W error run_bounded.py --stage profiles --output replay-nonlinear
```

Output directories must be new. The watchdog records sampled aggregate memory and exit status. Workers also record their operating-system lifetime peak RSS, library versions and numerical threadpool settings. Small tests are separate from the timed calibration and profile stages. All code and protocol hashes must be frozen before response execution. A repeated run should compare scientific arrays separately from timing, which is machine-dependent.

`preflight.py` allocates no configuration grid. It enumerates all 45 harmonic and 45 nonlinear requirements of the old protocol across the three geometries and stiffnesses, reconstructs the known congruent minima, and reports exact graph storage alongside explicit planning-workspace assumptions. Matrix storage alone is not a memory admission test. `setup-history.json` preserves the initial static JSON-serialization failure and the isolated-environment decision.

## Check the saved evidence and preview the later campaign

These commands read saved evidence or print an admission preview; they do not launch a convergence campaign:

```sh
.venv/bin/python verify_saved.py --output replay-saved-verification.json
.venv/bin/python supplementary_minima_preflight.py --output replay-supplementary-preflight.json
.venv/bin/python estimate_campaign.py --output replay-campaign-allocation.json
.venv/bin/python continue_campaign.py
```

The published `campaign-allocation-estimate.json` is an input to the frozen continuation protocol. `estimate_campaign.py` recomputes the estimate into a separate output. If an input or formula is changed, the continuation hash check deliberately rejects the mismatch. Use a copy of the experiment for alternative planning assumptions.

The separately recorded construction-only probe can be repeated into a new directory with `.venv/bin/python -W error run_capacity.py --output replay-capacity`. Its fixed limit is 120 seconds and 6 GB. It allocates the 4,358,144-state largest grid and forty additional vectors but performs no exponential propagation. Run it without other numerical work if its timings are to be compared.

The proposed later campaign requires an explicit execution flag:

```sh
.venv/bin/python -W error continue_campaign.py --output new-continuation-run --execute
```

That command can use nine hours and 6 GB. It runs all five corrected grids with both propagators, fixed stage/action limits, both spacing comparisons and both domain comparisons. Its physical admission is uncertain. The archived study has **not** run this command; the preview and source checks are the only completed continuation checks.

Results are interpreted only within the stated scope: the passing exact harmonic calibration is a completed numerical milestone. The two nonlinear profiles support a resource-backed later allocation, including verification, but cannot establish grid convergence, whole-plane coverage, a native-basin interpretation or protein accuracy.
