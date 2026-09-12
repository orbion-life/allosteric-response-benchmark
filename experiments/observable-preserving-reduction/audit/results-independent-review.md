# Independent review of the response-operator results

**Date:** 12 September 2026. **Status:** The recorded finite-model campaign and an independent full replay agree. The campaign does not pass all declared fidelity, cost, grid or domain gates.

The independent reviewer read the frozen protocol and implementation, recomputed every saved checkpoint decision from the numerical arrays, repeated the complete campaign in `replay-independent/`, and compared the development references with immutable files from public release `v0.3.0`. No scientific source, protocol, primary result, or public repository file was changed during this audit. The audit code is `audit/review_results.py`.

## Source and protocol binding

The replay used Python 3.12.14, NumPy 2.1.3 and SciPy 1.14.1. The numerical package versions match `requirements.txt`. `OPENBLAS_NUM_THREADS`, `OMP_NUM_THREADS` and `VECLIB_MAXIMUM_THREADS` were each set to `1`; `PYTHONWARNINGS=error` and Python's `-W error` were enabled. The unchanged watchdog enforced the recorded 1,800-second and 4,000,000,000-byte limits. The replay was a fresh directory, not a resumed run.

| File | SHA-256 |
|---|---|
| `preanalysis.json` | `530fede82261aec39165f06d0126e75e995ad131c0af0db170ce6f9c6880ada3` |
| `operator_study.py` | `13d80e28d8acaf9ca435938d8ef08668b3cf9d0e750b7b0970b6cf26f19434ef` |
| `bounded_run.py` | `9eb3b7c152cf2855224d0948794c9190848cb050ddf2c7c428215a65733a4e62` |

All 61 source, protocol and primary-result files recorded in `audit/independent-before-replay-hashes.json` retained their hashes after the replay. Both runs report the same source and protocol hashes. The seven source receipts, including the inherited model code and pre-execution tests, were checked before launch.

The reproducible launch, from the response-operator directory with the pinned environment active, was:

```sh
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 PYTHONWARNINGS=error python -W error bounded_run.py replay-independent
```

The replay directory must not already exist. Installation and the later independent audit are outside the campaign timing.

## Coverage and numerical repeatability

Each run contains all 27 declared finite problems: three geometries, three stiffnesses and two primary grids give 18 operator comparisons; nine larger-domain cases provide full-reference diagnostics only. Each problem has three query times. The study therefore contains 81 full-reference response matrices, not 27 independent operator reductions.

The replay reproduces all **27 archives and 1,065 arrays**. The maximum absolute difference is **1.1657341758564144 × 10⁻¹⁵**, below the declared **10⁻⁹** numerical comparison allowance. Shapes, data types, finite values and field sets agree. Integer and Boolean fields, where present, require exact agreement. Every selected rank and every bound-selection and finite-fidelity decision agrees. Timing is compared separately rather than subjected to a numerical-array tolerance.

`audit/independent-numeric-comparison.json` contains the per-field differences. `audit/independent-primary-review.json` and `audit/independent-replay-review.json` contain separately recomputed gate decisions and timings. The independent audit did not replace primary decisions with the replay's measurements.

## Rank selection, residual bounds and fidelity

The code constructs each basis and evaluates its bounds before it calls the full-response solver. Each basis uses the original equilibrium-weighted, centered quartic node observables, normalized by the fixed all-mode native harmonic standard deviations. All three internal coordinates remain in the full finite model. The basis is rebuilt for each geometry, stiffness and grid, using the same declared algorithm. Numerical basis vectors are not transferred from the development geometry.

The saved data contain **93 checkpoints**. At every checkpoint, the reported static, delayed and response errors and the calculated response bound match an independent calculation from the arrays. The final checkpoint is the first one to pass the prescribed bound, static-error and orthogonality conditions, unless the run reaches the rank cap without passing. No earlier successful checkpoint was bypassed. The largest observed static error among the final checkpoints is **2.275957200481571 × 10⁻¹⁴**; the largest orthogonality error over all checkpoints is **9.103828801926284 × 10⁻¹⁵**.

All 93 checkpoint bounds cover the observed response errors. The minimum unadjusted bound-minus-error slack is **1.50987229944027 × 10⁻¹⁴**: no coverage decision needs the 10⁻⁹ arithmetic allowance in these runs. This is empirical confirmation of the evaluated bounds. It is not a directed-rounding or interval-arithmetic proof.

The implementation uses the same clipped positive reduced operator in its propagation and residual. Its one-sided bound includes the seed projection defect. Its two-sided bound also includes the measured mass-matrix and parallel-residual terms, with coefficient norms. The measured static discrepancy enters the response bound. These details agree with the pre-execution mathematical audit in `../operator-mathematical-audit.md`; no unclipped residual is substituted for the clipped trajectory.

The results distinguish observed accuracy from acceptance by a calculated bound:

| Condition | Primary | Independent replay |
|---|---:|---:|
| Final observed response error ≤ 0.002 | 18/18 | 18/18 |
| Bound-backed finite-fidelity gate | 17/18 | 17/18 |
| Finite fidelity, dimension reduction and total-cost gate together | 16/18 | 16/18 |

For the development geometry at κ = 1 on the 33³ grid, rank 96 reaches the allowed cap. Its observed response error is **0.00023496999449096334**, but its bound is **0.008672054095020159**, above 0.002. The code correctly records `selected_by_bound=false`. The field `selected_rank=96` means the final attempted rank in this case; it must not be described as an accepted rank. This case fails conservative acceptance, although its measured finite response error passes.

Successful ranks range from 24 to 96. The fourfold dimension condition passes throughout, but the imposed cap already makes this condition easy for grids of 4,913 and 35,937 states. It is not independent evidence of useful total cost. All 12 compact/elongated primary-grid comparisons pass their finite-fidelity and practical conditions in both runs. This is transfer of the fixed construction across two synthetic geometries, not validation on an independent protein family.

## Cost accounting and its limits

The accounting includes construction of the native Hessian and normalizers, equilibrium quadrature, the full sparse finite operator, observable weighting, seed decomposition, basis growth, both reorthogonalization passes, all unsuccessful checkpoints, residual bounds, projected exponentials and reconstruction at all three times. The full-reference route uses the same common setup and three full sparse exponential actions. The operator route does not report only its final small-matrix solve, and no unreported amortization is applied.

For each route, total compute time equals common setup plus its corresponding solve. The independent audit checked this arithmetic and the ratio. Common setup appears once in each alternative's cold comparison; the overall experimental workflow constructs it once before running both alternatives. Per-route times exclude result serialization and most validation/reporting overhead. The campaign wall time includes execution and saved-output overhead; it is a different quantity from a single route's cost.

The development geometry at κ = 1 on the 17³ grid passes fidelity but misses the total-cost gate in both runs:

| Measurement | Primary | Independent replay |
|---|---:|---:|
| Operator route, including shared setup | 0.281434 s | 0.285575 s |
| Full-reference route, including shared setup | 0.267648 s | 0.282977 s |
| Operator/reference ratio | 1.051506 | 1.009180 |

This is a narrow empirical cost failure. Two same-host runs do not establish a stable timing margin or a hardware-independent speedup. The κ = 1, 33³ development case has a favorable time ratio but fails the bound condition; favorable timing cannot rescue that rejection.

The primary watchdog recorded **179.095080 s** and a sampled aggregate peak of **799,391,744 bytes**. The independent watchdog recorded **179.787357 s** and **808,501,248 bytes**. Both are within the declared limits. Sampling occurs every 0.25 seconds and can miss an instantaneous peak. Native process maximum-RSS values are also saved; these are cumulative process-lifetime peaks, not isolated per-case peaks. The watchdog wrapper's own memory is not included in the monitored worker process group. No claim is made that these measurements establish protein-scale storage or quantum-circuit cost.

## Grid and domain diagnostics

The independent audit computed the two diagnostics directly from the saved full responses, holding the declared physical normalization and times fixed. Grid refinement compares 17³ and 33³ on the same ±4 domain. Domain expansion compares 33³/±4 with 41³/±5, at the same spacing. Each gate requires the largest difference over all observable pairs and all three times to be at most 0.001.

| Geometry | κ | Maximum grid difference | Maximum domain difference |
|---|---:|---:|---:|
| Development | 1 | 0.01106825 | 0.04443822 |
| Development | 10 | 0.02284124 | 0.00946862 |
| Development | 100 | 0.02544532 | 0.00262070 |
| Compact | 1 | 0.00308797 | 0.01754573 |
| Compact | 10 | 0.01814884 | 0.03398653 |
| Compact | 100 | 0.02200244 | 0.00364401 |
| Elongated | 1 | 0.01924783 | 0.05694048 |
| Elongated | 10 | 0.04178030 | 0.03335059 |
| Elongated | 100 | 0.04663106 | 0.00649348 |

**Neither diagnostic passes for any of the nine geometry–stiffness combinations.** Some individual times pass the domain threshold, but acceptance was declared over all three times. The corresponding finite-grid projection can be accurate while the full finite reference still changes materially with the grid or domain. These data do not establish continuum accuracy, adequate basin coverage, or biological relevance. A universal pass would contradict the protocol's explicit failure rule.

## Parity with unchanged v0.3.0 references

The audit read the six development archives directly from the immutable `v0.3.0` Git object, rather than from a possibly changed working copy. The annotated tag object is `c58009f1574155a5f8053a66aac65e20e2d383f6`; its commit is `df69d29c0870ff636c56e7194de037d8693b3e19`.

For κ = 1, 10 and 100 and grids 17³ and 33³, the rebuilt full reference agrees with the saved response, equilibrium and delayed covariance arrays to at most **1.687538997430238 × 10⁻¹⁴**. Times and harmonic standard deviations agree exactly. The new implementation weights observables before contraction, whereas the earlier archives normalized after contraction; the measured differences are at floating-point scale. Individual archive hashes and differences are preserved in `audit/independent-v0.3.0-parity.json`.

## Scientifically supportable conclusion

The recorded experiment supports a reproducible, observable-seeded reduction of the complete nonlinear **finite** operator in most of this small synthetic test set, with bounds evaluated before access to the full response and with complete setup and query costs included. It also identifies a bound failure, a marginal cost failure, and unresolved discretization and domain sensitivity in every tested geometry–stiffness combination. Those failures remain part of the result. This study does not validate allosteric rankings, protein-scale operator access, continuum response, a quantum implementation, or quantum advantage.
