# Independent diagnosis of the Linux portability failure

**Date:** 12 September 2026. **Disposition:** The strict portable verification remains **FAIL**. No scientific source, protocol, tolerance, verifier, primary result or public file was changed in this diagnosis.

The Linux campaign completed within its recorded resource limits, and all 14 unit tests passed. Its observable responses agree with the primary run at the declared tolerance. However, 61 error-bound arrays fail that same cross-platform comparison. The reduced-operator spectra also differ, so the large bound changes cannot be dismissed as a pure rotation of an otherwise identical basis. These are different numerical reduced representations with closely agreeing queried responses. Their precise point of divergence has not been identified.

## Evidence and immutable binding

The workflow status and commit were read directly from GitHub: [run 34717448039](https://github.com/orbion-life/allosteric-response-benchmark/actions/runs/34717448039) completed with conclusion `failure` at commit `fa9d45f2531dc70c7dd086c38530c3f13a101082`. The downloaded verification receipt and numerical archives were then read directly. The reviewed receipt has SHA-256 `4e04a4aeacf7162636447c1be4f4842692e504e001c96a30595130b9f6814eb8`.

| Item | SHA-256 |
|---|---|
| Scientific implementation | `13d80e28d8acaf9ca435938d8ef08668b3cf9d0e750b7b0970b6cf26f19434ef` |
| Pre-execution protocol | `530fede82261aec39165f06d0126e75e995ad131c0af0db170ce6f9c6880ada3` |
| Strict portable verifier | `d5a2516f09ec463652b3f4eb490e6afe3e606138e14a3cd5aeb4b31f139c0fba` |

The Linux receipt records Python 3.12.14, NumPy 2.1.3 and SciPy 1.14.1, matching the named package versions of the original-host replay. Matching those versions does not establish identical BLAS/LAPACK kernels, compiler behaviour or floating-point execution. The watchdog recorded 348.150297 seconds and a sampled aggregate peak of 542,924,800 bytes, within the unchanged 1,800-second and 4,000,000,000-byte limits. The monitoring and native-RSS limitations documented in the primary audit continue to apply.

## What reproduces, and what does not

The full comparison covers 27 archives and all 1,065 recorded arrays. The portable checker compares 786 arrays numerically and checks 279 reduced-coordinate arrays through their shape, finite-value, positivity and physical-reconstruction conditions. The numerical failures are confined to 61 of the 372 bound arrays.

| Output group | Arrays | Largest absolute difference | Result at the existing allowance |
|---|---:|---:|---|
| Responses, static/delayed covariances, times and scales | 414 | 5.31902437025944 × 10⁻¹⁰ | All pass 10⁻⁹ |
| Combined response bounds | 93 | 8.038895370634548 × 10⁻⁶ | 10 arrays fail |
| One-sided bounds | 93 | 3.800293256233621 × 10⁻⁴ | 14 arrays fail |
| Two-sided bounds | 93 | 8.038895362318609 × 10⁻⁶ | 10 arrays fail |
| Uniform bounds | 93 | 16.71540424762327 | 27 arrays fail |

All 61 failure records have the context `semantic physical array`; none arises from a covariance/response mismatch, a changed selected rank, a changed fidelity decision, failed physical reconstruction from the stored reduced coordinates, or failed case/source/protocol identity. The static-covariance maximum difference is 1.432187701766452 × 10⁻¹⁴. The largest response difference is 5.319008566928574 × 10⁻¹⁰; the slightly larger maximum in the table belongs to delayed covariance.

This partition is a diagnostic summary of a failed full check. It is not a replacement acceptance policy or a new overall `PASS`. In particular, the bound-comparison tolerance remains 10⁻⁹, and every failed bound comparison stays in the evidence.

## Why a pure basis-gauge explanation is insufficient

For a fixed full operator and an orthogonal change of reduced basis, $V'=VO$, one has

$$
H_r'=O^{\mathsf T}H_rO,\qquad
\mathcal R'=\mathcal R O.
$$

Thus the eigenvalues of $H_r$ and $\mathcal R^{\mathsf T}\mathcal R$, and the spectral norm of $\mathcal R$, are invariant. Raw matrix entries may change under this transformation, but these spectra cannot.

Direct eigendecomposition of the saved arrays gives the following counterexample to a pure-gauge interpretation. For the elongated geometry at κ = 100, grid 33³ and rank 48:

| Quantity | Primary host | Linux |
|---|---:|---:|
| Residual spectral norm | 3,082.3106043087982 | 2,934.786088026689 |
| Maximum combined response bound | 0.0001445118591981439 | 0.00013647296382750934 |
| Observed error against the run's full reference | 3.2094751412081735 × 10⁻⁸ | 3.1605769557074836 × 10⁻⁸ |
| Uniform-bound entry at the largest cross-host difference | 349.24410576766576 | 332.5287015200425 |

The largest difference between sorted reduced-operator eigenvalues is **34.733287859430675**; for the residual-Gram eigenvalues it is **887,669.2789594736**. The two saved response matrices nevertheless differ by only **4.889881832781384 × 10⁻¹⁰**. The development geometry at the same stiffness, grid and rank also has a reduced-operator spectral difference of approximately 699.35. This is not explained by a sign convention for singular vectors.

The evidence is consistent with backend-sensitive Krylov directions that contribute little to the requested observables while affecting residual norms and bounds. This interpretation is an inference, not a located root cause. The existing archives do not retain every full basis vector, residual-block singular value or cutoff margin. They therefore cannot establish which block first diverged, whether a near-dependency decision or checkpoint truncation caused it, or how much comes from numerical differences in the full-model construction rather than subsequent basis growth. The closely matching full-reference responses constrain the discrepancy but do not prove equality of every full-operator entry.

## The mathematical bound checks remain informative

The residual argument does not require the primary and Linux reductions to be identical. It applies separately to each actual reduced trajectory under the stated positive-operator assumptions, with the recorded seed, static, mass-matrix and parallel-residual corrections. In particular, using the actual simulated projected operator in the residual remains necessary when tiny negative eigenvalues are clipped.

Both hosts have **93 of 93** checkpoint errors below their own calculated response bounds. The minimum raw bound-minus-error slack is **1.50987229944027 × 10⁻¹⁴** on the primary host and **1.5536277977954484 × 10⁻¹⁴** on Linux. These coverage checks pass even without adding the 10⁻⁹ comparison allowance. Every selected rank and bound-based admission decision agrees between hosts. This diagnosis identifies no counterexample to the finite-model inequality or defect in its existing source implementation. Empirical coverage does not convert floating-point evaluation into an interval-arithmetic certificate.

The numerical outcomes remain distinct from timing:

| Gate | Primary host | Linux |
|---|---:|---:|
| Bound-backed finite fidelity | 17/18 | 17/18 |
| Total-cost ratio alone ≤ 1 | 17/18 | 18/18 |
| Combined finite-fidelity/dimension/cost acceptance | 16/18 | 17/18 |
| Grid screen | 0/9 | 0/9 |
| Domain screen | 0/9 | 0/9 |

For the development κ = 1, 17³ case, the time ratio changes from 1.051506 on the primary host to 0.581689 on Linux. This is allowed by the existing comparison policy: measured performance can change across machines. It does not erase the original cost failure. The development κ = 1, 33³ case remains the unaccepted bound-screen case on both hosts, although its observed finite-response error is below 0.002.

## Recommended closure of this iteration

Preserve the failed workflow, its complete arrays and diagnostics, the original comparison policy, and the primary results. Report the narrower outcome explicitly: **the tested observable responses reproduce within 10⁻⁹, but the full portability check fails because the reduced representation's bounds do not reproduce to that tolerance.** Keep the response-only summary separate from the strict overall disposition. Do not relabel the result as a gauge-only mismatch, silently widen a bound tolerance, remove inconvenient bound fields, or substitute Linux's more favourable timings for the primary cost record.

A future stability investigation can be bounded and predeclared without changing the physical model. Three diagnostic instances could provide a bounded starting point: the elongated and development κ = 100, 33³ cases at rank 48, plus one coarser-grid control fixed before execution. Such a diagnostic should retain full bases, operator-action differences, residual-block singular values and cutoff margins at every attempted checkpoint. Principal angles and aligned residuals would help distinguish changes of subspace from a simple basis rotation and locate the first divergence. It can retain the existing one-worker, 30-minute and 4-GB ceilings. This investigation has not been run or promised a passing outcome.

Any change to the recurrence, rank-admission rule, arithmetic or acceptance policy would require a separately recorded implementation/protocol and fresh validation. No such change is needed to state the current limitation accurately, and none is made in this diagnosis. Grid/domain convergence, physical basin definition, protein-scale construction, biological validation and quantum usefulness remain unresolved independently of this portability issue.
