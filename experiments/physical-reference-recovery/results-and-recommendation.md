# The first numerical milestone is complete; nonlinear convergence remains open

13 September 2026. The revised harmonic calibration is executable and has passed. It retains the original quartic contact-energy observables, all three internal coordinates, the flat affine-coordinate measure, the prescribed mobility and times, and the full harmonic normalization. It does not establish a native protein basin or a nonlinear continuum reference.

The two independent continuum calculations agree within 6.67×10⁻¹⁶ in normalized static covariance, delayed covariance and response. The factorized finite-grid calculation agrees with an independently assembled 216-state full generator within 8.09×10⁻¹⁴. Eight small tests pass, including direct evaluation of the quartic observable, the complete small-case receipt, and an independent matrix-exponential check. Saved harmonic scales and times agree with the existing v10 reference within 3.11×10⁻¹⁴. [Evidence and source verification](saved-verification.json).

The complete harmonic stage took 0.715 seconds of outer elapsed time and reached a sampled aggregate peak of 135.7 MB. The worker's operating-system peak was 81.7 MB. One-thread environment controls were set; this macOS build exposed no threadpools to the introspection library, so an independently observed thread count is unavailable. [Harmonic watchdog](harmonic-primary/watchdog.json).

## Finer harmonic resolution passes without allocating the Cartesian product

The implementation solves three one-dimensional reversible generators and contracts exact degree-four polynomial moment matrices. The full harmonic generator separates, although the original quartic observables need not separate individually. This identity avoids the three-dimensional grid while preserving its finite-model answer. It does not make the nonlinear generator separable.

All five v11 harmonic rows remain recorded. At faces four and spacing 1/8, the largest normalized component error against the analytic reference is approximately 0.0115. At faces six and spacing 1/8, delayed covariance and response still exceed 0.001. Those failures have not been reclassified.

The finer schedule was fixed before outcomes. Both final grid refinements, 1/16→1/32 and 1/32→1/64 at faces six, pass 0.001 separately for every static, delayed and response component. Both fixed-spacing domain expansions, six→seven and seven→eight at spacing 1/64, also pass. At the final grid the largest component error against the exact harmonic reference is 2.30×10⁻⁵. The largest one-dimensional problem has 1,024 points; its billion-state Cartesian product was never allocated. [Complete harmonic results](harmonic-primary/calibration/receipt.json).

## The nonlinear profiles include an independent propagator

The two predeclared cost probes retain faces (4,4,17). They are full nonlinear finite models, not the tensorized harmonic approximation. Each fresh worker includes the native Hessian and exact harmonic normalizers, equilibrium weights, observables, operator construction, three full Taylor propagation queries, three independent Chebyshev queries, verification and serialization.

| Measured quantity | 17,408 states | 139,264 states |
| --- | ---: | ---: |
| Complete outer worker time | 1.829 s | 42.052 s |
| Three Taylor queries | 1.472 s | 40.768 s |
| Taylor vector-equivalent products, including norm estimation | 7,226 | 25,333 |
| Three independent Chebyshev queries | 0.059 s | 0.886 s |
| Chebyshev degrees at the three times | 15, 41, 126 | 27, 79, 249 |
| Maximum normalized propagator disagreement | 1.60×10⁻¹¹ | 1.87×10⁻¹¹ |
| Operating-system worker peak | 69.7 MB | 151.0 MB |

The combined outer profile took 43.882 seconds and reached a sampled aggregate peak of 206.4 MB. These measurements include both propagators. The independent truncation calculation accounts for observable-vector norms; its analytical bound is evaluated in floating point and is not an interval certificate. [Profile watchdog](nonlinear-primary/watchdog.json), [coarse receipt](nonlinear-primary/nonlinear-h0.5/receipt.json), [middle receipt](nonlinear-primary/nonlinear-h0.25/receipt.json).

The observed first refinement does **not** pass the physical-resolution screen. Its maximum normalized changes are 0.000278 for static covariance, 0.004138 for delayed covariance and 0.004285 for response. Agreement between propagators on a fixed generator therefore does not establish grid convergence. The profile also records floating-point underflow counts in remote equilibrium weights; no density floor or transition deletion was introduced. These remain finite-generator cost measurements, not complete thermal-coverage evidence.

The independent Chebyshev calculation is a useful classical baseline for subsequent comparisons. Its observed advantage over the current Taylor implementation concerns these two local finite models and the declared accuracy settings. It supplies no evidence of quantum advantage.

## The corrected coverage margin needs a larger planned grid

The independent reconstruction now includes three congruent rotated minima, rather than only the rotation about the normal axis. The two additional development configurations lie close to a face of the old cost-probe box. Applying the existing two-unit margin rule to every known minimum requires initial faces (6,5,17). This is a correction to the margin inventory, not evidence that every possible minimum or relevant tail has been found. The two completed cost probes remain unchanged. [Supplementary geometric and storage preflight](supplementary-minima-preflight.json).

| Later finite problem | Shape | States |
| --- | --- | ---: |
| Initial faces, spacing 1/2 | 24×20×68 | 32,640 |
| Initial faces, spacing 1/4 | 48×40×136 | 261,120 |
| Initial faces, spacing 1/8 | 96×80×272 | 2,088,960 |
| First domain expansion, spacing 1/8 | 112×96×288 | 3,096,576 |
| Second domain expansion, spacing 1/8 | 128×112×304 | 4,358,144 |

A separately recorded and reviewed capacity probe constructed the largest of these generators, checked its actions, and touched forty additional state vectors. It ran once, without exponential propagation or a response calculation. Outer elapsed time was 5.366 seconds; sampled aggregate peak memory was 3.394 GB. The allocation was 120 seconds and 6 GB, which must be distinguished from those measured values.

The largest grid has 30,332,416 stored CSR entries. The installed implementation uses 381.4 MB for CSR storage; the conservative 64-bit-index preflight is 520.2 MB. A three-column state block occupies 104.6 MB, and the forty-vector capacity scenario occupies 1.395 GB. Measured CSR actions took 0.05823–0.05850 seconds per three-column block. The slice-based stencil used less operator storage but took 0.20575–0.25635 seconds per action. Matrix-free storage alone is therefore not an established speed improvement, and propagation workspace remains substantial. [Capacity receipt](capacity-primary/case/receipt.json), [capacity watchdog](capacity-primary/watchdog.json).

## The next complete campaign is executable but remains proposed

The measured-input allocation model includes all three Taylor and all three independent Chebyshev queries on every corrected grid, preprocessing and serialization, a verification reserve, observed action overhead and a declared twofold engineering margin. It estimates 31,958 seconds, or 8.88 hours. The proposed allocation is **nine hours, 6 GB and one worker**. It is an extrapolation, not an executed runtime or a rigorous runtime bound. The memory admission scenario is 4.50 GB: 1.25 times the measured aggregate capacity peak plus 256 MB. A hard 6-GB watchdog is still required because unseen propagation allocations are not proven to remain below that scenario. [Cost formulas and per-case estimates](campaign-allocation-estimate.json).

The prospective protocol is runnable and uses fresh output directories. Its five stage limits are 60, 480, 6,900, 10,260 and 14,400 seconds, under the global 32,400-second limit. Per-query Taylor action limits are fixed from the recorded cost model. A source mismatch, finite-operator check failure, resource limit or action limit stops execution and preserves the diagnostic receipt. The continuation has **not** been launched. Only its source, admission preview and syntax have been checked.

Both consecutive spacing comparisons and both fixed-spacing domain comparisons remain controlling. Every comparison must satisfy a maximum absolute normalized difference of 0.001 for each of static covariance, delayed covariance and response, across all ordered nodes and prescribed times. The three-level trend is recorded as a diagnostic. Each finite case must retain independent-propagator agreement within 10⁻⁸. No error is divided by a Richardson factor, and no old failed screen is waived. A briefly drafted but unexecuted finest-pair-only variant was rejected, preserved in history and explicitly superseded. [Final continuation protocol](continuation-protocol.json).

If the complete campaign fails these screens, another refinement or numerical method needs a separately justified and costed protocol. No larger grid will run automatically. If all screens pass, whole-plane tail coverage and the native-basin interpretation still require independent justification. Neither a passing triangle nor a feasible four-million-state grid establishes protein-scale affordability. The completed harmonic calibration and measured capacity check are the defensible first milestone; nonlinear convergence is the next research question.
