# Proposed comparison of Gaussian and harmonic response with the original nonlinear triangle

**Status: proposed and not executed.** This study tests the adequacy of the physical surrogate on the existing three-site triangle. It does not test protein prediction, alter the completed Gaussian compression experiment or establish biological usefulness. The seven additional times below were already listed in the historical dynamics protocol; this document specifies how to use them for model comparison without changing the earlier results.

## The existing discrepancy concerns normalized response amplitudes

The saved Gaussian comparison has maximum absolute response discrepancy **0.40929799244468607** against the admitted finest finite nonlinear reference. This maximum is taken over all nine entries of each three-by-three response matrix at t/τ = 0.1, 1 and 10. It is expressed in the fixed harmonic-standard-deviation normalization, not as a percentage, a relative error or an error measured on a protein. The corresponding harmonic discrepancy is 0.7455994804843893.

The historical comparison verified identical geometry, edges, positive coordinate basis, stiffness, inverse temperature, mobility, original quartic residue observables, normalization and physical times. The Gaussian covariance was reused without a new fit. Its equilibrium law and dynamics differ from the nonlinear reference. Consequently, this discrepancy measures model adequacy for that triangle, rather than error from compressing the Gaussian response. The finite reference used box faces (10, 9, 21) and spacing 0.08; its whole-space dynamic error was not interval-certified.

The primary records are `representation/controls.json`, `representation/controls.npz`, `representation/compare_controls.py` and `dynamics/protocol.json` in the [published field-switch evidence archive](https://github.com/orbion-life/allosteric-response-benchmark/releases/download/v0.6.0/pulsar-field-switch-2026-09-14-public.zip). The [v0.6.0 asset manifest](../evidence-v0.6.0/ASSET-MANIFEST-v0.6.0.json) records its SHA256 as `6c6e8f490b67e5b0df4fdb0c3030ed0ebfc0dd1364686bc672fe3d87bc5ee1b0`. Its helper sources use the sibling [earlier execution archive](https://github.com/orbion-life/allosteric-response-benchmark/releases/download/v0.6.0/pulsar-execution-2026-09-13-public.zip), SHA256 `402597f6ecc300b92206633633e4a075ca18d4a9696a751bc3e651e5bcb878a1`. Extract both sibling roots and follow the [existing replay guide](../evidence-v0.6.0/README.md). These are existing records, not results from the proposed study.

## Fixed controls and additional response times

Retain the original nonlinear contact potential, all three residue observables, coordinates, model parameters and exact whole-harmonic standard deviations. Retain the fitted Gaussian covariance and native harmonic control unchanged for the primary comparison. Use the same physical time scale τ = 1/(μλmin), defined by the original positive harmonic mode, rather than rescaling each surrogate by its own relaxation time.

Evaluate both controls and the finite nonlinear reference at the seven additional times

**t/τ = 0.0075, 0.025, 0.075, 0.25, 0.75, 2.5 and 7.5.**

Report these together with the three original times, giving ten specified times and 90 response entries for each model. Preserve the original values and the full matrices, including diagonal entries. Do not select a subset of entries, receivers or times after seeing the discrepancies. No parameter refit, basis selection or proposed nonlinear remedy is part of the fixed-control evaluation.

## Qualify the finite nonlinear reference at the new times

Before interpreting surrogate error at a new time, repeat the historical domain and spacing checks at that time. The recorded hierarchy expands box faces from (8, 7, 19) to (9, 8, 20) and (10, 9, 21) at spacing 0.125, then refines spacing to 0.1 and 0.08 at the largest box. Require successive normalized static covariance, delayed covariance and response differences to be at most **0.001** for both domain expansion and spacing refinement. Preserve any failed comparison rather than silently promoting the finest calculation to an adequate reference.

Use the recorded independent propagation comparison on the same finite generator. The historical primary method is FP64 Chebyshev/Bessel propagation, with FP64 uniformization as its independent counterpart and an independent-propagation gate of **10⁻⁸**. Retain the existing stationarity, truncation-tail and small-generator action checks and their observable-norm allocations. Save method discrepancies at every additional time. A GPU kernel or propagation completion alone does not pass these gates.

The 0.001 domain and spacing screens are evidence for a reference error allowance, not a proof of a whole-space bound. Report each contribution separately and justify the combined reference allowance before making an accuracy claim. If that allowance cannot be supported, retain the result as a comparison with an explicitly identified finite operator. A failed reference screen prevents a broader original-physics fidelity claim.

## Compare amplitude and the ordering of two other sites

For each of the ten times, save every normalized response entry and the maximum absolute Gaussian–reference and harmonic–reference differences. Also save static covariance and delayed covariance comparisons so that cancellation in the response cannot hide an equilibrium mismatch.

For each receiver, compare the **two other sites** using the absolute magnitude of their response to that receiver. Report both signed responses, both magnitudes, their ordering and the reference gap between them. Preserve exact ties as ties and label a comparison unresolved when the combined error interval cannot distinguish the two amplitudes. A three-site triangle cannot support a top-five ranking, so no precision-at-five or five-site overlap metric belongs in this study. The amplitude matrix, rather than ranking agreement alone, remains the principal adequacy test.

## Qualification and the next decision

The unchanged Gaussian surrogate already fails the 0.002 amplitude criterion at the original times because its maximum discrepancy is 0.40929799244468607. Adding more query times cannot repair that result. The new comparison will identify the time dependence and type of disagreement and establish a qualified reference for evaluating a separately justified future model revision.

A revised model could support fidelity to the tested original nonlinear model only if its error meets **0.002 across every original and additional time**, with the reference and numerical allowances included inside that total. Concretely, the observed model–reference discrepancy plus the justified reference allowance and remaining numerical allowance must be no greater than 0.002. The 0.001 reference screens do not permit an extra 0.001 to be added outside this target. Finite sampling in time still does not establish fidelity at every untested time.

Any future revision requires its own derivation, recorded fitting or construction rule and comparison with these unchanged controls. This protocol does not assume that adding a basis, changing Gaussian closure or introducing another operator will succeed. Model adequacy and biological usefulness remain different tests: useful residue nominations require an independent, appropriately labelled biological evaluation, even if the nonlinear triangle comparison eventually passes.

The proposed work must fit within the existing total programme ceilings of **200 CPU-worker hours and 100 A100 hours**, shared with the other studies. These are planning limits, not a new allocation, spending authorization or measured runtime. Estimate the additional reference workload before execution, stop before exhausting the remaining allowance and retain incomplete or failed cases. No new cloud job or hardware execution is included with this document.

The deliverable is a reproducible table of all 90 response entries per model, reference-refinement and independent-propagation records, the three receivers' two-site comparisons, complete cost receipts and a stated decision on the fixed Gaussian and harmonic surrogates. This supplies evidence for choosing a physically adequate representation before making a larger quantum cost comparison.
