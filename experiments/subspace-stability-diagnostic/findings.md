# The reduced-basis recurrence amplifies small implementation differences in a selected finite model

**Project Pulsar · 13 September 2026 · Completed local diagnostic.** This report describes the two predeclared cases and retains all failed attempts. The term “recurrence amplification” is an interpretation of controlled comparisons, not attribution to a uniquely identified arithmetic instruction.

## Question and controlled comparisons

The earlier experiment produced response matrices that reproduced more closely than its calculated response-error bounds. Its archives did not contain the full reduced basis or the incoming frontier, so they could not distinguish a change of reduced coordinates from a change of subspace. This diagnostic saved the native finite generator H, observable matrix F, full thin basis V and every admitted residual block. It then compared native inputs, shared identical H/F, identical V, and one identical residual block W in two available local environments.

The scientific statements in the original source remained unchanged. The instrumented and uninstrumented native builds were bitwise equal in both environments; the corresponding control bases were also bitwise equal. The corrected logger required every anticipated event to occur. All six complete trace sets passed the saved-file checks in `inspection.json`.

## The native input differs before the basis is constructed

For the elongated κ = 100, 33³ case, the two native Hessians and CSR sparsity patterns were byte equal. The retained eigenvalues differed by at most 5.684 × 10⁻¹⁴, the retained frame by 4.441 × 10⁻¹⁶, and the harmonic scales by 2.220 × 10⁻¹⁶. The resulting generator entries differed by at most 4.529 × 10⁻¹⁰ and F by 9.766 × 10⁻¹⁶. These values establish that the native matrices were not identical. They do not establish a physically meaningful difference between those native models.

The shared-input comparison removed this ambiguity by supplying the same pinned H, F and times to both environments. At the final rank 48, the maximum response difference was 2.090 × 10⁻¹¹ and the maximum calculated-bound difference was 3.213 × 10⁻⁷. The full-subspace gap was 0.013446. The corresponding native-input gap was 0.154560 and its bound difference was 2.094 × 10⁻⁶. Input construction therefore contributes to the observed native discrepancy, while identical finite matrices do not eliminate it.

## The admitted subspace diverges without a marginal rank decision

In the elongated shared-input trace, the first measured complete-basis gap above 10⁻¹⁰ occurred at admission step 3, when the rank became 12; that gap was 2.532 × 10⁻¹⁰. The first complete-basis gap above 10⁻⁶ occurred at step 8, rank 27, and was 2.120 × 10⁻⁶. The final gap was 0.013446 at rank 48. The saved stage-by-stage comparison also identifies earlier threshold crossings in the raw and orthogonalized residual blocks. These crossings are observations at specified numerical thresholds, not a claim that all earlier arithmetic was identical.

Every incoming residual block in all six traces had three eligible columns, and every admission retained those three columns. There were no partial-block prefix cuts. The smallest residual singular-value-to-threshold ratio was approximately 1.288 × 10¹¹ in the elongated case and 3.532 × 10¹¹ in the control. The effective-rank diagnostic did not flag any compared span as rank deficient. Thus neither a singular value crossing the dependency threshold nor a checkpoint splitting a residual singular cluster explains these selected runs.

The development κ = 1, 17³ control had a much smaller final gap: 1.340 × 10⁻⁸ with native inputs and 1.260 × 10⁻⁹ with shared input. Its final response differences were below 2.5 × 10⁻¹⁵ and its final bound differences below 2.2 × 10⁻¹². These observations show that the consequential amplification is case dependent within the tested pair. They do not estimate its frequency across other geometries or grids.

## Downstream diagonalization does not account for the large bound discrepancy when the basis is fixed

For every pinned checkpoint, the diagnostic supplied the same H, F and V to each environment and evaluated the projected problem using three symmetric eigensolver drivers. Across all checkpoints and drivers, the largest elongated-case response difference was 1.221 × 10⁻¹⁵ and the largest bound difference was 7.136 × 10⁻¹². The corresponding control maxima were 9.437 × 10⁻¹⁶ and 9.721 × 10⁻¹². Each evaluation reconstructed its actual propagated operator and corresponding residual before calculating the bound.

These comparisons strongly narrow the explanation for the larger 10⁻⁶ to 10⁻⁷ discrepancies in the elongated run: those discrepancies are associated with different constructed subspaces and are not reproduced by downstream eigensolver or bound evaluation on the same basis. This inference remains conditional on the tested environments, drivers, checkpoints and precision.

## One isolated SVD is insufficient to identify a repair

The prespecified rule selected the first residual block whose shared-input eligible-span gap exceeded 10⁻¹⁰. The selected step was 3 in the elongated case and 10 in the control. When exactly the same saved W and threshold were supplied to `gesdd` and `gesvd` in each environment, all keep masks agreed. The maximum admitted-subspace gap was 1.666 × 10⁻¹⁴ for the elongated case and 1.803 × 10⁻¹⁵ for the control.

The isolated operation therefore did not reproduce the already amplified discrepancy between the two separately evolved residual blocks. This result is consistent with accumulated rounding in the sequence of sparse actions, projection subtractions and regenerated frontiers. It does not distinguish the contribution of each operation or prove that changing one driver would stabilize the full recurrence. No replacement recurrence, higher-precision model or forward continuation was run.

## What remains unresolved

The diagnostic reproduced the previously reported local elongated-case discrepancy and exposed the previously unavailable full-subspace evidence. It did not reproduce a Linux trace, so it cannot establish the exact cause of every historical Linux discrepancy. Both local environments used the same host and both reported `longdouble` precision equal to float64; this run provides no extended-precision comparison. Active thread-pool introspection found the system OpenBLAS pool at one thread, whereas the pinned environment reported no discoverable pool. All six thread-control environment variables were set to one before numerical imports.

A future implementation study could compare a deterministic polynomial/block recurrence or a reproducibly reduced action against the original recurrence on the same accepted H/F and response queries. That study should retain the original failure, keep the original scientific accuracy gates, and declare any changed arithmetic or construction rule as a new implementation. The current diagnostic provides the trace fixtures needed for such a test; it does not select or validate a remedy.

The separately identified domain and mesh-convergence failures remain outside this diagnostic. No accepted physical continuum limit, biological allostery prediction or useful quantum cost has been established by these runs. A matched quantum study should begin only after accepting the physical model and should compare the original accepted H against both a sparse classical reference and a validated reduced operator, including all preparation and reconstruction costs.

## Evidence ledger

| Claim | Direct evidence | Permitted interpretation |
|---|---|---|
| Native H/F differ between the two environments. | `results/comparison/*-input-comparison.json` and saved native `input.npz` files. | Native-input equality must not be assumed. |
| Identical H/F still produce different full subspaces. | Shared-input evaluation hashes, full-V comparisons, seed and admitted blocks. | Construction differs even when finite input matrices are fixed. |
| No threshold straddling or partial-block cut occurred. | All repaired residual SVD, keep-mask and admission files; `inspection.json`. | Those mechanisms do not explain these selected traces. |
| Fixed-V downstream differences are much smaller. | `results/*/*/fixed-V-drivers/` and the summary's driver comparisons. | The large observed discrepancy is associated with construction rather than identical-basis downstream evaluation. |
| The isolated same-W SVD does not reproduce the amplified gap. | Saved isolation input and all four driver/environment results per case. | One isolated operation does not identify a validated remedy. |
| Instrumentation preserves the recorded numerical outputs. | Build and control-basis parity receipts; 688 repaired prior-file parity checks. | The corrected logger changes trace coverage, not the recorded scientific results. |
| A separately recomputed reference can differ below 10⁻¹⁵. | Historical failed inspection and corrected reference-difference fields. | Bitwise reference equality was an unsupported auxiliary assumption and is not retained. |
