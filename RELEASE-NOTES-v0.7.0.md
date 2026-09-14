# Version 0.7.0

This release completes the GRB2 functional, expanded temporal and quantum cost experiments proposed in the preceding document revision. It also publishes the prior 14 September feasibility supplement and an updated nine-page proposal with its editable source.

## New findings

The GRB2 analysis uses 21 distal sites with at least ten eligible substitutions per site. Gaussian correlation with mean absolute inferred binding free-energy change is 0.3805, compared with 0.3351 for harmonic response and 0.4806 for degree. The Gaussian-minus-harmonic paired 95% interval is [-0.0686, 0.2090]; Gaussian-minus-degree is [-0.4569, 0.2445]. All four multiplicity-controlled 98.75% comparison intervals include zero. These data do not establish incremental usefulness.

The fixed four-protein projections pass their three specified response queries but fail broader step and pulse fidelity. A separately frozen addendum samples immediately after field removal and confirms failure for all twelve duration/target curves. Secondary distal-to-receiver errors also exceed the threshold. The tiny finite-box field-switch identity passes; it does not establish protein dynamics or physical scrambling.

The primary SELECT–SWAP experiment reduces calculated T counts by 43.7%, missing its 50% target. A separate phase-gradient follow-up reduces T counts by 69.2% against the strongest evaluated dense circuit, with preparation charged once per complete shot. The full 18-query ideal-sampling plan requires 30,383,928 sufficient readouts and about 1.0846e13 T gates. It uses 182 logical qubits and 1.0290e13 CX gates, 106.5 times the dense CX total. Low-order complete circuits and full-size components are exported. Maximum-order complete circuits exceed the export cap. No new QPU shots, physical fault-tolerance accounting, or classical advantage is claimed.

## Reproducibility and scope

The release includes frozen protocols, inputs, raw data, actual compiled circuits, independent checks, failures, source and replay instructions. All eight new full temporal kernels are retained. Twenty-eight redundant waveform archives are omitted with their hashes and a reconstruction command; a fresh local tree regenerated identical archive hashes and passed both saved-result audits.

The repository adds three scoped Linux workflows for GRB2 saved evidence, temporal protocol/tiny-model checks and quantum source semantics. Large GPU calculations and the complete quantum cost synthesis are separate replays described in the assets. Original workflows and failures are retained. Current CI status should be read by workflow and commit, not inferred from this release description.

The final proposal distinguishes a feasible Gaussian classical delivery baseline, a measured small quantum prototype, a scoped logical T improvement and the proposed Phase 2 biological/readout tests. Contributor commitments, formal eligibility and submission confirmations remain outstanding. Nothing has been submitted to the challenge.
