# Receiver indices correct the ancillary sampled-peak diagnostic

**Correction issued on 2026-09-14. All original peak diagnostics are withdrawn.** The sealed response calculations and original results remain available unchanged. This is an explicit correction to a diagnostic defect, not a revised physical model or a new temporal holdout.

In `analyze.py`, the expression `candidate[:, None] & receiver[None, :]` combined a boolean candidate mask with integer receiver indices. It selected positions according to integer parity and then used those positions as matrix columns. The intended operation selects each eligible candidate and each actual receiver index. The event-aligned analyzer imports the same function and inherited the defect.

The affected fields are the `peaks` arrays in all 28 derived waveform archives, all endpoint `peak` objects in the original and event-aligned `analysis.json` records, and their copied `sampled_peak` fields in both summary files. Every count, amplitude, index, timing difference and match fraction in those original peak fields must be treated as withdrawn. The original sealed replay intentionally reproduces them for provenance. It does not apply this correction silently.

`correct_peaks.py` implements explicit candidate-by-actual-receiver iteration without importing the original diagnostic function. A fixture with nonconsecutive odd/even receiver indices demonstrates the original defect and checks the corrected pair selection. An independent `np.take` extraction confirms the corrected peak indices for every endpoint. The complete corrected pair records and small corrected summary are in `results/`.

All 28 endpoints were independently rechecked. Every full-matrix response-error array and gate, signed-entry count, reference/reduced RMS ranking, candidate-to-receiver error array, reference peak and relative Frobenius value remains unchanged within the verification tolerance of 1e-12. The qualification audit already used actual receiver indexing: `full[:, candidate][:, :, receiver]`. For the event-aligned 3τ pulses, the candidate-to-receiver maximum errors remain 0.00885329679531905 (KRAS), 0.010276286229872065 (ABL), 0.0035158548491814145 (MYC) and 0.013978188139825138 (MYH7). All twelve event-aligned pulse curves still fail the full-matrix 0.002 criterion.

The corrected maxima remain **unique interior sampled maxima under the fixed 1e-12 tie rule**. The 0.002 amplitude floor and all sampling times are unchanged. These maxima are not robust timing estimates or continuous-time certificates; plateaus can produce large changes of the selected time under small amplitude errors. No biological delay or quantum scrambling is established.

From a complete temporal evidence tree containing the original derived waveforms, run:

```sh
python corrections/receiver-peak-indexing/correct_peaks.py --output NEW_CORRECTION_DIRECTORY
```

For the compact archive, first reconstruct the omitted historical waveforms, then apply the explicit correction:

```sh
python replay_waveforms.py --output NEW_PROTEIN_REPLAY_DIRECTORY
python corrections/receiver-peak-indexing/correct_peaks.py --root NEW_PROTEIN_REPLAY_DIRECTORY --output NEW_CORRECTION_DIRECTORY
```

Neither command refits a protein, changes a basis or calls a GPU launcher. `historical/` preserves the documentation and manifests from before this correction. The source hashes in both original seals and all original scientific arrays remain unchanged. The original saved-array audits checked responses, identities and provenance; they did not independently verify the erroneous peak selector. Their passing status must not be described as validation of the withdrawn peak metrics.
