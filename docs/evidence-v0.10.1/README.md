# Corrected receiver rankings and reproducible response calculations

Version 0.10.1 corrects the receiver-column selection used to score temporal responses. The fixed Gaussian models, full kernel arrays and measured matrix-error results are unchanged. The report and active analysis now use the actual integer receiver indices.

| Target | Fresh matrices within 0.002 | Informative fresh queries | Complete top-five order preserved |
|---|---:|---:|---:|
| KRAS | 72/72 | 47 | 47 |
| ABL | 72/72 | 47 | 46 |
| MYC/MAX | 72/72 | 36 | 36 |
| MYH7 | 72/72 | 46 | 46 |
| Total | 288/288 | 176 | 175 |

Informative means that the maximum eligible candidate receiver-RMS score exceeds the fixed 0.002 floor. One ABL query replaces residue 361 with 385 at the fifth position, where the reference fifth/sixth gap is 3.31 × 10⁻⁸. These queries are model-time evaluations, not independent biological replicates. KRAS still fails the original mass-orthogonality gate.

## Inspect the source and verification

- [Correction explanation and replay commands](../../experiments/response-evidence-2026-09-14/temporal/CORRECTION-2026-09-15.md).
- [Explicit mask and integer-index validation](../../experiments/response-evidence-2026-09-14/temporal/ranking.py) and [six regression tests](../../experiments/response-evidence-2026-09-14/temporal/test_receiver_ranking.py), including saved data for every target.
- [Independent verifier](../../experiments/response-evidence-2026-09-14/temporal/verify_correction.py) and [2,216-query verification receipt](../../experiments/response-evidence-2026-09-14/temporal/correction-verification/summary.json).
- [Scope check for 1,916 earlier selection queries](../../experiments/response-evidence-2026-09-14/temporal/selection-scope-verification.json).
- [Proposed centroid, coherent-loading and independent-family protocols](bounded-research-plan.md).

The verifier separately diagonalizes the saved reduced generator, selects receiver columns explicitly and checks the scores and above-floor shortlists. All current full-matrix errors agree with their historical values. The old analyzer and verifier, which shared the indexing error, remain in the history directory with their superseded results. A software-check pass does not change the preserved nonlinear, biological, hardware or bound-portability findings.

## Download and replay

The [v0.10.1 release](https://github.com/orbion-life/allosteric-response-benchmark/releases/tag/v0.10.1) contains the corrected report, editable source and `pulsar-temporal-receiver-correction-v0.10.1.zip`. The [asset manifest](ASSET-MANIFEST-v0.10.1.json) gives sizes and SHA256 hashes.

First extract the original [v0.9.0 temporal archives](https://github.com/orbion-life/allosteric-response-benchmark/releases/tag/v0.9.0) into one parent directory. MYH7 needs both its evidence and reference parts. Then extract the v0.10.1 correction archive into the same parent. It replaces active analysis and derived results while preserving the old versions; it does not duplicate the large numerical arrays.

From the resulting `temporal` directory, using Python 3.11 or later:

```bash
python -m pip install -r requirements.txt
python verify_package.py CORRECTION-MANIFEST.json
python -W error -m unittest discover -s . -p test_receiver_ranking.py
python replay_correction.py --output ../corrected-replay
python verify_correction.py --analysis-root ../corrected-replay --output ../corrected-verification
```

Run each extracted target's `PUBLIC-MANIFEST` check as described in the correction record. The saved-data replay constructs no new protein model and launches no GPU or QPU work. Output analyses must use a fresh destination.

The nonlinear and physical-factor circuit arrays remain in [v0.10.0](https://github.com/orbion-life/allosteric-response-benchmark/releases/tag/v0.10.0), and the biological and six-qubit records remain in [v0.9.0](https://github.com/orbion-life/allosteric-response-benchmark/releases/tag/v0.9.0). Original source licences and data acquisition requirements still apply.

## Reproduce the report

Extract `project-pulsar-editable-source.zip` and follow its README. Python and Tectonic rebuild the report without rerunning experiments. The archive includes its figure inputs, per-file manifest, format receipt and relocated-rebuild check. All nine pages were checked for identical text and rendered pixels in the recorded environment.
