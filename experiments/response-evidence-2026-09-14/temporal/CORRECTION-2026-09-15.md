# Receiver ranking correction, 15 September 2026

The earlier-time basis preserves the ordered top five in **175 of 176 fresh queries above the candidate signal floor**, across four fixed Gaussian protein models. It preserves membership in the same 175 queries. All **288 fresh full response matrices** still pass the 0.002 maximum entry-error threshold. The receiver correction changes the ranking analysis; it does not change the fitted models, saved kernels, reduced operators, construction protocol or raw numerical qualification.

| Target | Informative fresh queries | Ordered top five preserved | Top-five membership preserved |
|---|---:|---:|---:|
| KRAS | 47 | 47 | 47 |
| ABL | 47 | 46 | 46 |
| MYC/MAX | 36 | 36 | 36 |
| MYH7 | 46 | 46 | 46 |
| Total | 176 | 175 | 175 |

An informative query has a maximum eligible candidate receiver-RMS score strictly greater than 0.002. A receiver-RMS score is the root mean square of the signed response over the specified receiver columns. Eligible candidates are ordered by descending score, with ascending canonical residue number breaking exact ties. The 0.002 signal floor, candidate set and tie rule are unchanged.

The single informative mismatch occurs for ABL after a 0.7τ perturbation, at a post-removal lag of 0.006819311811870924τ. The reference shortlist is [367, 267, 303, 298, 361]; the reduced shortlist is [367, 267, 303, 298, 385]. The reference fifth/sixth gap is 3.306259362297059 × 10⁻⁸ and the maximum candidate-score error is 4.107838192937523 × 10⁻⁶. A matrix-error pass therefore does not guarantee the exact shortlist. These are model-time queries, not independent biological replicates or measured molecular delays.

## What failed and what was corrected

The original `analyze.py` passed an integer receiver-index array through `np.flatnonzero`. That operation returns positions of nonzero values, rather than the values themselves. For example, receiver indices [0, 3, 6] were interpreted as [1, 2]. The saved candidate array is a boolean mask and the saved receiver array is an integer index list. The corrected code declares and validates those encodings separately. It rejects boolean receivers, integer candidate masks, duplicate indices, invalid shapes, empty selections and out-of-range indices.

The previous independent verifier repeated the same receiver conversion. Its passing ranking check was not an independent validation of receiver semantics. The historical analysis, summary, README, packaging records, manifest files and replay receipt remain byte-for-byte under `history/receiver-ranking-before-2026-09-15/`. The original release archives remain historical evidence. Their claim of 148/148 informative rankings and absence of informative MYC/MAX queries is superseded.

The new `verify_correction.py` does not import the production analyzer or scoring helper. It independently diagonalizes the saved reduced generator, gathers the actual receiver columns with `np.take`, computes scores with a vector norm and sorts candidates with Python tuple ordering. It checks all 2,216 current query records across fresh, selected and historical grids. It verifies that every full-matrix error is unchanged and verifies exact shortlists above the signal floor; low-signal numerical ties are not used as shortlist evidence. `correction-verification/summary.json` records the checks and exact input hashes.

The preceding selection experiment used the integer receiver list directly for RMS scores. `verify_selection_scope.py` checked all 28 saved waveform archives, comprising 1,916 query records, and confirmed their complete candidate rankings. Its separately withdrawn ancillary peak fields remain withdrawn. The earlier basis's temporal failures and the historical bound-portability failure remain unchanged. KRAS still fails the 10⁻⁸ mass-orthogonality gate despite passing the sampled response criterion. No new GPU, QPU, model fitting, biological validation or quantum-advantage experiment was performed.

## Portable saved-data replay

Use Python 3.11 or later. Install `requirements.txt` for the released analysis environment. Extract the original v0.9.0 target evidence archives into one parent directory, then extract `pulsar-temporal-receiver-correction-v0.10.1.zip` into that same parent. All members begin with `temporal/`. Keep the original archive files for provenance. The correction replaces the active source and derived summaries and adds a historical copy of the superseded files. MYH7 still requires both its evidence and reference archives. The correction archive does not duplicate the large raw arrays.

From the resulting `temporal` directory:

```bash
python -m pip install -r requirements.txt
python verify_package.py CORRECTION-MANIFEST.json
python verify_package.py PUBLIC-MANIFEST-kras.json
python -W error -m unittest discover -s . -p test_receiver_ranking.py
python replay_correction.py --output ../corrected-replay
python verify_correction.py --analysis-root ../corrected-replay --output ../corrected-verification
```

The all-target replay needs all four target archives. The public manifest verifies the matching target plus common files; run it for each extracted target. The MYH7 reference manifest covers its separate part. `CORRECTION-MANIFEST.json` verifies the complete correction overlay without requiring any large raw array. Its input provenance lists the original metadata, operator and kernel hashes. Neither replay command writes into the source data tree. Existing output analyses are rejected by both the CLI and callable analysis function.

For one extracted target, use `python analyze.py kras --output ../kras-corrected-replay` and replace `kras` as appropriate. Separate data and reference trees are supported through `--data-root` and `--reference-root`. `PULSAR_REFERENCE_ROOT` remains supported by the analyzer.

The repository CI runs synthetic adversarial selector tests and actual receiver-score fixtures at τ for all four proteins. The fixtures contain signed reference/reduced matrices and original selector metadata; each corrected shortlist differs from the erroneous receiver interpretation. These bounded CI checks do not rerun all protein kernels. Full local replay and independent verification are recorded separately.
