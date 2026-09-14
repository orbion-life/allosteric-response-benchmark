# Earlier snapshots preserve sampled Gaussian protein responses

The fixed six-node basis uses times 0, 0.0001, 0.001, 0.01, 0.1 and 1 in units of the model relaxation time. All **288 fresh response matrices** across KRAS, ABL, MYC/MAX and MYH7 have maximum entry error below 0.002. The corrected receiver analysis preserves the ordered top five and membership in **175 of 176 queries above the defined candidate signal floor**. MYC/MAX contributes 36 informative queries, all with preserved shortlists. The one ABL mismatch is retained. See [the correction and full replay record](CORRECTION-2026-09-15.md).

| Target | Retained dimension | Largest fresh response error | Numerical qualification |
|---|---:|---:|---|
| KRAS | 810 | 0.000276 | Mass orthogonality fails |
| ABL | 1,261 | 0.000182 | Pass |
| MYC/MAX | 622 | 0.001149 | Pass |
| MYH7 | 2,972 | 0.000556 | Pass |

The successful construction and fresh-kernel campaign took 390.73 seconds on one A100, with 4.61 GB peak host resident memory. These recorded-host measurements do not include historical kernel generation. KRAS mass error is 9.25 × 10⁻⁷ against a fixed 10⁻⁸ gate. This failure remains even though its sampled response errors pass. Sampled Gaussian fidelity does not establish all-time preservation or original nonlinear fidelity.

## Reproduce response analysis without a GPU

First download the corresponding `pulsar-temporal-TARGET-evidence.zip` for each desired target from **v0.9.0** and verify its SHA256. Extract those original raw archives, then extract **`pulsar-temporal-receiver-correction-v0.10.1.zip` from v0.10.1 into the same parent directory**. The correction overlay supplies the corrected analyzer and superseding ranking results; apply it before running the commands below. Each original target package supplies the common source and protocol, its input model, complete operator, 72 fresh kernels and older kernels for historical comparisons. **MYH7 requires both `pulsar-temporal-myh7-evidence.zip` and `pulsar-temporal-myh7-reference.zip`** because its complete evidence exceeds a single release asset. Extract target archives into the same empty parent directory; identical common files may be replaced. From `temporal` run:

```bash
python -m pip install -r requirements.txt
python analyze.py kras --output ../kras-corrected-replay
```

Repeat with `abl`, `myc` or `myh7` after extracting its archive. Analysis uses `reference/` beside the script; `PULSAR_REFERENCE_ROOT` may point to an equivalent separately stored reference tree. It validates metadata, kernel and operator hashes before use and writes the target's derived `analysis.json` to the specified separate output. Existing output analyses are rejected. Timings are new measurements; scientific arrays and checks should be compared numerically.

The Git checkout provides source and compact results; the large input/operator/kernel arrays are release assets. `campaign.py` is the recorded GPU construction source and requires CuPy plus an explicitly provisioned CUDA GPU. It is never launched by saved-data analysis. No account-specific launcher or credentials are distributed. The initial infrastructure failure and analysis reporting correction are preserved under `history/` and in `analysis-reporting-corrections.json`.

`PACKAGING.json` records path adaptation, provider-ID removal and this ranking correction. Raw scientific inputs retain their original bytes. Superseded derived results are preserved under `history/receiver-ranking-before-2026-09-15/`. Download the v0.10.1 correction overlay after the original raw archives, as explained in the correction record. `DATA-AND-SOFTWARE-ATTRIBUTION.md` identifies inherited structure/model sources.

For a separately provisioned GPU replay, install CuPy 14.1.1 with the CUDA libraries appropriate to the machine, and use the recorded NumPy 2.2.6, SciPy 1.15.3 and Python 3.13.3 environment. From this directory, an explicit invocation such as `python -c "from campaign import run; run('kras', 'results-new/kras')"` constructs a new target record. It refuses an existing output directory. This optional command runs a GPU calculation; the saved-data commands above do not.
