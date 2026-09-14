# Earlier snapshots preserve sampled Gaussian protein responses

The fixed six-node basis uses times 0, 0.0001, 0.001, 0.01, 0.1 and 1 in units of the model relaxation time. All **288 fresh response matrices** across KRAS, ABL, MYC/MAX and MYH7 have maximum entry error below 0.002. All **148 queries above the defined signal floor** preserve the ordered top five. MYC/MAX has no fresh query above that floor.

| Target | Retained dimension | Largest fresh response error | Numerical qualification |
|---|---:|---:|---|
| KRAS | 810 | 0.000276 | Mass orthogonality fails |
| ABL | 1,261 | 0.000182 | Pass |
| MYC/MAX | 622 | 0.001149 | Pass |
| MYH7 | 2,972 | 0.000556 | Pass |

The successful construction and fresh-kernel campaign took 390.73 seconds on one A100, with 4.61 GB peak host resident memory. These recorded-host measurements do not include historical kernel generation. KRAS mass error is 9.25 × 10⁻⁷ against a fixed 10⁻⁸ gate. This failure remains even though its sampled response errors pass. Sampled Gaussian fidelity does not establish all-time preservation or original nonlinear fidelity.

## Reproduce response analysis without a GPU

Download the corresponding `pulsar-temporal-TARGET-evidence.zip` for each desired target from the release and verify its SHA256. Each target package supplies the common source and protocol, its input model, complete operator, 72 fresh kernels and older kernels for historical comparisons. **MYH7 requires both `pulsar-temporal-myh7-evidence.zip` and `pulsar-temporal-myh7-reference.zip`** because its complete evidence exceeds a single release asset. Extract target archives into the same empty parent directory; identical common files may be replaced. From `temporal` run:

```bash
python -m pip install -r requirements.txt
python analyze.py kras
```

Repeat with `abl`, `myc` or `myh7` after extracting its archive. Analysis uses `reference/` beside the script; `PULSAR_REFERENCE_ROOT` may point to an equivalent separately stored reference tree. It validates kernel and operator hashes before use and replaces only that target's derived `analysis.json`. Record a copy first if byte preservation is required. Timings are new measurements; scientific arrays and checks should be compared numerically.

The Git checkout provides source and compact results; the large input/operator/kernel arrays are release assets. `campaign.py` is the recorded GPU construction source and requires CuPy plus an explicitly provisioned CUDA GPU. It is never launched by saved-data analysis. No account-specific launcher or credentials are distributed. The initial infrastructure failure and analysis reporting correction are preserved under `history/` and in `analysis-reporting-corrections.json`.

`PACKAGING.json` records path adaptation and provider-ID removal. Scientific numeric inputs and outputs retain their original bytes. `DATA-AND-SOFTWARE-ATTRIBUTION.md` identifies inherited structure/model sources.

For a separately provisioned GPU replay, install CuPy 14.1.1 with the CUDA libraries appropriate to the machine, and use the recorded NumPy 2.2.6, SciPy 1.15.3 and Python 3.13.3 environment. From this directory, an explicit invocation such as `python -c "from campaign import run; run('kras', 'results-new/kras')"` constructs a new target record. It refuses an existing output directory. This optional command runs a GPU calculation; the saved-data commands above do not.
