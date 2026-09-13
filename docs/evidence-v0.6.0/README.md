# Reviewing the Gaussian, quantum and nonlinear response evidence

The v0.6.0 evidence assets extend the [preceding v0.5.0 release](https://github.com/orbion-life/allosteric-response-benchmark/releases/tag/v0.5.0). Download the [earlier scientific package](https://github.com/orbion-life/allosteric-response-benchmark/releases/download/v0.6.0/pulsar-execution-2026-09-13-public.zip) and the [field-switch and functional package](https://github.com/orbion-life/allosteric-response-benchmark/releases/download/v0.6.0/pulsar-field-switch-2026-09-14-public.zip), then check their hashes against [ASSET-MANIFEST-v0.6.0.json](ASSET-MANIFEST-v0.6.0.json). The full asset list also includes four new particle overlays and all six older stability source/input/trace archives. The proposal, applicant details and internal reviews are not public assets.

The scientific extraction preserves two sibling directories:

```
pulsar-phase1-execution-2026-09-13/
pulsar-field-switch-execution-2026-09-13/
```

Below, OLD means the first directory and NEW the second. The final release manifest, asset sizes and SHA-256 digests determine which files have actually been delivered. The 32 full SMC particle files are distributed in four separate ZIP assets; the selected higher-moment audit uses compact raw-moment files instead.

## Read the evidence before running verification

| Question | First files to inspect | Supported conclusion |
|---|---|---|
| What is computed at protein size? | OLD `gaussian/README.md`, `gaussian/four-target-compute-summary.json` and target result folders | Complete quartic-observable response matrices under the fitted Gaussian/OU surrogate; not original nonlinear protein fidelity. |
| Does the representation preserve response? | OLD `gaussian-galerkin/summary.json`, `gaussian-galerkin/temporal-holdout/summary.json`, `gaussian-galerkin-recovery/confirmation/summary.json` | Original three-time agreement and later early-time failures are separate results. |
| Was the quantum estimator measured? | OLD `quantum/gaussian-rank9/receipt.json`, `quantum/gaussian-rank9/actual-shots.json`, `quantum/portable/README.md` | Ideal simulator shots for one pair at three times; matrix-walk checks for 18 pair/time queries. No QPU or full protein quantum output. |
| Is the nonlinear reference admitted? | NEW `dynamics/reference-summary.json`, `static/summary.json`, `static/tail-analysis.json` | Finite-grid/domain and independent propagation checks pass; whole-space dynamics are not certified. |
| Does intervention agree with the formula? | NEW `dynamics/remote/direct-fields-finest/data/field/receipt.json` and neighboring raw arrays | Direct computational field-switch agreement in the triangle model, not a mutation or drug experiment. |
| Did the nonlinear representation succeed? | NEW `representation/summary.json`, `representation/rank-diagnostic/summary.json`, `representation/controls.json` | Rank-nine and full-span Hermite failures; separate Gaussian/harmonic model discrepancies. |
| Is original-density sampling admitted? | NEW `smc/precision-continuation/README.md`, `smc/precision-continuation/downloaded/production/summary.json`, `smc/precision-continuation/verification.json` | A separate finite-box precision continuation passes; the earlier failure and cross-kernel partition limitation remain. |
| Were higher moments independently checked? | NEW `generator-moment-audit/README.md`, `generator-moment-audit/summary.json`, `generator-moment-audit/comparison.csv` | Fifteen continuum differential/gradient quantities pass their separate diagnostic; not a finite-grid response-residual certificate. |
| Is there biological evidence? | NEW `biology/README.md`, `biology/PROTOCOL.md`, `biology/results/summary.json`, `biology/results/common-mask-comparison.csv` | A retrospective KRAS association, with competitive structural controls and unestablished incremental value. |

The preceding public v0.5.0 release contains earlier harmonic recovery, nonlinear cost probes and subspace-stability diagnostics. The new Gaussian, quantum, field-switch, sampling and functional packages should not be attributed to that tag. The v0.5.0 page records the earlier incomplete raw upload; v0.6.0 includes those previously unpublished files as separate, hashed assets while preserving the earlier status record.

## Verify saved evidence without repeating the campaigns

Use a **disposable extraction**. Several original scripts write fresh verification JSON or comparison tables beside their inputs. Keep downloaded archives and their manifest unchanged. The curated generator-moment and biological saved-result audits passed in a relocated copy with all original-project reads blocked. Other command targets and arguments were inspected; this preparation does not claim a fresh full campaign replay.

Run from the directory containing both sibling roots. Use Python with NumPy and SciPy; the biological audit also requires Biopython. Numerical environments differ across the recorded experiments: Python 3.9.6/NumPy 2.0.2/SciPy 1.13.1 and Python 3.13.5/NumPy 2.4.6/SciPy 1.16.0. A new environment is not a promise of identical cross-version output.

```sh
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
PULSAR_OLD=pulsar-phase1-execution-2026-09-13
PULSAR_NEW=pulsar-field-switch-execution-2026-09-13
PULSAR_REPLAY=replay-checks

python "$PULSAR_OLD/gaussian/check_independent.py" --output "$PULSAR_REPLAY/gaussian-moments"
python "$PULSAR_NEW/static/verify_static.py"
python "$PULSAR_NEW/dynamics/check_small.py" --output "$PULSAR_REPLAY/dynamics-small"
python "$PULSAR_NEW/representation/verify_and_summarize.py"
python "$PULSAR_NEW/generator-moment-audit/verify.py"
python "$PULSAR_NEW/biology/scripts/audit.py"
```

The Gaussian check tests moment identities and derivatives on bounded fixtures; it does not recompute four proteins. The dynamics check uses a small finite generator, not the 29.5-million-state grid. The representation verifier should preserve `FAIL_PRESCRIBED_HERMITE_BRANCH` while its implementation checks pass. The higher-moment verifier reads the 32 saved `raw-moments.npz`/`result.json` pairs. The biological audit reconstructs assay aggregation, receiver/readout, metrics and saved matched-null results; it does not create independent biological validation.

For quantum saved-data verification, add Qiskit and psutil; recorded versions include Qiskit 2.3.1 and psutil 7.0.0. The adapter creates a separate output directory and refuses an existing one:

```sh
python "$PULSAR_OLD/quantum/portable/replay.py" --mode audit --output "$PULSAR_REPLAY/quantum-audit"
```

It reconstructs saved counts and performs bounded component checks. It does not simulate millions of new shots or replay every complete circuit. `quantum/portable/README.md` documents optional full simulation modes, the Qiskit Aer dependency and resource caps.

If all full-particle assets are present:

```sh
python "$PULSAR_NEW/smc/precision-continuation/verify_downloads.py"
```

This reads approximately 7.27 GB of production archives and checks every saved final particle with independent contact geometry. The recorded verification took 43.28 seconds and about 563 MB peak process memory; these are host measurements, not runtime guarantees. It does not resample or allocate a GPU.

## Interpret each verdict separately

Check hashes and source/protocol versions first, then computation completion, numerical verification and scientific acceptance. A successful execution may correctly reproduce a scientific failure. Keep static covariance, delayed response, original-model discrepancy, rank uncertainty and biological value separate. Include construction, preparation, routing, shots and full readout in quantum/classical comparisons.

The response threshold is maximum absolute entry error 0.002 in fixed harmonic-standard-deviation units, not a relative percentage. Reference refinement has a separate 0.001 screen. The higher-moment diagnostic divides error, nominal uncertainty and reference discrepancy by `max(1, |reference|)` and compares with 0.01; it is not a universal relative-error bound. Neither SMC intervals nor finite-grid checks certify whole-space dynamics.

No new cloud calculation is needed to inspect these results. Full GPU replay is a separate resource-bearing operation that must follow the retained canary, source seal and bounded protocol.

## Public curation and third-party terms

Each scientific archive contains PUBLIC-CURATION.json and PUBLIC-FILE-MANIFEST.json. Historical execution manifests may name omitted internal/unused files or earlier verification receipts; they are retained as history. The functional archive explains its three omissions and updated verification-receipt hash in biology/RELEASE-CURATION-NOTE.md. Numerical sources and result arrays retained from the original compact bundles are unchanged. The data/software attribution files preserve PDB, UniProt, SIFTS, MaveDB, Lehner MIT and Ohm GPL terms. The recorded Ohm binary is macOS ARM; its corresponding source, license, patch and build records accompany it.
