# Gaussian and harmonic controls fail nonlinear response fidelity across ten fixed times

The completed comparison evaluates the unchanged Gaussian and harmonic controls against the original three-site nonlinear model at ten specified times. Five finite grids pass the recorded domain, spacing, stationarity and independent-propagation screens. Neither control meets the normalized 0.002 response criterion in any of its 90 entries. Gaussian maximum discrepancy remains 0.4092979924; harmonic maximum discrepancy is 0.7464356610.

The results distinguish reliable calculation of the chosen finite reference from adequacy of the Gaussian physical surrogate. They do not constitute protein prediction, biological validation or a whole-space dynamic certificate.

| Quantity | Gaussian | Harmonic |
|---|---:|---:|
| Largest normalized amplitude discrepancy | 0.409298 | 0.746436 |
| Entries within 0.002, of 90 | 0 | 0 |
| Raw two-other-site ordering agreements, of 30 | 16 | 16 |
| Raw ordering reversals, of 30 | 14 | 14 |
| Agreements after empirical error screening | 16 | 16 |
| Clear reversals after empirical error screening | 12 | 12 |
| Unresolved orderings after empirical error screening | 2 | 2 |

Each ordering compares the two sites other than a specified receiver. A three-site model cannot provide a top-five ranking. The screened counts are conditional on an empirical envelope formed from observed domain and spacing differences, independent-method disagreement and propagation tail allowances. This envelope is useful for sensitivity analysis but is not a certified bound on the infinite-domain solution. No fraction in this table is a biological success rate.

## Frozen experiment and reference qualification

`protocol.json` and `protocol-freeze.json` were written before new calculations. They preserve the original times 0.1, 1 and 10 tau and add 0.0075, 0.025, 0.075, 0.25, 0.75, 2.5 and 7.5 tau, all already listed in the historical dynamics protocol. Geometry, original quartic residue energies, beta, mobility, stiffness, Gaussian covariance and harmonic normalization remain fixed. Beta, mobility and stiffness are all one in this fixture; the control implementation must not be generalized silently to other values.

The five grids use faces (8,7,19), (9,8,20) and (10,9,21) at spacing 0.125, followed by spacings 0.1 and 0.08 on the largest box. Their maximum successive G0, K or C difference is 0.00076913, below 0.001. Across all 100 propagations, the largest Chebyshev–uniformization difference is 4.23e-10, below 1e-8. Each propagator's allocated covariance tail remains below 1e-9. Numerical qualification explicitly requires these gates and stationarity; a completed calculation alone does not qualify a reference.

The three historical nonlinear response matrices reproduce exactly. Independent polynomial Gaussian moments agree with the contracted control calculations within 2.00e-15. The local CPU pilot and A100 reference agree within 2.08e-13 in normalized response. These checks address computation; they do not repair the physical-model discrepancy.

## Files and replay

- `results/summary.json` records all headline results and evidence boundaries.
- `results/all-response-entries.csv` contains 180 rows: all 90 entries for each control, signed response, reference value, discrepancy and error accounting.
- `results/all-two-other-site-orderings.csv` retains all 60 control/receiver/time comparisons, signed amplitudes, gaps and unresolved cases.
- `results/reference-refinement.csv` and `results/independent-propagation.csv` retain every reference qualification result.
- `remote/all-five-grids-ten-times/data/case-*/response.npz` contains the full static and delayed three-by-three kernels and response matrices from both propagators. Each accompanying receipt records grid size, stationarity, polynomial degree, tails, timing and memory information.
- `source/` contains unchanged inherited reference and Gaussian source. `inputs/` supplies the frozen covariance, historical control records, static checks and protocol.
- `compute/resource-ledger.json` distinguishes actual instrumented duration, requested resources, delivered hardware and billing limitations.

For saved-data verification, install the numerical dependencies in a suitable Python environment and run:

```bash
python verify.py
```

This reads saved results, reconstructs all control responses independently, checks source/input hashes and reference gates, recomputes ordering counts, and replays the aggregate analysis. All required inputs are included in this directory; historical paths in provenance are labels, not replay dependencies. It launches no GPU, cloud or QPU task. Python 3.9.6 with NumPy 2.0.2 and SciPy 1.13.1 was used for this local verification. Dependency pins describe that environment; numerical tolerance and hashes have different roles when replaying on another platform.

To recompute derived tables, run `python analyze.py`; it replaces only the current derived CSV and summary files. `python controls.py` replaces the current control arrays and receipt, so use a fresh copy if original-file bytes must be preserved. The complete raw reference calculations are reproducible from `source/dynamics/core.py` with the cases and times in `protocol.json`. `run_reference(..., backend_name='numpy', times=..., independent=True)` provides CPU execution. Its command-line interface uses the historical three-time default; a Python caller must supply all ten frozen times.

`local_run.py` reruns the 288-state implementation check and first registered grid at 0.0075 tau, with an 8 GB resident-memory ceiling and a 300-second limit per worker. Its output directories must be absent; use a fresh copy. The full CUDA replay uses the explicit `modal_run.py` entrypoint and `job.json`, requires separately authorized resources, and creates a cloud job. It is not invoked by verification. The inherited FP64 CUDA stencil has no Metal backend.

## Resources and interpretation

The local first-grid pilot used 4.93 solver seconds and 2.01 GB peak resident memory. One remote invocation requested A100-40GB, eight CPU cores and 64 GiB host memory; the provider delivered an A100-80GB. Instrumented function time was 258.73 seconds, or 0.07187 A100 hours, within the one-hour ceiling. The full client interval including downloads was 285.64 seconds. Final account inventory contained no running containers. GPU pool reservations and host process RSS have separate meanings and are recorded as such. No billing invoice is inferred from execution times.

The remote image pins Python 3.12, NumPy 2.2.6, SciPy 1.15.3 and CuPy 14.1.1. NumPy and CuPy versions are also measured in the reference receipts; the exact remote Python patch version was not captured. The public copy replaces a provider task identifier and local path prefixes in metadata; measured resources, scientific values, source algorithms and failure records are retained.

The useful decision is now clearer: accurate reduction of the Gaussian model and adequacy of that model for nonlinear mechanics are separate requirements. Gaussian response remains a testable screening descriptor with its own independent biological evaluation, but this experiment does not support describing it as a faithful approximation to the tested original nonlinear dynamics. A future remedy needs a new physical derivation and must beat these frozen controls under the same amplitude and reference requirements.

## Public provenance and independent audits

`PUBLIC-MANIFEST.json` seals the public files. `PUBLIC-PACKAGING.json` maps original to public hashes for metadata adaptations and records changes to the replay wrappers. The original protocol and its freeze record are retained byte for byte. The freeze attests the original local input hashes; the public verifier checks the explicit mapping for redacted metadata. It does not claim that those modified metadata files retain their original seal. Numeric JSON values and all original raw NPZ arrays and CSV result tables are preserved.

`provenance/` retains the original local file manifest, verification receipt and verifier source for provenance, rather than as current replay entrypoints. `inputs/nonlinear-model-adequacy.md` is the historical proposal that preceded this completed study; its status and planning ceilings describe that earlier record. Only its two links have been made absolute for portability. The completed protocol, resource ledger and results govern this execution.

The additional independent audit evaluates the quartic observables directly by six-dimensional Gauss–Hermite quadrature, using 15,625 points at each model/time. It imports neither the contracted Gaussian kernel nor its polynomial-moment implementation. The resulting control moments agree within 1.96e-14. The second audit independently reconstructs every reference screen, response discrepancy and two-site ordering from saved arrays. To repeat both without modifying the package, choose a separate output directory:

```bash
python audits/verify_controls.py --output /tmp/pulsar-adequacy-audit
python audits/verify_nonlinear_results.py --output /tmp/pulsar-adequacy-audit --independent-controls /tmp/pulsar-adequacy-audit/independent-control-moments.npz
```

`audits/` contains the portable audit scripts and their packaging replay receipts. `audits/original/` preserves the original audit scripts, receipts and arrays. Replay receipts record their own wrapper hashes and timings; they are new local checks of the saved data, not additional nonlinear simulations. Code and input provenance are described in `SOURCES-AND-RIGHTS.md`.
