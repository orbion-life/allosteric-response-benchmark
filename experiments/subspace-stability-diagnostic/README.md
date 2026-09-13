# Project Pulsar: tracing numerical subspace divergence

This completed local diagnostic identifies a consequential difference in the construction of the reduced subspace. It does not fix the original portability failure. The most informative case retains an accurate finite-grid response while its constructed basis and calculated error bound differ between two local numerical environments. Holding the input matrices fixed reduces the discrepancy; holding the full basis fixed reduces the bound discrepancy to approximately 10⁻¹¹ or less. All original scientific settings and historical failures remain unchanged.

The [scientific findings](findings.md) explain what these measurements support. A [three-panel evidence figure](figures/computational-evidence.pdf) places the stability result beside the separately completed harmonic calibration and nonlinear timing profile. ![Finite-model evidence](figures/computational-evidence.png)

The machine-readable [inspection receipt](inspection.json), [original comparison summary](results/summary.json), and full traces allow the measurements to be checked directly. The historical summary predates the trace-completeness correction; the inspection receipt identifies the complete trace locations.

## Results at the final selected rank

The response differences below compare the same normalized response matrix across all three prescribed times. “Shared input” means that both environments receive exactly the pinned environment's saved H, F and time arrays. The subspace gap is the largest of the two orthogonal-complement projection norms; it is insensitive to signs or rotations within the same complete subspace.

| Case | Comparison | Final rank | Full-subspace gap | Maximum response difference | Maximum bound difference |
|---|---|---:|---:|---:|---:|
| Development, κ = 1, 17³ | Native inputs | 96 | 1.340 × 10⁻⁸ | 2.442 × 10⁻¹⁵ | 4.400 × 10⁻¹³ |
| Development, κ = 1, 17³ | Shared input | 96 | 1.260 × 10⁻⁹ | 1.776 × 10⁻¹⁵ | 2.116 × 10⁻¹² |
| Elongated, κ = 100, 33³ | Native inputs | 48 | 0.154560 | 1.614 × 10⁻¹¹ | 2.094 × 10⁻⁶ |
| Elongated, κ = 100, 33³ | Shared input | 48 | 0.013446 | 2.090 × 10⁻¹¹ | 3.213 × 10⁻⁷ |

Every residual block admitted all three available columns. No checkpoint cut a residual block, no compared span was diagnostically rank deficient, and the smallest residual singular value was more than 10¹¹ times its admission threshold. These observations rule out threshold straddling and partial-block truncation as explanations for these particular traces. They do not identify one uniquely causal floating-point operation.

## Frozen scientific contract

The experiment uses the two cases, in the order, specified by the earlier recovery protocol: development κ = 1 at 17³, then elongated κ = 100 at 33³. Both use extent 4, all three internal coordinates, the original biquadratic energy and quartic node observables, β = μ = 1, harmonic normalization, reflecting boundaries, and the original three observation times. The dependency threshold remains 10⁻¹², rank checkpoints remain 3, 6, 12, 24, 48 and 96, and the response-bound admission threshold remains 0.002. The historical 10⁻⁹ portability allowance is unchanged. No full reference response chooses the rank.

The scientific source is retained byte for byte in [immutable/operator_study.py](immutable/operator_study.py), SHA-256 `13d80e28d8acaf9ca435938d8ef08668b3cf9d0e750b7b0970b6cf26f19434ef`. The scientific preanalysis is [immutable/preanalysis.json](immutable/preanalysis.json), SHA-256 `530fede82261aec39165f06d0126e75e995ad131c0af0db170ce6f9c6880ada3`. The historical development κ = 1, 33³ rejection was executed in the original campaign and is preserved under `original-results/`. It was not rerun in this diagnostic or replaced by a favorable new case.

## Reproduce the diagnostic in a new directory

The launcher uses the Python standard library and accepts explicit interpreter paths. Both interpreters must already contain NumPy and SciPy. The executed environments were Python 3.12.14 / NumPy 2.1.3 / SciPy 1.14.1 and Python 3.9.6 / NumPy 2.0.2 / SciPy 1.13.1, on the same arm64 macOS host. These runs are not Linux replays. Interpreter versions and configuration were captured during execution; [interpreter binary hashes](environment-provenance.json) were added at closure, so they are not a pre-launch binary attestation. The names `pinned` and `system` identify the two supplied interpreters; the inventory records their actual versions and numerical libraries.

```sh
python3 bootstrap_replay.py \
  --output ./fresh-replay \
  --pinned-python /path/to/first/environment/bin/python \
  --system-python /path/to/second/environment/bin/python
```

This command only prepares a new directory and prints its freshly written protocol hash. Add `--run` to execute the ordered 23-job diagnostic under one worker, the aggregate 1,800-second allowance and the 4,000,000,000-byte sampled memory limit. The helper refuses to overwrite any existing output directory. It neither installs packages nor modifies the archived execution protocols. The complete fresh campaign has not been rerun after adding this launcher; its preparation, source hashes, refusal to overwrite, missing-interpreter behavior, and preservation of virtual-environment invocation paths passed the [public bootstrap tests](provenance/public-bootstrap-test.log).

To check the launcher's packaging behavior without numerical work:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 -W error test_bootstrap.py
```

The original execution protocols, local paths and host runtime labels are preserved byte for byte in the separate raw provenance archives. Public-readable summaries use neutral placeholders. The public replay-template.json is explicitly a template, and the bootstrap writes a new dated, hashed execution protocol using your supplied interpreter paths. Original and sanitized-copy hashes are recorded separately in provenance/public-copy-map.json.

## Trace contents and completeness

The instrumenter inserts logging calls into the parsed original source. Removing those calls recovers the original numerical syntax tree. Native input construction passed bitwise comparisons with the uninstrumented source in both environments, and the control basis passed the corresponding uninstrumented comparison in both environments.

The complete trace sets are:

- `results/pinned/development-k1-n17-e4/native/` and `results/pinned/elongated-k100-n33-e4/native/`.
- `trace-repair/system/development-k1-n17-e4/native/` and its `shared/` sibling.
- `trace-repair/system/elongated-k100-n33-e4/native/` and its `shared/` sibling.

Each set has a `trace-completeness.json` receipt. The development traces contain 31 admissions and six checkpoints; the elongated traces contain 15 admissions and five checkpoints. Each admission records the raw action, both orthogonalization passes, singular decomposition, eligibility mask, optional prefix, signed admitted block and orthogonality diagnostics. The seed and admitted blocks reconstruct the complete thin basis between checkpoints. Each checkpoint also saves the full basis, projected matrices, residual singular values, propagated response, and the signed-term decomposition of the two bounds.

Identical-basis tests use every pinned checkpoint in both environments with `evr`, `evd` and `ev`. Identical-residual tests use `gesdd` and `gesvd` in both environments. The saved same-residual selection is step 10 for the development case and step 3 for the elongated case. No alternative recurrence or forward-remedy campaign was executed.

## Retained failed attempts and amendments

1. The initial run stopped while collecting optional configuration metadata because a missing PyYAML package triggered a warning under warnings-as-errors. No model was built during that attempt. Its executable source, protocol, log and receipt are preserved in `history/attempt-1-metadata-warning/`.
2. The first continuation completed the numerical jobs, but the Python 3.9 trace omitted the per-step residual SVD file. A text match against the AST unparser's assignment format caused that logging omission. The original incomplete traces remain in `results/system/`, and the pre-correction code remains in `history/attempt-2-incomplete-system-svd-trace/`.
3. The second continuation replaced that text match with a structural AST match and added static hook counts and dynamic trace assertions. It reran only the four system traces. All 688 previously saved non-reference NPZ files matched bit for bit. The model, numerical statements and selection criteria did not change.
4. A saved-array inspection initially assumed that a separately recomputed sparse reference would also be bitwise equal. That assertion failed for the development native reference: its maximum response difference was 6.106 × 10⁻¹⁶, its delayed-moment difference was 5.829 × 10⁻¹⁶, and its equilibrium covariance was identical. The other three repaired references were bitwise equal. The failed assertion and log remain in `history/inspection-reference-bitwise-assertion/`. The corrected inspection reports these differences instead of assuming exact equality; it does not alter the original scientific tolerance.

The monitored execution receipts report 2.005, 59.739 and 22.331 seconds for the original attempt and two continuations. Their largest sampled aggregate resident memory was 746,192,896 bytes. The final continuation's conservative carry-forward total was 85.358 seconds before the two short saved-array inspections. The [resource receipt](resource-accounting.json) distinguishes measured work from a conservative allowance for manifest and packaging overhead. Documentation and coordination gaps are not reported as numerical execution time.

## Storage and publication scope

The complete evidence tree contains approximately 1.825 GB before final documentation and manifests. Approximately 629 MB belong to identical-basis driver traces, which repeat some large arrays. No raw arrays were deleted or replaced. The [storage note](storage-and-release.md) recommends a small source/summary package plus separately hashed trace assets; a future content-addressed layer must preserve byte-exact reconstruction.

This is a retrospective numerical diagnostic of two finite models. It establishes neither continuum convergence, protein predictive validity, a stable replacement recurrence nor quantum advantage. The original bound-portability failure and the separate physical-domain and mesh-convergence failures remain open.

## Raw evidence assets

The six companion ZIP files are listed with their exact sizes and SHA-256 hashes in [the archive index](provenance/archive-index.json). Each decompressed member was checked against the original source file. Download all six assets, then extract them into the same new directory to reconstruct the complete `stability/` evidence tree. The raw tree includes the original local paths and host runtime metadata; those records have not been sanitized because their original hashes are part of the evidence. The small source package contains compact JSON results and historical controls; large trace arrays are intentionally in the ZIP assets.

Release asset availability is recorded by GitHub; the hashes in the archive index allow independent verification. The source follows the existing Orbion MIT license; NumPy, SciPy and Matplotlib are external dependencies and their binaries are not bundled.

### v0.5.0 evidence assets

Asset availability is recorded on the [GitHub release](https://github.com/orbion-life/allosteric-response-benchmark/releases/tag/v0.5.0). Each link below identifies one checksum-bound archive; the index retains the exact byte size and SHA-256 hash.

| Raw evidence archive | Compressed size |
|---|---:|
| [source history](https://github.com/orbion-life/allosteric-response-benchmark/releases/download/v0.5.0/subspace-stability-2026-09-13-source-history.zip) | 1.13 MB |
| [inputs and comparisons](https://github.com/orbion-life/allosteric-response-benchmark/releases/download/v0.5.0/subspace-stability-2026-09-13-inputs-and-comparisons.zip) | 20.64 MB |
| [pinned traces](https://github.com/orbion-life/allosteric-response-benchmark/releases/download/v0.5.0/subspace-stability-2026-09-13-pinned-traces.zip) | 214.11 MB |
| [original system traces](https://github.com/orbion-life/allosteric-response-benchmark/releases/download/v0.5.0/subspace-stability-2026-09-13-original-system-traces.zip) | 396.08 MB |
| [repaired system traces](https://github.com/orbion-life/allosteric-response-benchmark/releases/download/v0.5.0/subspace-stability-2026-09-13-repaired-system-traces.zip) | 428.43 MB |
| [fixed V drivers](https://github.com/orbion-life/allosteric-response-benchmark/releases/download/v0.5.0/subspace-stability-2026-09-13-fixed-V-drivers.zip) | 496.56 MB |
