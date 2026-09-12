# Reproduce the cached-input pilot timing

Keep `replay_pilot_cost.py` and `pilot-cost-worker.py` together. The recorded run used Python 3.12.14, NumPy 2.1.3 and SciPy 1.14.1. The scripts additionally use only the Python standard library. RSS accounting supports macOS and Linux. Use the pilot's pinned environment; package installation is outside the measured workflow.

The supplied pilot directory must contain the frozen `prepare.py`, `run_pilot.py`, `preanalysis-protocol.json`, cached `raw/` inputs and the reference `model/` and `results/` arrays. The script checks the exact scientific source hashes used for the recorded timing. It refuses an incompatible source version rather than silently timing a different method.

```sh
python replay_pilot_cost.py --pilot /path/to/kras-pilot --work /path/to/new-cost-run
```

Use a fresh work directory outside the source pilot. The launcher copies the cached inputs and scripts there. It modifies only the copied `prepare.py`, using the exact timer-only transformation in `instrument_prepare()`: import `time`, start a monotonic timer immediately before contact-graph/Hessian construction, stop after eigendecomposition and eigenvector conventions, then write a separate timing record. No scientific statement, parameter or calculation is replaced. The expected hash of the instrumented copy is checked. `run_pilot.py` is copied unchanged.

The worker launches each phase once in a fresh process and measures wall time plus process peak RSS. It uses the existing per-dimension harmonic timers to identify the combined all-mode covariance, response and standard-deviation calculation. That time is not described as a pure normalizer-only benchmark. After the timed phases, it compares all 37 generated NPZ files with the supplied reference: 513 arrays in the recorded version. Integer, Boolean and ranking arrays must agree exactly; floating arrays must agree with absolute and relative tolerance 10⁻¹⁰. The original measured replay matched every array exactly.

The outputs are `full-pilot-cost.json` and an identical sanitized `full-pilot-cost-public.json` in the new work directory. `--output` can select another new JSON path. Existing output records are not overwritten. The record includes machine class, operating system, software versions, source/input hashes, timing scope, verification results and exclusions; it contains neither a host name nor the supplied absolute paths.

Timing and memory are observed measurements and will vary with hardware and system load. This replay includes structural preparation, two independent harmonic moment checks, eight harmonic dimensions, 24 grid cases, built-in physical diagnostics, graph/geometric/degree baselines and saved results. It excludes installation, cached-input copying, network fetches, external Ohm execution, reference-pocket/resampling evaluation, figures, PDF generation, quantum work and the post-run comparison. Shared preparation and harmonic computations are already within the total, so adding them again would double-count cost.

The portable launcher and worker were syntax checked, and the instrumentation transformation was verified against the executed copy. The scientific workload was not rerun when packaging these scripts.
