# Protein response circuits and complete estimation cost

**Choose the replay package first.** This Git checkout contains source and a small triangle fixture for the lightweight checks below. The complete KRAS input, saved results, circuits and verification records are in the [v0.8.0 quantum evidence ZIP](https://github.com/orbion-life/allosteric-response-benchmark/releases/download/v0.8.0/pulsar-protein-quantum-cost-evidence.zip). Verify the ZIP checksum against the [release asset manifest](../../docs/evidence-v0.8.0/ASSET-MANIFEST-v0.8.0.json), extract it, and run full verification from `pulsar-biological-quantum-resolution-2026-09-14/quantum/`. For regeneration, follow the fresh-copy instructions below. The archive manifest describes the sealed full evidence package; this repository README adds this checkout-specific instruction.


An exact spectral representation now produces actual protein response circuits for the accepted 347-dimensional KRAS Gaussian operator. Representative native circuits pass numerical statevector checks. An independently checked iterative amplitude-estimation study reduces the estimated complete T cost to **1.336% of the evaluated fixed sufficient sampling budget for the same spectral circuit**. The remaining cost is still too large to establish a practical full-protein quantum route, and the same classical preprocessing already yields the exact model response in milliseconds.

This is a new exact representation and estimator study. It does not change the previous biological, original-density, temporal holdout or physical-convergence findings. The four-target Gaussian development revision remains the source of the KRAS operator; no original nonlinear protein fidelity is inferred.

## What was implemented

The circuit retains every eigenmode of H = VΛVᵀ. It explicitly prepares the signed columns of VᵀB, implements exp(−tΛ) through a controlled flag rotation, and measures the complete response by a Hadamard test. This removes the factor representation's Γ-dependent polynomial degree, while charging the full classical eigendecomposition, transformed observables and data-dependent rotation tables. All 166 KRAS observable preparations and all three attenuation circuits were constructed and saved.[^2][^3]

The KRAS response circuit has **11 active qubits, 3,068 Ry rotations and 3,068 CX gates** for the fixed representative pair. Its actual Grover circuit has **20 qubits**, including nine clean reflection helpers, 6,136 Ry rotations, 6,191 CX gates and 72 T/T† gates in its exact reflection. The rotations in this native count remain arbitrary rotations; they are not free fault-tolerant gates. Six representative KRAS native probability checks, including an actual full Grover iteration, agree with the independent expected probabilities within 2.56×10⁻¹⁵. Triangle checks cover powers 0, 1, 2 and 4. These are noiseless numerical circuit executions, not QPU experiments.

The confidence method uses the monotonic-branch principle of iterative amplitude estimation.[^1] Our explicitly specified variant uses fresh independent batches and allocates failure probability 0.05/[Jr(r+1)] to output/round pairs. It estimates every requested unordered pair at the three original saved times: 18 triangle outputs and 41,583 KRAS outputs. The scheduler uses previous observations, not the known amplitude. Exact model probabilities generate explicitly labelled **classical binomial simulations** of the ideal measurements. No hardware or all-output fault-tolerant execution is claimed.

## Complete measured and estimated results

| Quantity | Triangle, rank 9 | KRAS, rank 347 |
|---|---:|---:|
| Requested outputs | 18 | 41,583 |
| Completed ideal interval estimates | 18 | 41,583 |
| True ideal amplitudes inside saved intervals | 18 | 41,583 |
| Readouts in realized IQAE simulation | 20,992 | 63,497,472 |
| Calls to A or its inverse | 1,133,056 | 10,568,331,008 |
| Maximum Grover power in a shot | 212 | 633 |
| Estimated complete T gates | 141,682,054,144 | 56,450,203,649,724,416 |
| Estimated complete CX gates | 118,141,952 | 32,712,522,454,784 |
| T ratio to the same spectral sampling control | 0.037284 | 0.013361 |
| Maximum CX gates in one coherent shot | 44,400 | 3,921,971 |
| Maximum T gates in one coherent shot | 53,170,034 | 7,642,625,052 |
| Maximum native serial depth upper estimate | 87,524 | 7,821,140 |

The sufficient k = 0 sampling control uses **1,687,979 shots per output** for the triangle and **18,983,467 per output** for KRAS, with the same output family, circuit and remaining error budget. The KRAS control totals 789,389,508,261 readouts and 4.22515×10¹⁸ constructive T gates. The reduction is relative to that evaluated control; it is not a comparison with an optimal quantum algorithm or the exact classical calculation.

All precision allowances remain inside 0.002 at the saved queries. The KRAS inherited Gaussian projection error is 0.000901653. Its maximum saved statistical response halfwidth is 0.001065304. Including the declared arithmetic and constructive encoding bounds gives **0.001968468**. This is a conditional simultaneous confidence calculation for ideal repeated circuits plus fixed encoding bias, not a physical-noise guarantee. The maximum point error in the ideal likelihood simulation is 0.000526322 relative to the reduced model.

## Fault-tolerant accounting and its limits

Every native Ry is assigned a finite modulo-4π word and decomposed into a fixed bank of actually synthesized Clifford+T Rz rotations. All words, bank circuits, synthesis records and counts are saved.[^4] Exact inverses reuse the emitted sequence. Thus the fixed approximate A has a bounded encoding bias, and its Grover iterations can target that same fixed unitary. The bound is explained in [CONFIDENCE-AND-ENCODING.md](CONFIDENCE-AND-ENCODING.md).

This binary compiler is constructive and conservative; its counts are **not minimum T costs**. In particular, the representative KRAS A costs 6,018,368 T gates and 15,154,252 operations in this recipe, exceeding the two-million-operation export cap. Its complete fault-tolerant circuit was not constructed. Actual bank circuits and every data word determine the count. The smaller triangle A and Q were exported and numerically executed: the emitted ten-qubit k = 1 circuit contains 946,908 operations and preserves the amplification probability relation to 8.72×10⁻¹³. No maximum-power adaptive circuit was built.

The earlier SELECT–SWAP/phase-gradient results remain in v0.7. This new spectral study does not claim that its binary compiler is better than every phase-gradient or joint-synthesis alternative. Extending such optimizations to the spectral tables could reduce T counts further, but would require its own workspace, reflection, precision and complete-output accounting. It was not used to replace the measured result here.

The adaptive cost records describe realized **ideal** sample schedules. A finite-precision circuit can slightly change the probabilities and therefore the realized schedule. Noise can invalidate the Grover likelihood entirely. Routing, logical-to-physical encoding, state factories, calibration, physical cycle time and available hardware are not inferred from the 11- or 20-qubit logical registers. The deepest KRAS shot alone contains nearly four million native CX gates. The evidence supports circuit construction and a resource reduction, while leaving practical execution unresolved.

## The decisive classical control

On the recorded host, KRAS diagonalization took 0.01185 seconds, transforming all B columns took 0.001081 seconds, and computing all three full response matrices afterward took **0.000808 seconds**. A separate matrix-exponential computation agrees within 1.19×10⁻¹³. These are one-host measured timings, not universal performance claims. The quantum representation uses this same diagonalization, so it has no demonstrated advantage over computing the weighted dot products directly.

The next useful quantum decision is therefore not a claim that a smaller register solves this protein faster. It is a comparison on a precisely defined problem where classical preprocessing and requested outputs remain expensive, and where the quantum circuit receives an explicitly costed data-access model. The current spectral construction is a checked control and compiler benchmark for that decision. A hardware study must separately show that a declared implementation can tolerate its compiled coherent depth and recover its stated output precision.

## Confidence and verification

The exact Clopper–Pearson finite-binomial check covers 1,001 probabilities at each of two failure levels. A separate set of 352 fixed-amplitude simulations yields 351 covered intervals under a nominal 99% individual procedure; that statistical miss is retained. The simultaneous claim follows the conditional interval/union-bound argument, not this finite diagnostic. Our variant does not inherit the published IQAE paper's exact asymptotic constant or termination theorem without a separate proof.[^1]

`independent-verification.json` reconstructs every saved adaptive interval, all readout/Grover counts, every binary rotation word, T counts from actual bank QPY files and both complete cost totals. The parent also performed independent native and accounting audits outside this branch. No hardware result or biological advantage is implied by these checks.

Two implementation corrections are preserved. The initial fixed-amplitude test exposed a counter that skipped the first doubling step of a repeated fresh batch; it was corrected before model IQAE results. The first KRAS cost run repeatedly decompressed immutable NPZ arrays and took 252.31 seconds. It completed naturally just before an attempted interrupt. The canonical source now caches those arrays and completes in 2.55 seconds with every scientific summary field unchanged. The original complete stdout and source are retained, but its raw NPZ was overwritten before the history copy. Therefore no original-versus-cached raw byte comparison is claimed; the current cached raw schedule is independently reconstructible.

## Portable replay

All required project inputs and own helper sources are included. Current entrypoints resolve them relative to this directory. The source uses no required original workspace path or older downloaded release. A normal installation of the pinned requirements makes the optional `PULSAR_SYNTH_RUNTIME` override unnecessary. No installed runtime, credentials, account inventory or full research article is redistributed.

Lightweight source-only CI requires `ci_smoke.py`, `iqae.py`, `spectral.py`, `vendor/multiplexed.py` and `requirements-ci.txt`. It writes no artifacts and does not need archived protein results:

```sh
python -m pip install -r requirements-ci.txt
python -B ci_smoke.py
```

To verify the retained full evidence in a working copy:

```sh
python -m pip install -r requirements.txt
python verify.py
```

To regenerate results, use a fresh copy containing the source, `inputs/` and `vendor/`, with no existing `results/` directory. Run sequentially; entrypoints deliberately refuse to overwrite model construction folders:

```sh
python test_iqae.py
python run_native.py triangle
python run_native.py kras
python fault_tolerant.py triangle
python fault_tolerant.py kras
python run_iqae.py triangle
python run_iqae.py kras
python check_compiled_triangle.py
python verify.py
```

Measured environment: Python 3.13.5, NumPy 2.4.6, SciPy 1.16.0, Qiskit 2.3.1 and Aer 0.17.2. The largest native test used 20 statevector qubits and stayed inside the 12 GiB worker cap. Native KRAS construction/checks took 19.72 seconds; no new GPU or QPU allocation was made. The package manifest distinguishes repository source/fixtures from large release evidence and excludes installed runtimes, caches and the internal process profiler dump.

## Sources

[^1]: D. Grinko, J. Gacon, C. Zoufal and S. Woerner, [*Iterative quantum amplitude estimation*, npj Quantum Information 7, 52 (2021); author version, Section III and Appendix C](https://arxiv.org/html/1912.05559v3). Used for the amplification probability and monotonic interval principle; the modified confidence allocation is derived separately in this evidence.
[^2]: G. Brassard, P. Høyer, M. Mosca and A. Tapp, [*Quantum amplitude amplification and estimation*](https://arxiv.org/abs/quant-ph/0005055). Original amplitude-amplification and estimation framework.
[^3]: M. Möttönen, J. J. Vartiainen, V. Bergholm and M. M. Salomaa, [*Transformation of quantum states using uniformly controlled rotations*](https://arxiv.org/abs/quant-ph/0407010). State-preparation construction underlying the signed real loaders.
[^4]: [Qiskit Algorithms IQAE implementation](https://github.com/qiskit-community/qiskit-algorithms/blob/main/qiskit_algorithms/amplitude_estimators/iae.py) and [pygridsynth](https://github.com/quantum-programming/pygridsynth). Primary implementation references; this branch implements its own fresh-batch confidence variant and includes own project helper code. Qiskit is Apache-2.0; pygridsynth is MIT-licensed and installed separately.
