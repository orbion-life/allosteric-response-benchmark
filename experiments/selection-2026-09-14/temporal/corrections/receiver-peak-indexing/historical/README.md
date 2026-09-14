# Temporal fidelity of fixed Gaussian protein representations

Completed on 2026-09-14. The source and protocol were sealed before the new response queries. This is a dense numerical characterization of previously examined models, not a fresh temporal or biological holdout. No model, covariance, basis, normalization, candidate mask or threshold was changed after evaluation.

**Result:** the separately frozen event-aligned addendum confirms that early errors recur after all three pulse durations; see [the addendum](event-aligned/README.md). The current representations fail the 0.002 response threshold somewhere on the declared time grid for every protein. The earlier three-time successes remain reproducible. The separate small original-model pulse experiment passes its intervention-identity check. These findings support an explicit temporal recovery milestone; they do not support a claim of unrestricted response preservation.

## Fixed-grid results

The test uses 74 observation times: 61 logarithmic times from 0.001 to 30τ, zero, and declared/historical comparison times. The three pulse durations are 0.03, 0.3 and 3τ. Evaluating all nonnegative shifted arguments requires 173 unique full Gaussian kernel times for each protein. τ remains an uncalibrated model time scale, not experimental seconds.

| Target | Fixed rank | Largest step error | Worst sampled t/τ | Passing step times | Largest pulse error, d/τ = .03 / .3 / 3 |
|---|---:|---:|---:|---:|---|
| KRAS | 347 | 0.0465222 | 0.0110831 | 46/74 | 0.0465222 / 0.0465222 / 0.0465222 |
| ABL | 506 | 0.0590525 | 0.00332913 | 50/74 | 0.0590525 / 0.0590525 / 0.0590525 |
| MYC | 172 | 0.0198672 | 0.001 | 62/74 | 0.0198672 / 0.0198672 / 0.0198672 |
| MYH7 | 934 | 0.0654595 | 0.001 | 55/74 | 0.0654595 / 0.0654595 / 0.0654595 |

Errors are maximum absolute signed matrix-entry differences in the fixed harmonic-standard-deviation normalization; they are not relative percentages. A pulse includes the initial on-phase, so its all-time maximum can equal the step maximum. After-removal values below answer a separate descriptive question and do not replace the frozen whole-curve gate. The observation grid is measured from pulse onset: the 0.3τ and 3τ pulses do not have a dense set of points immediately after removal. Their lower sampled errors therefore do not establish fidelity at short post-removal lags. An event-aligned follow-up would require an explicit new protocol; no points were silently added here.

| Target | Largest pulse error strictly after removal, d/τ = .03 / .3 / 3 |
|---|---|
| KRAS | 0.0369474 / 0.00957862 / 0.000381739 |
| ABL | 0.0486572 / 0.00100955 / 3.07193e-05 |
| MYC | 0.0184599 / 1.54194e-05 / 9.37973e-10 |
| MYH7 | 0.0645113 / 7.6033e-05 / 3.48252e-07 |

## Rankings, signs and timing

Complete signed responses, receiver RMS scores, top-five memberships/orders and sampled absolute-entry peak indices are in each `step.npz`, `pulse-*.npz` and `analysis.json`. The frozen `low_signal` flag uses the maximum receiver RMS score across all residues, including receivers, with a 0.002 floor. A separately labelled candidate-only qualification audit is provided because that global flag can classify distal candidate scores incorrectly. Fifth/sixth score gaps of at most 0.004 are marked uncertifiable under the allocated 0.002 entry-error budget. A matching observed shortlist does not repair a failed full response matrix.

The signed-entry comparison excludes reference magnitudes at or below 0.002 only from its sign-disagreement denominator; all entries remain in the primary error metric. A reference ranking that changes over time is a model prediction, whereas disagreement with the reference at the same time is representation error.

The peak diagnostic distinguishes low-amplitude, exactly flat/tied and boundary peaks using the frozen rules. These are unique interior sampled maxima under the fixed 1e-12 tie rule, not error-certified or continuous-time peaks. Nearly flat plateaus can produce large sampled peak-time differences under tiny amplitude changes. Do not interpret these timing figures as molecular signal velocities, biological delays or quantum scrambling.

## Independent pulse intervention check

The unchanged original quartic three-node model was propagated on a separate 64-state reflecting finite box. Fields were switched on until 0.3τ and then removed. Three shrinking signed amplitudes (.01, .005, .0025) were applied at each source and evaluated at six times. The finest-amplitude central response differs from the linear-response pulse identity by at most 6.5601604e-08; halving amplitude reduces the observed error by approximately four. Independent forward-generator and symmetric-eigensystem propagation differ by 3.0447866e-14 in the central response.

This verifies a finite-model intervention identity. It does not establish physical convergence, the validity of the Gaussian approximation for proteins or biological efficacy. An initial missing model-protocol file stopped the tiny test before any grid or scores; the setup failure and exact dependency correction are retained.

## Compute and numerical checks

One private A100-40GB worker completed the four-target reference campaign in 332.32 seconds (337.16 seconds client elapsed). The requested worker was limited to 3 CPU cores, 12 GiB host memory and 2,700 seconds, with no automatic retries. No QPU was used. Across the four targets, saved CPU versus GPU kernel parity is at most 3.0642155e-14. Each target also passes independent polynomial/Wick/derivative fixture checks before new protein times are evaluated.

The independent reduced matrix-exponential comparison is at most 2.6445697e-13. The complete offline verifier passed: True. This is an implementation/provenance verification result; the scientific temporal failures above remain failures. Full raw matrices, resource receipts, source/input hashes and query costs are retained.

## Portable replay

All numerical replay inputs and dependencies are included. The archived `prepare.py` records how the fixed inputs were assembled from preceding evidence; it is not required for replay. Machine-specific original source paths in its manifest are provenance labels, not runtime dependencies.

Install NumPy and SciPy. Run the compact protocol/source and saved tiny-pulse check without allocating GPU resources:

```sh
python verify.py --mode protocol --output protocol-audit.json
```

`COMPACT-CI.json` lists the small files needed for that CI check. A relocated copy containing only those files passes (`ci-relocation-receipt.json`). It does not replay the four proteins. With the complete inputs/results present, verify every full-kernel hash, fixed input, step/pulse identity and reported error:

```sh
python verify.py --mode saved --output saved-verification.json
```

Recompute a fixed-covariance full reference grid in a new output directory on an existing CPU or CUDA machine:

```sh
python full_kernel.py --target kras --backend numpy --output fresh-kras
python full_kernel.py --target myh7 --backend cupy --output fresh-myh7
```

The CUDA path additionally requires a compatible CuPy/CUDA runtime; the measured environment is recorded in the receipts. The CPU path uses the same float64 contractions. `analyze.py TARGET` reconstructs current results under `results/TARGET`; use a disposable package copy to avoid overwriting archived analyses. It never refits the model. `remote.py` is an explicitly invoked optional Modal launcher that allocates one GPU; no verifier calls it. Historical results are inspectable without the private volume.

## Files for review

- `protocol.json`, `source-seal.json`, `input-manifest.json`: fixed pre-execution decisions and identities.
- `summary.json`: compact target/endpoint metrics and scope.
- `results/TARGET/full-kernels.npy`: all 173 exact Gaussian reference matrices.
- `results/TARGET/step.npz`, `pulse-*.npz`: all 74 full/reduced response matrices, scores, orderings and sampled peaks.
- `results/TARGET/analysis.json`: complete frozen-metric results.
- `results/TARGET/receipt.json`, `saved-reference-canary.json`, `fixture-check.json`: timing, hardware and independent numerical admission.
- `tiny-field/raw.npz`, `tiny-field/receipt.json`: separately propagated on/off fields and linear-response comparison.
- `saved-verification.json`: complete offline audit, with scientific failure status retained.
- `figures/temporal-error.pdf`: static scientific figure, reproduced from the saved errors.
- `PUBLIC-FILE-MANIFEST.json`: publication-reviewed files and hashes. No biological outcome table or restricted third-party supplement is in this branch.

A successful recovery would require a separately frozen revised representation and new confirmation points after construction. This campaign supplies no basis change, fresh biological prediction or quantum advantage.

Figure: lines connect the finite sampled time points as guides; they do not certify behavior between the tested points. Values below the displayed 1e-9 plotting limit remain in the raw arrays.

## Distal responses and compact publication

`candidate-qualification-audit.json` is a secondary descriptive analysis using the predeclared candidate and receiver masks. It retains every frozen full-matrix gate and global-score flag, then reports candidate-specific signal floors, distal matrix errors, reference peaks, and relative errors. Relative Frobenius errors use the unweighted norm over the listed times and fixed candidate-by-receiver entries; they are not continuous-time integrals. A reference peak at or below 0.002 is marked indeterminate rather than treated as useful fidelity.

For the event-aligned 3τ pulses, the maximum distal errors are 0.0088533 (KRAS), 0.0102763 (ABL), 0.00351585 (MYC) and 0.0139782 (MYH7). Their reference peaks are 0.107566, 0.107461, 0.0710303 and 0.119742. The full response failures therefore cannot be attributed solely to irrelevant diagonal entries. These remain model-based numerical findings; they establish neither biological signal distortion nor the proposed quantum-regenesis mechanism.

`COMPACT-PUBLICATION-MANIFEST.json` is the authoritative publication selection after the addendum. It retains all eight immutable raw full-kernel arrays, fixed model/operator inputs, code, protocols, receipts, summaries and independent tiny-field measurements. Redundant derived step/pulse archives are explicitly omitted with original hashes and reconstruction commands. `FULL-LOCAL-FILE-MANIFEST.json` records the complete local evidence; the earlier `PUBLIC-FILE-MANIFEST.json` remains a historical pre-addendum review.

Run `python replay_waveforms.py --output NEW_REPLAY_DIRECTORY` to regenerate the omitted waveforms and perform both complete saved-array audits in a separate tree. The immutable raw files are hardlinked where possible or copied; source input files receive no writes. This workflow reconstructs derived response arrays and does not rerun the GPU kernels or fit models. The direct original saved-array audits and this regenerated-waveform audit are reported separately.

The compact reconstruction was executed in a fresh tree and passed both complete waveform audits. All 28 regenerated derived archives have identical SHA-256 hashes to their original local versions on the recorded environment (`portable-replay-verification.json`). The replay took 213.52 seconds and did not recompute a full reference kernel. This is separate from the two GPU reference campaigns.
