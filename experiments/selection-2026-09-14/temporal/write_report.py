"""Rebuild the historical absolute-grid index only; final README includes separate addendum qualifications."""
from pathlib import Path
import json,hashlib,datetime
ROOT=Path(__file__).resolve().parent
s=json.loads((ROOT/'summary.json').read_text());v=json.loads((ROOT/'saved-verification.json').read_text());cloud=json.loads((ROOT/'remote-receipt.json').read_text())
lines=['# Temporal fidelity of fixed Gaussian protein representations','',
'Completed on 2026-09-14. The source and protocol were sealed before the new response queries. This is a dense numerical characterization of previously examined models, not a fresh temporal or biological holdout. No model, covariance, basis, normalization, candidate mask or threshold was changed after evaluation.','',
'**Result:** the current representations fail the 0.002 response threshold somewhere on the declared time grid for every protein. The earlier three-time successes remain reproducible. The separate small original-model pulse experiment passes its intervention-identity check. These findings support an explicit temporal recovery milestone; they do not support a claim of unrestricted response preservation.','',
'## Fixed-grid results','',
'The test uses 74 observation times: 61 logarithmic times from 0.001 to 30τ, zero, and declared/historical comparison times. The three pulse durations are 0.03, 0.3 and 3τ. Evaluating all nonnegative shifted arguments requires 173 unique full Gaussian kernel times for each protein. τ remains an uncalibrated model time scale, not experimental seconds.','',
'| Target | Fixed rank | Largest step error | Worst sampled t/τ | Passing step times | Largest pulse error, d/τ = .03 / .3 / 3 |','|---|---:|---:|---:|---:|---|']
for r in s['rows']:
 a=r['step'];b=' / '.join(f"{x['max_absolute_error']:.6g}" for x in r['pulses'])
 lines.append(f"| {r['target'].upper()} | {r['rank']} | {a['max_absolute_error']:.6g} | {a['worst_time_over_tau']:.6g} | {a['passing_times']}/{a['total_times']} | {b} |")
lines+=['','Errors are maximum absolute signed matrix-entry differences in the fixed harmonic-standard-deviation normalization; they are not relative percentages. A pulse includes the initial on-phase, so its all-time maximum can equal the step maximum. After-removal values below answer a separate descriptive question and do not replace the frozen whole-curve gate.','',
'| Target | Largest pulse error strictly after removal, d/τ = .03 / .3 / 3 |','|---|---|']
for r in s['rows']:lines.append(f"| {r['target'].upper()} | "+' / '.join(f"{x['max_error_strictly_after_removal']:.6g}" for x in r['pulses'])+' |')
lines+=['','## Rankings, signs and timing','',
'Complete signed responses, receiver RMS scores, top-five memberships/orders and sampled absolute-entry peak indices are in each `step.npz`, `pulse-*.npz` and `analysis.json`. Reference scores at or below 0.002 are marked low-signal. Fifth/sixth score gaps of at most 0.004 are marked uncertifiable under the allocated 0.002 entry-error budget. A matching observed shortlist does not repair a failed full response matrix.','',
'The signed-entry comparison excludes reference magnitudes at or below 0.002 only from its sign-disagreement denominator; all entries remain in the primary error metric. A reference ranking that changes over time is a model prediction, whereas disagreement with the reference at the same time is representation error.','',
'The peak diagnostic distinguishes low-amplitude, exactly flat/tied and boundary peaks using the frozen rules. A unique interior maximum on this finite grid is not an error-certified or continuous-time peak. Nearly flat plateaus can produce large sampled peak-time differences under tiny amplitude changes. Do not interpret these timing figures as molecular signal velocities, biological delays or quantum scrambling.','',
'## Independent pulse intervention check','',
f"The unchanged original quartic three-node model was propagated on a separate 64-state reflecting finite box. Fields were switched on until 0.3τ and then removed. Three shrinking signed amplitudes (.01, .005, .0025) were applied at each source and evaluated at six times. The finest-amplitude central response differs from the linear-response pulse identity by at most {s['tiny_field']['max_identity_error_by_amplitude'][-1]:.8g}; halving amplitude reduces the observed error by approximately four. Independent forward-generator and symmetric-eigensystem propagation differ by {s['tiny_field']['independent_forward_central_max_abs']:.8g} in the central response.",
'',
'This verifies a finite-model intervention identity. It does not establish physical convergence, the validity of the Gaussian approximation for proteins or biological efficacy. An initial missing model-protocol file stopped the tiny test before any grid or scores; the setup failure and exact dependency correction are retained.','',
'## Compute and numerical checks','',
f"One private A100-40GB worker completed the four-target reference campaign in {cloud['wall_seconds']:.2f} seconds ({cloud['client_wall_seconds']:.2f} seconds client elapsed). The requested worker was limited to 3 CPU cores, 12 GiB host memory and 2,700 seconds, with no automatic retries. No QPU was used. Across the four targets, saved CPU versus GPU kernel parity is at most {max(x['GPU_saved_CPU_parity'] for x in s['rows']):.8g}. Each target also passes independent polynomial/Wick/derivative fixture checks before new protein times are evaluated.",
'',
f"The independent reduced matrix-exponential comparison is at most {max(x['independent_reduced_exponential_max_abs'] for x in s['rows']):.8g}. The complete offline verifier passed: {v['passed']}. This is an implementation/provenance verification result; the scientific temporal failures above remain failures. Full raw matrices, resource receipts, source/input hashes and query costs are retained.",
'',
'## Portable replay','',
'All numerical replay inputs and dependencies are included. The archived `prepare.py` records how the fixed inputs were assembled from preceding evidence; it is not required for replay. Machine-specific original source paths in its manifest are provenance labels, not runtime dependencies.','',
'Install NumPy and SciPy. Run the compact protocol/source and saved tiny-pulse check without allocating GPU resources:','',
'```sh','python verify.py --mode protocol --output protocol-audit.json','```','',
'`COMPACT-CI.json` lists the small files needed for that CI check. It does not replay the four proteins. With the complete inputs/results present, verify every full-kernel hash, fixed input, step/pulse identity and reported error:','',
'```sh','python verify.py --mode saved --output saved-verification.json','```','',
'Recompute a fixed-covariance full reference grid in a new output directory on an existing CPU or CUDA machine:','',
'```sh','python full_kernel.py --target kras --backend numpy --output fresh-kras','python full_kernel.py --target myh7 --backend cupy --output fresh-myh7','```','',
'The CUDA path additionally requires a compatible CuPy/CUDA runtime; the measured environment is recorded in the receipts. The CPU path uses the same float64 contractions. `analyze.py TARGET` reconstructs current results under `results/TARGET`; use a disposable package copy to avoid overwriting archived analyses. It never refits the model. `remote.py` is an explicitly invoked optional Modal launcher that allocates one GPU; no verifier calls it. Historical results are inspectable without the private volume.','',
'## Files for review','',
'- `protocol.json`, `source-seal.json`, `input-manifest.json`: fixed pre-execution decisions and identities.',
'- `summary.json`: compact target/endpoint metrics and scope.',
'- `results/TARGET/full-kernels.npy`: all 173 exact Gaussian reference matrices.',
'- `results/TARGET/step.npz`, `pulse-*.npz`: all 74 full/reduced response matrices, scores, orderings and sampled peaks.',
'- `results/TARGET/analysis.json`: complete frozen-metric results.',
'- `results/TARGET/receipt.json`, `saved-reference-canary.json`, `fixture-check.json`: timing, hardware and independent numerical admission.',
'- `tiny-field/raw.npz`, `tiny-field/receipt.json`: separately propagated on/off fields and linear-response comparison.',
'- `saved-verification.json`: complete offline audit, with scientific failure status retained.',
'- `figures/temporal-error.pdf`: static scientific figure, reproduced from the saved errors.',
'- `PUBLIC-FILE-MANIFEST.json`: publication-reviewed files and hashes. No biological outcome table or restricted third-party supplement is in this branch.','',
'A successful recovery would require a separately frozen revised representation and new confirmation points after construction. This campaign supplies no basis change, fresh biological prediction or quantum advantage.']
(ROOT/'absolute-grid-index.md').write_text('\n'.join(lines)+'\n')
print(ROOT/'absolute-grid-index.md')
