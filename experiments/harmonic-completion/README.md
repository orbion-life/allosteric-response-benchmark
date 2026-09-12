# Does restoring omitted harmonic response repair nonlinear compression?

**The tested additive correction does not pass the declared accuracy criterion.** It restores the harmonic response algebraically, but this is insufficient for the nonlinear response. Every fixed result is retained.

For common full-harmonic observable scales and time, the tested expression is

\[\widetilde C_d(t)=C^{H}_{all}(t)+C^{NL}_{d,grid}(t)-C^{H}_{d,grid}(t).\]

The three-node triangle has three positive Cartesian Hessian modes. All three are included in the nonlinear finite-grid reference. This is a fully retained *coordinate chart*, not an exact global quotient by rigid rotations or a certified continuum equilibrium. The local chart can contain distant rigidly equivalent minima; see the independently computed [domain audit](../gaussian-closure-design/coordinate-domain-check.json). No response or population claim follows from that audit.

| Stiffness κ | Restricted two-mode error | Error after completion |
|---:|---:|---:|
| 1 | 0.062145 | 0.600675 |
| 10 | 0.427594 | 0.055961 |
| 100 | 0.454687 | 0.008208 |

These are maximum absolute normalized response differences across all residue pairs at t = τ, using 65 points per coordinate. The proposed allocation is 0.002. It is an engineering choice, not a biological or literature-derived tolerance. Grid and domain differences also fail some fixed conditions, so these values are finite-grid diagnostics, not certified errors against a continuum solution. Neither an exact harmonic limit nor the full-coordinate same-grid identity demonstrates nonlinear fidelity.

The protocol fixes κ = 1, 10 and 100; t/τ = 0.1, 1 and 10; retained dimensions 1, 2 and 3; and 17, 33 and 65 points per coordinate. The domain extends four slow-mode standard deviations; 81 points at five deviations hold spacing fixed. The 72 finite-grid archives each contain all three times. Nine analytical harmonic archives complete the reference record. All 27 dimension/stiffness/time joint decisions and all 81 grid comparison rows are published in `analysis/`.

## Reproduce and check

Use Python 3.12 and an unused output directory, from this folder:

```sh
python3.12 -m venv .venv
.venv/bin/pip install -r requirements.txt
.venv/bin/python verify.py
VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python run.py --output replay --workers 4
.venv/bin/python verify.py --repeat replay --output checks/replay.json
.venv/bin/python analyze.py --results replay --output replay-analysis
.venv/bin/python protein_diagnostic.py --output replay-protein
```

`verify.py` also calculates the derivative independently by evolving the unperturbed equilibrium under positive and negative interventions on a smaller 125-state grid. It checks the linear-response identity, generator stationarity, the zero-time response, common normalization, harmonic restoration and saved-array repeatability. This separate intervention test is algebra/implementation evidence, not proof of continuum accuracy.

A complete local run can take tens of minutes. The original 20-minute/4-GB allocation underestimated the three-coordinate grids. `execution-amendment.json` records a 60-minute/8-GB local allocation and four-worker continuation, after unfavorable results were observed. No geometry, stiffness, time, grid, reference, response gate or failed result was changed. Source history and checkpoint hashes retain the original execution. The first run including the pause took about 31 minutes. Its original memory receipt measured the parent only; the current runner and clean replay report a conservative parent-plus-worker bound. Machine timings and memory vary; an exceeded allocation produces an unresolved run, not a scientific pass. The clean replay completed in 1,243 seconds with a conservative concurrent memory bound of 5.79 GB. All 567 arrays across 81 archives agreed within 5.99 × 10⁻¹⁵; see `checks/full-replay-verification.json` and `checks/clean-run-receipt.json`. No paid compute was used.

## KRAS diagnostic and interpretation

After the synthetic outputs were frozen, the same formula was applied to the original two-mode KRAS 33² and 65² arrays. Both yield residues 34, 85, 122, 123 and 81. Only residue 34 is within 5 Å of MOV. This is a diagnostic, not a replacement primary prediction. There is no full nonlinear KRAS reference, so this result cannot certify protein fidelity. The corrected missing-reference labels are used; the original released predictions are unchanged.

The source manifest binds the reused scientific functions. Original code and vendored Orbion functions use the MIT licence. The synthetic geometry is explicitly defined and requires no download. The optional KRAS diagnostic uses cached public data from `../kras/pilot` and the evaluation amendment. No molecular trajectories, fitting, remote services or ligand labels enter the synthetic experiment.
