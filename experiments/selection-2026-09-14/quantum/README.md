# Quantum execution cost experiment — 14 September 2026

The initial SELECT–SWAP circuit reduced the estimated complete T cost by 43.66%, but failed the frozen requirement of at least 50%. A separately frozen phase-gradient follow-up passed that requirement: **69.25% fewer T gates than the strongest dense implementation evaluated**, for the same rank-nine operator, response precision and full output schedule. Both outcomes remain in the evidence.

| Implementation | Logical qubits | Compiled T gates per controlled walk | Expected T gates for the complete ideal sampling plan | Ratio to strongest evaluated dense | Half-cost gate |
|---|---:|---:|---:|---:|---|
| Dense, norm-optimal normalization | 6 | 383,234 | 35,265,567,413,911 | 1.0000 | Comparator |
| Dense, factor normalization | 6 | 388,652 | 58,613,303,217,224 | 1.6621 | Comparator |
| Primary SELECT–SWAP with synthesized bit rotations | 147 | 121,512 | 19,869,883,844,711 | 0.5634 | FAIL |
| Separate SELECT–SWAP with phase-gradient rotations | 182 | 66,312 | 10,845,685,943,780 | 0.3075 | PASS within evaluated comparison |

The T totals combine actual synthesized component counts with the exact coefficient-weighted repetitions of the specified estimator. They include observable preparation and inverse preparation, lookup and uncomputation, controlled reflections, signed rotations, degree, fresh gradient preparation where used, and all 30,383,928 planned readouts. They are **conditional ideal resource estimates, not shots executed on a QPU**. They omit routing, error correction, distillation and hardware runtime. The phase-gradient circuit uses more logical qubits and more CX gates; a T reduction alone does not establish lower physical cost or quantum advantage over classical computation.

## Fixed problem and precision

The input is the previously accepted, fixed-centroid Gaussian triangle operator, rank nine with three residue observables, padded to dimension 16. Its SHA256 is `34ed3e91ee2a4d98569be2b047359c1cecfda6ad1040bec4a63e0912b4fa7ee8`. All candidates use the same real upper Cholesky factor of the **dense whitened operator**. This is not an assumed sparse protein oracle. The factor normalization is Γ = 44.023640376039786 and α = Γ/2; the strong dense normalization has α = 8.336073961228562. The degrees at 0.1, 1 and 10τ are 10, 25 and 76 for the factor, and 7, 17 and 47 for the strong dense control.

There are 18 unordered residue-pair/time outputs, response tolerance 0.002 and family failure probability 0.05. The sufficient common plan has 1,687,996 shots per output. Angle rounding and Clifford+T synthesis each receive a response error allowance of 10⁻⁶. Signed angles are rounded modulo 4π; conditional 2π phases are retained. Numerical projection, polynomial tail and arithmetic allowances are recorded separately in the protocols and receipts.

## What was checked

`protocol.json` and `source-seal.json` precede primary model outcomes. `results/semantic-tests.json` checks clean lookup, arbitrary dirty workspace restoration, relative phases, signed rotation and uncomputation. Actual deterministic Ross–Selinger gate strings are retained and independently checked at 70 decimal digits. Both six-qubit dense Clifford+T circuits were replayed as complete component operators. The primary candidate's large workspace was checked by exact small semantic tests and composition of the actual synthesized rotation matrices; no 147-qubit statevector was run. Its maximum Gaussian response error was 2.09724×10⁻⁵.

`independent-verification.json` reconstructs the primary cost, response, synthesis and protein accounting and verifies the **failed** half-cost decision. The primary walk spent 92.25% of its T gates on bit rotations. Even making all other T gates free would leave the fixed rotation subtotal at 0.5233 times the dense cost. That measured bottleneck motivated the separately recorded follow-up; the original `comparison.json` preserves the proposal as it stood before the follow-up ran.

See [the phase-gradient evidence](phase-gradient-followup/README.md) for the later protocol, emitted full-size components, low-order complete circuits, exact tiny tests, 69-qubit numerical component probes and independent cost/error audit. Its Gaussian response error upper bound is 2.10246×10⁻⁵. The result does not certify the original nonlinear protein model or biological predictions.

## Protein accounting is a separate estimate

The four copied common-revision Gaussian operators have ranks 347, 506, 172 and 934. Actual angle tables, rounded words and full bit-rotation banks were produced. Protein lookup networks, complete walks, high-order circuits and quantum responses were **not** compiled or executed. The following costs are literal implementation upper counts before compiler cancellation, not minimum bounds.

| Model | Residues | Rank | Angle bits | Clean lookup workspace | Controlled-walk T upper count | Sufficient full-matrix readouts, all three times | Complete polynomial admission |
|---|---:|---:|---:|---:|---:|---:|---|
| KRAS | 166 | 347 | 39 | 5,002 | 1,277,332 | 789,403,937,562 | Degrees 233, 737, 2329 |
| ABL | 252 | 506 | 40 | 5,130 | 1,305,176 | 628,084,806,426 | Third time exceeds degree 4096 cap |
| MYC/MAX mechanics, MYC candidates | 162 | 172 | 39 | 2,505 | 769,092 | 123,873,859,953 | Second and third times exceed cap |
| MYH7 | 694 | 934 | 40 | 10,251 | 2,299,280 | 4,642,268,045,760 | Second and third times exceed cap |

KRAS's primary estimated complete T upper total is 1.54535×10²⁰. Other complete totals are deliberately absent because not all times meet the polynomial cap. Bit precision for capped cases is conditional through degree 4096, not admission of a missing polynomial. System registers require another 20, 20, 18 and 22 qubits respectively. The small phase-gradient improvement has **not** been extended to these protein circuits.

## Replay

All current entrypoints locate dependencies relative to this `quantum/` directory. Own prior inputs, dense controls and required helpers are copied into `inputs/`, `baseline/` and `vendor/`; provenance hashes are in `inputs/provenance.json`. Historical paths in receipts are provenance, not a required original workspace. For these quantum experiments the complete curated branch is sufficient; reproducing how the earlier Gaussian operators were derived requires the separate earlier evidence.

Source-only CI needs Python 3.13 and only two pinned dependencies:

```sh
python -m pip install -r requirements-ci.txt
python ci_smoke.py
```

This test does not download assets or test the complete cost claim. Full archived verification requires the release overlay and `requirements.txt`:

```sh
python -m pip install -r requirements.txt
python verify_results.py
python phase-gradient-followup/verify.py
```

Run scientific regeneration in a **copy** of the archive because entrypoints write their declared result files. Keep the archived evidence immutable. The sequence is:

```sh
python test_semantics.py
python bounded.py dense_norm_optimal
python bounded.py dense_factor_normalization
python bounded-candidate.py select_swap
python protein_accounting.py
python verify_results.py
python phase-gradient-followup/test_semantics.py
python phase-gradient-followup/run.py
python phase-gradient-followup/verify.py
python phase-gradient-followup/full_bank_probe.py
python phase-gradient-followup/low_word_probe.py
```

The two dense runs each took about 250 seconds locally; all measured workers stayed below 0.7 GB. The guards enforce 12 GiB per worker and 600 seconds per compilation or operator check. Run sequentially to keep resource use bounded. The maximum-order complete circuits exceed the one-million-operation cap and are explicitly accounted from cached components instead of constructed. The follow-up's order-76 circuit would contain 11,321,459 operations. No uncapped full statevector or hardware job is part of these commands.

`experiment.py` preserves the original candidate indexing failure. `candidate_experiment.py` contains the documented empty-control-address array fix and is the correct candidate entrypoint. `history/` also retains early missing-dependency failures and the corrected literal-versus-compiled count assertion. Neither fix changes H, normalization, precision, candidate selection or acceptance thresholds.

## Next Phase 2 decision

The concrete engineering result supports advancing the phase-gradient loader as a **candidate**, with explicit cost limits. First compile the same method on the actual KRAS operator and compare it with dense controls that also receive specialized constant-addition and joint-rotation optimization. Charge fault-tolerant workspace, routing and state factories, retain the full error budget, and stop if the advantage disappears. Separately, the 30.4 million readouts for three observables and up to trillions for proteins make output estimation a major unresolved cost. Freeze a biologically justified selected receiver/candidate query set and implement an estimator for those exact outputs; compare against both the current sampling method and equally scoped classical computation. Amplitude estimation or reduced output counts remain proposed until implemented and checked. The protein Γ/degree penalty also needs admission before any end-to-end advantage claim.

## Publication split

`publication-manifest.json` provides exact SHA256 and byte size for every public file, separated into repository source/fixtures and release evidence. `requirements-ci.txt`, `ci_smoke.py` and `select_swap.py` are the minimum CI footprint. The complete source and protocols should also live in Git; large QPY circuits, numerical arrays and raw logs belong in a release asset that overlays this same directory. `runtime/`, caches and installed libraries are excluded. No full copied research article is included. See [NOTICE.md](NOTICE.md) for algorithm and dependency attribution.
