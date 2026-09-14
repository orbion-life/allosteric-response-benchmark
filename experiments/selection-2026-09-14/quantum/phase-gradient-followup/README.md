# Separately frozen phase-gradient follow-up

This follow-up passes the predeclared 50% T-cost reduction test against the strongest dense implementation evaluated: **0.307543 times its estimated complete ideal T cost**, or a 69.25% reduction. The earlier SELECT–SWAP bit-rotation experiment failed at 0.563436 and remains unchanged in the parent directory. This is a scoped logical circuit improvement, not evidence of practical quantum advantage or original nonlinear protein accuracy.

`protocol.json` was frozen at 09:55:34 UTC on 14 September 2026; `source-seal.json` precedes full-size outcomes. The same rank-nine H, Cholesky factor, Γ, 33-bit angles, degrees, outputs, precision and sampling schedule are retained. Only the implementation of the data-dependent rotation changes.

## Construction and fair control

An exact fixed-width CDKM ripple adder implements phase kickback from a Fourier gradient state. Explicit controls, relative-phase Toffoli compute/uncompute, signed Ry basis changes and two clean helper qubits are included. Tiny exhaustive full-unitary checks verify addition, both target states, all word values, restored helper/gradient states and the conditional 2π sign. This follows the phase-gradient alternative in [Low, Kliuchnikov and Schaeffer, Appendix D.1.2](https://arxiv.org/html/1812.00954v2), with an explicit Qiskit CDKM implementation rather than an asymptotic oracle count.

The gradient is prepared freshly once per complete estimator shot, reused through that shot's preparations and walk powers, and discarded after readout. Its 2,099 T preparation cost is charged every shot. Its independently reconstructed preparation-state error gives a 5.21930×10⁻⁸ response error bound for the entire shot, including possible entanglement; perfect catalytic reuse is not assumed. No unpreparation is needed when the register is discarded. The probes apply its inverse only to check restoration.

The same implemented phase-gradient bank was offered to every known-angle dense rotation. At their own whole-estimator precision, the norm-optimal and matched-normalization controls require 40 and 41 bits. Their compiled banks cost 4,312 and 4,420 T gates per rotation, while every directly synthesized dense Rz costs at most 126 and 132 T respectively. Both controls therefore retain their cheaper direct synthesis. Specialized constant adders, joint optimization across rotations and globally optimal dense circuits remain unassessed; the result is not a universal minimum comparison.

## Measured evidence and remaining estimates

`results/receipt.json` records actual compiled components: a 69-qubit rotation bank with 3,556 T gates, a 182-qubit controlled walk with 66,312 T gates, and complete order-0 and order-1 estimator circuits with 30,795 and 97,107 T gates. Exact cached composition plus coefficient-weighted repetitions gives 10,845,685,943,780 expected T gates for all 18 outputs and 30,383,928 ideal planned readouts. Maximum order 76 would require 11,321,459 operations, exceeding the frozen one-million-operation construction cap. That full circuit was not built or run.

`independent-verification.json` independently reconstructs preparation error from emitted Clifford+T QPY gates at 70 digits, saved model responses, QPY gate counts, all full-plan totals and the dense choices. The Gaussian response error upper bound is 2.1024541×10⁻⁵, below 0.002. The large-workspace response is established by checked semantic composition and preparation error propagation; there was no 182-qubit full estimator statevector.

`results/component-probe-summary.json` indexes deterministic MPS checks of the **actual compiled 69-qubit bank**. Four fixed inputs test zero, positive and negative π/2, and coherent conditional 2π phase. A separate fifth input exercises mixed high and low carry bits. All pass; observed maximum bond dimensions were 4 and 7, below the configured cap of 256, with truncation threshold 10⁻²⁰. These are numerical component checks with no noise or QPU execution, not a full high-order estimator simulation. The low-word raw receipt inherited the phrase “four declared inputs” from its template; its actual `cases` list and separately frozen protocol contain exactly one input. The compact summary records this reporting correction without altering raw evidence.

## Replay and publication

Retain the parent `quantum/` directory; this follow-up imports its unchanged source, inputs, dense results and synthesized primary components. Install the parent's pinned requirements and run the commands in its README. The source-only CI smoke test covers the primary exact lookup/phase logic; the full follow-up tests and cost verifier require the archived artifacts and full dependencies. Large QPY files and raw MPS logs belong in the release overlay, while source, protocols and this explanation belong in Git. No phase-gradient protein compilation was performed.
