# Audit of the proposed direction-register walk

**Date: 12 September 2026. Status: source audit and small-model algebra complete; the first native-operation gate failed, so full readout and matched compilation remain incomplete.**

The proposed direction register can encode the same discriminant as the existing two-address walk. The 21-versus-29 arithmetic is correct for the stated packed-address encoding and degree-42 overlap plan, before additional workspace. It is not a measured reduction in compiled qubits or total cost. The proposal correctly restricts this construction to the original sparse grid; a projected response operator generally does not inherit the required nonnegative stochastic transition structure.

## What the current implementation establishes

I read the released circuit and costing source directly and compared their bytes with tag `v0.3.0`, which resolves to commit `df69d29c0870ff636c56e7194de037d8693b3e19`. The current circuit prepares each complete transition row on a second address register. It applies row preparation, register SWAP, inverse row preparation and the reflection `2|0_Y><0_Y|−I`. Thus its mathematical walk is `W = R_Y V† S V`, in right-to-left operator order. Its double-controlled implementation keeps V and V† unconditional, and controls only SWAP and the reflection. This retains both loading calls on every walk invocation. Its full overlap controls all five preparation, selection and inverse components by the readout bit. The reflection's absolute global minus sign matters once the operation is controlled.

The primary costing record has 1,089 grid states, 11 bits per packed address, degree K=42, six coefficient bits, one readout bit and 29 abstract qubits. Its existing binary SELECT executes 63 walk calls and 126 preparation/inverse calls. The cheaper readout-basis plan also retains K=42. These are planning counts; the complete primary circuit is not compiled. The separate fixed-error plans use K=21 or 25 and five coefficient bits, so their two-address count is 28 rather than 29.

Source hashes:

| Source, relative to the public repository | SHA-256 |
|---|---|
| `experiments/kras/quantum/circuit.py` | `33a8711b458480cf8f566f52e2920716c348fa21790d68ff6af759b42a4b27f2` |
| `experiments/kras/quantum/model.py` | `f46e14722b674b4c4cf49bd72752aecb1ab74853b00d5929c607a910f14b6879` |
| `experiments/kras/same-model-cost/analyze.py` | `4df046f7f63ff58cda9160a7da5176f87dd6aebf049532fc3a27db2fb6fe632e` |
| `experiments/kras/same-model-cost/results/results.json` | `6d3e22af67e54a75e046969e113b3183940d5092260e1f70ac826e4e50327729` |

## Reversible direction map and boundaries

For a two-dimensional n×n grid, use labels 0=stay, 1=+first coordinate, 2=−first coordinate, 3=+second coordinate and 4=−second coordinate. Labels 5–7 are reserved. Pack address x=ni+j. A valid first-coordinate move changes x by ±n; a valid second-coordinate move changes x by ±1, with separate checks on i and j. Merely adding ±1 without the second-coordinate boundary check would incorrectly join adjacent rows.

Define S by exchanging each valid directed edge `(x,+axis)` with its reverse `(y,−axis)`. Fix every stay label, every outward-pointing invalid boundary label, every reserved label and every label attached to a padded address. Each nontrivial orbit is a two-cycle, so this full-domain permutation is bijective and satisfies S²=I. Do not flip an invalid outward label into an inward label: the inward label may already belong to a valid edge and the resulting rule need not be bijective.

The self-loop probability must be the residual `P_xx=1−Σ_{y≠x}P_xy` after summing valid transitions. Missing boundary moves receive zero direction amplitude. They must not be assigned an additional invented reflected transition or counted twice in the self-loop. Padded addresses have P_xx=1, and their direction preparation is the stay state. Physical observables have zero padded-address amplitude.

For the present logistic generator, `P_xy = [2/(νδ²)]/[1+exp(β(U_y−U_x))]` on valid neighbours. The current ν=4d/δ² bounds the exit rate and makes the residual probability nonnegative. Transition preparation therefore does not require the global partition function. This observation does not remove the equilibrium-weighted observable preparation cost.

## Exact compression identity

Let A be a unitary extension of the controlled direction preparation

`A|x,0> = |x> Σ_d sqrt(p_d(x))|d>`.

Its workspace must be clean at this interface. With J|x>=|x,0>, define T=A†SA. T is a Hermitian involution. For a valid grid edge, the signal block satisfies

`<y|J†TJ|x> = sqrt(P_xy P_yx)`;

its diagonal is P_xx. This follows by matching the unique edge label at x with its reverse at y. Fixed zero-amplitude labels contribute nothing. Thus `J†TJ=D`, where `D_xy=sqrt(P_xyP_yx)`. On the reversible physical block, detailed balance gives `D=diag(sqrt(pi)) P diag(1/sqrt(pi))=I−H/ν`. The padded block is defined as identity directly, without dividing by a zero physical padding weight.

The sandwiched operator T alone is not the walk. The required reflection is `R_0=2JJ†−I`, giving `W=R_0T`. Its projected powers obey `J†W^kJ=T_k(D)`. This is why the existing nonnegative Bessel/Chebyshev coefficient preparation can be reused for the same finite grid and time. The discriminant/reflection framework is related to Szegedy (2004); the register substitution and its resource count are derived here, not a reported result of that paper. [Szegedy, 2004](https://arxiv.org/abs/quant-ph/0401053).

## Controls, inverses and workspace

The current control factorization carries over: apply A unconditionally, a doubly controlled S, A† unconditionally and a doubly controlled R_0. If either control is false, the sequence is exactly A†A=I, including phase. If both controls are true, it is W. The inverse walk is `W†=TR_0`; reversing only the direction shift is insufficient because S is already its own inverse. In circuit order, the inverse applies R_0, A, S and A†. All five readout controls in the existing overlap estimator must remain.

A future arithmetic loader must compute rates and rotation angles, prepare the direction state and uncompute temporary values before the shift. If it leaves a record of the old address or direction, the subsequent A† at the shifted address need not erase that record. That changes the signal block. Tests must therefore inspect clean-workspace probability and false-control coherence, rather than just output probabilities on a few input addresses.

The proposed primary count is `11 address + 3 direction + 6 coefficient + 1 readout = 21`. Using two six-bit coordinates gives 22. Both omit reversible division/modulo 33 or an alternative packed-index implementation, boundary flags, carries, energy and rate fixed-point values, rotation synthesis, multi-control decomposition workspace and observable/coefficient loading workspace. A table-based controlled preparation can avoid arithmetic ancillas at the cost of gate count and classical table construction. It is an implementation proof, not evidence of efficient primary-model access. The retained SELECT still invokes 63 walks for K=42, and equivalent overlaps retain the same shot requirements under the same estimator and error allocation.

## Frozen executable canary

The next local test uses a synthetic 3×3 asymmetric quartic grid, padded to 16 addresses. It retains the logistic generator, original reflection order, signed observables, three-bit direction register and all five controlled estimator components. A lookup-based direction loader and a Gray-code X/MCX implementation of the edge permutation will be compared with the unchanged v0.3.0 two-address construction on exactly the same P, coefficient array and observables. No whole-walk matrix or matrix exponential will be inserted into either compiled circuit. Dense matrices are reference calculations only.

The complete specification, including the exact energy and observables, is frozen in `direction-register/protocol.json`, SHA-256 `4f83a3a77e90648331a1b20d9be777773be0a60fc9ea1cf1509648530e9146f6`. Its local cap is 30 minutes and 4 GB, with no paid execution or 1,089-state circuit. It requires full permutation, padding, boundary, inverse, inactive-control, absolute-phase, projected-Chebyshev and off-diagonal sign checks. Compiled all-to-all and line circuits will use the same abstract u/cx basis. Saved QPY circuits must be reloaded and checked against the classical polynomial reference within 10⁻¹⁰. The difference from the dense exponential will be recorded separately from circuit error. Failures and partial results will be retained without changing the frozen workload.

Passing this canary would establish implementation on that small lookup-loaded model. It would not establish the proposed 21-qubit primary construction, efficient reversible energy arithmetic, protein response fidelity, hardware robustness, reduced shot cost or quantum advantage.

## Executed canary outcome

The fixed 3×3 canary completed its recorded algebra checks below the 10⁻¹⁰ threshold. The full shift is bijective and involutive on 128 basis states, with 104 fixed stay, boundary, reserved or padded states. The largest recorded algebra error is 7.30 × 10⁻¹² for controlled inverse composition. Projected Chebyshev powers, absolute-phase control factorization and false-control identity checks also lie below the threshold.

The first direction-register overlap then exceeded the frozen cap of 1,000,000 total compiled operations in the abstract all-to-all u/cx basis. The guard executes before native QPY and exact-count persistence. Accordingly, only the threshold exceedance is available; the exact count and depth cannot be recovered from the saved receipt. The ten-qubit logical QPY was saved and structurally reloaded, but no native circuit or complete readout result was saved. The two-address compilation, line comparison, sign-flip execution and full persisted-statevector replay were not reached. No matched compiled saving is established.

The original attempt was stopped after 972.56 seconds because exhaustive composite-operator construction was slow. A recorded verification-only amendment enabled optimized NumPy contraction with the same Qiskit index mapping. Eight independent complex parity cases passed, with maximum difference 1.78 × 10⁻¹⁵; RuntimeWarnings were errors and finiteness was checked. The unchanged continuation ran for 771.99 seconds before the operation cap failed. The aggregate execution time was 1,744.55 seconds, plus a 1.02-second structural reload, within the shared 1,800-second budget. Peak polled RSS was 1,637,793,792 bytes. Both attempts, source snapshots and failure records are retained. No gate, model, instance, synthesis setting or criterion was changed.

This outcome establishes the tested small-model algebra and exposes an implementation resource failure. It does not refute the algebraic encoding, establish primary-model hardware feasibility or justify the 21-qubit estimate as a measured saving. The next implementation should preserve failing counts before raising its cap error; no extra campaign was run in this task.
