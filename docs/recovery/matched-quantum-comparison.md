# Project Pulsar: a matched quantum and classical response comparison

Date: 13 September 2026. Status: proposed protocol, not an executed experiment. This protocol extends section 5 of [the stability recovery note](subspace-stability-recovery.md). It changes no scientific code, historical results or existing acceptance rules.

The question is whether a complete quantum estimator can calculate the same signed thermal response as the strongest available classical solvers at a useful cost. The quantum target will be the **accepted original finite relaxation operator H**. A classical projection Hᵣ will be a comparator. The projected matrix will not be substituted into the direction-register circuit without a separately derived encoding.

## The operator and observable contract

Before any comparative cost calculation, record one model that has independently passed its declared physical-domain, equilibrium-measure, boundary, mesh and observable-convergence requirements. A small convenient circuit fixture cannot substitute for that accepted model. If no model qualifies, retain this protocol as proposed and stop before making a useful-quantum-cost claim.

The input package must contain the model definition, units, coordinate and residue order, domain and boundary rules, equilibrium weights π, original node observables Eᵢ, positive reference scales sᵢ, query times, sparse generator L and symmetric positive-semidefinite operator H. Hash the actual arrays and source files. For the row-generator convention, verify L1 = 0, nonnegative off-diagonal rates and detailed balance, and define H = −diag(√π)Ldiag(1/√π) on the physical support. Record the derivation of the rates if the physical model has changed.

Use exactly the same centered, normalized columns Fᵢ(x) = √πₓ[Eᵢ(x) − ⟨Eᵢ⟩π]/sᵢ in every route. The dimensionless response is

\[
C_{ij}(t)=F_i^{\mathsf T}e^{-tH}F_j-F_i^{\mathsf T}F_j.
\]

Under this convention, the unnormalized linear response is Rⱼᵢ = βsᵢsⱼCᵢⱼ. Preserve the sign and the equilibrium subtraction. A zero observable column has a known zero response and must not be passed to normalized state preparation.

Choose the first off-diagonal pair by the first two nonzero observable columns in the recorded canonical order, and use the earliest recorded positive time. Do not choose the pair using its response, signal magnitude or apparent quantum advantage. After that query passes, calculate the complete symmetric response matrix at every recorded time. Record all unique pairs before execution. Exact reduction of the observable columns is permissible only with a verified reconstruction and its amplification of errors, and both classical and quantum routes may use it.

## Three routes answer the same query

| Route | Computation and required comparison |
|---|---|
| Sparse classical propagation | Compute exp(−tH)F directly and reconstruct C. Batch columns and reuse work where the implementation permits it. Record the supplied trace, numerical method and accuracy checks. |
| Classical response-preserving projection | Construct V from the same H and F, calculate Hᵣ = VᵀHV and B = VᵀF, and reconstruct Bᵀexp(−tHᵣ)B − BᵀB. Select rank by the recorded residual-bound rule before revealing reference responses. Include unsuccessful ranks and construction costs. |
| Quantum original-operator estimator | Prepare the signed normalized F columns, approximate exp(−tH) through a valid controlled walk polynomial, estimate the overlaps, restore all norms and subtract the same equilibrium covariance. Include loading, inverses and reconstruction. |

Sparse exponential action is a substantive comparator because it calculates exp(tA)B without constructing the full dense exponential (Al-Mohy and Higham, 2011, [primary paper and implementation](https://eprints.maths.manchester.ac.uk/1591/)). SciPy exposes this operation through `expm_multiply`, including block inputs and a supplied trace ([SciPy documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.expm_multiply.html)). Its successful return is not an interval certificate. Validate the finite reference through an independent numerical method appropriate to its size and record an uncertainty allowance; a semigroup consistency check alone is insufficient.

For the walk route, derive ν such that P = I + L/ν is stochastic, and verify that M = I − H/ν has the required discriminant representation. Specify the reversible neighbor map, self-loops, padding and the unitary extension of preparation. A direction encoding requires a reversible local transition construction. Szegedy's walk framework supplies relevant spectral structure, but its search result is not evidence of a speedup for this response task (Szegedy, 2004, [primary paper](https://arxiv.org/abs/quant-ph/0401053)). A projected Hᵣ need not have nonnegative transition rates or retain local sparsity.

## Precision and statistical decisions

The launch record must assign a numerical solver error εC in the displayed C units, an independent reference allowance εref, and a family failure probability η before quantum estimates are examined. The reference allowance must be smaller than εC. The physical-model error remains a separate allowance. The numerical values must follow the accepted model's intended decision; this proposed protocol does not assign a new biological sufficiency threshold while that model remains unresolved.

For every query, require the measured discrepancy from the finite reference plus εref to be within εC. Also require a defensible deterministic-error bound and a simultaneous confidence interval at the recorded η. Allocate polynomial, preparation, arithmetic, equilibrium-subtraction, sampling and, if applicable, device-bias errors so that their total plus εref is at most εC. The classical projection's approximation and numerical errors must meet the same overall numerical budget. Report each allowance in the physical response units after reconstruction.

For direct observable loading, let ρᵢ = ‖Fᵢ‖ and let s be the retained positive polynomial-coefficient sum. An overlap half-width h contributes sρᵢρⱼh to the response interval. A conservative sufficient count for independent ±1 outcomes is

\[
n_{ij}=\left\lceil 2\log(2m/\eta)/h_{ij}^{2}\right\rceil,
\]

where m includes all distinct overlap estimates across the complete query family. This is the Hoeffding and union-bound calculation used by the inspected local estimator, not a lower bound or an optimal shot allocation. Include the norm and reconstruction factors when selecting h. If sampling is adaptive, record an error-spending rule in advance and use fresh samples when the implemented polynomial changes. Exact reference responses may evaluate the result, but may not choose polynomial degree, sample count or a favorable pair retrospectively.

A first-query interval that includes zero does not establish its sign. The complete matrix must retain signed values and intervals even if a rank calculation uses squared responses. A later protein shortlist requires intervals that separate the chosen residues from outsiders; passing a matrix tolerance alone does not establish that separation.

## Complete cost and implementation evidence

Report separate cold-start totals for the single query and the complete matrix. Charge shared model construction once to each alternative. Provide a separate, explicitly stated reuse scenario only when the same input and number of future queries justify it. For classical timings, retain at least three repetitions if they fit the declared budget, with identical thread limits, timing boundaries and all observations. Record fewer repetitions as a limitation.

For all routes, count input parsing, Hessian and scale construction, finite-state enumeration or operator access, equilibrium calculation, observable centering, preprocessing, queries and reconstruction. The projection total includes basis construction, every attempted bound and reduced solve. The quantum total additionally includes reversible rate arithmetic or lookup-table construction, state and coefficient loading, every controlled walk and inverse, arithmetic precision, workspace, compilation, routing, shots, resets, readout, calibration and any mitigation. Preserve complete native operation counts and depths **before** applying an operation cap, so a failed compile still produces a cost receipt where compilation has completed.

Save both logical and native circuits, compiler version and settings, hardware connectivity, gate set, all workspace qubits and measurement assumptions. Verify inactive-control identity with phase, controlled inverse, boundary and padded states, off-diagonal reconstruction and a sign-flipped observable test. The sign-flip check is an implementation test, not a second physical model. Compare the compiled circuit with its logical circuit and the common classical finite reference.

A simulator's runtime measures simulator cost, not quantum processor runtime. If no device is used, report native circuit resources and sufficient shot requirements as estimates. Do not divide a hypothetical hardware duration by a measured classical runtime and call the ratio an observed speedup. The smaller classical Hᵣ must remain in the comparison even when it makes quantum computation unattractive.

## Execution limits and outcomes

Before launch, record one aggregate wall-time limit, memory limit, worker count, native-operation ceiling and total-shot ceiling, including failed attempts. These are implementation limits and do not establish scientific sufficiency. No paid or hardware execution is authorized by this proposed protocol. An execution record with unset limits or unset precision values is incomplete and must not be launched.

Stop if the physical model has not been accepted, the transition/preparation construction is invalid, required complete circuits cannot be compiled within the recorded limit, error intervals cannot meet the response requirement, or the compute or shot budget is exhausted. Preserve partial results and the exact last completed stage. The existing direction-register canary already records successful small-model algebra followed by a native-operation-cap failure, with no complete readout or matched compiled comparison; that failure remains historical evidence.

The final disposition must distinguish: implementation incomplete; implementation correct at the required finite-response precision; complete cost measured or estimated; and useful cost relative to **both** classical routes. Only an actual matched execution can support a measured advantage claim. If the quantum route has no credible complete-cost regime against classical projection, retain the classical result and report that quantum usefulness is unestablished.

## Evidence record

This note directly rechecked `rebuild-v10/source/quantum-methods.tex`, `v10-preparation/direction-register/results/outcome.json`, and the public `experiments/kras/quantum/{circuit.py,model.py}` and `experiments/kras/same-model-cost/analyze.py`. These sources establish the current estimator, normalization, shot formula and retained canary failure. The future accepted model, its precision contract and execution limits remain to be supplied before execution. No new numerical or circuit result is claimed here.
