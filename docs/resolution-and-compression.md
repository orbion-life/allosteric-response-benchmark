# Choosing resolution and compression from the response

Status: research and proposed protocol, 2026-09-12. No new experiments were run for this note. The numerical allocations below are proposed engineering tolerances; no cited paper validates 0.01, or any fixed coordinate count, for protein allostery. Primary-source support and applicability limits are distinguished below.

## Recommendation

Choose the least expensive tested configuration that meets an output-based error ledger and resolves the shortlist relative to a named larger reference. Keep four distinct objects: the physical contact model; the restriction `r = r0 + B_d q`; the finite domain and grid; and the numerical or quantum evaluation of that grid. Accurate evaluation of the last object does not validate the preceding three.

Use one harmonic reference for `s_i^H`, one set of physical times, the same fixed observables `E_i`, receiver `F`, candidate universe and energy parameters throughout each comparison. In particular, do not recompute the normalization or primary horizon separately at each dimension and interpret their normalized difference as coordinate convergence. Export raw response alongside normalized response. This is a mathematical comparison requirement, not a literature-derived biological threshold.

## Three primary sources and exactly what they support

1. **Celik et al. (2008)**. *Procedure for estimation and reporting of uncertainty due to discretization in CFD applications*. Journal of Fluids Engineering 130, 078001. [DOI](https://doi.org/10.1115/1.2960953); [author institutional record](https://pure.kfupm.edu.sa/en/publications/procedure-for-estimation-and-reporting-of-uncertainty-due-to-disc/). The original paper's recommended procedure, pp. 1–2, uses three grids, the actual grid-spacing ratios, observed order, Richardson extrapolation and a grid-convergence index. Equation 7 uses the safety factor 1.25. This is a discretization-uncertainty estimate; it excludes modelling error. See also the [ASME numerical-accuracy guidance](https://www.asme.org/wwwasmeorg/media/resourcefiles/shop/journals/jfenumaccuracy.pdf) and [NASA grid-convergence worked example](https://www.grc.nasa.gov/www/wind/valid/tutorial/spatconv.html). **Scope:** transfer the refinement procedure to the diffusion output, not a CFD tolerance, a 95% coverage claim, or a rigorous continuum certificate.

2. **Grepl and Patera (2005)**. *A posteriori error bounds for reduced-basis approximations of parametrized parabolic partial differential equations*. ESAIM M2AN 39, 157–181. [Publisher PDF](https://www.numdam.org/article/M2AN_2005__39_1_157_0.pdf), [DOI](https://doi.org/10.1051/m2an:2005006). Directly read pp. 157–160, 170–171 and §5.1: reduced-basis Galerkin approximation, residual/output bounds and greedy enrichment until a prescribed error tolerance is met. Their reference discretization must itself be sufficiently rich. Their stated coercivity and other operator assumptions must be checked. **Scope:** supports output-guided adaptation within a common state space. It does not certify a restriction of the protein's physical coordinates or establish continuum accuracy of the reference.

3. **Lelièvre and Zhang (2019)**. *Pathwise estimates for effective dynamics: the case of nonlinear vectorial reaction coordinates*. Multiscale Modeling & Simulation 17, 1019–1051. [DOI](https://doi.org/10.1137/18M1186034); [author preprint](https://arxiv.org/pdf/1805.01928). Directly read preprint pp. 2–4 and 6–8: effective coefficients use conditional expectations, and pathwise estimates require regularity, integrability and conditional Poincaré assumptions. The latter quantifies mixing along omitted directions. **Scope:** closely related reversible diffusion theory explains what a physical reduction argument needs. The present slice `U(r0+B_d q)` does not integrate over omitted coordinates; these theorems cannot be attached to it without a new argument.

## Concrete selection protocol

### 1. Declare the reference and accuracy contract first

For a tractable larger coordinate model `B_D`, name its domain, reference discretization, observable normalization and numerical uncertainty. Define the initial `0.01` absolute entrywise response target relative to this declared reference, not to the unknown protein response. Use a provisional additive ledger:

| Contribution in normalized C units | Allocation | Selection rule |
|---|---:|---|
| Restricting D coordinates to d | 0.002 | Resolved pilot comparisons must give a discrepancy, including reference-comparison uncertainty, no larger than this allocation for the declared output family. |
| Production grid | 0.001 | The finest accepted refinement uncertainty estimate must fit this allocation. |
| Finite domain | 0.001 | The domain test below must fit this allocation; report whether it is a stabilization indicator or an actual bound. |
| Propagation, preparation, means, normalization, reconstruction and algebraic approximation combined | 0.001 | Their propagated output bounds or measured errors must sum below this allocation. |
| Measurement | 0.005 | Simultaneous 95% half-width over the prespecified family; an individual 95% interval is insufficient. |
| **Total** | **0.010** | Add systematic contributions; do not combine them in quadrature without a justified probabilistic model. |

The split is a design choice to be recorded before held-out outcomes. It makes resolution and compression decisions testable; it is not an optimal allocation established by data. If the larger reference is inaccessible or unresolved, label the coordinate term unresolved and state any achieved 0.01 numerical target **conditional on the fixed d-coordinate model**. Never replace the missing coordinate term by zero. Biology beyond `B_D` and contact-law bias remain outside this ledger. The 95% statement belongs to the measurement component only; an empirical grid or domain indicator does not turn the combined interval into a certified 95% error bound.

### 2. Choose the grid independently of dimension and domain

At fixed `d`, physical domain and model, halve spacing on genuinely nested grids, for example 17, 33, 65 and 129 endpoint-inclusive nodes per axis. These are an illustrative refinement sequence, not a promised production size. Compute every entry needed for the claimed response accuracy and every score used for selection, across all compared energy laws and declared times. A pilot subset can triage cost but cannot establish a whole-matrix claim.

For each scalar output `c`, three equal-ratio grids give observed order `p = log2(|c_coarse-c_medium|/|c_medium-c_fine|)` only when successive differences have the same sign and are above arithmetic/solver error. Do not assume second order. Estimate the finest error with `1.25 |c_fine-c_medium|/(2^p-1)` when a positive-order asymptotic pattern is supported. Use the fourth grid to check whether the extrapolated limits from the two triples agree within the 0.001 allocation; require both the estimated error and extrapolated-limit discrepancy to fit that allocation. If differences oscillate, the order is unstable, or the fitted limit fails this check, refine or record unresolved grid error. A disappearing difference near roundoff is not automatically convergence.

Choose the smallest grid in the validated sequence meeting the rule for the entire declared family, then check it again at the final expanded domain. This is a conservative proposed operational use of Celik et al., not a universal error theorem. The current toy's 4, 8, 16, 32, 64 endpoint-inclusive node counts are not nested: spacing is proportional to `1/(n-1)`. Their observed order requires the actual unequal ratios; `log2` cannot be used literally on those counts.

### 3. Choose the domain at fixed spacing

Increase the physical box while holding spacing fixed, adding grid nodes rather than stretching the same grid. Start from the declared harmonic scale and expand in one harmonic-standard-deviation increments. Accept provisionally only when two consecutive expansions change every declared normalized response by at most 0.001 and do not reverse a resolved shortlist decision. Recheck grid resolution in the accepted box. The two-expansion rule is a prespecified diagnostic, not a theorem bounding the remaining infinite-domain error.

Record equilibrium boundary-shell probability and relevant observable moments as warning diagnostics. Shell probability alone cannot bound the response: the observables grow with displacement, and a reflecting boundary changes dynamics. A rigorous infinite-domain claim would require a separate bound on the relevant weighted moments and boundary-induced dynamical error. Without it, call this **domain stabilization**, not certified truncation error.

### 4. Choose d by output discrepancy, with a clearly bounded claim

Construct nested nonrigid `B_d` from the apo/reference geometry, in an order fixed before biological labels are exposed. Start with the proposed small circuits, then add the next mode or prescribed block; keep near-degenerate eigenspaces together to avoid a basis-dependent selection. For each candidate d, first resolve its own grid/domain. Resolve a named larger D reference separately. Compare both using the D-reference harmonic normalization and identical physical times.

For a resolved pilot output pair `c_d*`, `c_D*` with validated error bounds `u_d*`, `u_D*`, the triangle inequality gives a conservative coordinate discrepancy bound relative to D of `|c_d*-c_D*| + u_d* + u_D*`. If the `u` values are refinement estimates instead, label the result an estimated discrepancy. Choose the smallest d for which this quantity is at most 0.002 for all declared outputs and shortlist intervals remain separated. Add another reference block, when tractable, to test whether the D reference itself changes the conclusion. Do not average a large local error away with a Frobenius norm.

This criterion may reject every affordable d; that is a legitimate feasibility result. The current 2–4-coordinate tests do not yet pass it. Modal variance retained, matrix rank, or top-five overlap alone does not bound response error. A residual-driven reduced-basis alternative can strengthen numerical compression on a fixed operator, but it requires its own implementation and assumptions; it is not already part of the verified toy.

### 5. Distinguish physical restriction from compression on a fixed operator

For the fixed reversible discrete operator `H >= 0`, let `u_j(t)=exp(-tH)e_j`, and let `u_tilde_j` be an embedded approximation in the **same** space. With residual `r_j = d u_tilde_j/dt + H u_tilde_j`, semigroup contraction and the variation-of-constants identity give

`||u_j(t)-u_tilde_j(t)|| <= ||e_j-u_tilde_j(0)|| + integral_0^t ||r_j(s)|| ds`.

If the equal-time covariance term and observable vectors are exact, the corresponding normalized response error is at most this bound times `||e_i||/(s_i^H s_j^H)`; beta cancels from C. Errors in those vectors, means and normalization need additional terms. This is a direct mathematical derivation, not a theorem quoted verbatim from Grepl and Patera. Residual quadrature must be bounded to make the result a bound; sampled residuals alone produce an indicator. The derivation is valid for finite times without a spectral gap. Equilibrium needs a separate stationary calculation. It does not bound omitted physical coordinates that change the operator and equilibrium distribution.

This supplies a concrete adaptive choice for any later numerical basis compression: enrich using the largest output-residual contribution, recompute the bound, and stop only below the assigned ledger component. If an exact monomial representation is used without truncation, verify its algebraic identity and retain the numerical reconstruction error; its small rank says nothing about discarded physical coordinates.

### 6. Make the rank decision follow the error

For `S_i = ||C_Fi||_2/sqrt(|F|)` and simultaneous entrywise bounds or covered half-widths `epsilon_ji`, the reverse triangle inequality gives `eta_i = ||epsilon_Fi||_2/sqrt(|F|)` as a score uncertainty. Use intervals `[max(0,S_i-eta_i), S_i+eta_i]`. A selected five-residue set is resolved only when its smallest lower bound exceeds every excluded candidate's upper bound. With uniform 0.01 entrywise bounds, a point-score fifth/sixth gap greater than 0.02 is sufficient; a smaller gap is not automatically evidence of a wrong ranking, but it remains unresolved. If a new error source is only an empirical sensitivity estimate, the rank statement has that same limitation. Allocate additional precision to candidates whose intervals overlap, while accounting for sequential selection in measurement coverage.


## What the present data show

At κ=100, changing n=32 to n=64 changes C43 by 0.0046207 in the harmonic model and 0.0027175 in the biquadratic model. At fixed spacing, expanding the harmonic domain from four to five slow-mode standard deviations changes C43 by 0.0067626, despite an outer-shell probability near 4.18×10⁻⁵ at extent four. The d=4 harmonic reference still changes by about 0.0157394 between n=12 and n=16 under common d=4 normalization. These observations motivate separate checks; they do not establish successful selection or compression validity. All values can be recovered from `reference/sensitivity-results.json`.
