# A Gaussian variational closure is mathematically well defined, but its response accuracy must be tested

The candidate has a concrete rationale: retain all internal coordinates and replace the exponentially large equilibrium grid by a full covariance matrix. For the fixed-centroid biquadratic network, the covariance optimization is strictly convex. This removes arbitrary coordinate truncation from this particular approximation. It replaces that error with a Gaussian-distribution assumption and an additional dynamical assumption. Neither has been validated for Project Pulsar.

## The fixed-centroid Gaussian objective is exact within its trial family

Use **Σ** for displacement covariance to distinguish it from the report's normalized response matrix **C**. Let q ∈ ℝᴰ contain all internal coordinates in an orthonormal Cartesian basis B. Rigid translations and rotations are removed; any additional zero mode must be reported, not inverted. The existence result below assumes a positive-definite native harmonic Hessian K on this space.

For edge e, let aₑ be the native separation vector, ℓₑ = ‖aₑ‖ > 0, Tₑ the map from q to the relative displacement of the edge's endpoints, and κₑ > 0 its stiffness. Then δₑ = Tₑq and the original potential is

\[
U(q)=\sum_e\frac{\kappa_e}{8\ell_e^2}
\left(2a_e^\mathsf{T}\delta_e+\delta_e^\mathsf{T}\delta_e\right)^2.
\]

The native harmonic Hessian is

\[
K=\sum_e\frac{\kappa_e}{\ell_e^2}
T_e^\mathsf{T}a_ea_e^\mathsf{T}T_e.
\]

Choose q ∼ N(0,Σ), with Σ ≻ 0, and define Sₑ = TₑΣTₑᵀ. The cubic term has zero expectation because this Gaussian is centered. Gaussian fourth moments give

\[
\mathbb{E}_{\Sigma}U
=\tfrac12\operatorname{tr}(K\Sigma)
+\sum_e\frac{\kappa_e}{8\ell_e^2}
\left[2\operatorname{tr}(S_e^2)+\operatorname{tr}(S_e)^2\right].
\]

This expectation is exact for the original polynomial potential within the fixed linear coordinate space and the chosen Gaussian family. It is not a quartic-term truncation or a stochastic estimate. The mean is constrained to remain at the native position. The cubic potential still exists; its mean vanishes under the trial distribution, not under every possible physical equilibrium.

For inverse temperature β > 0, omit only covariance-independent entropy constants and minimize

\[
\mathcal F(\Sigma)=
\tfrac12\operatorname{tr}(K\Sigma)
+\sum_e\frac{\kappa_e}{8\ell_e^2}
\left[2\operatorname{tr}(S_e^2)+\operatorname{tr}(S_e)^2\right]
-\frac{1}{2\beta}\log\det\Sigma.
\]

Restoring the constant −D[1 + log(2π)]/(2β) gives the configurational Gaussian trial free energy. Its difference from the exact configurational free energy is β⁻¹KL(qΣ‖p), where p is the equilibrium density proportional to exp(−βU). It is therefore an upper bound when both distributions are defined on the same unbounded internal-coordinate space and the partition function is finite. A reflected finite-box numerical reference has a different domain; its free energy must not be compared with this unbounded-Gaussian bound without matching domains or controlling truncation.

This use of a Gaussian trial distribution follows the variational approach underlying SCHA. The published methods generally optimize centroids as well as fluctuations and can include nuclear quantum effects; the present candidate is a narrower **classical, fixed-centroid** construction. It is not the full SSCHA method (Monacelli et al., 2021; [primary paper](https://arxiv.org/html/2103.03973v1), sections II–III).

## Strict convexity gives a unique covariance, not a guarantee of physical accuracy

For a nonzero symmetric covariance variation Δ, define ΔSₑ = TₑΔTₑᵀ. The second directional derivative is

\[
\begin{aligned}
\mathrm d^2\mathcal F[\Delta,\Delta]
={}&\sum_e\frac{\kappa_e}{8\ell_e^2}
\left[4\operatorname{tr}(\Delta S_e^2)
+2\operatorname{tr}(\Delta S_e)^2\right]\\
&+\frac{1}{2\beta}
\operatorname{tr}(\Sigma^{-1}\Delta\Sigma^{-1}\Delta)>0.
\end{aligned}
\]

Each edge term is nonnegative because ΔSₑ is symmetric. The log-determinant term equals a positive multiple of the squared Frobenius norm of Σ⁻¹ᐟ²ΔΣ⁻¹ᐟ² and is strictly positive for Δ ≠ 0. Thus the objective is strictly convex on the positive-definite cone. With K ≻ 0, the linear energy controls large eigenvalues, and the negative log determinant excludes singular covariance. The minimum exists and is unique. Additional unconstrained zero modes would invalidate this simple existence argument and must be removed or separately treated.

The gradient equation is

\[
0=\tfrac12K
+\sum_e\frac{\kappa_e}{8\ell_e^2}
T_e^\mathsf{T}\left[4S_e+2\operatorname{tr}(S_e)I_3\right]T_e
-\frac{1}{2\beta}\Sigma^{-1}.
\]

Equivalently, the optimal covariance obeys

\[
H_{\rm eff}\equiv(\beta\Sigma_*)^{-1}
=K+\sum_e\frac{\kappa_e}{\ell_e^2}
T_e^\mathsf{T}\left[S_e+\tfrac12\operatorname{tr}(S_e)I_3\right]T_e
=\mathbb E_{\Sigma_*}\nabla_q^2U.
\]

This stationary equation is self-consistent. Naively iterating its right-hand side is not guaranteed to converge; use a convex optimizer with positive-definite line search and an explicit stationarity tolerance. All edge corrections are positive semidefinite. Consequently, **Hₑff ≽ K and Σ* ≼ (βK)⁻¹**. At a fixed native centroid, this candidate can only stiffen the Gaussian matrix and shrink its covariance relative to the native harmonic model. That limitation is a consequence of the stated objective, not an empirical benefit.

The problem is a **convex log-determinant optimization with quadratic covariance terms**. It is not a linear semidefinite program. Introducing auxiliary variables can express quadratic epigraphs with conic constraints, but the log determinant remains a distinct nonlinear convex term. General determinant maximization and linear SDP formulations are distinguished in the primary optimization literature (Vandenberghe et al., 1998; [author-hosted paper](https://web.stanford.edu/~boyd/papers/maxdet.html); Vandenberghe and Boyd, 1996; [author-hosted paper](https://web.stanford.edu/~boyd/papers/sdp.html)).

## Holding the centroid fixed omits a measurable force residual

The expectation of the original force need not vanish at the fixed-centroid optimum. Its negative is the excluded centroid gradient,

\[
g_m=\mathbb E_{\Sigma_*}\nabla_qU
=\sum_e\frac{\kappa_e}{2\ell_e^2}
T_e^\mathsf{T}\left[\operatorname{tr}(S_e)a_e+2S_ea_e\right].
\]

The cubic energy produces this term. A nonzero gₘ means the optimal Gaussian would still lower its trial free energy by moving its centroid if that direction were allowed. The fixed-centroid covariance is therefore not a fully relaxed Gaussian approximation. It can miss thermal mean shifts, expansion, asymmetric fluctuations and state changes. Allowing centroids to vary is a different candidate; it must not be introduced silently when the fixed-centroid test fails, and the covariance-only convexity argument does not automatically extend to joint centroid optimization.

The primary SCHA literature supplies precedent for variational treatment of strong anharmonicity in materials. Errea et al. (2014) optimized harmonic trial Hamiltonians, including positions and vibrational parameters, and applied the method to platinum and palladium hydrides. That establishes an empirical precedent for the general family, not this fixed-centroid closure or protein allostery ([primary paper](https://arxiv.org/pdf/1311.3083)).

## The proposed OU response is an additional approximation

After selecting Σ*, define the auxiliary overdamped process

\[
\mathrm dq=-\mu H_{\rm eff}q\,\mathrm dt
+\sqrt{2\mu/\beta}\,\mathrm dW_t.
\]

It has stationary covariance Σ*, and

\[
\mathbb E[q(t)q(0)^\mathsf T]
=e^{-\mu H_{\rm eff}t}\Sigma_*.
\]

Keep every original quartic residue observable Eᵢ, the original all-mode native-harmonic standard deviations σᵢ, and the original time horizon. Do not normalize by new closure-specific deviations or choose new times to improve agreement. Gaussian joint moments then give the candidate normalized response,

\[
C^{\rm G}_{ij}(t)=
\frac{\operatorname{Cov}_{\rm OU}(E_i(t),E_j(0))
-\operatorname{Cov}_{\Sigma_*}(E_i,E_j)}{\sigma_i\sigma_j}.
\]

The factor β cancels because the report defines C = R/(βσᵢσⱼ). This is the normalized stiffness response of an auxiliary OU equilibrium potential perturbed by hEⱼ **with its fitted Gaussian parameters held fixed**. It is not automatically the derivative obtained by reoptimizing Σ for U + hEⱼ. The chosen estimand must remain explicit.

A lower variational free-energy upper bound does not bound errors in quartic covariances or delayed correlations. The trial stiffness matrix is an auxiliary parameter, not a validated relaxation spectrum of the original anharmonic process. Bianco et al. (2017) explicitly distinguish the SCHA auxiliary matrix from the free-energy curvature and discuss a separate dynamical extension. This supports treating dynamics as a separate validation problem; it does not validate the OU ansatz proposed here ([primary paper](https://arxiv.org/pdf/1703.03212), sections VIII–IX and conclusion).

## The two-week experiment has an explicit stopping rule

The following is a **proposed protocol**. It must be frozen with code and input hashes before computing any new nonlinear-response comparison. The existing triangle is a development case. Two additional geometries are reserved for evaluation after the implementation is frozen. They are method-development holdouts, not independent biological systems.

| Role | Native coordinates in Å | Edges and internal coordinates |
|---|---|---|
| Development triangle | (0,0,0), (4,0,0), (1,3.5,0) | All three edges; all three positive harmonic modes. |
| Held-out compact triangle | (0,0,0), (3,0,0), (1.2,2.4,0) | All three edges; all three positive harmonic modes. |
| Held-out elongated triangle | (0,0,0), (5,0,0), (0.8,1.5,0) | All three edges; all three positive harmonic modes. |

Use κ = 1,10,100, β = μ = 1 and t/τ = 0.1,1,10, with τ = 1/(μλₘᵢₙ(K)). This defines nine geometry–stiffness conditions and 27 delayed-response matrices. Fit one covariance per geometry–stiffness condition; no time-dependent refitting is allowed. Use all residue pairs for the error statistic. There are no target-pocket labels, receiver selection or favorable-pair selection.

**Days 1–2 establish algebra and solver correctness.** Implement analytic objective, gradient and Hessian action. Verify the quartic expectation by independent Gaussian quadrature and derivatives by finite differences. The purely harmonic limit must recover Σ = (βK)⁻¹ and the existing full-mode Gaussian response. Require relative stationarity residual below 10⁻⁸, independent-start covariance agreement below 10⁻⁸ and harmonic-response agreement below 10⁻¹⁰ in these small fixtures. These are proposed numerical engineering tolerances, not literature-derived biological cutoffs. Check positive definiteness, basis-rotation invariance and the predicted Loewner covariance contraction.

**Days 3–5 establish a usable nonlinear reference on the development geometry.** First audit the coordinate domain. Infinitesimal removal of rigid modes does not globally remove finite rotations: this triangle's full linear internal span contains a 180°-rotated configuration at ‖q‖ = 8.2057 Å with unchanged edge lengths and zero energy to numerical precision. The four-standard-deviation half-width at κ = 1 is only 3.5494 Å. Local box refinement therefore cannot alone certify integration over the full unbounded coordinate model. The check and its explicit scope are recorded in `coordinate-domain-check.json`.

The first screening endpoint must be labeled as the **specified finite-box nonlinear generator with all three coordinates**, matching the existing implementation. Its response comparison is an engineering test; the unbounded-Gaussian free-energy bound does not apply to that conditional box distribution. A claim about the full unbounded nonlinear equilibrium additionally requires a domain search that accounts for remote minima, or a separately justified nonlinear rigid-motion quotient. Neither is supplied by a 4σ/5σ box check. If the desired endpoint is physical single-basin equilibrium, its basin and alignment must be defined before results, with a separate treatment of the Gaussian's support outside that basin. Do not silently change the target distribution to obtain a pass.

Reuse already frozen reference arrays when their provenance and numerical requirements match the declared finite-box screening endpoint. Use the full three-coordinate original quartic generator, never a truncated-coordinate reference. Refine the reference grid from 33³ to 65³ at four slowest-mode standard deviations and check a five-standard-deviation domain at 81³ with matching spacing. Existing evidence already shows that this refinement can fail the response tolerance. A single predeclared escalation to 129³ and the corresponding 161³ domain check is allowed only if its estimated memory fits the fixed cap. Do not continue unbounded refinement until a favorable result appears.

For reference admission, require maximum changes below 0.001 separately for grid and domain refinement, with stationary-distribution and harmonic-generator checks. Report these as observed numerical convergence tests, not a mathematically certified error bound. If the reference fails or exceeds the compute cap, mark the comparison **numerically unresolved** and do not claim closure accuracy. A failed reference is not a passing model.

**Days 6–8 evaluate the frozen closure on both held-out geometries.** Apply the same stiffnesses, times, solver settings, common normalization, reference procedure and escalation rule. Report the maximum response error across every pair for each case. Also report mean displacement, equilibrium quartic covariance, delayed quartic covariance, force residual and contact-strain distributions so that static and dynamical error cannot cancel unnoticed. The native all-mode harmonic model is the mandatory control. The two failed representation remedies remain historical context, not baselines to optimize against.

**Days 9–10 repeat and decide.** Rebuild the calculation in a fresh pinned environment and compare numeric arrays within declared tolerances. Publish all converged, failed and unresolved cases with inputs, seed-free solver histories, elapsed time, peak memory and source hashes. No geometry, stiffness, time, normalization or success threshold may be removed after results are seen.

The primary response gate is **max|Cᴳ − Cᴺᴸ| ≤ 0.002 in every admitted geometry–stiffness–time case**, using the existing Project Pulsar representation allowance. To prevent cancellation from hiding poor static or dynamic approximations, also require each normalized equilibrium and delayed covariance discrepancy to be at most 0.002. These are engineering acceptance criteria inherited or explicitly extended from the current design; the SCHA papers do not establish them for proteins. Compare the closure and native harmonic error per case, and require no deterioration larger than 0.001; this is an explicit proposed safeguard, not a claim of existing superiority.

A standardized mean-displacement diagnostic, ‖Σᴴ⁻¹ᐟ²mᴺᴸ‖/√D, is reported alongside each error. A value above 0.1 flags a material centroid restriction under a declared engineering convention. It is not a validated biological boundary. Such a flag prevents describing the closure as a relaxed equilibrium model even if a response threshold happens to pass.

If any admitted case fails the response or component gates, reject this fixed-centroid OU closure as the general representation remedy for the tested range. If a reference is unresolved, the stage is incomplete and protein-scale response claims remain deferred. If all admitted cases pass and all required cases have usable references, advance only to a broader nonlinear-network validation. Passing the finite-box screen is not a pass for unbounded equilibrium or for a physical basin whose domain has not been justified. Passing three small geometries would not establish protein fidelity, prediction accuracy or quantum advantage.

## The cost is dominated by the reference, not the 3 × 3 fit

The exploratory algebra check took approximately 0.06 seconds on the current local environment. That timing covers six tiny optimizations and algebra checks, not a protein calculation. It should not be extrapolated to full nonlinear response validation.

The current full-grid experiment recorded approximately 365 seconds and about 1.3 GB resident memory for an 81³ nonlinear case. These are planning measurements from a related implementation, not guarantees for every new geometry or machine. A 129³ or 161³ calculation contains approximately four or eight times as many states; time can grow faster than the state count because grid refinement changes the generator spectrum. Estimate allocation before running, and retain the escalation failure if the cap cannot be met.

Reserve at most **32 local core-hours, 8 GB peak aggregate memory and two concurrent reference workers** for this proposed two-week stage. The reference budget is an allocation cap, not a promise of completing every refined case. Track cumulative core time across retries. No paid compute, protein fitting or contact-label analysis is included. If the projected cost exceeds the cap, report an unresolved reference rather than weaken its gate.

At D = 492, a full symmetric covariance has 121,278 independent entries. Storing it densely is modest, but building a dense Hessian over those entries would require roughly 118 GB for one float64 matrix. A future implementation must use Hessian actions or another structured convex solver, and measure actual cost. The all-coordinate Gaussian model still has a classical polynomial-moment route to its response; this candidate provides **no quantum advantage**. If it becomes useful, it is initially a classical approximation and comparator. A separate, validated non-Gaussian task would be needed to motivate scaled quantum computation.

## What has actually been checked

The completed calculation is confined to the development triangle at κ = β = 1. Three positive-definite starts converge to the same optimum for both the harmonic and quartic objectives. The harmonic covariance error is 7.55 × 10⁻¹⁸; the quartic Gaussian energy agrees with independent tensor Gauss–Hermite quadrature to 4.44 × 10⁻¹⁶; the stationarity residual is 8.23 × 10⁻¹⁶. The numerical gradient and Hessian checks agree within 5.8 × 10⁻¹⁰ and 5.9 × 10⁻⁹. The excluded centroid-gradient norm is 1.177 in this uncalibrated model. These checks validate the written formula and small optimizer. **No delayed-response benefit has been measured.**
