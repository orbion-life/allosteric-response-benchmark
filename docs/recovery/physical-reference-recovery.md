# Project Pulsar: proposed physical-reference recovery protocol

**Status: design only; 13 September 2026. No recovery campaign has been run.** This protocol preserves the published finite-model results and failed convergence checks. It defines a new numerical experiment, including a new boundary discretization. Its thresholds and resource limits are proposed engineering decisions, not biologically validated tolerances.

The immediate objective is to obtain a defensible reference for the full affine internal-coordinate model. It is not yet possible to identify that model with a native protein basin or an exact quotient by finite rotations. Those physical interpretations remain open decisions. A converged numerical answer for the affine model would be useful for testing solvers and response approximations; it would not, by itself, establish protein realism.

## 1. The reference object must be fixed before numerical recovery

For each frozen three-node geometry, retain all three positive harmonic modes in the orthonormal matrix \(B_0\). Define

$$
r(q)=r_0+B_0q,\qquad q\in\mathbb R^3,\qquad B_0^{\mathsf T}B_0=I.
$$

The native configuration, contact graph, rest lengths and coordinate frame remain fixed. The potential and observable remain

$$
u_{ij}(q)=\frac{\kappa}{8\ell_{ij}^2}
\left(\|r_i(q)-r_j(q)\|^2-\ell_{ij}^2\right)^2,
\quad U=\sum_{(i,j)}u_{ij},
\quad E_i=\frac{1}{d_i}\sum_{j\in\mathcal N_i}u_{ij}.
$$

The proposed whole-plane equilibrium measure is

$$
\pi(dq)=Z^{-1}e^{-\beta U(q)}\,dq.
$$

It uses flat Lebesgue measure in the retained affine coordinates. For this connected, translation-free finite network, the leading quartic term is positive in every nonzero coordinate direction: vanishing edge displacements would require a common translation, which the retained subspace excludes. Compactness of the unit sphere therefore gives a positive quartic lower coefficient. This establishes confinement at sufficiently large \(\|q\|\), rather than a native-basin interpretation.

The proposed dynamics are overdamped reversible diffusion with **constant mobility \(\mu I\) in \(q\)**:

$$
dq_t=-\mu\nabla U(q_t)\,dt+\sqrt{2\mu/\beta}\,dW_t,
\qquad \mathcal Lf=-\mu\nabla U\cdot\nabla f+(\mu/\beta)\Delta f.
$$

Keep \(\beta=\mu=1\), the three prescribed stiffnesses \(\kappa\in\{1,10,100\}\), and \(t/\tau\in\{0.1,1,10\}\), where \(\tau=(\mu\lambda_1)^{-1}\) uses the native harmonic Hessian. These are model settings, not calibrated physiological units. Recovery begins with the development geometry at \(\kappa=1\); the other cases are conditional confirmation work.

Let \(s_i^{\mathrm H}>0\) be the standard deviation of the **same quartic \(E_i\)** in the full, unbounded native harmonic model. Freeze these scales across all grids, domains and nonlinear comparisons for a given geometry and stiffness. Define three separately reported quantities:

$$
G_{ij}^{0}=\frac{\operatorname{Cov}_\pi(E_i,E_j)}{s_i^{\mathrm H}s_j^{\mathrm H}},
\qquad K_{ij}(t)=\frac{\operatorname{Cov}_\pi(E_i(q_0),E_j(q_t))}{s_i^{\mathrm H}s_j^{\mathrm H}},
\qquad C_{ij}(t)=K_{ij}(t)-G_{ij}^{0}.
$$

All nine ordered node pairs remain in the tests. There is no top-five interpretation for a three-node system. Under the frozen-domain perturbation \(U+hE_i\), equilibrium initialization and reversible dynamics, this \(C\) is the normalized step-field response after dividing the derivative by \(\beta s_i^{\mathrm H}s_j^{\mathrm H}\). A field-dependent domain, coordinate basis or observable would define a different response and is excluded.

## 2. The affine plane is not a native basin or an exact rotation quotient

Removing infinitesimal rigid motions at the native structure does not remove every finite rotated configuration from an affine normal-mode plane. Eckart's rotating-frame analysis is a local vibration framework, not a justification for treating an arbitrary finite affine displacement as an exact rotation quotient ([Eckart, 1935](https://journals.aps.org/pr/abstract/10.1103/PhysRev.47.552)).

The saved assessment constructs congruent, approximately zero-energy, 180-degree in-plane rotations inside each retained three-coordinate plane. At \(\kappa=1\), their required half-widths in the **old slowest-mode-SD box convention** are 9.2475, 6.8586 and 7.5034 for the development, compact and elongated geometries. Both old half-widths, four and five, exclude them. These are directly recorded constructions; their contribution to the observed convergence failures has not been established. They nevertheless rule out interpreting agreement between two smaller boxes as adequate evidence of whole-plane coverage.

Two physical alternatives require new derivations before use:

| Alternative | Required definition and consequence |
| --- | --- |
| A native attraction basin | Define the set of coordinates whose deterministic gradient descent under the unperturbed potential reaches the named native minimum. Basin classification by descent has a statistical-mechanical precedent ([Stillinger and Weber, 1982](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.25.978)). Conditioning equilibrium on this basin and imposing reflection gives a trapped, restricted model. It does not reproduce unrestricted escape dynamics. The basin, metric, boundary and classification tolerances must be fixed before response evaluation. |
| A quotient by finite rotations | Specify the equivalence relation, a coordinate chart or atlas, the induced measure and the induced mobility. Alignment or masking alone does not justify flat measure and unchanged diffusion. Coordinate constraints can introduce metric factors in equilibrium averages ([Fixman, 1974](https://pmc.ncbi.nlm.nih.gov/articles/PMC388618/)). No corresponding Pulsar quotient measure or generator has been derived. |

A stiff-constraint interpretation needs an additional distinction: its equilibrium weights can depend on the shape of the constraining potential, so metric factors alone need not give the stiff-spring limit ([Waszkiewicz and Lisicki, 2025](https://doi.org/10.1063/5.0265585)). This is a warning against importing a constraint formula without deriving the intended limit; it does not supply the missing Pulsar rotation quotient.

For a future coordinate density \(\rho(q)=J(q)e^{-\beta U(q)}\) and mobility \(M(q)\), the reversible generator would have the divergence form \(\mathcal Lf=\beta^{-1}\rho^{-1}\nabla\cdot(\rho M\nabla f)\). Neither \(J\) nor \(M\) may be dropped without derivation. The existing Cartesian-grid walk does not automatically implement either alternative.

**Decision for this recovery:** test the explicitly defined whole affine plane by expanding reflecting Cartesian boxes. Keep the physical interpretation open. If adequate coverage is unaffordable, report that limitation and return to the model-choice decision; do not choose a smaller basin because it happens to pass.

## 3. New cell-centred grids will separate spacing from boundary location

The published implementation uses equally weighted endpoint nodes and suppresses outside edges. Its finite Markov chain is well defined. A cell-volume interpretation places the reflecting faces half a spacing beyond the endpoint centres, however, so fixed endpoint positions do not define exactly fixed faces as spacing changes. Recovery will use explicit cell-centred grids with fixed outer faces. This is a disclosed discretization change, not a reinterpretation of the old outputs.

For faces \([-L_k,L_k]\), take \(n_k\) cells, \(\delta_k=2L_k/n_k\), and centres \(-L_k+(m+1/2)\delta_k\). Cell volumes are equal within each grid and cancel in the normalized Boltzmann weights. Nearest-neighbour rates along axis \(k\) are

$$
L_{xy}=\frac{2\mu}{\beta\delta_k^2\,[1+e^{\beta(U_y-U_x)}]},
\qquad L_{xx}=-\sum_{y\ne x}L_{xy}.
$$

Outside edges are omitted to impose zero flux at the declared faces. The rate ratio is \(e^{-\beta(U_y-U_x)}\), giving detailed balance. A Taylor expansion gives the stated interior drift and diffusion as spacing tends to zero. Structure-preserving finite-volume discretization of reversible diffusions provides the relevant numerical principle; it does not certify this particular implementation or its achieved error ([Latorre et al., 2011](https://publications.imp.fu-berlin.de/896/)).

A useful, precisely equivalent coordinate rescaling is \(x_k=\sqrt{\beta\lambda_k}\,q_k\). Its Jacobian is constant, but the diffusion becomes anisotropic. With \(V(x)=\beta U(q(x))\),

$$
\mathcal L_xf=\sum_k\mu\lambda_k\left[\partial_{kk}f-(\partial_kV)\partial_kf\right],
\qquad L_{xy}=\frac{2\mu\lambda_k}{h_k^2[1+e^{V_y-V_x}]}.
$$

The proposed recovery implementation uses this rescaling to resolve each native harmonic width. The \(\mu\lambda_k\) factors are mandatory; isotropic rates in \(x\) would change the dynamics and time scale. These are still Cartesian boxes, with explicitly mapped faces in \(q\).

## 4. Harmonic calibration comes before nonlinear propagation

Replace only the sampling potential by \(U_{\mathrm H}=q^{\mathsf T}\Lambda q/2\). Retain the original quartic observables. The unbounded harmonic reference is Gaussian with

$$
\Sigma=(\beta\Lambda)^{-1},\qquad
\operatorname{Cov}(q_0,q_t)=\Sigma e^{-\mu\Lambda t}.
$$

This is an Ornstein–Uhlenbeck calibration, whose equilibrium and relaxation can be calculated analytically ([Uhlenbeck and Ornstein, 1930](https://journals.aps.org/pr/abstract/10.1103/PhysRev.36.823)). Expand each observable in probabilists' Hermite polynomials of \(x=\Sigma^{-1/2}q\): \(E_i=\sum_{|\alpha|\le4}c_{i\alpha}\operatorname{He}_\alpha(x)\). The exact delayed covariance is

$$
\Gamma^{\mathrm H}_{ij}(t)=
\sum_{1\le|\alpha|\le4}\alpha!\,c_{i\alpha}c_{j\alpha}
\exp\!\left[-\mu t\sum_k\alpha_k\lambda_k\right].
$$

An independent joint-Gaussian moment contraction, through total degree eight, must reproduce this result and the normalizers. The required moment factorization is standard Gaussian algebra ([Isserlis, 1918](https://academic.oup.com/biomet/article-abstract/12/1-2/134/193428)). Odd Hermite degrees must not be discarded: the contact-energy expansion contains cubic terms, even though the potential is expressed through a squared distance mismatch.

First compare exact moments and full finite-grid propagation on the same harmonic problem. A factorized one-dimensional harmonic generator may accelerate this calibration, but its result must match the full three-dimensional generator on a small grid. No observable projection is used. Harmonic boundary probability alone is insufficient: products of quartic observables involve eighth-order moments, whose tails can matter after most probability mass has been included.

## 5. Predeclared minimum experiment and resource limits

The following is a proposed bounded first campaign. Its caps are deliberate local resource limits, not predictions that convergence will be achieved.

1. **Freeze inputs and algebra.** Reuse the development triangle and \(\kappa=1\), with all coordinates retained. Record source, protocol and input hashes, exact normalizers and known congruent-minimum coordinates. Verify the small-grid harmonic generator against direct symmetric eigendecomposition. Do not use reduced responses to choose numerical settings.
2. **Calibrate the harmonic implementation.** In \(x\) coordinates, use faces \(A=(4,4,4)\) and spacings \(h=1/2,1/4,1/8\). Then use faces \((5,5,5)\) and \((6,6,6)\) at \(h=1/8\). These are five finite problems. They have 4,096, 32,768, 262,144, 512,000 and 884,736 states. Failure of the harmonic criteria below stops nonlinear propagation.
3. **Check whole-plane coverage before nonlinear allocation.** Transform every already known congruent minimum into \(x\) coordinates. Set each initial face distance to the smallest integer \(A_k\ge\max(4,\max_m|x_k^{(m)}|+2)\). The two-unit margin is a declared planning margin, not a tail bound. Refuse a whole-plane coverage claim if a known minimum is outside the largest tested box. Do not discard a known minimum by changing the native coordinate frame after inspecting results.
4. **Run the full nonlinear reference without reduction.** For that frozen rectangular box, run the same three spacings. Then enlarge every face distance by one and by two at the finest spacing. These are another five finite problems. Apply the resource preflight to every allocation; an oversized required case is an unresolved result, not permission to reduce the domain. Report all completed cases, including failures.

The first campaign is capped at **two hours of wall time, one worker, one BLAS thread, 4 GB aggregate resident memory and one million grid states per finite problem**, excluding environment installation. This cap covers preparation, calibration, nonlinear propagation and verification. A hard watchdog must enforce time and memory limits. No larger grid, extra stiffness or alternative basin is added after inspecting responses under this protocol. Any extension requires a new dated protocol and an explicit resource decision.

This is at most ten full finite problems; the harmonic factorization can reduce their actual cost. It can also terminate during preflight. Two successive grid refinements in three dimensions multiply node count by eight per refinement. A sparse six-neighbour generator has at most \(7G\) stored entries. With 64-bit values, column indices and row pointers, its CSR storage alone is at most approximately \(120G\) bytes; observable arrays, rates, geometry, propagation workspace and construction temporaries are additional. Preflight must include those arrays and a measured memory margin rather than equating sparse-matrix size with peak memory.

Existing Linux measurements provide a warning about propagation cost: for the development \(\kappa=1\) reference, the 4,913-state problem took 0.514 s for full propagation, the 35,937-state problem took 12.862 s, and the 68,921-state problem at the same spacing as the second took 24.640 s. These are measurements of the old endpoint-grid implementation, not forecasts for the new campaign. At fixed physical time, halving spacing increases the largest generator rates approximately fourfold as well as increasing state count eightfold. This explains why a state-count-only cost forecast is inadequate.

Only a passing development campaign permits a new, separately budgeted confirmation protocol for compact and elongated geometries and \(\kappa=10,100\). They are no longer unseen geometries after v10; they must not be relabelled as fresh held-out evidence. No protein or quantum scaling claim follows from this first campaign.

## 6. Acceptance requires component checks and coverage evidence

For every quantity below, use the maximum absolute difference over all ordered pairs and, where applicable, all three times. Keep the fixed harmonic scales; do not normalize each matrix by its own largest value.

| Check | Proposed decision rule |
| --- | --- |
| Harmonic analytic calibration | The finest-grid \(G^0\), every \(K(t)\), and every \(C(t)\) must each agree with their unbounded analytic values within 0.001. The separate grid and domain requirements below must also pass. |
| Grid refinement | On the same fixed faces, both successive spacing refinements must change each of \(G^0\), \(K\) and \(C\) by at most 0.001. |
| Domain enlargement | At the finest tested spacing, both successive face expansions must change each of \(G^0\), \(K\) and \(C\) by at most 0.001. Passing \(C\) cannot compensate for cancelling errors in \(G^0\) and \(K\). |
| Finite-operator propagation | Independent propagation on the same finite generator must agree within 0.0001 in every normalized component. On small grids use dense symmetric eigendecomposition; on larger grids use a separately implemented polynomial method with a recorded truncation bound, if affordable. If this check cannot be completed within the cap, mark independent propagation verification incomplete. It is not an independent continuum calculation. |
| Structural numerical identities | Require finite values, positive normalizers, nonnegative rates and weights, mass conservation, detailed balance, symmetry of the transformed operator, stationary \(\sqrt\pi\), and \(C(0)=0\). Apply a 10⁻¹⁰ relative operator-scale tolerance to matrix identities and a 10⁻¹⁰ absolute tolerance to normalized observable identities. Record raw residuals and denominators. No density floor or transition deletion may be introduced to hide underflow or disconnected numerical components. |
| Whole-plane coverage | All known congruent minima must lie inside the largest box. Record occupied regions, connections, boundary mass, and boundary contributions to \(E_i\) and \(E_iE_j\). A positive coverage decision additionally requires a justified tail or independent whole-plane comparison; two agreeing boxes alone do not rule out an unvisited distant region. |

Passing two refinements is **numerical convergence evidence**, not a rigorous error bound to the unbounded diffusion. In particular, small equilibrium boundary mass does not by itself bound reflected-versus-unbounded delayed covariance for unbounded polynomial observables. If numerical refinement passes but whole-plane coverage remains unjustified, report “finite-box numerical convergence; whole-plane reference unresolved”. Do not report a physically converged reference.

A primary rectangular-grid run uses no configuration mask. A future native-basin variant must classify cells as inside, outside or unresolved against a frozen unperturbed descent rule. Unresolved cells cannot be silently dropped. It must test boundary classification under refinement and compare inner and outer classifications, with a separately allocated error budget for each normalized component. The same frozen domain must be used under perturbation; an adaptive \(h\)-dependent basin would change the estimand.

## 7. The stopping decision preserves useful outcomes

There are three distinct possible outcomes. A failed harmonic calibration identifies a numerical implementation or resolution problem. A calibrated implementation with unresolved nonlinear domain coverage identifies a physical-reference problem. A converged affine-plane reference permits renewed testing of response reduction, with its setup and total query costs included, but still does not establish a native protein model.

If the specified whole-plane reference cannot be justified within the cap, stop protein fitting and keep the physical-domain decision open. Retain the exact full harmonic results and the already verified finite-model calculations as bounded comparators. The proposed Gaussian approximation remains a separate classical comparator whose nonlinear response performance has not been established. Neither fallback is evidence of a quantum advantage or validated allosteric prediction.

## Evidence binding and report wording check

This design was checked against the v10 assessment, its saved `math/verify_math.py` and `math/independent-math-results.json`, and the frozen response-operator experiment. The operator source SHA-256 is `13d80e28d8acaf9ca435938d8ef08668b3cf9d0e750b7b0970b6cf26f19434ef`; its protocol SHA-256 is `530fede82261aec39165f06d0126e75e995ad131c0af0db170ce6f9c6880ada3`. Recorded timings come from the archived Linux replay of run 34717448039. They do not replace the strict failed portability receipt.

The concise report states the proposed coordinate model and its unresolved physical interpretation. The detailed recovery remains a proposed experiment. Finite-box agreement must not be promoted to whole-plane coverage or protein validation.
