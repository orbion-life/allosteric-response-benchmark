# The thermal response, its quantum estimator and its error budget

This methods specification defines the quantity to be compared across numerical solvers. The sibling pilot provides the implemented response and circuit fixtures. The finite, reduced model is an approximation to protein mechanics; exact agreement between its solvers does not validate the omitted motion or biological predictions.

## Static contact mechanics defines the intervention

Each retained residue supplies an alpha-carbon position. A contact connects residues within 10 Å; consecutive mapped sequence neighbours are also retained, and sequence gaps receive no invented backbone link. The contact cutoff has 8/12 Å sensitivities. With reference position r₀, retained-coordinate vector q and harmonic-mode matrix B, r=r₀+Bq. The retained dimension is d. Rigid translations and rotations are removed. Near-degenerate modes are retained together. Omitted coordinates are held fixed rather than thermally integrated out.

For contact (i,j) with reference length ℓᵢⱼ, the energy is

\[
u_{ij}(r)=\frac{\kappa_{ij}}{8\ell_{ij}^{2}}
\left(\lVert r_i-r_j\rVert^2-\ell_{ij}^2\right)^2,
\qquad U(r)=\sum_{(i,j)}u_{ij}(r).
\]

The local observable averages the adjacent contact energies:

\[
E_i(r)=\frac{1}{d_i}\sum_{j\in\mathcal N_i}u_{ij}(r),
\qquad U_h(r)=U(r)+hE_i(r).
\]

Here dᵢ is the number of contact neighbours, distinct from the retained dimension d. A small positive h increases the weight of those contacts and remains applied after time zero. Averaging is a chosen per-residue intervention; it does not equalize every residue's geometric response or establish how a ligand acts. An unaveraged contact-sum intervention is a sensitivity.

The harmonic control is the quadratic expansion of U with the same Hessian. The distance-Hookean control uses κᵢⱼ(‖rᵢ−rⱼ‖−ℓᵢⱼ)²/2. Each control retains the same quartic Eᵢ observable, receiver set, comparison time and named normalization. The membrane spring-law precedent is [Delingette (2008)](https://doi.org/10.1109/TVCG.2007.70431); transfer to protein prediction is a hypothesis. The protein pilot uses κ=β=μ=1 and angstrom coordinates. Energy, temperature and time are not physiologically calibrated. Solvent, ligand, cofactor, DNA and modification nodes do not enter this mechanical model.

## Reversible thermal relaxation defines the response

The equilibrium distribution is π(q)∝exp[−βU(q)], where β is inverse thermal energy. On a uniform d-dimensional configuration grid with spacing δ and diffusion coefficient D=μ/β, neighbouring configurations have rates

\[
L_{xy}=\frac{2D}{\delta^2[1+\exp\{\beta(U_y-U_x)\}]},
\qquad L_{xx}=-\sum_{y\ne x}L_{xy}.
\]

Removing outward edges gives reflecting boundaries. This row generator conserves probability and satisfies detailed balance with the finite-grid Boltzmann distribution. The primary time is τ=1/(μλ₁), where λ₁ is the smallest positive harmonic eigenvalue. The protocol reports 0.1τ, 10τ and the equilibrium limit as sensitivities.

Starting in unperturbed equilibrium, the first-order response to the held perturbation is

\[
R_{ji}(t)=\left.\partial_h\langle E_j(t)\rangle_h\right|_{h=0^+}
=-\beta\left[\operatorname{Cov}_{\pi}(E_i,E_j)
-\operatorname{Cov}_{\pi}(E_i(0),E_j(t))\right].
\]

The response is a change in mean contact energy. The two covariances measure simultaneous and delayed fluctuations used to calculate it. Linear response refers to the first order in h, while U and the observables remain nonlinear. A negative R denotes lower mean receiver energy, not demonstrated inhibition. Reversible reciprocity does not demonstrate a directed biological pathway.

A named harmonic reference supplies sᵢᴴ, the standard deviation of the same quartic Eᵢ. Define

\[
C_{ji}=R_{ji}/(\beta s_i^{\mathrm H}s_j^{\mathrm H}),
\qquad S_i=\sqrt{|F|^{-1}\sum_{j\in F}C_{ji}^2}.
\]

F is the fixed functional receiver. Squaring prevents opposite response signs from cancelling. The matrix retains signed responses; S ranks the five unique eligible residues. The candidate masks and contact-label tests are specified separately in the validation protocol. A zero fluctuation scale or zero-norm vector requires explicit handling; it is not silently divided through.

## The quantum walk estimates the same finite-grid covariance

Let Π=diag(π) and eᵢ(x)=√πₓ[Eᵢ(x)−〈Eᵢ〉π]. Set

\[
H=-\Pi^{1/2}L\Pi^{-1/2},\quad \nu=4dD/\delta^2,
\quad P=I+L/\nu,\quad M=I-H/\nu.
\]

Detailed balance makes H symmetric and positive semidefinite. P is stochastic, and the spectrum of M lies in [−1,1]. The response becomes −β〈eᵢ|(I−exp(−tH))|eⱼ〉. Exact small-matrix calculations and batched sparse exponential action are classical references for this identical discrete model; the all-mode analytic harmonic calculation is a different, continuum reference.

A unitary extension of V|x,0〉=|x〉Σᵧ√Pₓᵧ|y〉 gives Uₘ=V†SWAP V. Its zero-ancilla block is M. With Q=I⊗|0〉〈0|, the walk W=(2Q−I)Uₘ satisfies the projected Chebyshev relation for Tₖ(M). For z=νt,

\[
e^{-tH}=\sum_{k\ge0}a_kT_k(M),\quad
a_0=e^{-z}I_0(z),\quad a_{k>0}=2e^{-z}I_k(z).
\]

The modified Bessel coefficients are nonnegative and sum to one. Degree K truncates the expansion. Its omitted tail ε bounds the operator error and gives a raw-response bound β‖eᵢ‖‖eⱼ‖ε. This is the positive-coefficient walk construction associated with [Marteau (2023)](https://arxiv.org/abs/2306.14993); it is used here for a classical thermal-relaxation calculation, not protein quantum coherence.

PREP preserves signs and phases of the normalized observable vectors. COEF prepares √(aₖ/s), where s=Σₖ₌₀ᴷaₖ. SELECT applies Wᵏ. Controlled preparation, selected powers and inverse preparation define the full branch overlap. A final Hadamard test assigns +1 to computational-basis outcome 0 and −1 to outcome 1. If its sample mean is χ̂ᵢⱼ, reconstruct

\[
\widehat R_{ji}^{(K)}=-\beta\left[\operatorname{Cov}_{\pi}(E_i,E_j)
-s\lVert e_i\rVert\lVert e_j\rVert\widehat\chi_{ij}\right].
\]

The estimator uses the complete unitary without postselection. The degree and Y registers may be nonzero after inverse preparation; temporary arithmetic workspace must be counted and coherently uncomputed. Zero-norm observables bypass preparation. The sibling circuit fixture tests a particular small lookup implementation; it does not establish an efficient generic protein state-preparation oracle.

## Reconstruction and sampling must preserve the ranking

Quartic contact observables share at most b=choose(d+4,4)−1−d monomials. With centered monomials φₐ, weighted vectors fₐ=√πφₐ and coefficients Aᵢₐ, the residue response is R=A Rφ Aᵀ. An exact within-grid basis reduction is accepted only with a reconstruction residual. It does not recover omitted physical coordinates.

For m independent overlap estimates, family failure probability η and overlap half-width εₐ, a sufficient shot count per overlap is ceil[2 log(2m/η)/εₐ²]. Let γᵢ=Σₐ|Aᵢₐ|‖fₐ‖/sᵢᴴ. Then the sampling contribution obeys |ΔCⱼᵢ|≤s εₐ γᵢγⱼ. This reconstruction amplification prevents overlap error from being reported as though it were automatically the residue-response error.

The provisional absolute C budget is 0.01: coordinates 0.002, grid 0.001, domain 0.001, polynomial and other deterministic errors 0.001, sampling 0.004 and device bias 0.001. These are engineering caps, not bounds on biological model error. Coordinate comparisons use common time and full-reference scales. Grid refinement fixes the domain; domain expansion fixes spacing. Boundary mass, contact strain and noncontact overlaps are reported. A failed coordinate reference prevents a faithful-compression claim even when the shortlist appears stable.

Ranking confidence additionally requires the fifth selected candidate's lower score bound to exceed every excluded candidate's upper bound. Device bias is separate and cannot be removed merely by taking more shots. Complete costs include grid enumeration, vector and coefficient preparation, every controlled walk power, inverses, routing, shots, mitigation and full-matrix reconstruction. No computational advantage follows from counting only the state qubits.
