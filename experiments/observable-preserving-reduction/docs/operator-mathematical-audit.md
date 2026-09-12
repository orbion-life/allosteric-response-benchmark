# Mathematical audit of observable-preserving reduction

**Date:** 12 September 2026. This public derivation preserves the finite-model mathematics and pre-execution implementation audit. It omits the earlier internal planning/source-inventory passages; the executed protocol is `../preanalysis.json`. The completed result and replay audit is `../audit/results-independent-review.md`. A pre-execution statement below is not a later result claim.

The equal-time identity and residual inequalities require the assumptions stated below. They are derived explicitly here. Krylov approximation of an exponential action has established numerical precedent ([Saad, 1992](https://doi.org/10.1137/0729014)); this citation does not establish the specific biological model, its continuum accuracy or its cost advantage.

## 1. Define the preserved mathematical problem

Let there be \(G\) configurations with a fixed equilibrium probability vector \(\pi\). For residue observable \(f_i\), define

$$
e_i(x)=\sqrt{\pi_x}\,[f_i(x)-\langle f_i\rangle_\pi],\qquad
E=[e_1,\ldots,e_p].
$$

The discrete inner product, domain, quadrature convention, means and observables must be the same in both calculations. Let \(H=H^{\mathsf T}\succeq0\) be the symmetric finite relaxation operator. Then

$$
G_0=E^{\mathsf T}E,\qquad
K(t)=E^{\mathsf T}e^{-tH}E.
$$

Here \(G_0\) denotes the Gram matrix, not the number of configurations. To distinguish response and residual notation, write the projection residual as \(\mathcal R\). Let

$$
V^{\mathsf T}V=I_m,\quad P=VV^{\mathsf T},\quad
A=V^{\mathsf T}HV,\quad b_i=V^{\mathsf T}e_i,\quad
\mathcal R=HV-VA=(I-P)HV.
$$

The reduced delayed covariance is

$$
\widetilde K_{ij}(t)=b_i^{\mathsf T}e^{-tA}b_j.
$$

Compression preserves symmetry and positive semidefiniteness: \(A\succeq0\). It need not preserve Markov off-diagonal signs, row sums, a positive stationary vector or sparse grid adjacency. An implementation for transition probabilities on the original grid therefore cannot be attached automatically to \(A\).

## 2. Exact equal-time preservation and its numerical limitation

If \(PE=E\), then \(\widetilde K(0)=E^{\mathsf T}PE=G_0\). Conversely, equality of all Gram-matrix diagonal entries implies \(PE=E\), because

$$
G_0-E^{\mathsf T}PE=E_\perp^{\mathsf T}E_\perp,\qquad E_\perp=(I-P)E.
$$

Thus dropping a numerically dependent direction is not automatically an exact preservation claim. Its actual projection residual and normalized Gram error must be recorded. A fixed singular-value tolerance is an algorithmic rule, not a proof that discarded components are zero.

Column scaling by the original positive harmonic standard deviations is useful for numerical seed construction: \(F=E\operatorname{diag}(\sigma_i^{-1})\) has the same exact span as \(E\), while its components correspond directly to the normalized response. If the calculation uses \(F\) throughout, do not divide the final covariances by the harmonic scales a second time. Undefined or zero scales require an explicit failure or a justified constant-observable convention.

When \(PE=E\), the first delayed-covariance derivative also matches:

$$
-E^{\mathsf T}HE=-E^{\mathsf T}VAV^{\mathsf T}E.
$$

This is a useful algebra check; it is not a proof of agreement at finite times. In exact arithmetic, a block Krylov space containing \(E,HE,\ldots,H^{k-1}E\) matches the bilinear moments through degree \(2k-1\). Numerical dependency removal and incomplete blocks require checking the actual retained span before asserting that property.

## 3. Delayed-covariance error: the stated bound

Assume both \(e_i\) and \(e_j\) lie in the range of \(V\). Define

$$
x_j(t)=e^{-tH}e_j,\qquad y_j(t)=Ve^{-tA}b_j,\qquad
r_j(s)=\mathcal R e^{-sA}b_j.
$$

The approximation has defect \(y_j'+Hy_j=r_j\). Since its initial value is exact, variation of constants gives

$$
x_j(t)-y_j(t)=-\int_0^t e^{-(t-s)H}r_j(s)\,ds.
$$

Both full and reduced semigroups are contractions. Consequently,

$$
|K_{ij}(t)-\widetilde K_{ij}(t)|
\le \|e_i\|\int_0^t\|r_j(s)\|\,ds
\le t\,\|e_i\|\,\|\mathcal R\|_2\,\|e_j\|.
$$

This proves the stated inequality. It is a finite-model inequality for a specified operator and input space. The dimensionful product \(t\|\mathcal R\|\) is appropriate because \(H\) and \(\mathcal R\) have inverse-time units.

## 4. Stronger bounds from the reduced eigensystem

Define the integrated squared residual

$$
I_j(t)=\int_0^t\|r_j(s)\|^2\,ds.
$$

Cauchy–Schwarz gives

$$
|K_{ij}-\widetilde K_{ij}|
\le \min\{\|e_i\|\sqrt{tI_j},\;\|e_j\|\sqrt{tI_i}\}.
$$

For \(A=U\operatorname{diag}(\lambda_a)U^{\mathsf T}\), define \(c_j=U^{\mathsf T}b_j\) and \(B=U^{\mathsf T}\mathcal R^{\mathsf T}\mathcal R U\). Then

$$
I_j(t)=\sum_{a,b}c_{j,a}c_{j,b}B_{ab}\,\phi_t(\lambda_a+\lambda_b),
\quad
\phi_t(z)=\begin{cases}(1-e^{-tz})/z,&z>0,\\t,&z=0.\end{cases}
$$

This requires the full-space residual Gram matrix, but no additional full-space exponential. `expm1` avoids cancellation for small \(tz\). The integral is nonnegative in exact arithmetic. Materially negative numerical values indicate a failed calculation; they must not be silently clipped into apparently favorable bounds.

A further bound uses the orthogonality \(V^{\mathsf T}\mathcal R=0\). Substitute the same error identity for the left observable to obtain

$$
K_{ij}(t)-\widetilde K_{ij}(t)
=\int_{u,s\ge0,\;u+s\le t}
r_i(u)^{\mathsf T}e^{-(t-u-s)H}r_j(s)\,du\,ds.
$$

Therefore, with \(J_j(t)=\int_0^t(t-s)\|r_j(s)\|^2\,ds\),

$$
|K_{ij}-\widetilde K_{ij}|\le\sqrt{J_iJ_j}
\le \frac{t^2}{2}\|\mathcal R\|_2^2\|e_i\|\|e_j\|.
$$

The same reduced sum evaluates \(J_j\), replacing \(\phi_t\) with

$$
\psi_t(z)=\begin{cases}[tz+e^{-tz}-1]/z^2,&z>0,\\t^2/2,&z=0.\end{cases}
$$

For small \(x=tz\), use \(\psi_t(z)=t^2(1/2-x/6+x^2/24-\cdots)\). A declared implementation may take the minimum of these independently valid bounds. It must account for numerical evaluation error in each bound. The two-sided bound needs exact inclusion of both observables and orthogonality of the projection residual; it is not a generic estimate for an arbitrary fitted reduced matrix.

## 5. Projection, access and floating-point error cannot be omitted

For imperfect seed inclusion, write \(e_i=p_i+d_i\), where \(p_i=Pe_i\), \(d_i=(I-P)e_i\), and \(\delta_i=\|d_i\|\). A safe delayed-covariance bound adds

$$
\delta_i\|p_j\|+\|p_i\|\delta_j+\delta_i\delta_j
$$

to the projected-subspace propagation bound. The static Gram error is exactly \(d_i^{\mathsf T}d_j\), bounded by \(\delta_i\delta_j\). Adding these two bounds gives a conservative response bound. If the original static covariance is retained exactly while delayed covariance is projected, its zero-time response need not vanish; that choice must be disclosed and tested.

Finite-precision orthogonalization requires an orthogonality check and a measured seed reconstruction error. Symmetrizing or clipping the reduced matrix changes the numerical approximation and needs its own error accounting. In particular, replacing \(A\) with a different matrix destroys the exact identity \(V^{\mathsf T}\mathcal R=0\) unless the difference is handled separately.

Approximate equilibrium weights or observable moments add errors in \(E\). Approximate operator actions add errors in \(H\). For two symmetric positive-semidefinite operators, their exponential difference is bounded by \(t\|H-\widehat H\|\), but access to such a norm bound is itself a requirement. A residual computed with the same inaccurate operator used to construct the basis cannot certify the intended exact operator.

The spectral residual norm can be obtained from its small Gram matrix; the Frobenius norm supplies a conservative upper bound in exact arithmetic. Ordinary floating-point evaluations are numerical estimates of these mathematical bounds. A strict computer-certified statement requires outward error bounds or an equivalent validated error ledger. Report this distinction explicitly rather than calling a raw residual estimate a complete certificate.

## 6. Normalized response and ranking

With the fixed harmonic scales used by Pulsar,

$$
C_{ij}(t)=\frac{K_{ij}(t)-(G_0)_{ij}}{\sigma_i\sigma_j}.
$$

The factor \(\beta\) cancels because the physical response is \(\beta[K(t)-G_0]\) and normalization divides by \(\beta\sigma_i\sigma_j\). The finite-equilibrium observable norm \(\|e_i\|\) is generally not the original harmonic scale \(\sigma_i\). Replacing one by the other without evidence would invalidate the numerical bound.

For a fixed receiver set \(\mathcal F\), suppose each response error is bounded by \(\varepsilon_{ij}\). The reverse triangle inequality gives

$$
|\widetilde S_i-S_i|\le
\left(\frac1{|\mathcal F|}\sum_{j\in\mathcal F}\varepsilon_{ij}^2\right)^{1/2}
\equiv\eta_i.
$$

For a measured top-five set \(T\), certify separation only if

$$
\min_{i\in T}(\widetilde S_i-\eta_i)
>\max_{j\notin T}(\widetilde S_j+\eta_j).
$$

The same principle applies to any shortlist size. A uniform response bound \(\varepsilon\) implies uniform score error at most \(\varepsilon\), so a true fifth–sixth gap exceeding \(2\varepsilon\) is sufficient. The known KRAS gap can inform a retrospective planning study but must not be supplied to the operational stopping rule. Numerical representation, measurement and any other admitted error contributions must all be included. A three-node toy has no meaningful five-residue ranking test.

## 7. Pre-execution audit of the implemented operator canary

The implementation was directly read before the primary study. The reviewed `operator_study.py` has SHA-256 `13d80e28d8acaf9ca435938d8ef08668b3cf9d0e750b7b0970b6cf26f19434ef`; `preanalysis.json` has SHA-256 `530fede82261aec39165f06d0126e75e995ad131c0af0db170ce6f9c6880ada3`. This review did not execute the response study or its tests.

An initial version evolved with clipped reduced eigenvalues but formed its residual from the unclipped reduced matrix. It also added a parallel-residual correction without a mass-matrix correction for nonorthogonal numerical basis vectors. Both issues were identified before execution and repaired. The revised code consistently uses the simulated positive-semidefinite reduced matrix \(A_*\) and residual \(\mathcal R_*=HV-VA_*\).

For general numerical \(V\), write \(b_i=V^{\mathsf T}F_i\), \(p_i=Vb_i\), \(d_i=F_i-p_i\), \(M=V^{\mathsf T}V\). Here \(F\) already includes harmonic normalization. Define \(I_i,J_i\) using \(A_*\) and \(\mathcal R_*\), exactly as above. A one-sided delayed-covariance bound that does not require exact basis orthogonality is

$$
B^{(1)}_{ij}=
\min\{\|F_i\|(\|d_j\|+\sqrt{tI_j}),
       \|F_j\|(\|d_i\|+\sqrt{tI_i})\}.
$$

A corresponding two-sided bound is

$$
B^{(2)}_{ij}=a_{ij}+\sqrt{J_iJ_j}
+[\|M-I\|_2+t\|V^{\mathsf T}\mathcal R_*\|_2]\,\|b_i\|\,\|b_j\|,
$$

where

$$
a_{ij}=\min\{\|d_i\|\|F_j\|+\|p_i\|\|d_j\|,
\|d_j\|\|F_i\|+\|p_j\|\|d_i\|\}.
$$

To derive it, first compare \(F_i,F_j\) with \(p_i,p_j\). Duhamel applied to the projected vectors gives the two-residual triangular integral, a parallel term involving \(V^{\mathsf T}\mathcal R_*\), and the readout discrepancy \(b_i^{\mathsf T}(M-I)e^{-tA_*}b_j\). Contractivity then gives the displayed bound. The mass and parallel terms use coefficient norms \(\|b_i\|\), not reconstructed-vector norms \(\|p_i\|\). Adding the measured static Gram discrepancy to either delayed bound bounds the normalized response error. The revised implementation follows these formulas and takes their minimum, with an additional explicitly numerical allowance. No separate unsupported clipping correction is needed.

The grid construction retains the original three internal coordinates, quartic node energies, equilibrium weights, harmonic normalization and reversible rates. Its normalized-observable arithmetic is ordered differently from the original reference; the supplied parity test is intended to quantify that difference. Rank selection uses only calculated bounds and declared numerical gates before the full response is computed. Common setup, basis construction, residual evaluation, reduced queries and reconstruction are counted in the operator cost.

The final protocol now describes the development geometry as previously studied, rather than algebra-only. It explicitly declares the `1e-9` bound-comparison arithmetic allowance and the rule that may discard residual-block columns at a checkpoint. That rule is a truncated block expansion; no complete block-moment-matching claim follows when directions are discarded. The watchdog now records monitoring failures and terminates the owned child process group instead of waiting for an unmonitored calculation.

No remaining mathematical defect was identified in the revised projection/bound core or the directly inspected finite-model construction. This preflight disposition is conditional on the independent tests, source/protocol freeze and resource checks. It does not report a passing canary or replace the required all-case fidelity, discretization, cost and repeatability evaluation.
