# Quantum cost of the evaluated KRAS response model

This experiment connects the quantum resource calculation to the **same KRAS model that produced the evaluated shortlist**. It uses the frozen 4OBE model with two harmonic coordinates, a 33 × 33 grid, biquadratic contact energy, 18 receiver residues and 126 eligible candidates. All responses retain the original all-mode harmonic standard deviations. No structure, candidate rule, energy parameter, receiver or validation label is changed.

The primary response is reproduced classically to a maximum absolute difference of **2.60 × 10⁻¹⁸**. The numerical readout rank is 12. A conservative ranking-aware sampling plan for the existing monomial-SVD basis allocates **121,876,716 shots**, before accounting for an adequate device-bias bound. A better-conditioned, also compiled monomial readout variant lowers this sufficient planning total to **37,190,400 shots**. The existing generic loading constructor also exceeds the declared local memory limit on individual transition rows. These results do **not** demonstrate primary-model quantum execution, biological advantage, physical coordinate convergence or quantum speedup.

The smaller, previously executed four-state quantum fixture is a separate implementation check. Its qubit, gate and shot counts are not substituted for this model.

## What was measured and what was derived

| Quantity | Result | Scope |
|---|---:|---|
| Valid configurations | 1,089 | Frozen primary model |
| Configurations after flat binary padding | 2,048 | Explicit encoding; 959 invariant unused states |
| Weighted monomial / physical-observable rank | 12 / 12 | Numerical SVD with recorded rank threshold |
| Shared symmetric overlaps | 78 | The evaluated full-matrix strategies; not a universal lower bound |
| Fifth score, residue 65 | 0.0000762282607077 | Classical fixed-model reference |
| Sixth score, residue 64 | 0.0000733437207896 | Classical fixed-model reference |
| Fifth–sixth gap | 0.00000288453991812 | Used only in explicitly labelled retrospective planning |
| Fixed-budget polynomial degree | 25 | Conservative tail bound and 0.0005 response allocation |
| Fixed-budget sampling shots | 12,480 | Sufficient simultaneous ideal bound; fails shortlist separation |
| Oracle-calibrated degree / shots | 42 / 121,876,716 | Sufficient ideal planning estimate for the current SVD basis |
| Abstract register count at degree 42 | 29 | Two 11-qubit configuration registers, six coefficients, one readout; extra compiler workspace excluded |
| Full complex128 statevector at 29 qubits | 8 GiB | Derived memory requirement; no full statevector constructed |
| One controlled observable preparation | 16,310 CX | Each of all 12 actual basis states compiled; all-to-all `u/cx` |
| One controlled coefficient preparation, degree 42 | 468 CX | Actual compiled primitive |
| Transition-row preparation | Stopped at memory limit | Three actual rows attempted; no invented full SELECT cost |

The top five remain 60, 69, 62, 61 and 65. The largest absolute response is 0.00307063, so a zero response would satisfy a loose 0.004 sampling error criterion. That criterion alone is therefore insufficient for this shortlist.

## Exact observable reconstruction

Let the weighted centred monomial matrix be

\[
 F_{x\alpha}=\sqrt{\pi_x}\,[\phi_\alpha(q_x)-\langle\phi_\alpha\rangle_\pi].
\]

The physical weighted observable matrix is \(E=FA^T\). A thin SVD of \(F\) defines the primary orthonormal basis \(U\), and \(D=E^TU\) gives \(E=UD^T\). Each basis vector retains its signs. Its largest-magnitude component is made positive to fix the column-sign convention. The rank cutoff is \(\max(F.shape)\,\epsilon_{64}\,\sigma_1\), approximately \(1.10\times10^{-10}\); the smallest retained singular value is 3.309. A separate SVD of the physical observable matrix confirms rank 12. The maximum reconstruction residual is below \(8\times10^{-18}\).

With the frozen harmonic standard deviation \(s_i^{H}\), define

\[
 \rho_i=\|D_i\|_2/s_i^{H},\qquad
 \gamma_i=\|D_i\|_1/s_i^{H}.
\]

The primary basis gives \(\rho_{\max}^2=0.00324214\) and \(\gamma_{\max}^2=0.0129230\). The first controls the physical response error from the polynomial tail. The second controls error amplification when shared overlaps are estimated separately. Exact algebraic factorization within a fixed grid does not recover the 490 omitted physical harmonic coordinates.

The three evaluated representations all use 78 symmetric overlaps, but their conditioning differs. At the same retrospectively chosen degree 42 and the same confidence rule, the sufficient shot totals are 121,876,716 for the monomial-SVD basis, 61,988,238 for the physical-observable SVD basis, and 37,190,400 for separately normalized monomials. Thus orthogonalizing the readout does not automatically reduce sampling cost. These alternatives use known observable coefficients, not pocket labels. Their gap-dependent shot totals remain retrospective estimates. All 12 separately normalized monomial preparations were also compiled and verified; the physical-observable SVD preparations were not compiled. No optimality claim is made.

For normalized monomials, \(U_{\rm raw}=F\operatorname{diag}(\|F_\alpha\|^{-1})\) and \(D_{\rm raw}=A\operatorname{diag}(\|F_\alpha\|)\). The columns are normalized but need not be orthogonal. The tail bound therefore retains the physical norm \(\rho_i=\|E_i\|_2/s_i^H\); it must not be replaced by the Euclidean norm of nonorthogonal expansion coefficients.

## Polynomial degree is set by a response error bound

The symmetric generator \(H\) is reconstructed from the frozen sparse arrays. With \(\nu=4d/\delta^2=9.73347\), \(M=I-H/\nu\) and \(z=\nu\tau=128\), the coefficients are

\[
 e^{-\tau H}=e^{-z}e^{zM},\qquad
 a_0=e^{-z}I_0(z),\quad a_k=2e^{-z}I_k(z)\ (k>0).
\]

The [NIST Bessel generating function](https://dlmf.nist.gov/10.35) gives the moment-generating function \(\exp[z(\cosh\theta-1)]\). Applying Markov's inequality to both tails and setting \(\theta=\operatorname{asinh}((K+1)/z)\) yields the following conservative bound, with \(k=K+1\):

\[
 b_K=\min\left\{1,2\exp\left[-k\operatorname{asinh}(k/z)+z\left(\sqrt{1+(k/z)^2}-1\right)\right]\right\}.
\]

Therefore \(\|e^{-\tau H}-\sum_{k=0}^{K}a_kT_k(M)\|\le b_K\), and the maximum normalized response error is at most \(b_K\rho_{\max}^2\). The script selects the smallest degree satisfying this criterion. The bound is evaluated in float64, not directed interval arithmetic; it is not described as a machine-certified enclosure. Independent numerical tail sums and the classical response error are checked against it.

The existing 0.001 polynomial/other deterministic allocation is divided equally here: 0.0005 for the explicit tail, with 0.0005 reserved for other numerical errors. This is the earlier fixture's engineering split, not a newly fitted biological tolerance. Degree 25 has a retained coefficient sum of 0.975748, a conservative normalized-response tail bound of 0.000466603, and an observed classical polynomial response error of 0.0000272626. Neither the grid/domain budget nor the failed physical-coordinate budget is included in this fixed-model estimator bound.

## Confidence intervals determine when to stop

Let \(s_K=\sum_{k=0}^{K}a_k\), and let \(\widehat\chi_{ab}\) estimate the normalized shared-basis polynomial overlap. Each Hadamard-test result is ±1. A simultaneous overlap half-width \(\epsilon\) for \(m=78\) pairs has the sufficient per-pair shot count

\[
 n=\left\lceil\frac{2\log(2m/\eta)}{\epsilon^2}\right\rceil.
\]

For receiver-RMS score \(S_i\), define \(\rho_F=\sqrt{|F|^{-1}\sum_{j\in F}\rho_j^2}\) and similarly \(\gamma_F\). The score half-width on the fixed model is bounded by

\[
 h_i=b_K\rho_i\rho_F+s_K\epsilon\gamma_i\gamma_F+d_i,
\]

where \(d_i\) is an independently justified device-bias score bound. The reverse triangle inequality for the RMS norm supplies this score bound. An implementation may stop only when the **smallest lower bound among its measured top five exceeds the largest upper bound among all other eligible candidates**. The decision function receives measured means and known reconstruction constants; it has no exact-ranking argument.

The proposed adaptive schedule uses \(0.0005/2^\ell\) as its tail allocation and \(0.004/2^\ell\) as its sampling allocation. It assigns \(\eta_\ell=0.05/[(\ell+1)(\ell+2)]\), whose sum is 0.05. All pairs and all inspected stages are therefore covered by the union bound. A new polynomial changes the experiment, so fresh shots are charged at every stage. There is no silent reuse of measurements from another degree.

For planning only, this artifact evaluates those intervals using exact classical polynomial means as their hypothetical centres. Ideal separation first occurs at stage 6, degree 41, after 106,348,788 cumulative sufficient shots. This is **not** an executed adaptive experiment or a guarantee that a sampled run would stop at that stage. The separate oracle-calibrated degree-42 estimate uses the exact reference ranking and limits deterministic pairwise uncertainty to 10% of each true top-five/outside gap; it then uses half the remaining admissible sampling half-width. It is clearly separated from the operational decision rule.

A uniform 0.001 device response-bias allowance also bounds each RMS score error by 0.001 and prevents separation here. An adequate, smaller measured or justified device bound is needed. No noise result from the four-state fixture is transferred to this primary model.

The declared local limits are 2,000,000 total shots, a 512 MiB full statevector, and at most 20 seconds, 1 GiB process-tree RSS and 200,000 accepted native operations per isolated compilation. The current full statevector already exceeds its cap at degree 25. Even if that limit were overcome, the shot schedule would stop before stage 4: cumulative planned shots rise from 1,492,920 at stage 3 to 6,255,444 at stage 4 while intervals still overlap. The correct output is **unresolved within the local budget**, not a certified shortlist.

## Loading and complete-circuit limitations

The reversible stochastic matrix is reconstructed as \(P=\operatorname{diag}(\sqrt\pi)^{-1}M\operatorname{diag}(\sqrt\pi)\). Row sums, nonnegativity and detailed balance are checked. Its valid block has 5,313 nonzero entries and at most five entries per row. Invariant padded states give 6,272 nonzeros. A dense 2,048 × 2,048 float64 table consumes 32 MiB; a separate dense amplitude table consumes another 32 MiB. The actual sparse valid matrix arrays occupy 68,116 bytes. These are array-size counts, not a compiler memory bound.

The existing constructor creates one row-controlled amplitude preparation for each of 2,048 rows in \(V\). Its controlled walk uses both \(V\) and \(V^\dagger\); the binary SELECT implementation invokes \(2^a-1\) walks. Thus degree 25 uses 126,976 logical row-PREP calls per overlap, and degree 42 uses 258,048. These counts include all row loading in the current construction; they are not CX lower bounds and do not establish an optimal method.

All 12 actual signed SVD basis preparations were compiled in Qiskit 2.5.2 into all-to-all `u/cx` operations. Each requires 16,310 CX, 20,403–20,405 U gates and depth 30,584–30,585. The largest phase-sensitive prepared-state error is \(6.92\times10^{-11}\), against a declared \(10^{-9}\) check; coherent inverse preparation is also checked. The 12 separately normalized monomial preparations have the same gate-count and depth ranges, with maximum phase-sensitive state error \(1.74\times10^{-10}\). Coefficient preparations at degrees 25 and 42 require 218 and 468 CX respectively. No topology routing or hardware calibration is included.

For both measured readout families, the sum of separately compiled observable/coefficient preparation and inverse modules is 33,056 CX per degree-25 overlap or 33,556 CX per degree-42 overlap. Each family's shot budget is joined only to its own measured primitive receipts. This modular subtotal excludes SELECT, readout and routing, and is not a whole-circuit gate count or a lower bound after global optimization. Attempts to compile actual central, corner and padded transition rows stopped when sampled process-tree RSS crossed 1 GiB. The receipt preserves the modest monitoring overshoot. No full transition loading circuit, native SELECT count, full circuit depth, noisy primary circuit or hardware cost has been measured.

A possible research direction is a reversible sparse transition oracle that evaluates the fixed quartic energy difference at the grid's neighbours. It could avoid the current exhaustive row constructor. However, reversible arithmetic, logistic/square-root rotations, valid-grid encoding, precision, workspace and uncomputation remain unimplemented. Preparing signed Gibbs-weighted observable states and their centring/norms remains a separate problem. A structured transition oracle alone would not establish scalable state preparation or a quantum advantage.

## Classical cost includes shared preparation

The preserved earlier cached-input receipt reports 67.7007 seconds and a peak worker RSS of 256,999,424 bytes for the complete pilot workflow. That includes preparation, 24 finite-grid cases, eight harmonic dimensions, built-in physical checks, baseline calculations and saved arrays. Shared preparation took 0.489824 seconds including process startup; its contact/Hessian/eigenbasis subsection took 0.0495524 seconds. The combined full-492-mode covariance, delayed covariance, normalizer extraction and save took 5.99114 seconds. These costs are already in the total and supply both computational routes; they must not be added twice. A pure normalizer-only timing was not measured.

The original primary grid case took 0.405528 seconds including its model diagnostics. This artifact separately times a loaded-array 12-vector classical exponential action and response reconstruction; the observed value is stored in `results.json`. It excludes the upstream structure/Hessian/normalizer and is not an end-to-end timing. Current analysis timings are single local observations with uncontrolled system load, not benchmark averages. The retained full-pilot receipt records exact inclusion/exclusion lists and input/source hashes. It excludes external reference evaluation, resampling, plots, quantum compilation and installation.

## Reproduce the analysis and bounded compilations

The frozen inputs and their SHA-256 hashes are included. The experiment performs no network calls. Python 3.12 was used. Start in this directory and use separate environments for the pinned scientific and quantum dependencies:

```sh
python3.12 -m venv .venv-analysis
.venv-analysis/bin/pip install -r requirements-analysis.txt
.venv-analysis/bin/python analyze.py --output results-replay
.venv-analysis/bin/python verify.py --repeat results-replay

python3.12 -m venv .venv-quantum
.venv-quantum/bin/pip install -r requirements-quantum.txt
.venv-quantum/bin/python compile_primitives.py --output results-replay/primitives
.venv-quantum/bin/python compile_primitives.py --raw-only --output results-replay/raw-primitives
```

`analyze.py` validates input hashes and uses the stored sparse generator and fixed normalization. Its matrix exponential is a classical reference, never a circuit component. `verify.py` compares every saved analysis array, requires exact integer/ranking identities, checks the strict response tolerance, rejects a planted meaningful perturbation, checks the Bessel bounds, and tests that the interval decision rejects insufficient precision and excessive device uncertainty. Compilation timing and whether a local cap is reached can differ by machine. The bounded compilation script refuses to overwrite an existing primitive receipt.

`results/verification.json` records a successful independent process replay on the same host. It does not claim cross-platform CI execution for this new experiment. `results/primitives/` and `results/raw-primitives/` retain successful and stopped attempts. `summarize_primitives.py` reproduces the modular cost arithmetic from the preserved receipts in `results/`. No full primary quantum shots were executed.

The upstream implementation and frozen pilot provenance belong to [Orbion's allosteric-response benchmark](https://github.com/orbion-life/allosteric-response-benchmark). This folder is prepared for public integration; its presence here does not claim that it has already been published. Code is MIT licensed. The source model derives from public protein structures; upstream source and data terms remain applicable.
