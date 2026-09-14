# Computing Gaussian residue responses from local ABL force factors

We reconstructed the Gaussian drift of the 252-residue ABL model from its 2,146 physical contact blocks, implemented the complete polynomial observable representation on a small ABL fragment, and compiled a quantum circuit that calculates one normalized residue response. The local force construction agrees with the fitted Gaussian model, and the small circuit reproduces its reference response. At full ABL scale, the implemented generic loading route exceeds the declared resource budget, while a cached classical calculation returns a complete response matrix in 1.29–1.79 ms after a 0.237 s eigensystem setup.

This experiment concerns the fixed Gaussian approximation and its quantum implementation. It does not establish fidelity to the original nonlinear dynamics, biological usefulness or practical quantum advantage. Its useful outcome is a constructive route from contact forces to a measured small response circuit, together with a specific access bottleneck for the full protein.

## What was measured

| Check | Result | Interpretation |
|---|---:|---|
| ABL contact factor | 2,146 local blocks; 38,628 sparse Cartesian entries; 489,292 bytes | Sparse access is available before the dense coordinate restriction. |
| Coordinate factor identity | Relative Frobenius error 2.8040 × 10⁻¹⁵ | The assembled force factor reproduces the contact stiffness. |
| Gaussian stationarity | Relative Frobenius residual 9.9994 × 10⁻¹⁴ | The fitted covariance is consistent with that force drift. |
| Full ABL Gaussian response | Maximum reduction error 2.3801 × 10⁻⁵ at four fixed times | The saved rank-1261 response agrees with the full Gaussian reference within the 0.002 allocation at these times. |
| Factor-consistent covariance | Maximum normalized response change 2.5180 × 10⁻¹³ | Enforcing this force identity changes the tested Gaussian response negligibly. |
| Small physical circuit at 0.01τ | 17 qubits, 120,840 operations, 60,250 CX gates; error 6.5940 × 10⁻⁹ | A complete ideal calculation passed the 10⁻⁶ response tolerance. |
| Small physical circuit at 0.1τ | 18 qubits, 258,802 operations | Compilation exceeded the frozen 200,000-operation cap, so this circuit was not simulated. |
| Explicit full Hermite embedding | 13,360,196,875 rows; 134,777,666,075,000 bytes | This dense embedding construction fails the 8 GiB gate. This is not a memory lower bound for every representation. |

All response errors use the original harmonic standard deviations, with β = μ = κ = 1. The four ABL times are t/τ = 0.001, 0.03, 0.1 and 1, where τ = 39.804673596275215 in the model's time units. Some of these times coincide with snapshot nodes; this is a fixed diagnostic panel, not an all-time accuracy guarantee. The selected scalar queries use canonical residues 242–245 against receiver residues 248 and 249. Full output contains every 252 × 252 entry at each time.

The measured records are [classical-summary.json](results/classical-summary.json), [quantum-summary.json](results/quantum-summary.json), [complete-cost-accounting.json](results/complete-cost-accounting.json) and [precision-checks.json](results/precision-checks.json). Raw matrices, coefficient rows, circuits and timing repetitions accompany them.

## Local contact forces determine the fitted Gaussian drift

For contact e, write its reference displacement as aₑ, squared length as ℓₑ², and fluctuating relative displacement as δₑ = Tₑq. The quartic contact energy is

```math
u_e=\frac{\kappa}{8\ell_e^2}(2a_e^T\delta_e+\delta_e^T\delta_e)^2.
```

For a centred Gaussian covariance Σ, the local displacement covariance is Sₑ = TₑΣTₑᵀ. Taking the Gaussian expectation of the contact Hessian gives the positive semidefinite 3 × 3 block

```math
M_e=\frac{\kappa}{\ell_e^2}\left[a_ea_e^T+S_e+\tfrac12\operatorname{tr}(S_e)I_3\right].
```

We took only these 3 × 3 matrix roots and combined them with the Cartesian contact incidence matrix. Each row of the resulting Cartesian factor L has at most six entries. If B₀ is the supplied dense matrix that removes rigid coordinates, then F = √μ LB₀ and Γ = FᵀF = μK. The equality βK = Σ⁻¹ is the Gaussian variational stationarity condition; we measured its residual instead of assuming that the numerical fit satisfies it exactly.

The smallest local block eigenvalue was 0.0069651. The factor identity has maximum absolute error 8.7041 × 10⁻¹⁴, and βK − Σ⁻¹ has maximum absolute entry 4.9297 × 10⁻¹². The dense coordinate projection F occupies 38,628,000 bytes; B₀ occupies 4,536,000 bytes. These transformations are counted as classical preparation and are not treated as free quantum oracles.

We also set Σₚ = (βK)⁻¹, evaluated the Gaussian response again with the same original harmonic normalization, and compared it with the response obtained from the frozen fitted Σ. This changes both the covariance and its associated Gaussian drift. The maximum full-response difference was 2.5180 × 10⁻¹³. A matrix-free exponential action through L and B₀ independently reproduced the first covariance column at all four times; full factor-consistent response matrices used a separately charged covariance eigensystem and Wick contractions.

## Quartic residue observables require polynomial sectors one through four

Diagonalize the coordinate covariance as Σ = V diag(σ)Vᵀ and write q = V diag(√σ)y, with y a standard normal vector. This covariance eigensystem is a dense classical input. It is not the final rank-1261 response eigensystem, but its cost still belongs in the comparison.

For each contact, define Aₑ = TₑV diag(√σ), lₑ = Aₑᵀaₑ, Qₑ = AₑᵀAₑ and cₑ = κ/(8ℓₑ²). The centred quartic observable can be written in normalized Hermite functions through degree four. With `sym` denoting the arithmetic average over the three distinct placements or pairings, its Wick-ordered tensors are

```math
\begin{aligned}
C_1&=4c[\operatorname{tr}(Q)l+2Ql],\\
C_2&=c[4ll^T+2\operatorname{tr}(Q)Q+4Q^2],\\
C_3&=4c\operatorname{sym}(l\otimes Q),\\
C_4&=c\operatorname{sym}(Q\otimes Q).
\end{aligned}
```

For a degree-k occupation α, the coefficient of ψα = Heα/√α! is k! Cα/√α!. Incident contact coefficients are averaged by the residue degree and divided by that residue's frozen harmonic standard deviation. The first sector is necessary because the cubic part of the contact energy has a nonzero first Hermite component.

The Gaussian relaxation rate is ωα = Σⱼ αⱼμ/(βσⱼ). If hαi is the coefficient of residue i, the exact Gaussian covariance is Σα hαi exp(−tωα)hαj. This is an explicit constructive formula for an individual coefficient; it does not, by itself, prepare a coherent superposition of all coefficients or their row norms.

For saved snapshot nodes sₐ and whitening matrix W, the embedding into the reduced space is

```math
\Psi_{\alpha r}=\sum_{a,i} e^{-s_a\omega_\alpha}h_{\alpha i}W_{(a,i),r},
\qquad F^{\rm reduced}_{\alpha r}=\sqrt{\omega_\alpha}\Psi_{\alpha r}.
```

In exact arithmetic, ΨᵀΨ = I and the reduced generator is Ψᵀdiag(ω)Ψ. We implemented 64 prespecified ABL rows across all four sectors and saved the coefficients, embedding rows and factor rows. Independent conditional Gaussian quadrature agrees with their observable coefficients to 3.9465 × 10⁻¹⁶ and their weighted reduced-factor rows to 3.2882 × 10⁻¹³. Selected rows verify that the formula can be evaluated; they do not verify the full ABL embedding identity or construct its coherent row-norm preparation. A separate complete small-system check verifies the finite polynomial construction against direct Gaussian kernels.

For d = 750 coordinates, the complete centred polynomial space has binomial(754,4) − 1 = 13,360,196,875 dimensions. An explicit float64 table for Ψ with 1,261 columns occupies 134.78 TB in decimal units. The corresponding weighted reduced factor table has the same shape and storage. The underlying physical lowering factor is a different object. No argument here requires every possible algorithm to materialize either table; a compact coherent representation could change the cost, but it has not been constructed in this experiment.

The implemented generic binary-tree table loader would use 46 qubits and 35,184,372,088,832 rotation slots for this full reduced factor. It fails the frozen 20-qubit simulation and 8 GiB memory gates before allocation. A coefficient routine that computes one row classically is insufficient to remove these limits: coherent access to coefficient values, row norms, observable amplitudes and the dense whitening transform remains to be designed and costed.

## A physical ABL fragment completes the quantum response calculation

The bounded circuit uses the first three actual ABL Cα positions, canonical residues 242–244, their inherited contacts and the first two positive harmonic modes. We fitted its own Gaussian covariance and retained its own harmonic normalization. The complete observable basis has 14 Hermite functions across degrees one through four.

We constructed the polynomial factor directly from the physical coordinate factor using the normalized Hermite derivative identity. Its entries are

```math
D_{(r,\beta),\alpha}=(FV)_{ra}\sqrt{\alpha_a}
\quad\text{when }\beta=\alpha-e_a.
```

The factor has 90 rows and 14 columns. Its Gram matrix agrees with the covariance-derived polynomial generator to 4.3219 × 10⁻¹² in maximum absolute entry. All 14 columns of the compiled factor block and all 14 projected walk columns were independently simulated without constructing a dense unitary. Their maximum errors were 1.6098 × 10⁻¹⁵ and 4.3770 × 10⁻¹⁵. The loader used 12 qubits and 4,093 operations; its controlled walk used 13 qubits and 17,246 operations.

The loader normalization was the actual factor Frobenius square γ = 73.70669891726268, which is 7.0265 times the generator spectral norm. With M = I − 2DᵀD/γ and a = γt/2, a positive linear combination of Chebyshev walk powers implements

```math
e^{-tD^TD}=e^{-a}\left[I_0(a)I+2\sum_{k\ge1}I_k(a)T_k(M)\right].
```

We truncated the sum only when its norm-restored tail bound was at most 2.5 × 10⁻⁷, subject to degree ≤32. The circuit prepared the actual normalized residue coefficient vectors. If αLCU is the retained coefficient sum, L = ‖hᵢ‖‖hⱼ‖ and p₊ is the Hadamard X = +1 probability, the restored response is

```math
C_{ij}(t)=L\alpha_{\rm LCU}(2p_+-1)-h_i^Th_j.
```

At 0.01τ, the circuit used degree six and seven controlled walk calls. For residues 242 and 243, L = 0.34250273876727916 and the static covariance was 0.27223529951917896. The ideal probability was 0.8777057029954194, yielding a restored response of −0.01350483324977414, compared with −0.013504839843795358 from the fitted Gaussian reference. The 17-qubit circuit contained 120,840 native `u`/`cx` operations, including 60,250 CX gates, and its ideal CPU simulation took 15.48 s. These gate counts are compilation results in the declared abstract basis, not a routed device schedule or an estimate of quantum wall time.

At 0.1τ, degree thirteen required fifteen controlled walk calls. Compilation produced 18 qubits and 258,802 operations, including 129,046 CX gates, exceeding the frozen 200,000-operation limit. The circuit was retained as evidence but not simulated. No hardware or noisy execution was performed in this experiment.

## Classical calculations provide the matched performance reference

Every method below uses the same frozen ABL observables, original harmonic normalization and four query times. The reduced methods calculate the same rank-1261 Gaussian response. The full analytic Gaussian calculation supplies its reference. Cheap timings exclude one warm call and report the median of three timed repetitions; expensive full kernels and exponential actions were measured once. These are technical timing repetitions on one CPU host, not independent scientific replicates.

| Calculation | One scalar query | Eight selected entries | Full 252 × 252 matrix |
|---|---:|---:|---:|
| Cached reduced spectral sum | 10.9–12.8 µs | 69.5–77.0 µs | 1.29–1.79 ms |
| Matrix-free reduced exponential action | 0.041–5.91 s | 0.039–6.33 s | 0.645–23.50 s |
| Full analytic Gaussian contractions | 1.75–2.14 ms | 14.1–16.5 ms | 5.17–5.51 s |

The eight-entry exponential action shares two receiver right-hand sides. The spectral setup costs 0.2372 s for a fresh eigendecomposition and projection of all observable columns. The analytic Gaussian setup costs 0.06645 s; its shared full zero-time kernel costs 5.1488 s. The timed scalar and eight-entry analytic routines recompute their zero-time scalar references, which can also be cached. Exact repeated queries may be cached on both classical and quantum routes.

The matrix-free reduced calculation agrees with the spectral calculation to 8.7319 × 10⁻¹⁴. “Matrix-free” here means exponential action through a supplied dense reduced matrix without forming its exponential; it does not mean that the reduced input matrix is sparse. The physical force action separately uses the sparse Cartesian contact factor and the dense B₀ restriction.

Fresh common preprocessing measured 0.0281 s for harmonic stiffness assembly, 0.0878 s for its dense eigendecomposition and 2.7705 s for the Gaussian covariance refit. The refit reproduces the frozen covariance within 2.4869 × 10⁻¹⁴. Snapshot-kernel construction and whitening were precomputed common inputs and were not retimed on this CPU. Historical A100 construction receipts are included separately, so these heterogeneous measurements must not be added into a single same-machine end-to-end benchmark.

## Full quantum cost remains an access and readout problem

An alternative explicit fallback is to form a dense Cholesky factor of the saved reduced generator. We measured this step at 0.04684 s with 12,720,968 bytes of factor storage and a Gram error of 4.0856 × 10⁻¹⁴. This is a valid dense factor, but it performs the reduced classical construction first and therefore does not supply the desired inexpensive physical oracle.

The generic table loader for that dense fallback needs 23 qubits and 4,194,304 rotation slots per loading unitary before optimization. Its Frobenius normalization is 12,663.5261, or 422.43 times the reduced spectral norm. Both the qubit and operation counts exceed the frozen simulation gates, so we did not construct this full circuit. The coefficient table for one normalized reduced observable has 2,048 padded entries and 2,047 binary-tree rotation slots; coherent loading and norm restoration remain part of every uncached scalar query.

For an explicitly labelled planning allocation of 0.0005 to semigroup approximation, this normalization gives Chebyshev degrees 55, 300, 548 and 1,732 at the four ABL times. In the named generic binary-control implementation, the loading rotations alone total approximately 5.28 × 10⁸ to 1.72 × 10¹⁰ per scalar-query circuit before additional controls, reflections, observable preparation, routing, fault tolerance and repetition. These are formula-based extrapolations for this implementation, not experimentally compiled full-protein costs or lower bounds on alternative algorithms.

A separate planning budget allocates an absolute 0.001 response error to statistics. Using the largest measured observable norm product and an ideal independent Bernoulli Hoeffding bound gives a sufficient 5,839,223 shots for one scalar query at 95% confidence. Simultaneous 95% coverage of the 127,512 symmetry-distinct entries in all four response matrices gives a sufficient 24,448,045 shots per entry, or approximately 3.117 × 10¹² shots overall. These conservative bounds assume no hardware bias and are not minimal costs. Adaptive amplitude estimation would require a separate coherent implementation and coverage analysis.

Classical rounding of the dense loading factor and reduced observable columns to float32 changes the four response matrices by at most 1.822 × 10⁻⁸. Float16 changes them by at most 1.460 × 10⁻⁴, although its conservative perturbation bound becomes too loose to certify the long-time allocation. This coefficient test does not validate hardware angle precision or noise.

The next constructive step is to encode contact factors, polynomial coefficients, row norms and observable amplitudes without expanding the full Hermite table or forming the dense reduced generator first. Such a circuit must be compared with the cached spectral calculation and the direct Gaussian contractions under the same query workload. The present results establish a small implementation and a full-protein diagnostic; they do not support a practical quantum speedup claim.

## Reproduction and provenance

For reproduction, download and extract `pulsar-physical-factor-access-2026-09-14.zip` from the [v0.10.0 release](https://github.com/orbion-life/allosteric-response-benchmark/releases/tag/v0.10.0). The complete archive includes the raw input arrays and QPY circuits. The Git tree contains the source, metadata, receipts and manifest; its large numerical and circuit files are supplied through the release asset.

The package contains `source/`, the exact inherited `vendor/` files, four frozen ABL input arrays, raw results and compiled QPY circuits. Source paths resolve relative to the experiment directory. Python 3.13.5, NumPy 2.4.6, SciPy 1.16.0, Qiskit 2.3.1, Qiskit Aer 0.17.2 and threadpoolctl 3.6.0 were used on an arm64 macOS CPU. Each runner limits native numerical libraries to four threads. Peak resident memory was 568 MB for the classical diagnostic and 444 MB for the quantum diagnostic. No GPU, cloud or QPU execution was used in these runs.

Run the following from an extracted copy to preserve the archived evidence:

```sh
python source/run_classical.py
python source/run_quantum.py
python source/run_costs.py
python source/check_precision.py
python checks/verify.py
python checks/verify.py --quantum
python checks/verify.py --abl
```

The local protocol and four stage-specific amendments were frozen before their corresponding measurements. The public export removes private absolute paths from metadata while preserving every scientific input, source byte and result. `provenance/public-export.json` records the original local protocol hash and the exported protocol hash; this is an explicit metadata transformation, not a claim that their bytes are identical. `MANIFEST.sha256` checks the exported files. The original input hashes remain unchanged.

The portable verifier loads unchanged independent audit programs and overrides only their data and output directory variables. It rechecks physical dimensions, normalization, response arithmetic and resource formulas; `--quantum` additionally performs direct Gaussian quadrature and replays the saved small QPY circuits. The `--abl` option checks the complete physical-force blocks and all 64 selected ABL coefficient rows by independent Gaussian quadrature. New receipts go into `checks/replay/`, preserving the archived results. The accompanying three-coordinate triangle files are provenance inputs for a separate complete 34-term Hermite identity audit, not the two-coordinate ABL circuit fixture.

`vendor/multiplexed.py` retains a legacy helper that materializes an operator matrix; the new runners never call it. All quantum checks use individual statevectors and projected amplitudes. The retained later-time QPY is a compiled failure record and is not evidence of completed simulation.
