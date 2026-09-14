# Spectral encoding and simultaneous confidence

For the saved positive Gaussian reduced operator, write H = VΛVᵀ and Z = VᵀB. All eigenmodes are retained. A residue vector is ψᵢ = Zᵢ/‖Bᵢ‖. A multiplexed flag rotation has its zero-flag block Fₜ = diag(exp(−tλ)). The explicitly constructed Hadamard circuit applies the controlled preparation of ψⱼ, controlled flag rotation and inverse controlled preparation of ψᵢ. Its readout-zero probability is

\[
a_{ijt}=\tfrac12\left(1+\frac{B_i^T e^{-tH}B_j}{\|B_i\|\|B_j\|}\right),\qquad
C_{ijt}=\|B_i\|\|B_j\|(2a_{ijt}-1)-B_i^TB_j.
\]

These identities are independently checked against matrix exponentiation and actual representative native circuits. The representation needs the complete classical eigendecomposition and transformed vectors. Those same arrays already permit exact classical responses through weighted dot products; this construction establishes no quantum advantage over that comparator.

## Adaptive confidence allocation

Amplitude amplification gives success probability sin²((2k+1)θ), with a = sin²θ. The monotonic-branch principle is described by [Grinko et al., *Iterative quantum amplitude estimation* (2021), Section III](https://arxiv.org/html/1912.05559v3). Our implementation is a separately specified variant: it uses fresh independent batches, a summable error budget and a finite round cap. We do **not** transplant that paper's particular asymptotic constant or termination theorem to this modified implementation.

Let J be the total number of requested outputs and δ = 0.05/J. At adaptive round r, the previously observed history determines the multiplier m = 2k+1, a monotonic branch, and the new batch size. The scheduler does not receive the true amplitude. Conditional on that history, fresh successes follow Binomial(n, sin²(mθ)). A two-sided Clopper–Pearson interval has conditional noncoverage at most δᵣ = δ/[r(r+1)]. Inverting it on the branch containing the previous θ interval and intersecting with that interval preserves the true θ whenever all previous intervals and the new binomial interval cover.

For each output, the probability of a first failed interval is at most Σᵣδᵣ ≤ δ. A union bound across J outputs is at most 0.05; neither independence between different outputs nor a fixed stopping time is required. Fresh batches are essential to this simple conditional argument. Adaptive cumulative data are not treated as a new fixed-sample interval with the same confidence budget. An empty intersection is recorded as failure, and an unfinished run at the finite cap is not returned as an admitted estimate.

The algorithm stops when the amplitude interval halfwidth is at most εC/(2s), where s = maxᵢ‖Bᵢ‖² and εC is the remaining response budget after inherited projection and declared arithmetic/encoding allowances. The midpoint then has response error at most εC simultaneously on the coverage event. The receipts also retain each actual, usually smaller, interval halfwidth. The 352 fixed-amplitude trials are empirical diagnostics, not a proof of simultaneous coverage; their one uncovered interval is retained under the nominal 99% individual procedure.

## Fixed approximate circuits

Every emitted Ry rotation is formed from a fixed binary Rz bank, with its exact emitted gate-sequence inverse reused in A†. Rounding angles modulo 4π preserves conditional 2π signs. If at most R rotations occur in A, b is the angle precision, and each selected bank primitive has operator error at most η, a telescoping operator bound gives

\[
\|\widetilde A-A\|\le R(\pi/2^b+b\eta).
\]

A projector probability changes by at most twice this norm. Mapping probability back to response therefore changes C by at most 4sR(π/2ᵇ+bη). The separate rounding and synthesis allowances are each at most 10⁻⁶. Actual primitive synthesis records and all bit words are saved.

IQAE applied to the exact repeated unitary Ã and its exact inverse estimates the fixed probability ã. The encoding bias is then added once to the final response error. This does not excuse gate noise: random physical errors, inaccurate reflections, drift or using a different approximate inverse can invalidate the ideal Grover likelihood. No available device fidelity is inferred from the qubit count.

The full-output schedules in this study use classical binomial draws from the ideal model probabilities, supported by representative actual native circuit checks. They are not QPU samples or an execution of every fault-tolerant circuit. They determine realized ideal resource scenarios; the slightly different compiled probability can change an adaptive schedule. The finite algorithm caps bound possible execution, but no practical worst-case runtime guarantee is inferred from the reported scenarios.

## Fair comparisons and hardware interpretation

The matched sampling control uses the same spectral preparations, output set and remaining precision: n = ceil(log(2J/0.05)/(2εa²)) independent k = 0 shots per output. IQAE costs include (2k+1) calls to A or A† per shot, k zero/objective reflections, all preparation and uncomputation, all outputs and every readout. Native CX/depth and constructive T totals are reported separately. T estimates use a fully specified binary compiler and are upper construction counts, not a globally optimized synthesis claim.

The old factor polynomial remains a historical control, with its normalization and failed or costly routes preserved. Estimating every retained term separately would require 684 amplitude estimates for the triangle factor polynomial and 45,769,022 for the KRAS factor polynomial at the original three times, before repetitions. This is an enumeration of retained terms, not an IQAE query-complexity bound. The spectral circuit encodes the complete exponential directly and has no discarded polynomial terms.

A 1% campaign failure allowance divided by a total logical gate count is a sufficient per-location requirement under a union bound, not a measured noise model or a necessary hardware threshold. Mapping logical gates into a physical device requires a declared error-correction architecture, routing, state factories, cycle times and calibrated errors. No physical runtime, device availability or fault-tolerance overhead is assumed here.

The uniformly controlled real preparation follows the established state-preparation construction of [Möttönen et al., *Transformation of quantum states using uniformly controlled rotations*](https://arxiv.org/abs/quant-ph/0407010). The amplification identity traces to [Brassard et al., *Quantum amplitude amplification and estimation*](https://arxiv.org/abs/quant-ph/0005055). Software versions, our exact circuit choices and measured limits accompany the saved evidence rather than being inferred from these general results.
