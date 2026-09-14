# Small target counts limit inference on mean precision@5 improvement

Date: 2026-09-14. Status: statistical design audit, with proposed contract language. No experimental labels were audited, no biological analysis was run, and no measured performance result is asserted here. The possibility that only two targets have usable labels is a supplied scenario, not a verified dataset fact.

The [v0.7.0 proposal](https://github.com/orbion-life/allosteric-response-benchmark/releases/tag/v0.7.0), pages 1 and 6, was reviewed before the present revision. It specified five candidates for each of four targets, a mean precision@5 improvement of at least 0.10 over Ohm, a paired 95% interval above zero and structural comparisons. The present version separates four-target delivery from broader inference.

The statistical issue is specific. With two or four independent, nonzero target differences, a conventional exact paired sign or sign-flip test cannot attain a 5% rejection threshold, even if every difference favors Pulsar. A confidence interval obtained by inverting the corresponding exact tests cannot exclude zero at 95% confidence. This is a limitation of those procedures and assumptions. It is not a theorem that every conceivable interval method must include zero.

## Precision@5 identifies the effect and its resolution

Let each method return exactly five distinct eligible candidates for target t. Let h(P,t) and h(B,t) be the numbers with positive labels among the Pulsar and baseline lists. When every selected candidate has an interpretable binary label for the recorded endpoint,

\[
P_t=h(P,t)/5,\qquad B_t=h(B,t)/5,\qquad D_t=P_t-B_t,
\qquad \overline D_n=\frac{1}{n}\sum_{t=1}^{n}D_t.
\]

Each precision is one of 0, 0.2, 0.4, 0.6, 0.8 or 1. Each paired difference lies on the 0.2 grid from −1 to 1. The equally weighted mean difference changes in steps of 0.2/n. These are arithmetic identities for binary labels and five-item lists, not assumptions about biological performance.

| Fully labelled targets | Selected candidates per method | Mean-difference step | Net additional positive hits corresponding to mean improvement 0.10 |
|---|---:|---:|---:|
| 2 | 10 | 0.10 | 1 |
| 4 | 20 | 0.05 | 2 |

“Net additional” means the sum of Pulsar hit counts minus the sum of baseline hit counts. The 0.10 requirement is therefore attainable as a point estimate. It does not require every target to improve. For example, differences (0.2, 0, 0, 0.2) have mean 0.10, while only two differences are nonzero. This vector is illustrative, not observed.

There are two different possible estimands. A fixed four-target average describes those targets. A population mean, μ = E[D], describes a specified population and selection process for future targets. Generalization to that population requires a defensible sampling or modelling argument. Twenty selected candidates do not constitute twenty independent observations of cross-target performance. Repeated structures, timepoints, random seeds, simulated measurements or resampled controls also do not increase the independent target count.

## Exact paired inference has a small-sample probability floor

For the sign test, define n* as the number of independent nonzero differences and W as the number that favor Pulsar. Under the directional null, positive and negative signs must each have probability 1/2 conditional on a nonzero difference, independently across units. Then W has the Binomial(n*, 1/2) distribution. The ordinary greater-tail probability is

\[
p_+=2^{-n^*}\sum_{j=W}^{n^*}\binom{n^*}{j}.
\]

The sign test concerns direction or a median under appropriate additional conditions. It is not a general test of μ = 0: a distribution can have many small positive differences but a negative mean because of fewer large negative differences. NIST documents the binomial construction and exclusion of ties in the [paired sign test](https://www.itl.nist.gov/div898/software/dataplot/refman1/auxillar/signtest.htm).

For a paired sign-flip test using the sum or mean difference, there are 2^n* equally weighted sign assignments under the required null invariance. The observed assignment must be included in the tail. With all differences positive, exactly one assignment has the maximum positive sum; the all-negative assignment gives its equally extreme opposite. Under the usual inclusive tail convention, the smallest possible one-sided probability is 2^−n*, and the smallest two-sided probability is 2^(1−n*). SciPy specifies this [paired sample-swapping construction and tail convention](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.permutation_test.html).

| Independent nonzero target differences n* | Sign assignments | Minimum one-sided p | Minimum two-sided p |
|---:|---:|---:|---:|
| 2 | 4 | 0.25 | 0.50 |
| 4 | 16 | 0.0625 | 0.125 |
| 5 | 32 | 0.03125 | 0.0625 |
| 6 | 64 | 0.015625 | 0.03125 |
| 7 | 128 | 0.0078125 | 0.015625 |

These floors are exact conditional on the stated sign model or sign-flip invariance, nonrandomized testing and inclusive tails. Zero differences reduce n*. Tied statistic values must remain in the tail and can make probabilities larger. Counting more Monte Carlo draws cannot create additional independent sign assignments. A mid-p convention or randomized rejection rule changes the test and its error properties; it does not remove the information limitation.

For the sum/mean sign-flip test, the invariance assumption is stronger than merely E[D] = 0. It can follow from a valid paired randomization design or an explicitly justified symmetric difference model. Deterministic algorithm names were not randomly assigned to targets merely because the two algorithms are compared. If related targets are dependent, sign flips must preserve the independent blocks. Calling an analysis “family-stratified” does not itself establish exchangeability or independence.

Inverting an exact two-sided test with a minimum p of 0.125 cannot produce a 95% confidence set excluding zero at n* = 4. A one-sided 95% lower bound also cannot exclude zero because 0.0625 exceeds 0.05. A median interval from a sign test must not be labelled a confidence interval for mean improvement. A symmetric location-model interval requires that model to be stated.

The proposal must name its interval method. A paired t interval can numerically exclude zero at n = 4. For example, the hypothetical differences (0.2, 0.2, 0.2, 0.4) give mean 0.25, sample standard deviation 0.10 and a t interval of approximately [0.09088, 0.40912]. Exact normal-theory coverage requires its distributional assumptions; four discrete observations cannot establish them. NIST describes the [mean t interval and its small-sample normality limitation](https://www.itl.nist.gov/div898/software/dataplot/refman1/auxillar/conflimi.htm). Likewise, an ordinary percentile bootstrap can put its entire interval above zero when every observed difference is positive, even with only two targets. That numerical output does not establish reliable population coverage.

## Multiplicity depends on which claims can be declared successful

The proposal names a baseline precision contrast and two score contrasts, while structural-control comparisons may introduce additional tests. The exact family of confirmatory hypotheses, sidedness, endpoints, methods, subgroup claims and timepoint choices must be recorded before the external outcomes are inspected. The text does not establish the final number of hypotheses.

If individual positive findings may be reported from a family of m hypotheses, Holm adjustment is a defensible default. Sort valid raw probabilities and compare p(i) sequentially with 0.05/(m−i+1). The first rejection must clear 0.05/m, regardless of dependence among the tests. This follows [Holm (1979), A simple sequentially rejective multiple test procedure](https://www.ime.usp.br/~abe/lista/pdf4R8xPVzCnX.pdf).

| Hypotheses m | First Holm threshold | Smallest n* that could pass, one-sided | Smallest n* that could pass, two-sided |
|---:|---:|---:|---:|
| 1 | 0.05 | 5 | 6 |
| 2 | 0.025 | 6 | 7 |
| 3 | 0.016666… | 6 | 7 |
| 4 | 0.0125 | 7 | 8 |

The table applies the exact probability floors and strict inequality. These are feasibility minima in the most favorable possible outcome, not sample sizes that provide adequate power. A target-level endpoint can also have fewer usable units than the precision endpoint.

There is an important logical distinction. If the only permitted success claim is a single conjunction requiring every endpoint to succeed, testing each component at level 0.05 controls the error of that joint claim: whenever the joint null is true, at least one component null is true, and rejecting every component implies rejecting that true component. The joint false-positive event therefore has probability at most 0.05. This is a direct set-inclusion argument. It does not provide adjusted individual claims or simultaneous 95% intervals. The existing proposal explicitly promises multiplicity adjustment; a joint-claim formulation would require an explicit contract revision rather than a silent reinterpretation.

A comparison with many matched residues can have many valid conditional permutations within a target, if its matching and exchangeability conditions are justified. Such a result addresses that within-target comparison. It does not increase n* for the cross-target precision@5 endpoint. The random-background contrast and the supported-negative contrast also answer different questions; neither automatically establishes superiority over the baseline or structural controls.

## Unknown labels create identification limits before sampling uncertainty

A missing functional label is unknown. Absence from a positive-reference list is not, by itself, a supported negative. If a five-item list contains p known positives and u unknowns, its true binary precision is bounded by

\[
P@5\in[p/5,(p+u)/5].
\]

These are identification bounds, not 95% confidence intervals. Reporting p/5 alone as “supported-hit yield” is possible if that different endpoint is explicitly defined. Restricting a ranking to labelled candidates changes the candidate universe; dividing only by labelled selections changes the denominator. Neither silently estimates precision for the original top five.

Paired bounds can account exactly for overlap. Let c be the number of known positives selected only by Pulsar minus the number selected only by the baseline. Let uP and uB count unknown candidates selected only by each respective method. Under common, consistent labels for shared candidates,

\[
D_t\in[(c-u_B)/5,(c+u_P)/5].
\]

Shared candidates cancel, including shared unknowns. These bounds are sharp given otherwise unconstrained binary unknown labels. Average the lower and upper endpoints to bound a fixed target-set mean. Sampling uncertainty, label uncertainty and uncertainty about transfer to other targets are separate quantities and should be reported separately.

If two of the required four targets are entirely unassessable, their differences can each range from −1 to 1. Given the mean difference d of the two assessed targets, the fixed four-target mean is bounded by

\[
\overline D_4\in[0.5d-0.5,\;0.5d+0.5].
\]

For the illustrative d = 0.10, the interval is [−0.45, 0.55]. Even d = 1 gives a lower bound of zero. Thus two wholly unknown target outcomes cannot establish positive mean improvement across the fixed four by worst-case identification alone. Reporting the two-target complete-case mean is legitimate with its denominator and scope stated; extending it to four targets requires additional assumptions or labels.

## A development and external-validation contract preserves the four-target delivery

The following contract is proposed statistical language. It does not claim that external data, resources or sample sizes are available.

1. **Development delivery.** Retain all four required targets, their matrices, five-item lists and maps. Report baseline and structural-control rankings, candidate overlap, label coverage, per-target effects or identification bounds, failures and unassessable endpoints. A missing label does not excuse a missing prediction deliverable. Describe analysis on targets used to choose the method as development or retrospective evaluation.

2. **External cohort and population.** Before inspecting external outcomes, record the population of intended use, target eligibility, family relationship rules, development exposure and target weights. Use independent targets or independent family blocks that were not used to select the method. One block contributes one independent effect for inference when dependence requires that aggregation. A held-out label subset from a development target does not by itself establish transfer to new targets.

3. **Endpoint definition and coverage.** Record one functional endpoint definition, positive and supported-negative rules, candidate eligibility, exactly five selections per method, tie handling and baseline versions before outcome comparison. Apply the same labels to overlapping candidates. Record how unknown outcomes will be assessed or bounded; do not select assessable targets after viewing comparative performance.

4. **Primary success criterion.** For the prespecified mean effect, require both an observed improvement of at least 0.10 and a lower 95% confidence limit above zero, using a stated mean-valid procedure and the chosen multiplicity policy. The first condition is an observed effect threshold. It is not evidence that the population improvement is at least 0.10; that stronger claim requires a lower confidence limit above 0.10. Record the structural-control and supported-negative decisions separately and state which are required for the overall success claim.

5. **Sample-size feasibility.** Determine external sample size from the independent unit, expected effect variation, ties, coverage, dependence, multiplicity and desired power or precision. Do not use five or six independent targets as a powered-study recommendation merely because an exact p-value first becomes attainable there. Do not increase sample size or replace targets in response to disappointing interim effects unless the sequential rule was specified in advance.

6. **Decision status.** If the external cohort, labels or required controls are unavailable, deliver the development package and record external validation as unassessable or pending. If the external assessment is completed and a valid success criterion is not met, report the estimate and uncertainty without claiming superiority. Do not pool development targets into the external test to obtain a smaller probability after the result is seen.

A conservative, fully specified mean-interval option is available when independent target differences lie in [−1,1]. Applying Hoeffding’s inequality to X = (D+1)/2 gives the 95% confidence interval

\[
\left[\max(-1,\overline D-\sqrt{2\log(40)/n}),\;
\min(1,\overline D+\sqrt{2\log(40)/n})\right].
\]

Coverage follows for the average expected difference; a common target population adds the sampling interpretation μ = E[D]. The half-width is 1.92065 at n = 2 and 1.35810 at n = 4. This particular conservative bound first has half-width below 0.10 at n = 738. That is a property of this bound, not an inherent minimum sample size for the task. For m simultaneous mean intervals, replace log(40) with log(40m). More efficient alternatives require their own justified assumptions and validation. Source: [Hoeffding (1963), theorem 1, equation 2.3 and range rescaling, printed page 15](https://www.cs.rpi.edu/academics/courses/spring06/random/hoefding.pdf).

## Numerical verification and permissible claims

Exact integer enumeration checked all 121 possible two-target and all 14,641 possible four-target difference vectors on the precision@5 lattice. Each vector was evaluated over every sign assignment, retaining equal tail values. The observed minima were 1/4 and 1/2 for one- and two-sided tests at n = 2, and 1/16 and 1/8 at n = 4. This verifies the numerical statements, not the applicability of the sign-flip assumptions to a future dataset.

The minimal calculation can be repeated without experimental data:

```python
from itertools import product
from fractions import Fraction

for n in (2, 4):
    masks = list(product((-1, 1), repeat=n))
    best_one = best_two = Fraction(1)
    for d in product(range(-5, 6), repeat=n):
        observed = sum(d)  # Units of one net hit; dividing by 5*n is immaterial.
        null = [sum(s*x for s, x in zip(mask, d)) for mask in masks]
        p_one = Fraction(sum(x >= observed for x in null), len(null))
        p_two = Fraction(sum(abs(x) >= abs(observed) for x in null), len(null))
        best_one, best_two = min(best_one, p_one), min(best_two, p_two)
    print(n, best_one, best_two)
```

The proposal can promise a reproducible four-target development deliverable and a prespecified external validation protocol. It can state the 0.10 threshold as a desired outcome. It cannot promise that a two- or four-target exact paired test will satisfy the stated 95% inferential success criterion. That prior text needed an explicit estimand, interval method, independent unit, multiplicity family and missing-label policy before that criterion becomes statistically reviewable.
