# KRAS evaluation amendment: missing reference residues remain unknown

This amendment corrects reference-label handling in the retrospective KRAS pilot. It preserves the released model, response matrices, candidate mask, scores and prediction ranks. KRAS residues 105–107 are present in the input model but have no deposited protein heavy atoms in reference structure 6OIM. They therefore remain prediction candidates with **unknown reference labels**. The amended evaluation uses 123 observable candidates; the prediction still ranks all 126 candidates.

The amendment is based on public release [v0.2.0, commit 527712c](https://github.com/orbion-life/allosteric-response-benchmark/tree/527712c281f1a366b9fb1524342a8daf31a06d7a). It adds an evaluation record and does not replace historical results. The result remains a retrospective engineering comparison, with no demonstrated nonlinear advantage or validated protein compression.

## Reproduce the amendment

Place this directory at `experiments/kras/evaluation-amendment/` beside the unchanged `pilot/` and `protocol/` folders. The default command reads those siblings and writes only its own `results/` folder. Python 3.9.6, NumPy 2.0.2 and SciPy 1.13.1 produced the supplied results.

```sh
python3 -m pip install -r experiments/kras/evaluation-amendment/requirements.txt
python3 experiments/kras/evaluation-amendment/amend_evaluation.py
python3 -m unittest discover -s experiments/kras/evaluation-amendment -p 'test_*.py' -v
```

For an amendment directory stored elsewhere, supply `--kras-root /path/to/experiments/kras` and optionally `--out /path/to/new-results`. The script verifies every required input against `source-manifest.json` and verifies all 59 files in the original prediction freeze. It refuses an output directory inside the frozen `pilot/` or `protocol/` folders. A changed source file causes an error before evaluation.

`verification.json` records software versions, source and output checksums, historical draw parity, and the correction's limits. `matched-configurations.npz` preserves every sampled residue set. Unit tests cover explicit missing labels, immutable prediction eligibility, unknown top-five predictions, exclusion from corrected null pools, exact finite-space cases, fixed-cut sensitivity and historical reproduction. The package contains no network call or protein-model fit.

## Define prediction and evaluation separately

The original evaluator computes `distance <= 5`. For a missing distance, NumPy returns `False`, which silently turns an unavailable contact label into an observed noncontact. The amended code uses three states:

| State | Meaning in this deposited structure | Contact field |
|---|---|---|
| `contact` | A deposited protein heavy atom lies within 5 Å of a deposited MOV heavy atom. | `true` |
| `observed_noncontact` | At least one protein heavy atom is deposited, but none of those atoms lies within 5 Å of MOV. | `false` |
| `unknown` | No protein heavy atom is deposited at that mapped input position. | `null` |

An observed noncontact is **not** a verified non-functional pocket or a biological negative. A residue with some deposited atoms can also have an unresolved side chain, which can hide a contact. This amendment resolves the completely absent-residue error; it does not impute missing atoms or claim complete chemical observability for every other residue.

The evaluation rule is `evaluation_eligible = prediction_candidate AND reference_observed`. The input candidate field remains unchanged. All 126 positions retain their original prediction ranks. An unknown residue appearing in the top five would retain its place: precision at five would be reported as unavailable with lower and upper bounds. No lower-ranked residue would be substituted to improve apparent performance.

In this KRAS dataset, every top-five residue is observable for each of the ten scored methods. The primary prediction remains **60, 69, 62, 61, 65**, with four reference contacts and precision at five of **0.8**. The same-grid harmonic and distance-Hookean controls also retain four contacts. Equal hit counts do not establish equivalent predictors, but they do not demonstrate a nonlinear benefit.

## Preserve both historical matching schemes

The score statistic is the mean of a method's fixed residue scores across the 17 eligible MOV-contact residues. The null compares that mean with means from residue sets having the same positive-label counts within each matching stratum. Reference labels do not enter the score calculation.

1. The original pilot matches joint tertiles of contact degree and receiver distance. It makes 10,000 draws, sampling without replacement within each draw; configurations can repeat between draws.
2. The later protocol also matches the exposure class RSA ≥ 0.20. It enumerates all configurations when the space contains at most 100,000 configurations. Otherwise, it samples 10,000 distinct configurations, excluding the observed positive set. The supplied KRAS cases all use sampling. The existing RSA values are reused; the exposure cutoff is a declared protocol choice, not a fitted biological boundary. The normalization of RSA follows Tien et al. (2013).

Both schemes use PCG64 with seed 20260912. For sampled tests, the reported raw probability is `(1 + number of null statistics ≥ observed statistic) / 10001`. Its smallest reported value is 1/10001; observing that floor does not determine the exact tail probability. For exact enumeration, the tail count is divided by the full number of configurations. The nominal finest exact resolution is 1/M, and ties can make the achievable tail probability larger. Conditional random sets express a geometric exchangeability assumption. They are not biological replicates or independently evidenced non-functional pockets.

Each scheme is run three ways. `historical126` reproduces the historical universe and cutpoints. `evaluable123` removes the three unknowns from evaluation and recomputes tertiles by the inherited protocol rule. `evaluable123_fixed_historical_cuts` removes the unknowns while retaining the historical cutpoints, isolating exclusion from rebinning. The original pilot and the later protocol computed receiver distances with slightly different parser precision; the amendment retains each scheme's original numerical source. Their maximum input-distance difference is 2.79 × 10⁻⁶ Å, with identical prediction eligibility.

| Conditional space | Historical 126 | Evaluable 123 with new tertiles | Evaluable 123 with historical cuts |
|---|---:|---:|---:|
| Degree and receiver distance | 2,872,175,981,875,200 | 3,764,901,254,400,000 | 2,872,175,981,875,200 |
| Degree, receiver distance and exposure | 26,676,557,107,200 | 11,997,455,097,600 | 26,676,557,107,200 |

The three missing residues occupy a historical stratum containing **zero positive contacts**. Consequently, none was selected in the historical sampled configurations. Excluding them while keeping those cuts produces identical configurations and statistics. Differences in the amended re-tertiled null arise from moving observable residues between strata. The correction does not show that missing residues directly diluted the historical sampled null scores.

## Report corrected enrichment and paired controls

The following values use the later degree–distance–exposure scheme. Complete results for both schemes and all ten methods are in `results/null-comparisons.csv`.

| Fixed score | Historical raw probability | Evaluable 123 raw probability |
|---|---:|---:|
| Biquadratic, two coordinates | 0.00019998 | 0.00009999 |
| Harmonic, same two-coordinate grid | 0.00009999 | 0.00009999 |
| Distance-Hookean, same two-coordinate grid | 0.00019998 | 0.00009999 |
| Biquadratic equilibrium, two coordinates | 0.00019998 | 0.00009999 |
| Harmonic, all 492 internal coordinates | 0.33306669 | 0.40405959 |

Enrichment against a matched set distribution does not compare the biquadratic method directly with a control. The paired statistic therefore first converts each method's scores to average-tie percentile ranks over the **original 126 prediction candidates**, with higher scores receiving higher ranks. For each residue, it subtracts the control percentile rank from the biquadratic percentile rank. The mean difference across the 17 contact residues is compared with the same matched sets. Keeping the 126-residue rank denominator prevents reference missingness from changing the predictor being compared.

| Paired comparison | Mean contact percentile difference | Historical raw probability | Evaluable 123 raw probability |
|---|---:|---:|---:|
| Biquadratic minus same-grid harmonic | −0.00233427 | 0.91570843 | 0.85081492 |
| Biquadratic minus same-grid distance-Hookean | 0 | 0.28337166 | 0.34706529 |

These one-sided comparisons provide no evidence of a positive rank advantage. The raw probabilities are retrospective sensitivities. They do not complete the reserved six-test comparative family across three targets, and they do not establish equivalence or absence of every possible nonlinear effect. Failed coordinate-compression and noise guarantees remain separate findings in the original pilot.

As a descriptive baseline, five distinct uniformly chosen residues would have at least one contact with probability **52.1632%** in the historical 126-residue universe and **53.0984%** in the 123-residue evaluable universe. The corrected expected number of hits is 0.6911; the probability of at least four hits is 0.0011962. The observed precision of 0.8 is 5.7882 times the corrected contact fraction, 17/123. These calculations do not match degree, receiver distance or exposure and must not be presented as the primary matched test.

## Inspect the records

| File | Purpose |
|---|---|
| `source-manifest.json` | Pins all directly used public inputs and code by SHA-256. |
| `historical/` | Preserves the original evaluator, matched-null code, evaluation and prediction-freeze record. |
| `results/reference-labels.json` | Gives explicit contact, observed-noncontact and unknown states for all 166 modeled residues. |
| `results/unchanged-rankings.csv` | Preserves all ten complete candidate rankings and adds evaluation eligibility. |
| `results/null-plans.json` | Gives cutpoints, every stratum, positive counts and exact configuration-space counts. |
| `results/matched-configurations.npz` | Stores the six sampled configuration arrays. |
| `results/null-comparisons.csv` | Reports all methods under both matching schemes and three universe/cut policies. |
| `results/paired-comparisons.csv` | Reports the paired control statistics, with the fixed rank denominator. |
| `results/uniform-hit-baselines.json` | Separately records descriptive uniform-sampling hit probabilities. |
| `results/verification.json` | Records frozen-input verification, historical parity and all output checksums. |

The deposited [6OIM structure](https://www.rcsb.org/structure/6OIM) and its cached PDB `REMARK 465` records independently identify the three absent residues. The [PDBe SIFTS mapping](https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/6oim) maps the evaluated chain to UniProt P01116. The original reference study describes the inhibitor-bound structure (Canon et al., 2019; [doi:10.1038/s41586-019-1694-1](https://doi.org/10.1038/s41586-019-1694-1)). Exposure normalization uses theoretical maximum accessible areas (Tien et al., 2013; [doi:10.1371/journal.pone.0080635](https://doi.org/10.1371/journal.pone.0080635)). The unchanged public protocol records data attribution, input mutations, surface calculations and their limits in [DATA-AND-SOFTWARE-ATTRIBUTION.md](../protocol/DATA-AND-SOFTWARE-ATTRIBUTION.md).

The amendment code is provided under the repository's MIT license. Existing data and third-party software retain their original terms. No new protein coordinates, negative-pocket evidence or experimental measurements are introduced.
