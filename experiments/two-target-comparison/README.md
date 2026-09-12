# Project Pulsar: a joint retrospective statistical audit

The corrected KRAS result retains conditional contact-score enrichment after the specified three-slot Holm adjustment. The locked ABL result does not transfer. **No paired comparison with harmonic or distance-Hookean controls passes the six-slot Holm correction at 0.05.** The small positive ABL rank effects are retained below; they are not presented as demonstrated advantage.

![Paired input structures and reference labels](figures/protein-comparison.png)

The figure shows actual frozen input Cα traces and the five ordered reference labels for each method. A filled marker denotes a resolved heavy-atom contact within 5 Å of the reference ligand; an open marker denotes an observed residue farther away; a question mark denotes an unresolved residue. The reference ligands are 6OIM/MOV for KRAS and 5MO4/AY7 for ABL. The first four rows use two-mode, 33²-grid quantities. All-mode harmonic uses 492 modes for KRAS and 750 for ABL. The traces have the same Å-to-point display scale. Contact recovery is not functional validation.

## Families and source selection

`audit.py` reads hash-bound copies of the corrected KRAS and ABL evaluation outputs. It selects the KRAS `degree_distance_exposure / evaluable123` records, not the historical 126-member null. ABL uses its corrected 175-member evaluable universe. Both analyses keep the original prediction ranks fixed: 126 candidates for KRAS and 193 for ABL. Unknown reference residues are excluded from evaluation and random-label pools.

The primary family has three slots: KRAS, ABL and MYH7. Each measured primary test compares the mean biquadratic score at known reference contacts with a degree-, receiver-distance- and exposure-matched label null. The paired family has six slots: each of those targets against its same-grid harmonic and distance-Hookean controls. Paired effects are the mean primary-minus-control percentile rank among eligible known contacts.

MYH7 has not been run. Its one primary and two paired slots use p = 1 only for multiplicity correction. They have no measured effect and are explicitly labelled as unrun. Their inclusion does not imply that an MYH7 experiment was completed. These are two separate correction families, as specified for this retrospective audit.

## Primary enrichment

| Target | Status | Raw conditional p | Holm-adjusted p, three slots |
|---|---|---:|---:|
| KRAS | Retrospective measurement | 0.0000999900 | 0.0002999700 |
| ABL | Retrospective measurement | 1.000000 | 1.000000 |
| MYH7 | Unrun correction slot | 1 by convention | 1 by convention |

KRAS enrichment concerns the score under the stated conditional-label null. It does not show that nonlinear dynamics improve ranking over the matched physical controls. ABL's primary shortlist contains zero known contacts among four reference-observed members and one unresolved member, so its precision is bounded by 0–1/5. Its primary p value is 1.0. The ABL data and complete replay are available in [the companion pilot](../abl/pilot/README.md).

## Paired comparisons

Positive effects mean that known contacts receive higher percentile ranks under the biquadratic model. One percentage point equals 0.01 on the reported 0–1 percentile scale.

| Target | Control | Mean percentile difference | Difference in percentage points | Raw p | Holm-adjusted p, six slots |
|---|---|---:|---:|---:|---:|
| KRAS | Harmonic | −0.00233427 | −0.233427 | 0.850815 | 1.000000 |
| KRAS | Distance-Hookean | 0.00000000 | 0.000000 | 0.347065 | 1.000000 |
| ABL | Harmonic | 0.01709845 | 1.709845 | 0.012599 | 0.075592 |
| ABL | Distance-Hookean | 0.00025907 | 0.025907 | 0.539046 | 1.000000 |
| MYH7 | Harmonic | Unmeasured | Unmeasured | 1 by convention | 1 by convention |
| MYH7 | Distance-Hookean | Unmeasured | Unmeasured | 1 by convention | 1 by convention |

The smallest adjusted paired p value is 0.075592 for ABL versus harmonic dynamics. It exceeds 0.05. The ABL difference against distance-Hookean dynamics is approximately 0.026 percentage points, and its adjusted p value is 1.0. These values do not support a claim of demonstrated nonlinear advantage. The equilibrium comparison remains in the ABL source output, but it is not silently added to or substituted for either specified control slot.

## Reproduce and verify

The statistical audit uses only the Python standard library. The bundled source JSON files make it independent of local report paths or an internet connection.

```bash
python3 audit.py
python3 verify.py
```

The audit verifies every input hash and writes `results/audit.json`, `results/primary.csv` and `results/paired.csv`. The verifier checks Holm boundary and tie cases, rejects nonfinite p values, confirms the family sizes and unrun slots, and repeats the audit in a fresh process and output directory. All three outputs reproduce byte for byte. The [verification receipt](verification-receipt.json) records code and result hashes.

The optional figure can be rebuilt with NumPy and Matplotlib:

```bash
python3 figures/protein-comparison.py
```

Its portable coordinate/label data, native-size PDF, SVG, PNG and provenance receipt are in `figures/`. Its native size is 493 × 216 PDF points, with 10.4-point text. Source-file hashes accompany the exported data.

## Limits of this audit

The references were known retrospectively. Conditional random-label tests depend on exchangeability after matching; these strata do not reproduce the spatial contiguity of a physical pocket. Holm arithmetic does not repair that modelling assumption or make the result prospective. The stated adjustment does not cover every exploratory analysis in the larger project.

Both protein pilots use deliberately reduced mechanics, and their two-mode harmonic comparisons fail the locked full-mode representation criterion. Stable grids and repeatable computations do not establish physical adequacy. This audit therefore preserves the distinction between statistical enrichment, improvement over controls, approximation accuracy and functional biological validation. It establishes no finite-time, biological or quantum advantage.

Original audit code and text are MIT-licensed. The source evaluations derive from public scientific records whose attribution remains in the companion pilot and corrected KRAS package. No original v0.2.0 artifact is replaced by this new folder.

The adjustment follows Holm (1979), [A Simple Sequentially Rejective Multiple Test Procedure](https://www.jstor.org/stable/4615733), *Scandinavian Journal of Statistics* 6, 65–70. The implementation sorts raw p values, applies the remaining family size at each step, takes the running maximum and caps adjusted values at 1.
