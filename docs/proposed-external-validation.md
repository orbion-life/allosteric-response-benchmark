# Proposed external validation of Project Pulsar

Prepared 12 September 2026. This is an analysis specification under development, not a registered or executed benchmark. No external functional labels have been scored in this revision.

## Decision and endpoint

The intended decision is which five eligible residues should receive experimental follow-up. The official ligand-contact endpoint and a separate functional endpoint answer different questions. Precision@5 (P@5) is the number of predeclared positive labels among all five predictions divided by five. Unknown reference positions are retained through lower and upper bounds. Missing predictions or failed calculations do not disappear from the cohort; report coverage and a prespecified conservative failure policy.

A mean paired improvement of 0.10 is a proposed practical target, corresponding to half an additional hit per five-residue list across the panel. It has not been empirically established as an economically or biologically sufficient effect. Before registration, the accountable biological reviewer must assess the target against assay cost, baseline yield and the practical value of another correct residue. The statistical criterion is a paired 95% confidence interval excluding zero plus a point estimate of at least 0.10. This supports a positive effect if assumptions hold; demonstrating that the true improvement exceeds 0.10 requires its lower confidence limit to exceed0.10.

## Independent cohort and comparator

KRAS and ABL are development families and must remain outside external confirmation. MYH7 and MYC retain their official target roles; another structure or conformer from one family does not create another independent family. Select the external families from a documented source universe using criteria fixed before accessing their evaluation labels. Record relatedness, receiver availability, construct, structure quality, ligand state, unresolved atoms, assay provenance and comparator training overlap. An assessor can sequester labels while the prediction team works from permitted inputs. Public historical structures remain a retrospective benchmark even when operationally held out.

Use equal family weights. If a family contributes several structures or states, first aggregate its paired method difference using the fixed rule; then average across families. Report the full family-level paired differences, coverage and effect interval. Do not treat residues, conformers, technical repeats, noise draws or random null sets as independent biological replicates.

One primary comparator must be nominated before label access using independent development evidence and feasibility. AlloPred, Ohm, harmonic, distance-Hookean and equilibrium methods form the supporting comparison set. Several confirmatory comparator claims require a declared multiplicity procedure; selecting whichever comparator is weakest after observing the external result is invalid. Check and document each baseline's actual software, required inputs, calibration, failed cases and reproducibility.

The proposed 12–20-family range from the strategy is a resource-planning range only. Before locking the cohort, quantify attainable precision from an independent pilot or explicit variance assumptions, and report sensitivity to those assumptions. Use simulation of family-level paired outcomes to assess interval coverage, false-positive behaviour and power at 0, 0.05, 0.10 and 0.20 gains. Such simulations assess the statistical design; they do not supply biological observations. Fix the final cohort, grouping, primary test, confidence method, exclusions and failure policy before unblinding. If precision is inadequate within the sprint, label the result an estimation pilot and report its uncertainty without a generalization claim.

## Contact labels, functional labels and controls

Keep the public KRAS/ABL predictions and existing hypothesis families unchanged. Ligand contacts use observed reference atoms and the declared 5 Å cutoff; unobserved positions remain unknown. The historical enrichment family contains KRAS, ABL and MYH7, and the historical component-value family contains six target-by-control slots. Their unrun slots are not reallocated after results. The external family comparison is a new hypothesis family.

Weng et al. (2024) report mutation-induced RAF1-binding effects in four KRAS surface pockets. These pockets cannot be assumed to be inactive alternatives to the sotorasib pocket. The study supplies a possible functional development analysis, not an independent protein-family validation panel. Before any functional score comparison, fix dataset/version, sequence mapping, nucleotide state, binding partner, direct-interface exclusion, mutation aggregation, abundance/folding filters, assay uncertainty and coverage. The existing GDP receiver and a proposed RAF1-interface receiver are separate analyses. Local contact stiffening is not an amino-acid substitution; mechanical-response sign does not establish a mutation's or ligand's effect sign.

Use the label “experimentally non-functional” only for the named assay and perturbation, with a justified negligible-effect range and enough precision to support it. Nonsignificance alone does not show inactivity. Unknown pockets are structural decoys. A non-contact to one ligand is not a functional negative. Random-residue controls, matched structural decoys and experimentally inactive pockets must be reported separately.

For a newly specified pocket-enrichment sensitivity, define connected patches and matching on size, exposure, degree and receiver distance before evaluation. Report matching failures and the number of admissible sets. Exact enumeration over M equally weighted sets, including the observed set, has minimum nonzero probability 1/M. A plus-one Monte Carlo estimator using B draws has floor 1/(B+1). More draws refine numerical inference; they do not add proteins. Use consistent assignments of unknown labels across paired methods.

## Four-target completion and unresolved facts

KRAS 4OBE→6OIM and ABL 1OPL→5MO4 have public development outputs. MYH7 5TBY and MYC 1NKP matrices and shortlists remain proposed. The prescribed 6C1H reference is myosin-Ib, so its intended cardiac validation role requires authoritative clarification. MYC's named assembly, receiver and structural sensitivities should be fixed before prediction; consensus or docking cannot substitute for functional labels. No organizer contact, profile commitment, eligibility determination or submission has been performed here.

## Sources and interpretation

- Weng, C., Faure, A.J., Escobedo, A. and Lehner, B. (2024). The energetic and allosteric landscape for KRAS inhibition. Nature 626, 643–652. https://doi.org/10.1038/s41586-023-06954-0. Primary paper directly checked; this source supports the KRAS assay interpretation, not this proposed cohort size or success threshold.
- Official enterprise statement and Phase 1 submission guidelines, cached primary PDFs checked 12 September 2026. The proposed 0.10 threshold is absent from the official rules. Current programme: https://quantumai.thequantuminsider.com/program/.
- Statistical rules above are the proposed analysis logic, not claimed experimental findings. Exact null-probability floors follow directly from the stated counting or plus-one estimators.
