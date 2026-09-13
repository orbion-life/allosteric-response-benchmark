# Biological comparison and delivery planning

Read [the specified comparison](biological-validation-design.md) before using these calculations. The source inventory is empty, the general comparator adapter is unimplemented and contributor allocations remain unconfirmed. The files define an auditable route to evaluation; they are not completed biological evidence.

The precision scenarios use assumed normal paired differences. Their interval widths and success probabilities are planning calculations, not empirical uncertainty or power guarantees for a selected protein cohort. The independent unit is a protein family. Unknown reference labels and failed method outputs must retain their uncertainty.

## Reproduce the planning calculations

Copy this folder to a fresh location first, because the precision script writes its CSV and JSON beside its source. Use Python 3.12 with the pinned NumPy 2.1.3 and SciPy 1.14.1 environment specified by `../../experiments/physical-reference-recovery/requirements.txt`. Then run:

```sh
python -W error precision_scenarios.py
python -W error paired_label_bounds.py
python -W error audit_cohort.py cohort-template.json
```

The first command recomputes 96 assumption scenarios and six internal checks. The second verifies shared-label bounds by exhaustive enumeration, including failed-method intervals and fixed-five validity. The last command intentionally reports that the empty inventory cannot launch. Populate and freeze an eligible independent-family inventory, label provenance and the declared comparator before changing that status.

`commitments.csv` records proposed responsibilities and missing author input. It must not be interpreted as consent, current affiliation or committed effort. Published research supports subject expertise, not availability.
