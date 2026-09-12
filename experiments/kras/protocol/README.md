# Four-target allosteric-response validation protocol

This package contains reproducible input preparation, exact residue mappings, fixed evaluation rules and the executed KRAS Ohm comparator. It supports a development protocol, not completed four-protein validation or a prospectively blinded study. [validation-protocol.md](validation-protocol.md) explains every biological scope decision; [validation-contract.json](validation-contract.json) records the machine-readable settings. [mechanical-and-quantum-methods.md](mechanical-and-quantum-methods.md) defines the numerical response and error accounting.

## Reproduce the structural preparation and finite null spaces

The recorded environment was Python 3.9.6, NumPy 2.0.2 and Biopython 1.85. Work in a fresh copy if preserving the supplied receipts byte for byte.

```bash
python3 -m pip install -r requirements.txt
python3 build_validation_inputs.py
python3 extract_reference_labels.py
python3 matched_null.py --target KRAS --out evaluation/KRAS-null-plan.json
python3 matched_null.py --target ABL --out evaluation/ABL-null-plan.json
```

Input preparation reads only the selected input structures and sequence records. It writes an input-freeze receipt before reference-label extraction. This execution order makes future edits traceable; it does not erase the retrospective development history. [preparation-replay-verification.json](preparation-replay-verification.json) records a fresh-directory replay in which all 19 regenerated mapping, label and finite-null files were bitwise identical. This tests relocation within the stated host environment; it does not claim cross-platform bitwise identity. The four outputs contain 166 KRAS, 252 ABL, 694 MYH7 and 162 MYC/MAX nodes. Functional receivers contain 18, 19, 8 and 16 target residues respectively. Each record preserves author identifiers and current sequence mappings. Source URLs, retrieval date and SHA-256 values are in [sources/download-manifest.json](sources/download-manifest.json).

The null command without `--scores` calculates the finite sample space; it does not generate a prediction. With `--scores scores.csv --score-column S`, the file must provide a unique `canonical` position and finite numeric score for every eligible candidate. The routine returns a raw conditional-null p-value. Apply `holm([p_KRAS,p_ABL,1.0])` from `matched_null.py` to the full three-slot primary family while MYH7 remains unresolved. MYC has no designated ligand-contact endpoint. An absent negative-pocket annotation is not a true-negative pocket.

[ohm/README.md](ohm/README.md) documents the pinned public comparator, patch, score differences, parameters and independent reconstruction. Its portable build source is [../pilot/replay_ohm.py](../pilot/replay_ohm.py), documented in [../pilot/README.md](../pilot/README.md). The sibling pilot also contains the thermal-response calculations and circuit fixtures. This folder does not silently substitute Ohm scores for the primary full-law response test.

Input preparation, two reference-label extractions, finite-null planning and KRAS Ohm scoring have been executed. Remaining protein predictions, full-dimensional accuracy gates, additional comparators and eligible negative-pocket assessment remain proposed work. The response matrices and shortlist visuals for all four targets remain required outputs even when an evaluation endpoint is unresolved.

[DATA-AND-SOFTWARE-ATTRIBUTION.md](DATA-AND-SOFTWARE-ATTRIBUTION.md) identifies public data sources, original structure citations and third-party terms. [artifact-manifest.json](artifact-manifest.json) binds the included files. The package contains scientific material only; team profiles, administrative records and private proposal documents are excluded.
