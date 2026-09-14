# Reproduce the GRB2 comparison

Use Python 3.13 and the pinned packages in `requirements.txt`. The measured environment was Python 3.13.5 on macOS, with one BLAS thread. Scripts resolve their inputs relative to this folder. Original receipts can contain the original working directory as provenance; no script needs that directory.

## Saved-result audit

This route needs no assay download, compiler, Ohm binary, fit or sampling run:

```sh
python -m pip install -r requirements.txt
python scripts/verify_saved.py
python scripts/verify_inputs.py
python scripts/verify_math.py
```

The first audit checks the protocol and prediction hashes, contact graph, mask, covariance restoration, upstream Ohm seed mean, direct SciPy correlations, 500 deterministic bootstrap rechecks, all reported difference intervals, top-five sets and the public file manifest. It uses saved own site aggregates. The mathematical audit compares independent polynomial/Wick calculations, objective derivatives and the two harmonic moment implementations on the actual GRB2 model. It does not compare the surrogate with original nonlinear dynamics.

## Fresh score replay

Make a separate copy of the packet before rerunning, so executed originals and failures remain intact. Install GCC 15, make, patch and zlib development files through your normal environment process. Then:

```sh
python scripts/build_ohm.py --compiler g++-15
python scripts/predict.py
```

The build uses the included exact upstream archive and recorded two-line C++23 portability patch. It installs only into this folder. The comparator's numerical formulas are unchanged. The prediction script reads no assay outcomes. It regenerates primary Gaussian, harmonic, degree, distance and three-seed Ohm scores, followed by all 20 NMR sensitivity models. Every fit and raw response matrix is saved.

Fresh execution changes timing/path receipts and can change floating-point bytes, so the original `PUBLIC-MANIFEST.json` remains an archival integrity check, not a promise that freshly generated files retain every byte. Compare response arrays and statistics with the archived copy. Different standard libraries can also alter Ohm random streams; record compiler/platform and the three seed outcomes rather than assuming bitwise portability.

## Raw source replay

Independently acquire Faure Supplementary Table 7 from the primary link in `SOURCES-AND-RIGHTS.md`, subject to the applicable source terms. The exact SHA-256 is pinned. The workbook is deliberately absent from this packet.

```sh
python scripts/ingest_faure.py --source-xlsx /your/download/41586_2022_4586_MOESM10_ESM.xlsx
python scripts/evaluate.py
python scripts/verify_source.py
python scripts/plot_results.py
```

`ingest_faure.py` validates the workbook hash and creates `local-only/`. Do not add that folder to a public release. `verify_source.py` reads the workbook independently of the production CSV and uses decimal accumulation to check all site aggregates. The full evaluation fixes all decisions in `PROTOCOL.md`; changing them creates a new exploratory analysis and must not replace this result.
