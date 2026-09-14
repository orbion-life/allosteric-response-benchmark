# Input provenance and reuse terms

Original analysis code and documentation use the repository MIT license. This does not relicense upstream experimental data.

The input CSVs contain computed per-site aggregates and fixed model predictions. `input-manifest.json` gives their SHA256 and original study names. No publisher workbook, mutation-level assay rows or literature PDF is included.

- KRAS aggregates come from the earlier published KRAS–RAF1 study using MaveDB 00000115-a-6 and 00000115-a-7. The retained source metadata in that release state CC0 1.0. Their model and aggregation provenance remain in the v0.6.0 evidence.
- GRB2 aggregates derive from Faure et al., *Nature* (2022), Supplementary Table 7. The source workbook SHA256 is `58abcbd1a08ad9fed8778e3b3e8672bbca84ac5a576f1b0530ddd9f5164fbce2`. Unrestricted redistribution rights for that workbook were not established, so obtain it directly under the publisher's applicable terms if reconstructing the original aggregates. This release starts from the computed site aggregates, following the same boundary as the previous public packet.
- Contact geometry, receiver masks, Gaussian and harmonic scores, and Ohm comparison values are inherited unchanged. This reanalysis distributes no Ohm executable or new third-party implementation.

[Faure et al. primary article](https://www.nature.com/articles/s41586-022-04586-4) · [prior GRB2 source and rights record](https://github.com/orbion-life/allosteric-response-benchmark/blob/v0.8.0/experiments/selection-2026-09-14/biology/SOURCES-AND-RIGHTS.md) · [earlier release license record](https://github.com/orbion-life/allosteric-response-benchmark/blob/v0.8.0/LICENSES.md).
