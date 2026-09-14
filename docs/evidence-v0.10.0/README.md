# ABL contact factors and nonlinear response fidelity

Version 0.10.0 adds two completed experiments to the Gaussian protein response benchmark. One constructs the fitted ABL drift from local contact factors and carries the corresponding polynomial construction through a complete small quantum circuit. The other tests the unchanged Gaussian and harmonic approximations against the original nonlinear reference across ten fixed times.

| Experiment | Result | Meaning |
|---|---|---|
| Full ABL contact factor | Relative force-Gram error 2.804 × 10⁻¹⁵; construction takes 0.040 s given the fitted covariance. | Local contacts reproduce the fitted Gaussian drift. |
| Complete small ABL circuit | Seventeen qubits reproduce one Gaussian response within 6.594 × 10⁻⁹ in ideal simulation, with 60,250 CX gates. | The physical-factor construction is implemented on three ABL residues with two retained coordinates and fourteen Hermite functions. |
| Full protein access | An explicit Hermite-to-reduced table would occupy 134.8 TB; only 64 selected rows were constructed. | The implemented table route fails the resource gate. Compact coherent loading remains unresolved; the table size is not a universal lower bound. |
| Matched ABL classical calculation | Cached full matrices take 1.29–1.79 ms after a 0.237 s eigensystem setup. | These query timings exclude prior model construction and provide a demanding comparator. |
| Nonlinear response amplitudes | Neither fixed approximation meets 0.002 for any of its 90 entries across ten times. Maximum errors are 0.4093 for Gaussian and 0.7464 for harmonic. | Accuracy within the Gaussian model does not establish fidelity to the tested nonlinear model. |
| Nonlinear site ordering | Each approximation agrees in 16 comparisons, reverses 12 and leaves two unresolved under the stated empirical allowance. | These are repeated comparisons from one three-site fixture, not independent protein predictions. |

The later-time small ABL circuit was compiled but exceeded the frozen operation cap and was not simulated. No new QPU execution, biological experiment or practical quantum speedup is reported.

## Read and reproduce the experiments

- [ABL physical factors, circuits and complete cost accounting](../../experiments/physical-factor-access-2026-09-14/README.md).
- [Ten-time nonlinear reference comparison](../../experiments/model-adequacy-2026-09-14/README.md).
- [Release files](https://github.com/orbion-life/allosteric-response-benchmark/releases/tag/v0.10.0), including the full experimental archives, report and editable source.
- [Asset sizes and SHA256 hashes](ASSET-MANIFEST-v0.10.0.json).

Download and extract the complete physical-factor archive before running its saved-data or circuit replay: its large input arrays and QPY files are release assets. The repository contains the source, protocols, summaries and audit programs. The nonlinear archive is self-contained and includes all five grids, both propagators, complete result tables and independent quadrature checks. Each README specifies its environment and separates saved-data verification from a fresh computational run.

Independent checks cover the local force algebra, all 64 selected ABL coefficient rows, the complete small polynomial model, all fourteen factor and walk columns, the seventeen-qubit response, every nonlinear response entry and the ordering calculation. Portable packages preserve scientific source and raw arrays; explicit hash mappings identify metadata redactions.

## Earlier results remain available

The [v0.9.0 guide](../evidence-v0.9.0/README.md) and [release](https://github.com/orbion-life/allosteric-response-benchmark/releases/tag/v0.9.0) retain the GRB2 and KRAS analyses, 288-matrix temporal study and compiled six-qubit diagnostic. The [v0.9.1 literature note](../evidence-v0.9.1/README.md) explains the mechanical screening and quantum estimation rationale. The new nonlinear comparison expands the original three-time test; it does not replace its arrays or erase earlier failures.

The report distinguishes physical model adequacy, response reduction, biological utility and quantum implementation. Those are separate tests, even when they share a Gaussian calculation. Data-source attribution and third-party terms remain in each experiment and the repository's licensing record.
