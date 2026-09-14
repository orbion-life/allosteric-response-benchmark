# Computing ABL responses from contact factors and testing nonlinear fidelity

This release adds a physical-factor implementation of the fitted ABL Gaussian model and a ten-time comparison with the original nonlinear reference.

The ABL experiment reconstructs the coordinate drift from 2,146 local contacts and implements a complete response circuit for a three-residue fragment. Its seventeen-qubit ideal replay agrees with the fitted Gaussian response within 6.594 × 10⁻⁹. The archive contains the actual circuits, coefficient arrays, normalization, independent checks and classical timings. Full-protein coherent loading remains unresolved; the implemented explicit-table route exceeds its resource budget.

The nonlinear experiment retains the fixed Gaussian and harmonic controls, five reference grids and two independent propagators. Neither approximation meets the 0.002 response threshold across any of its 90 entries at ten times. Both preserve 16 site-order comparisons, reverse 12 and leave two unresolved under the stated empirical reference allowance. The complete arrays and failures are included.

The updated scientific report connects these results to the earlier protein, biological and six-qubit studies. Its editable archive rebuilds all nine pages with identical text and rendered pixels in the recorded environment.

Start with the [experiment guide](https://github.com/orbion-life/allosteric-response-benchmark/blob/v0.10.0/docs/evidence-v0.10.0/README.md). Download the full experimental archives below for raw inputs and saved circuit replay; the [asset manifest](https://github.com/orbion-life/allosteric-response-benchmark/blob/v0.10.0/docs/evidence-v0.10.0/ASSET-MANIFEST-v0.10.0.json) records file sizes and SHA256 hashes. Earlier releases preserve the biological and temporal data, original experiments and failure history.
