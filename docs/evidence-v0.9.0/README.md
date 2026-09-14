# Project Pulsar response evidence in version 0.9.0

Version 0.9.0 adds an exploratory comparison of Gaussian response with structural predictors, a four-protein transient preservation study, and compiled six-qubit response circuits. The [experiment guide](../../experiments/response-evidence-2026-09-14/README.md) links each method, result and replay command.

GRB2 shows a 9.92% reduction in held-out binding-percentile mean squared error when Gaussian response is added to degree and distance. KRAS shows 1.02%. These are within-protein exploratory results, not independent biological validation or precision-at-five gains. All 288 fresh Gaussian transient matrices pass the 0.002 response threshold, while KRAS still fails numerical orthogonality. The six-qubit study remains local compilation and simulation; its proposed hardware budgets and the separate shallow-amplification study are not hardware executions.

The [quantum rationale](quantum-response-rationale.md) explains the mathematical motivation, direct operator-access requirement and strongest classical comparison. The [release asset manifest](ASSET-MANIFEST-v0.9.0.json) records every evidence archive. Historical model and hardware failures remain available in earlier releases and the included records.

The [proposed nonlinear model comparison](nonlinear-model-adequacy.md) fixes seven additional response times and separates reference qualification, amplitude fidelity and the ordering of the triangle's two other sites. It has not been executed.

The [scientific report](https://github.com/orbion-life/allosteric-response-benchmark/releases/download/v0.9.0/project-pulsar-scientific-report.pdf) and [editable source](https://github.com/orbion-life/allosteric-response-benchmark/releases/download/v0.9.0/project-pulsar-editable-source.zip) connect the completed results and proposed studies. The editable source reproduces all nine pages with identical text and rendered pixels.
