# CI grid-comparison portability, 15 September 2026

Linux workflow [34904837140](https://github.com/orbion-life/allosteric-response-benchmark/actions/runs/34904837140) failed the historical `log_grid_contained` check. Its source seals and all saved tiny-pulse checks passed. The original failed receipt is preserved in `initial-linux-failure.json`.

The historical verifier regenerates a NumPy geometric grid and requires exact floating-point membership in saved timestamps. The event-aligned verifier similarly requires exact list equality for regenerated lags. These regenerated floating-point values can differ in their last bits across environments. The original verifiers remain unchanged.

The CI-only `ci_protocol_check.py` adapter retains every historical check and records its exact result. For regenerated grid comparisons, it also requires a unique nearest stored node, at most 32 binary64 units in the last place per value, and distinct increasing node lists. It rejects duplicate, missing, ambiguous or materially shifted nodes. The event-aligned lag list must also have the exact expected length. Each receipt reports the maximum absolute, relative and ULP difference and the matched stored indices.

The adapter does not rewrite a stored timestamp, kernel mapping, protocol, source seal, response array or scientific tolerance. Every other original check remains binding. This change does not repair or relabel the historical protein temporal failures or bound-portability failures.

Three adversarial tests cover one-ULP rounding, values beyond the allowance, and missing, duplicated, shifted or ambiguous grid nodes. The existing source/protocol and tiny-pulse CI check passes locally with the adapter. A separate Linux workflow job runs all six receiver-correction tests independently, so failure of the historical temporal check cannot suppress those tests. Both jobs must pass for the workflow to pass. The subsequent Linux result is recorded by the workflow itself; this local note does not claim that rerun has completed.
