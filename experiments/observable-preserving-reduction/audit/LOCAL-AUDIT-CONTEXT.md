# Context of the archived local audit

`review_results.py` is the original read-only audit used before packaging. Its project-relative lookup of the historical Git checkout belongs to the original workspace layout. It is preserved as audit history, not supplied as the portable reproduction command. Its immutable-release parity receipt identifies all six historical Git blobs and their hashes.

Use the experiment-root `verify.py --replay PATH --output RECEIPT.json` for a portable check. It does not require a Git checkout elsewhere or a private workspace path. The full public mathematical derivation is `../docs/operator-mathematical-audit.md`, relative to this audit directory.

The original same-host audit compared all 1,065 arrays, including basis-dependent coordinates. The portable cross-platform policy compares the 786 physical covariance, response, bound, time and scale arrays at the same absolute 1e-9 allowance. It also checks all 279 reduced-coordinate arrays for finite values, dimensions, symmetry/positivity where applicable, and reconstruction of the physical results. It does not compare their raw entries across basis gauges.

This is an explicit comparison-policy distinction made during packaging, before the new Linux campaign. It changes neither the experiment's scientific source, numerical parameters, original arithmetic allowance nor any frozen result. Cost decisions are recomputed for the replay but need not match another machine's timings. A software-verification pass retains every failed scientific gate.
