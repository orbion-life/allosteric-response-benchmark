# Version 0.4.1: corrected evidence descriptions and proposed recovery protocols

This is a documentation-only update. Scientific code, numerical data, original protocol settings, verification thresholds and all historical results are unchanged from v0.4.0. The existing strict Linux bound-portability failure remains unresolved.

The update makes four evidence descriptions precise:

- Four physics-control comparisons were completed within a six-slot Holm correction family; two MYH7 slots remain unrun.
- The four-state circuit's 1.64 × 10⁻⁷ normalized response error is sampled in the ideal simulator. Its deterministic same-grid error is 2.95 × 10⁻⁷.
- The original and v0.4.0-tag Linux replays each fail 61 of 372 bound arrays; the v0.4.0 main-branch replay fails 62. All 414 physical-response/covariance/time/scale comparisons pass in each run. These historical run counts do not predict the outcome of a future replay.
- The observable-preserving README now identifies the completed v0.4.0 release instead of describing its tag as pending. Same-host repeatability is explicitly qualified.

New proposed protocols define an explicit affine-coordinate physical reference, analytic harmonic calibration and separate grid/domain checks; a bounded identical-input subspace diagnostic; and a matched quantum comparison on the accepted original operator against sparse and reduced classical solvers. A recovery index includes the dependency diagram and stopping rules. Two small existing-array receipts and a prior known-minimum receipt support the design. No recovery campaign has been executed.

The content manifest checks every old entry. Only README descriptions and citation metadata change among prior files. Immutable prior release notes, source code, test code, workflow definitions, protocols, results and failure receipts remain unchanged. Current release checks validate documentation and source identity; they do not convert the known failed scientific or numerical acceptance checks into passes. Automatic numerical replays, if triggered, retain their unchanged verifier and report their actual outcomes separately.

This release contains public scientific methods and evidence. It excludes the manuscript, team profiles, submission details and private company records. Existing code and third-party licences continue to apply.
