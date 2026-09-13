**Archived planning note.** The following storage recommendation preceded packaging. The current completed archive list and verification receipts are in [provenance/archive-index.json](provenance/archive-index.json).

# Evidence storage and a portable release

The saved-array inspection counted 1,824,081,384 bytes in 3,021 files before its final receipt and documentation were added. The first result tree accounted for 1,343,195,818 bytes, repaired traces for 479,355,186 bytes, and identical-basis driver files for 629,417,854 bytes. These sizes are uncompressed file sizes, not a measured ZIP size.

Keep the scientific source, instrumenter, coordinator, portable bootstrap, protocol history, compact findings, environment summaries, complete file manifest and small test receipts in the source package. Large raw traces should be separate, checksum-bound release assets. Do not put approximately 1.8 GB of repeated arrays into ordinary Git history.

The most useful raw-asset division is: native input bundles; pinned complete traces; original incomplete system traces; corrected system traces; and fixed-basis driver traces plus isolated residuals. Every asset needs a manifest of its member paths, byte sizes and SHA-256 hashes. Keep the original failed execution records in the history asset. No raw trace is disposable merely because a corrected replay exists.

A future content-addressed representation can remove repeated array payloads while preserving each original NPZ's fields, dtypes, shapes and array bytes. It should additionally preserve or reconstruct the original ZIP metadata if it promises byte-identical NPZ files. Semantic reconstruction alone must not be described as preserving the original file hash. The seed and append-only admitted blocks already permit full-basis reconstruction between checkpoints, but the full checkpoint V files should remain in this evidence release until a tested reconstruction receipt exists.

The existing receipts preserve actual local executable and work paths. A public package should use separate sanitized display copies, with explicit links to the hashes of their original receipts. Do not alter an archived protocol while retaining its SHA-256. The new `bootstrap_replay.py` accepts interpreter paths as arguments, writes a new dated protocol, and uses only relative source paths inside the generated replay directory.

No archive compression, deduplication, public upload or fresh numerical campaign was performed as part of this packaging note. The parent task owns any later public release and its license/provenance review.
