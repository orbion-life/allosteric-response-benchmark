# Verify the temporal experiment from the repository

**Peak correction:** all original ancillary `peaks`, `peak` and `sampled_peak` fields are withdrawn. See [the receiver-indexing erratum](corrections/receiver-peak-indexing/README.md) and its corrected summary. The original sealed replay preserves the erroneous historical peak fields; the separate correction command applies the fix explicitly. All response-error gates and candidate-to-receiver errors were independently checked and remain unchanged.

This repository tree contains the frozen source and protocols, saved numerical summaries, and the complete small original-model pulse evidence. The companion `pulsar-temporal-response-evidence.zip` contains the fixed protein inputs and all eight immutable full-kernel arrays. `REPOSITORY-SELECTION.json` records the files deliberately kept outside Git; the compact publication manifest lists the full archive contents.

Install Python 3.13 and `requirements-ci.txt`, then run from this directory:

```sh
python -W error ci_check.py --output NEW_CHECK_DIRECTORY
```

The check verifies every copied public file hash and both frozen protocols. It audits the saved tiny-pulse arrays, constructs a fresh source tree without protein inputs, reruns the 64-state piecewise-field experiment, checks independent forward propagation, and compares every fresh tiny array with its saved counterpart. It also runs the corrected receiver-selector fixture with nonconsecutive odd/even indices. It does not allocate a GPU, replay any full protein response or corrected protein peak metric, or establish protein accuracy. All twelve event-aligned protein pulse curves retain their recorded failures of the 0.002 criterion. The scoped GitHub workflow runs this command on Linux and uploads the verification outputs.

For the complete protein waveform replay, extract the companion archive into `experiments/selection-2026-09-14/` from the repository root. Its members start with `temporal/`. This supplies the missing raw kernels and fixed inputs while preserving the identical archived source. From this temporal directory, run:

```sh
python replay_waveforms.py --output NEW_PROTEIN_REPLAY_DIRECTORY
python corrections/receiver-peak-indexing/correct_peaks.py --root NEW_PROTEIN_REPLAY_DIRECTORY --output NEW_PEAK_CORRECTION_DIRECTORY
```

The first command reconstructs the 28 deliberately omitted, algebraically derived waveform archives in a separate tree and runs both full saved-array audits. The second applies the explicitly documented receiver-index correction and independently rechecks the unaffected response and ranking metrics. Neither refits models or recomputes full Gaussian reference kernels. The original full audit, reconstructed-waveform audit, correction audit and source/tiny CI check have different scopes; none turns the retained scientific failures into successes.
