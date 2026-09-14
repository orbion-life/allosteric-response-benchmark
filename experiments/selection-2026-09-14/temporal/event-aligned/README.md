# Early response errors recur after finite perturbations are removed

**Correction, 2026-09-14:** the original ancillary peak fields used incorrect receiver-column indexing and are withdrawn. All original `peaks`, `peak` and `sampled_peak` fields are historical, including those recreated by the sealed replay. Use the [receiver-indexing erratum](../corrections/receiver-peak-indexing/README.md) and its separately corrected summary. Independent rechecking confirms that all response-error gates, signed-entry counts, RMS rankings and candidate-to-receiver error values are unchanged.

The initial absolute-time grid under-sampled the interval immediately after long pulse removals. This separate addendum was declared and sealed before the added kernel values or response scores were calculated. It preserves the original protocol, time grid, measurements and failures. No covariance, model, representation, basis or response threshold was changed.

Each duration (0.03, 0.3 and 3τ) now has 61 observations: exact removal and the original fixed logarithmic positive lags after removal, restricted to absolute t≤30τ. These use 180 new kernel times and 64 previous kernel times verified by SHA-256. The prior reference arrays remain unchanged.

| Target | d=0.03τ: max error / passing points | d=0.3τ: max error / passing points | d=3τ: max error / passing points |
|---|---|---|---|
| KRAS | 0.0377267 / 35/61 | 0.046568 / 36/61 | 0.0465222 / 36/61 |
| ABL | 0.0585032 / 39/61 | 0.0590528 / 39/61 | 0.0590525 / 39/61 |
| MYC | 0.0197916 / 49/61 | 0.0198672 / 49/61 | 0.0198672 / 49/61 |
| MYH7 | 0.0652329 / 43/61 | 0.0654595 / 43/61 | 0.0654595 / 43/61 |

These are maximum absolute signed matrix-entry errors in the fixed normalized response units, compared with the unchanged 0.002 threshold. The observation set starts at removal, so it directly measures the previously under-sampled post-removal response. The complete error arrays, sign counts and low-signal/gap-qualified rankings remain in `results/TARGET/analysis.json`.

The signed pulse equals K(d+s)−K(s) in linear response. The two terms are evaluated at their fixed declared arguments; the stored time-closure roundoff is also reported. Newly sampled errors can therefore expose deficiencies in the short-lag K(s) term even when a sustained response at the later absolute time appears accurate. This is numerical representation error within a fixed Gaussian model, not evidence of biological scrambling.

A source setup attempt failed before any new kernels because the direct CuPy eigensolver import did not locate libcusolver. Running the unchanged independent GPU fixture checks first, as in the successful original campaign, resolved the initialization issue. Both setup sources and the correction receipt are retained; the time/model protocol stayed unchanged.

The original peak fields are withdrawn because the shared historical analyzer selected incorrect receiver columns. Separately corrected fields mean unique interior sampled maxima under the fixed 1e-12 tie rule. They are not error-certified or continuous-time peaks; low-amplitude, flat/tied and boundary cases remain separate. No behavior between sampling points is certified.

Verification: `python verify_addendum.py --mode protocol --output protocol-audit.json` checks seals and exact grid construction without GPU work. With all original and addendum raw kernels and reconstructed waveforms present, `python verify_addendum.py --mode saved --output saved-verification.json` checks input hashes, full pulse identities, sign counts, rank arrays and error summaries. A passing implementation audit preserves the scientific failures.

The compact public bundle omits redundant derived `pulse-*.npz` files with explicit hashes. From the temporal package root, `python replay_waveforms.py --output NEW_DIRECTORY` reconstructs the original and event-aligned waveforms in a separate tree and performs both saved-array audits. This does not refit proteins or repeat the GPU kernel campaign.
