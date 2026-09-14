# Gaussian response adds exploratory information in the GRB2 assay

Adding Gaussian response to contact degree and receiver distance reduced held-out binding-percentile mean squared error by **9.92% for GRB2** and **1.02% for KRAS**. The GRB2 increment remained **4.56%** after harmonic response was also included. These are supervised exploratory combinations on previously inspected assays, with five consecutive sequence sections held out in turn. They are not gains in precision at five, independent protein validation or a demonstrated improvement in the original unsupervised shortlist.

The fixed protocol retains both outcomes and all comparator results. Conditional sequence-block bootstrap intervals resample saved out-of-fold losses without refitting models; they do not account for training instability or fully resolve spatial dependence. Gaussian and harmonic response still nominate the same five GRB2 residues.

## Reproduce the calculation

Download `pulsar-incremental-biology-evidence.zip`, check its SHA256 against the release asset manifest, and extract it. From its `biology` directory, run:

```bash
python -m pip install -r requirements.txt
python verify.py
python analyze.py replay-fresh
```

The saved-result verification independently reconstructs partial correlations, every held-out prediction and the first 100 block bootstrap draws, then checks the original same-host replay. A fresh analysis performs the complete 20,000-draw site and block bootstraps and writes a new output directory; keep the recorded results intact. Floating-point versions and platform can affect bytes, so evaluate numerical differences separately from file hashes.

The repository includes compact inputs and results. The archive additionally includes all bootstrap draws and the original replay arrays. Input CSVs are our derived site aggregates, not the publisher's mutation-level data or workbook. See `SOURCES-AND-RIGHTS.md`, `PROTOCOL.md`, `PACKAGING.json` and `verification.json`.
