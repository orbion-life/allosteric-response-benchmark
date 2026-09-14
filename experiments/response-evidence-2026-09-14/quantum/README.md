# Six qubit circuits estimate eighteen Gaussian response entries

This local study compiles the nine-mode, three-site Gaussian fixture to six qubits on a documented IQM Garnet connectivity subset. Eighteen data circuits use 120–164 CZ gates, 313–454 PRx rotations and depth 279–415. Their largest ideal response error is 1.77 × 10⁻⁷. All circuits, native mappings, assumed-noise calculations and statistical records are included. **No quantum hardware was executed.**

The direct-readout engineering plan uses 376,832 shots over 25 tasks, including seven controls. Its calculated on-demand QPU charge is US$553.91 at the 14 September 2026 AWS tariff; it targets a separate 0.05 engineering tolerance only if circuit bias is qualified. The conservative 0.002 plan requires 54,443,954 shots and about US$79,767 under an assumed separately established physical bias of at most 0.0005. These are sufficient plans for this sampler, not minimum costs for all algorithms. Simulated noise channels are assumptions, not current device calibrations.

## Reproduce the saved evidence

The repository contains the small fixture, circuits and compact results. The release asset `pulsar-six-qubit-diagnostic-evidence.zip` also includes the first numerical-allowance failure and its original circuit records. Verify its SHA256 before extraction. From `quantum` run:

```bash
python -m pip install -r requirements.txt
python verify.py
```

`verify.py` independently applies the exported PRx/CZ instructions, reconstructs all probability intervals and checks shot and monetary accounting. It does not contact a provider. For complete regeneration, make a fresh copy, rename the original `results` to `results-recorded`, and run `python run_study.py`; the entrypoint refuses to overwrite an existing results directory. That regeneration uses only local classical simulation.

`PROTOCOL.md` states the original plan and `AMENDMENT.md` preserves the numerical allowance correction. The proposed shallow-amplification experiment in `PROPOSED-AMPLIFICATION.md` has not been run. Counts from the earlier rank 347 KRAS operator do not apply to the newer temporal representations.
