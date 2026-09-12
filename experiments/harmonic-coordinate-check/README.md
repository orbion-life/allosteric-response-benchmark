# An independent harmonic-coordinate check

This calculation uses the original synthetic five-node geometry at stiffness 100, inverse thermal energy 1, mobility 1 and one slow-mode relaxation time. It computes the same quartic contact-energy observable with two independent Gaussian moment/Hermite routes and a direct-contact quadrature check.

Run `python verify_math.py` with NumPy and SciPy (the KRAS pilot pinned environment is sufficient). The script writes `math-calculations-independent.json`. No input download is required.

All retained dimensions use common nine-mode harmonic standard deviations. The receiver-4/sender-3 value in zero-based indexing is −0.0850991019 with four modes and −0.0281769023 with all nine. The absolute difference is 0.0569221996, or 28.46 times the provisional 0.002 coordinate allowance. The JSON uses the equivalent one-based label R54/C54.

This rejects those four coordinates as a validated approximation under this criterion. It does not test full nonlinear protein behavior. Shot arithmetic is a sufficient overlap-estimation illustration, not a lower bound or a hardware execution.
