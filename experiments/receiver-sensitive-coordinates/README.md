# Can receiver-strain coordinates preserve the KRAS response?

**The only tested cutoff that avoids a tied singular subspace, 64 coordinates, fails the response criterion.** Its maximum candidate-to-receiver harmonic response error is 0.198519, compared with 0.129861 for 64 low-frequency modes and the 0.002 allocation. These are harmonic-reference comparisons, not measured full nonlinear protein errors.

The selector uses the frozen input receiver and native contact geometry only. Let J contain the linear bond-strain gradients of receiver-incident contacts, B the complete positive harmonic eigenbasis and Λ its stiffnesses. Leading right singular vectors of J B Λ⁻½ select thermal receiver-strain directions. The resulting Cartesian span is orthonormalized and its projected Hessian diagonalized. Exact Gaussian moments of the original quartic residue observables give the response, with unchanged full-harmonic scales and time. No validation pocket labels or response-matrix fitting enter selection.

The preanalysis fixes dimensions 2, 4, 8, 16, 32 and 64. The first 54 singular values are tied, so the smaller cuts choose arbitrary parts of a degenerate subspace. Their recorded outcomes are diagnostics for one numerical orientation and do not establish that every orientation fails. The later admission audit explicitly marks those cuts inadmissible; no replacement dimension was selected after seeing results. The 64-coordinate cutoff has a nonzero gap and also fails to preserve the reference shortlist.

| Dimension | Maximum C error | Cutoff admissible? |
|---:|---:|:---|
| 2 | 0.218220 | No |
| 4 | 0.217605 | No |
| 8 | 0.217111 | No |
| 16 | 0.214181 | No |
| 32 | 0.209810 | No |
| 64 | 0.198519 | Yes |

## Reproduce

Use Python 3.12 from this directory and the sibling frozen KRAS pilot:

```sh
python3.12 -m venv .venv
.venv/bin/pip install -r requirements.txt
VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python run.py --output replay
.venv/bin/python verify.py --repeat replay
```

The locked cap is ten minutes and 2 GB; the local runs took under one minute. The verifier compares the unique 64-coordinate response and ranking, and checks finite values and orthogonality at tied cuts rather than demanding an arbitrary singular-vector orientation across libraries. Input and imported-function hashes are recorded. The 0.002 threshold is an engineering allocation, not a biologically validated cutoff.

`preanalysis.json` and its hash preserve the original specification. `admission-audit.json` records the later degeneracy correction. `results/` retains every tested dimension. `checks/verification.json` records a fresh same-host replay. Successful verification means repeatable computation, not successful compression. Original code is MIT licensed; cached KRAS structures retain their upstream data terms.
