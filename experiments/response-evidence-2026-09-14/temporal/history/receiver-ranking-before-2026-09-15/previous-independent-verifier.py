"""Independently reconstruct fresh query errors from saved matrices and kernels."""
from pathlib import Path
import json, sys, hashlib
import numpy as np
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]
T = ROOT / 'temporal'
P = json.loads((T / 'protocol.json').read_text())
def digest(p):
    with p.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()

with threadpool_limits(limits=3):
    for target in sys.argv[1:]:
        directory = T / 'results' / target
        z = np.load(directory / 'operator.npz')
        summary = json.loads((directory / 'analysis.json').read_text())
        lam, U = np.linalg.eigh(z['Hr'])
        V = U.T @ z['B']
        tau = float(z['tau'])
        kernels = np.load(directory / 'fresh-kernels.npy', mmap_mode='r')
        old_root = Path(P['original_root'])
        op = json.loads((old_root / 'protocol.json').read_text())
        old = np.load(old_root / 'results' / target / 'full-kernels.npy', mmap_mode='r')
        K0 = old[op['all_exact_kernel_times_over_tau'].index(0.0)]
        lookup = {float(t): i for i, t in enumerate(P['fresh_kernel_times_over_tau'])}
        meta = np.load(T / 'inputs' / target / 'ranking-metadata.npz')
        candidate = np.flatnonzero(meta['candidate'])
        receiver = np.flatnonzero(meta['receiver'])
        ids = meta['canonical']
        def full(t):
            return kernels[lookup[float(t)]] - K0
        def reduced(t):
            return V.T @ (np.expm1(-t * tau * lam)[:, None] * V)
        groups = [('fresh_step', None, summary['results']['fresh_step'])]
        groups += [('fresh_removal', float(s['duration_over_tau']), s)
                   for s in summary['results']['fresh_removal']]
        rows = []
        for name, duration, reported in groups:
            errors = []
            useful = identical = 0
            for t in P['fresh_step_times_over_tau']:
                if duration is None:
                    f, r = full(t), reduced(t)
                else:
                    f = full(duration + t) - full(t)
                    r = reduced(duration + t) - reduced(t)
                errors.append(float(np.max(abs(f-r))))
                fs = np.linalg.norm(f[:, receiver], axis=1) / np.sqrt(len(receiver))
                rs = np.linalg.norm(r[:, receiver], axis=1) / np.sqrt(len(receiver))
                if max(fs[candidate]) > .002:
                    useful += 1
                    fi = sorted(candidate, key=lambda k: (-fs[k], ids[k]))[:5]
                    ri = sorted(candidate, key=lambda k: (-rs[k], ids[k]))[:5]
                    identical += fi == ri
            err = max(errors)
            assert abs(err-reported['max_entry_error']) < 1e-8
            assert useful == reported['non_low_signal_query_count']
            assert identical == reported['non_low_signal_identical_top5_order_count']
            rows.append(dict(group=name, duration=duration, max_entry_error=err,
                             informative_queries=useful, identical_top5_order=identical))
        mass = float(np.max(abs(z['W'].T @ z['M'] @ z['W'] - np.eye(len(lam)))))
        static = float(np.max(abs(z['B'].T @ z['B'] - z['reference_G0'])))
        reported_mass = summary['checks']['mass_error']
        assert (mass <= 1e-8) == (reported_mass <= 1e-8)
        result = dict(status='PASS', target=target, queries_reconstructed=72,
                      independent_eigenvalue_min=float(min(lam)),
                      mass_error=mass, remote_mass_error=reported_mass,
                      mass_residual_platform_difference=abs(mass-reported_mass),
                      static_error=static, mass_gate_pass=mass <= 1e-8,
                      groups=rows, operator_sha256=digest(directory/'operator.npz'),
                      kernel_sha256=digest(directory/'fresh-kernels.npy'))
        out = ROOT/'review'/f'parent-temporal-{target}-verification.json'
        out.write_text(json.dumps(result, indent=2)+'\n')
        print(json.dumps(result), flush=True)
