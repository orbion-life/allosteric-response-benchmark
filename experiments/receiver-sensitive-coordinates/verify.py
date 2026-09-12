"""Check the recorded selector and optional fresh replay without choosing a basis in a tie."""
from pathlib import Path
import argparse, hashlib, json
import numpy as np
ROOT=Path(__file__).resolve().parent
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def run(pilot, results, repeat, output):
    receipt=json.loads((results/'summary.json').read_text())
    audit=json.loads((ROOT/'admission-audit.json').read_text())
    s=np.array(receipt['singular_values']); checked=[]
    for row in receipt['results']:
        d=row['d'];z=np.load(results/f'receiver-sensitive-d{d}.npz')
        gap=(s[d-1]-s[d])/s[0]
        admissible=gap>1e-8
        assert np.isfinite(z['C']).all()
        B=z['basis'];assert np.max(np.abs(np.einsum('ij,ik->jk',B,B,optimize=False)-np.eye(d)))<1e-10
        assert row['max_C_error']>.002
        r=dict(d=d,relative_singular_gap=float(gap),admissible_cut=bool(admissible),accuracy_pass=False)
        if repeat is not None:
            fresh=np.load(repeat/f'receiver-sensitive-d{d}.npz')
            # A cut inside a tied singular subspace has no unique selected span.
            if admissible:
                r['replay_C_error']=float(np.max(np.abs(z['C']-fresh['C'])))
                assert np.allclose(z['C'],fresh['C'],atol=1e-10,rtol=1e-10)
                assert np.array_equal(z['order'],fresh['order'])
            else:
                r['replay_check']='Finite arrays and orthogonality only; tied cut is scientifically inadmissible.'
                bf=fresh['basis'];assert np.isfinite(fresh['C']).all()
                assert np.max(np.abs(np.einsum('ij,ik->jk',bf,bf,optimize=False)-np.eye(d)))<1e-10
        checked.append(r)
    assert [r['d'] for r in checked if r['admissible_cut']]==[64]
    inputs=['model/static-model.npz','results/analytic-harmonic-d492.npz']+[f'results/analytic-harmonic-d{d}.npz' for d in [2,4,8,16,32,64]]
    result=dict(status='PASS: verification succeeds; representation fails',checks=checked,
                source_sha256={p:sha(pilot/p) for p in inputs},
                vendor_sha256=sha(ROOT.parent/'harmonic-completion/vendor/run_pilot.py'),
                admission_audit_sha256=sha(ROOT/'admission-audit.json'),
                note='The cutoff-admission audit was added after inspecting results; no failed dimension was replaced.')
    output.parent.mkdir(parents=True,exist_ok=True);output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--pilot',type=Path,default=ROOT.parent/'kras/pilot');p.add_argument('--results',type=Path,default=ROOT/'results');p.add_argument('--repeat',type=Path);p.add_argument('--output',type=Path,default=ROOT/'checks/verification.json');a=p.parse_args();run(a.pilot,a.results,a.repeat,a.output)
