"""Verify compact saved evidence without running any grid or propagator."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parent


def run(root,legacy_reference=None):
    checks=[]
    for manifest in ['execution-source-hashes.json','capacity-source-hashes.json']:
        hashes=json.loads((root/manifest).read_text())['sha256']
        for name,expected in hashes.items():
            actual=hashlib.sha256((root/name).read_bytes()).hexdigest()
            checks.append(dict(kind='source_hash',file=name,manifest=manifest,pass_check=actual==expected))
    arrays=0
    for folder in ['harmonic-primary/calibration','nonlinear-primary/nonlinear-h0.5','nonlinear-primary/nonlinear-h0.25']:
        for f in sorted((root/folder).glob('*.npz')):
            with np.load(f) as z:
                for key in z.files:
                    arrays+=1
                    if not np.isfinite(z[key]).all():
                        raise AssertionError(f'Nonfinite saved array: {f.name}/{key}')
    analytic=np.load(root/'harmonic-primary/calibration/analytic.npz')
    legacy=None
    if legacy_reference:
        with np.load(legacy_reference) as z:
            differences={key:float(np.max(np.abs(analytic[new]-z[old]))) for key,new,old in [('harmonic_sd','sd','harmonic_sd'),('times','times','times')]}
        legacy=dict(artifact_name=legacy_reference.name,sha256=hashlib.sha256(legacy_reference.read_bytes()).hexdigest(),differences=differences,pass_check=all(v<1e-12 for v in differences.values()))
        checks.append(dict(kind='legacy_normalization_and_times',pass_check=legacy['pass_check']))
    with np.load(root/'nonlinear-primary/nonlinear-h0.5/response.npz') as a,np.load(root/'nonlinear-primary/nonlinear-h0.25/response.npz') as b:
        refinement={key:float(np.max(np.abs(a[key]-b[key]))) for key in ('G0','K','C')}
    solver=[]
    for folder in ['nonlinear-primary/nonlinear-h0.5','nonlinear-primary/nonlinear-h0.25']:
        with np.load(root/folder/'response.npz') as z:
            error=max(float(np.max(np.abs(z['K']-z['K_chebyshev']))),float(np.max(np.abs(z['C']-z['C_chebyshev']))))
            relation=float(np.max(np.abs(z['C']-(z['K']-z['G0']))))
        solver.append(dict(case=folder,maximum_independent_difference=error,response_identity_error=relation,pass_check=error<=1e-8 and relation<=1e-14))
        checks.append(dict(kind='solver_agreement',case=folder,pass_check=solver[-1]['pass_check']))
    continuation=json.loads((root/'continuation-protocol.json').read_text())
    checks.append(dict(kind='all_four_controlling_pairs',pass_check=len(continuation['numerical_screen_pairs'])==4 and all(p['controlling'] for p in continuation['numerical_screen_pairs'])))
    out=dict(status='PASS' if all(c['pass_check'] for c in checks) else 'FAIL',checks=checks,finite_arrays_checked=arrays,
             legacy_reference=legacy,old_box_coarse_refinement=refinement,
             old_box_refinement_screen_pass=all(v<=.001 for v in refinement.values()),solver_checks=solver,
             scope='Saved-array and source-hash checks only. No new grid, propagation or campaign. Scientific failed screens do not fail software/evidence verification.')
    return out


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--root',type=Path,default=ROOT);parser.add_argument('--legacy-reference',type=Path)
    parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    receipt=run(args.root,args.legacy_reference);args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({k:receipt[k] for k in ['status','finite_arrays_checked','legacy_reference','old_box_coarse_refinement']},indent=2))
    if receipt['status']!='PASS':raise SystemExit(1)
