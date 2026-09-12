"""Apply the locked additive correction to frozen KRAS arrays, without refitting.

This is a post-prediction diagnostic, not a replacement primary prediction or
an accuracy test against a full nonlinear protein reference (none exists).
"""
from pathlib import Path
import argparse, csv, hashlib, json
import numpy as np
ROOT = Path(__file__).resolve().parent

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def run(pilot, labels, output):
    output.mkdir(parents=True, exist_ok=False)
    protocol = json.loads((ROOT/'preanalysis.json').read_text())
    assert sha(ROOT/'preanalysis.json') == (ROOT/'preanalysis.sha256').read_text().strip()
    assert 'protein_diagnostic' in protocol
    model = np.load(pilot/'model/static-model.npz')
    full = np.load(pilot/'results/analytic-harmonic-d492.npz')
    canon, receiver = model['canonical'], model['receiver']
    candidates = np.flatnonzero(model['candidate'])
    with labels.open() as stream:
        reference = {int(row['canonical']): row for row in csv.DictReader(stream)}
    rows, matrices, sources = [], [], [pilot/'model/static-model.npz', pilot/'results/analytic-harmonic-d492.npz', labels]
    for n in (33, 65):
        files = [pilot/f'results/{kind}-d2-n{n}-e4.npz' for kind in ('biquadratic', 'harmonic')]
        nl, hg = [np.load(p) for p in files]
        C = full['C'] + nl['C'] - hg['C']
        score = np.sqrt(np.mean(C[:, receiver]**2, axis=1))
        order = candidates[np.lexsort((canon[candidates], -score[candidates]))]
        assert np.isfinite(C).all() and C.shape == (166, 166)
        np.savez_compressed(output/f'completed-kras-d2-n{n}.npz', C=C, score=score, order=order, canonical=canon)
        top = []
        for idx in order[:5]:
            row = reference[int(canon[idx])]
            top.append(dict(canonical=int(canon[idx]), score=float(score[idx]),
                            label_state=row['label_state'],
                            MOV_distance_A=float(row['MOV_distance_A']) if row['MOV_distance_A'] else None))
        rows.append(dict(n=n, top5=top, score_gap_5_6=float(score[order[4]]-score[order[5]])))
        sources += files
        matrices.append(C)
    receipt = dict(status='COMPLETED DIAGNOSTIC; original predictions unchanged',
                   scope='No full nonlinear KRAS reference exists; contact recovery cannot certify the correction.',
                   protocol_sha256=sha(ROOT/'preanalysis.json'), script_sha256=sha(__file__),
                   source_sha256={str(p.relative_to(pilot)) if p.is_relative_to(pilot) else p.name:sha(p) for p in sources},
                   grid_max_C_difference=float(np.max(np.abs(matrices[1]-matrices[0]))), cases=rows)
    (output/'diagnostic.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps(receipt, indent=2))

if __name__ == '__main__':
    p=argparse.ArgumentParser(); p.add_argument('--pilot',type=Path,default=ROOT.parent/'kras/pilot')
    p.add_argument('--labels',type=Path,default=ROOT.parent/'kras/evaluation-amendment/results/reference-labels.csv')
    p.add_argument('--output',type=Path,default=ROOT/'protein-diagnostic')
    a=p.parse_args(); run(a.pilot,a.labels,a.output)
