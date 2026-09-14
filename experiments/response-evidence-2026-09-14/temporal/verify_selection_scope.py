"""Check historical selection RMS scores from the saved waveform arrays.

The separate historical ancillary peak fields are intentionally not validated.
"""
from pathlib import Path
import argparse, hashlib, json
import numpy as np


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def verify(root):
    rows = []
    for target in ['kras', 'abl', 'myc', 'myh7']:
        with np.load(root/'inputs'/target/'ranking-metadata.npz') as meta:
            receiver=meta['receiver']; candidate=meta['candidate']; ids=meta['canonical']
        assert receiver.dtype.kind in 'iu' and candidate.dtype.kind == 'b'
        candidates=[i for i,included in enumerate(candidate) if included]
        for prefix in ['', 'event-aligned']:
            for path in sorted((root/prefix/'results'/target).glob('*.npz')):
                if path.name != 'step.npz' and not path.name.startswith('pulse-'):
                    continue
                with np.load(path) as saved:
                    discrepancies=[]
                    for matrix_key, score_key, order_key in [('full','reference_score','reference_order'),('reduced','score','order')]:
                        matrices=saved[matrix_key]
                        score=np.linalg.norm(np.take(matrices,receiver,axis=2),axis=2)/np.sqrt(len(receiver))
                        discrepancies.append(float(np.max(abs(score-saved[score_key]))))
                        assert np.allclose(score,saved[score_key],atol=1e-12,rtol=0)
                        # Compare complete candidate ordering using the original RMS
                        # accumulation, which avoids changing order at numerical ties.
                        exact=np.sqrt(np.mean(np.take(matrices,receiver,axis=2)**2,axis=2))
                        order=np.array([sorted(candidates,key=lambda i:(-r[i],ids[i])) for r in exact])
                        assert np.array_equal(order,saved[order_key])
                rows.append(dict(target=target,file=str(path.relative_to(root)),sha256=digest(path),queries=len(score),maximum_score_discrepancy=max(discrepancies),full_candidate_orders_match=True))
    assert len(rows)==28
    return dict(status='PASS',scope='Historical selection receiver RMS and complete candidate rankings, checked from all 28 saved waveform archives. The separately withdrawn ancillary peak fields remain withdrawn.',waveforms=len(rows),queries=sum(r['queries'] for r in rows),rows=rows)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--root',type=Path,required=True);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    if args.output.exists():parser.error('Output exists')
    result=verify(args.root);args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(result,indent=2)+'\n');print(result['status'],result['waveforms'],'waveforms',result['queries'],'queries')
