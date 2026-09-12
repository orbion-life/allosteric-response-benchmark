"""Input-only thermal receiver-strain basis: a harmonic compression test."""
from pathlib import Path
import argparse, datetime, hashlib, json, resource, sys, time
import numpy as np
from scipy.linalg import eigh, svd, qr
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent/'harmonic-completion/vendor'))
from run_pilot import harmonic_covariances

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dot(a,b):return np.einsum('ij,jk->ik',a,b,optimize=False)

def run(pilot,out):
    protocol=json.loads((ROOT/'preanalysis.json').read_text())
    assert sha(ROOT/'preanalysis.json')==(ROOT/'preanalysis.sha256').read_text().strip()
    out.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    m=np.load(pilot/'model/static-model.npz');full=np.load(pilot/'results/analytic-harmonic-d492.npz')
    r0,edges,B,lam=m['r0'],m['edges'],m['B'],m['eigenvalues']
    tau=1/lam[0];sd=full['sd'];ref=full['C'];candidate=m['candidate'];receiver=m['receiver'];canonical=m['canonical']
    rec_edges=edges[np.any(np.isin(edges,receiver),axis=1)]
    J=np.zeros((len(rec_edges),r0.size))
    for k,(i,j) in enumerate(rec_edges):
        u=r0[i]-r0[j];u/=np.linalg.norm(u);J[k,3*i:3*i+3]=u;J[k,3*j:3*j+3]=-u
    whitened=B/np.sqrt(lam)
    thermal_J=dot(J,whitened)
    _,sing,Vh=svd(thermal_J,full_matrices=False,lapack_driver='gesvd')
    threshold=np.finfo(float).eps*max(thermal_J.shape)*sing[0]
    rank=int(np.sum(sing>threshold));assert rank>=max(protocol['dimensions'])
    rows=[]
    for d in protocol['dimensions']:
        if time.monotonic()-started>protocol['cap_seconds']:raise TimeoutError('locked compute cap')
        before=time.monotonic()
        cart=dot(whitened,Vh[:d].T);Q,_=qr(cart,mode='economic')
        Kr=dot(dot(Q.T,m['Hessian']),Q);lr,vr=eigh(Kr);Br=dot(Q,vr)
        for column in Br.T:
            if column[np.argmax(np.abs(column))]<0:column*=-1
        a,b=harmonic_covariances(r0,edges,Br,lr,tau)
        C=-(a-b)/np.outer(sd,sd)
        score=np.sqrt(np.mean(C[:,receiver]**2,axis=1));idx=np.flatnonzero(candidate)
        order=idx[np.lexsort((canonical[idx],-score[idx]))]
        block=np.ix_(candidate,receiver);error=float(np.max(np.abs((C-ref)[block])))
        raw=np.load(pilot/f'results/analytic-harmonic-d{d}.npz')['C']
        row=dict(d=d,seconds=time.monotonic()-before,receiver_strain_variance_retained=float(np.sum(sing[:d]**2)/np.sum(sing**2)),max_C_error=error,relative_Frobenius_error=float(np.linalg.norm((C-ref)[block])/np.linalg.norm(ref[block])),lowest_frequency_error=float(np.max(np.abs((raw-ref)[block]))),top5=canonical[order[:5]].tolist(),same_top5_set=bool(set(order[:5])==set(full['order'][:5])),coordinate_gate=error<=protocol['primary_gate']['max_C_error'],orthogonality_error=float(np.max(np.abs(dot(Br.T,Br)-np.eye(d)))),singular_gap_at_cutoff=float(sing[d-1]-sing[d]))
        np.savez_compressed(out/f'receiver-sensitive-d{d}.npz',C=C,covariance=a,time_covariance=b,score=score,order=order,canonical=canonical,basis=Br,eigenvalues=lr,common_harmonic_sd=sd)
        rows.append(row);print(json.dumps(row),flush=True)
    peak=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024)
    assert peak<=protocol['cap_peak_bytes']
    receipt=dict(status='COMPLETE',created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),protocol_sha256=sha(ROOT/'preanalysis.json'),script_sha256=sha(__file__),input_sha256={f:sha(pilot/f) for f in ['model/static-model.npz','results/analytic-harmonic-d492.npz']},receiver_contact_count=len(rec_edges),thermal_J_shape=list(thermal_J.shape),rank=rank,singular_values=sing.tolist(),seconds=time.monotonic()-started,peak_bytes=peak,results=rows,interpretation='Harmonic reference test only; no reference ligand labels read, no nonlinear fidelity or quantum cost conclusion')
    (out/'summary.json').write_text(json.dumps(receipt,indent=2)+'\n')
    (out/'freeze.json').write_text(json.dumps({p.name:sha(p) for p in sorted(out.iterdir()) if p.is_file()},indent=2)+'\n')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--pilot',type=Path,default=ROOT.parent/'kras/pilot');p.add_argument('--output',type=Path,default=ROOT/'results');a=p.parse_args();run(a.pilot,a.output)
