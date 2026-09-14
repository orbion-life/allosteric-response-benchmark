"""Portable saved-result audit; no downloads, model fits or private assay table required."""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import json,hashlib,sys
import numpy as np,pandas as pd
from scipy.stats import spearmanr
P=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 s=json.loads((P/'results/summary.json').read_text());assert sha(P/'PROTOCOL.md')==s['protocol_sha256'];assert sha(P/'results/predictions-before-label-join.csv')==s['prediction_sha256']
 m=np.load(P/'results/primary-model/static-model.npz');hm=np.load(P/'results/primary-model/harmonic-response.npz');raw=np.load(P/'results/primary-model/gaussian/raw-moments.npz');g=np.load(P/'results/primary-model/gaussian/response-matrices.npz');p=pd.read_csv(P/'results/predictions-before-label-join.csv');a=pd.read_csv(P/'results/site-aggregates-and-predictions.csv');e=a[a.eligible].copy().reset_index(drop=True);n=len(e);assert n==21
 assert np.array_equal(m['canonical'],np.arange(159,215));assert np.array_equal(np.flatnonzero(p.receiver),m['receiver'])
 dist=np.linalg.norm(m['r0'][:,None]-m['r0'][None],axis=-1);edges=np.array([(i,j) for i in range(56) for j in range(i+1,56) if dist[i,j]<=10 or j-i==1]);assert np.array_equal(edges,m['edges']);assert np.array_equal(np.bincount(edges.ravel(),minlength=56),p.degree)
 assert np.array_equal((dist[:,m['receiver']].min(axis=1)>=6)&~p.receiver,p.structural_candidate)
 den=np.outer(raw['harmonic_sd'],raw['harmonic_sd']);gc=(raw['delayed_covariance'][1]-raw['covariance'])/den;hc=(hm['delayed_covariance']-hm['covariance'])/den
 restore=max(float(np.max(abs(gc-g['C'][1]))),float(np.max(abs(hc-hm['C']))));assert restore<1e-11
 for name,C in [('gaussian',gc),('harmonic',hc)]:assert np.allclose(np.sqrt(np.mean(C[m['receiver']]**2,axis=0)),p[name],atol=1e-12,rtol=0)
 ohm=np.column_stack([np.loadtxt(P/'ohm'/f'seed{i}.nodes') for i in [11,29,47]]);assert np.allclose(ohm.mean(axis=1),p.ohm,atol=1e-12,rtol=0);assert (P/'ohm/seed11.nodes').read_bytes()==(P/'ohm/repeat11.nodes').read_bytes()
 names=['gaussian','harmonic','ohm','degree','negative_distance'];y=e.binding_magnitude.to_numpy();rho=np.array([spearmanr(e[name],y).statistic for name in names]);assert np.allclose(rho,[r['rho'] for r in s['metrics']],atol=1e-14,rtol=0)
 b=np.load(P/'results/paired-site-bootstrap.npz');assert b['site_indices'].shape==(20000,n);assert np.array_equal(b['canonical'],e.canonical)
 maxerr=0.
 # Independently call SciPy's direct paired routine for500 deterministic bootstrap indices.
 for k in np.linspace(0,19999,500,dtype=int):
  ix=b['site_indices'][k];check=np.array([spearmanr(e[name].to_numpy()[ix],y[ix]).statistic for name in names]);maxerr=max(maxerr,float(np.max(abs(check-b['rho'][k]))))
 assert maxerr<1e-13;assert np.allclose(b['gaussian_minus_controls'],b['rho'][:,[0]]-b['rho'][:,1:],atol=1e-14,rtol=0)
 for j,r in enumerate(s['gaussian_minus_controls']):
  v=b['gaussian_minus_controls'][:,j];assert np.allclose(np.quantile(v,[.025,.975]),r['paired_site_bootstrap95']);assert np.allclose(np.quantile(v,[.00625,.99375]),r['paired_site_bootstrap98_75'])
 for r in s['metrics']:
  top=e.sort_values([r['method'],'canonical'],ascending=[False,True]).head(5);assert top.canonical.tolist()==r['top5_canonical'];assert abs(top.binding_magnitude.mean()-r['top5_mean_binding_magnitude'])<1e-14
 assert not s['strong_within_assay_incremental_signal_gate'];assert len(pd.read_csv(P/'results/structure-sensitivity-metrics.csv'))==20
 if (P/'PUBLIC-MANIFEST.json').exists():
  manifest=json.loads((P/'PUBLIC-MANIFEST.json').read_text());assert all(sha(P/r['path'])==r['sha256'] for r in manifest['files'])
 print(json.dumps({'status':'PASS','canonical_sites':n,'full_response_restore_max_abs':restore,'independent_scipy_bootstrap_max_abs':maxerr,'bootstrap_resampling_unit':'residue, paired across all methods','raw_source_required':False,'scope':'Saved predictions, raw moments, masks, upstream Ohm means, reported statistics and intervals; does not independently establish experimental truth.'},indent=2))
if __name__=='__main__':main()
