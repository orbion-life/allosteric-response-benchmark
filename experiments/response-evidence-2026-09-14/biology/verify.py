"""Independent matrix-identity, ridge-fit and no-held-out-label-leakage checks."""
from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
from scipy.stats import rankdata,spearmanr
B=Path(__file__).resolve().parent
s=json.loads((B/'results/summary.json').read_text());checks={}
for target,r in s['assays'].items():
 d=pd.read_csv(B/'results'/f'{target}-per-site.csv');x=d[['gaussian','outcome','degree','negative_distance']].rank(method='average').to_numpy();precision=np.linalg.inv(np.corrcoef(x.T));partial=-precision[0,1]/np.sqrt(precision[0,0]*precision[1,1]);assert abs(partial-r['partial_spearman']['gaussian_degree_distance'])<1e-12
 maxerr=0.;loss=[]
 # Independent scalar ECDF and augmented least-squares ridge calculation.
 for f in range(5):
  tr=d[d.fold!=f];te=d[d.fold==f]
  def pct(t,q):return np.array([np.mean(t<v)+.5*np.mean(t==v) for v in q])
  for name,features in [('geometry',['degree','negative_distance']),('geometry_gaussian',['degree','negative_distance','gaussian']),('geometry_harmonic',['degree','negative_distance','harmonic']),('geometry_harmonic_gaussian',['degree','negative_distance','harmonic','gaussian'])]:
   X=np.column_stack([pct(tr[c].to_numpy(),tr[c].to_numpy()) for c in features]);V=np.column_stack([pct(tr[c].to_numpy(),te[c].to_numpy()) for c in features]);mu=X.mean(0);sd=X.std(0);sd[sd<1e-12]=1.;Z=np.column_stack([np.ones(len(X)),(X-mu)/sd]);T=np.column_stack([np.ones(len(V)),(V-mu)/sd]);Y=pct(tr.outcome.to_numpy(),tr.outcome.to_numpy());reg=np.diag([0]+[1]*len(features));fit=np.linalg.lstsq(np.vstack([Z,reg]),np.r_[Y,np.zeros(len(features)+1)],rcond=None)[0];pred=T@fit;maxerr=max(maxerr,float(np.max(abs(pred-te[name+'_prediction']))));assert np.allclose(pred,te[name+'_prediction'],atol=1e-12,rtol=0)
  assert set(tr.canonical).isdisjoint(te.canonical)
 # Directly verify saved bootstrap first 100 samples with statsmodels-style rank residuals.
 z=np.load(B/'results'/f'{target}-sequence_block10-bootstrap.npz')
 for i in range(100):
  ix=z['site_indices'][i,:z['sample_lengths'][i]];b=d.iloc[ix];a=b[['degree','negative_distance']].rank().to_numpy();X=np.column_stack([np.ones(len(b)),a]);v=b[['gaussian','outcome']].rank().to_numpy();res=v-X@np.linalg.pinv(X)@v;den=np.linalg.norm(res[:,0])*np.linalg.norm(res[:,1]);rho=(res[:,0]@res[:,1])/den if den>1e-10 else np.nan;assert np.allclose(rho,z['partial'][i,0],atol=1e-10,equal_nan=True)
 checks[target]={'partial_inverse_correlation_identity':float(partial),'all_fold_predictions_independently_reproduced_max_abs_error':maxerr,'fold_disjointness':True,'bootstrap_first100_directly_reproduced':True}
# Use implementation but mutate held-out outcomes to detect accidental exposure in transformations/fits.
import analyze
for target in s['assays']:
 original=pd.read_csv(B/'inputs'/f'{target}.csv');allpos=np.sort(original.canonical.to_numpy());d=pd.read_csv(B/'results'/f'{target}-per-site.csv')
 for f in range(5):
  changed=d.copy();mask=changed.fold==f;changed.loc[mask,'outcome']=np.arange(mask.sum())*12345.+9876.;pred,_,_=analyze.oof(changed,allpos)
  for name in analyze.MODEL_FEATURES:
   assert np.allclose(pred.loc[mask,name+'_prediction'],d.loc[mask,name+'_prediction'],atol=1e-12,rtol=0)
 checks[target]['heldout_outcome_mutation_leaves_predictions_unchanged']=True
# Replay arrays and numerical summary match; timestamps and runtime are naturally new.
for target in s['assays']:
 for suffix in ['per-site.csv','fold-models.json']:
  p=f'{target}-{suffix}';assert (B/'results'/p).read_bytes()==(B/'replay'/p).read_bytes()
 for kind in ['site','sequence_block10']:
  p=f'{target}-{kind}-bootstrap.npz';a=np.load(B/'results'/p);b=np.load(B/'replay'/p);assert a.files==b.files
  for k in a.files:
   if a[k].dtype.kind in 'fci':assert np.array_equal(a[k],b[k],equal_nan=True)
   else:assert np.array_equal(a[k],b[k])
r=json.loads((B/'replay/summary.json').read_text());assert s['assays']==r['assays'];out={'status':'PASS','checks':checks,'replay_all_csv_json_and_bootstrap_arrays_equal':True,'summary_sha256':hashlib.sha256((B/'results/summary.json').read_bytes()).hexdigest()};(B/'verification.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
