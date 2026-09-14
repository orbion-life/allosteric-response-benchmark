"""Prespecified GRB2 functional join and paired site-level uncertainty."""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import json,hashlib,time,datetime,resource,sys
import numpy as np,pandas as pd
from scipy.stats import rankdata
P=Path(__file__).resolve().parents[1]
METHODS=['gaussian','harmonic','ohm','degree','negative_distance']
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def corr(x,y):
 x=rankdata(x,axis=-1);y=rankdata(y,axis=-1);x=x-x.mean(axis=-1,keepdims=True);y=y-y.mean(axis=-1,keepdims=True);den=np.sqrt(np.sum(x*x,axis=-1)*np.sum(y*y,axis=-1))
 return np.divide(np.sum(x*y,axis=-1),den,out=np.full(np.shape(den),np.nan),where=den>0)
def partial(x,y,c):
 X=np.column_stack([np.ones(len(y))]+[rankdata(a) for a in c]);xr=rankdata(x);yr=rankdata(y);xr=xr-X@np.linalg.lstsq(X,xr,rcond=None)[0];yr=yr-X@np.linalg.lstsq(X,yr,rcond=None)[0];den=np.linalg.norm(xr)*np.linalg.norm(yr)
 return float(xr@yr/den) if den>1e-12 else None

def main():
 start=time.monotonic();freeze=json.loads((P/'protocol-freeze.json').read_text());assert sha(P/'PROTOCOL.md')==freeze['protocol_sha256'];receipt=json.loads((P/'results/prediction-receipt.json').read_text());assert sha(P/'results/predictions-before-label-join.csv')==receipt['prediction_sha256']
 predictions=pd.read_csv(P/'results/predictions-before-label-join.csv');source=pd.read_csv(P/'local-only/faure-table7.csv');df=source[(source.protein=='GRB2-SH3')&(source.mut_order==1)].copy()
 seq='TYVQALFDFDPQEDGELGFRRGDFIHVMDNSDPNWWKGACHGQTGMFPRNYVTPVN';df['canonical']=df.Pos.astype(int)+158
 mapping=np.array([1<=int(r.Pos)<=56 and seq[int(r.Pos)-1]==r.WT_AA and r.Mut in 'ACDEFGHIKLMNPQRSTVWY' and r.Mut!=r.WT_AA for r in df.itertuples()]);fields=['b_ddg_pred','f_ddg_pred','b_ddg_pred_sd','f_ddg_pred_sd'];df[fields]=df[fields].apply(pd.to_numeric,errors='coerce');finite=np.isfinite(df[fields]).all(axis=1);positive=(df.b_ddg_pred_sd>0)&(df.f_ddg_pred_sd>0);df['valid']=mapping&finite&positive;df['exclusion']=np.select([~mapping,~finite,~positive],['invalid_variant_mapping','missing_nonfinite_estimate_or_uncertainty','nonpositive_uncertainty'],default='');df.to_csv(P/'local-only/variants-with-exclusions.csv',index=False)
 good=df[df.valid].copy();rows=[]
 for canonical,g in good.groupby('canonical'):
  b=g.b_ddg_pred.to_numpy();f=g.f_ddg_pred.to_numpy();bw=1/g.b_ddg_pred_sd.to_numpy()**2;fw=1/g.f_ddg_pred_sd.to_numpy()**2
  rows.append(dict(canonical=int(canonical),mutation_count=len(g),binding_magnitude=float(abs(b).mean()),folding_magnitude=float(abs(f).mean()),weighted_binding_magnitude=float(np.average(abs(b),weights=bw)),weighted_folding_magnitude=float(np.average(abs(f),weights=fw))))
 site=pd.DataFrame(rows);d=predictions.merge(site,on='canonical',how='left',validate='one_to_one');d['eligible']=d.structural_candidate & (d.mutation_count>=10);d['exclusion']=np.select([d.receiver,~d.structural_candidate,~(d.mutation_count>=10)],['receiver','less_than_6A_from_receiver','fewer_than_10_common_valid_substitutions'],default='');d.to_csv(P/'results/site-aggregates-and-predictions.csv',index=False)
 e=d[d.eligible].copy().reset_index(drop=True);n=len(e);assert n>=5,'Fewer than5 eligible sites: top-five endpoint unavailable.'
 y=e.binding_magnitude.to_numpy();X=e[METHODS].to_numpy().T;rho=corr(X,y);rng=np.random.default_rng(20260914);idx=rng.integers(0,n,size=(20000,n));yr=y[idx];boots=np.stack([corr(x[idx],yr) for x in X],axis=1);diff=boots[:,[0]]-boots[:,1:];np.savez_compressed(P/'results/paired-site-bootstrap.npz',site_indices=idx,methods=np.array(METHODS),rho=boots,gaussian_minus_controls=diff,canonical=e.canonical.to_numpy())
 metrics=[];tops=[]
 for j,name in enumerate(METHODS):
  order=e.sort_values([name,'canonical'],ascending=[False,True]).index.to_numpy();top=order[:5];valid=boots[:,j][np.isfinite(boots[:,j])];r=dict(method=name,n=n,rho=float(rho[j]),rho_site_bootstrap95=np.quantile(valid,[.025,.975]).tolist(),undefined_bootstrap_draws=int(20000-len(valid)),top5_canonical=e.loc[top,'canonical'].astype(int).tolist(),top5_mean_binding_magnitude=float(y[top].mean()),top5_mean_folding_magnitude=float(e.loc[top,'folding_magnitude'].mean()),folding_partial_rho=partial(e[name],y,[e.folding_magnitude]),weighted_aggregation_rho=float(corr(e[name],e.weighted_binding_magnitude)))
  if name in ['gaussian','harmonic','ohm']:r['folding_degree_distance_partial_rho']=partial(e[name],y,[e.folding_magnitude,e.degree,e.receiver_CA_distance_A])
  metrics.append(r);t=e.loc[top,['canonical','binding_magnitude','folding_magnitude','mutation_count',name]].copy();t=t.rename(columns={name:'prediction_score'});t.insert(0,'method',name);t.insert(1,'rank',np.arange(1,6));tops.append(t)
 comparisons=[]
 for j,name in enumerate(METHODS[1:]):
  v=diff[:,j];v=v[np.isfinite(v)];comparisons.append(dict(comparator=name,delta_rho=float(rho[0]-rho[j+1]),paired_site_bootstrap95=np.quantile(v,[.025,.975]).tolist(),paired_site_bootstrap98_75=np.quantile(v,[.00625,.99375]).tolist(),undefined_bootstrap_draws=int(20000-len(v)),fraction_draws_above_zero=float(np.mean(v>0))))
 pd.concat(tops).to_csv(P/'results/top-five.csv',index=False)
 # Propagate reported marginal fit spreads; no independent-replicate confidence claim.
 eligiblevariants=good[good.canonical.isin(e.canonical)].copy();b=eligiblevariants.b_ddg_pred.to_numpy();bs=eligiblevariants.b_ddg_pred_sd.to_numpy();f=eligiblevariants.f_ddg_pred.to_numpy();fs=eligiblevariants.f_ddg_pred_sd.to_numpy();rgen=np.random.default_rng(2026091401);bdraw=rgen.normal(b,bs,size=(1000,len(b)));fdraw=rgen.normal(f,fs,size=(1000,len(f)));by=np.column_stack([abs(bdraw[:,eligiblevariants.canonical.to_numpy()==c]).mean(axis=1) for c in e.canonical]);fy=np.column_stack([abs(fdraw[:,eligiblevariants.canonical.to_numpy()==c]).mean(axis=1) for c in e.canonical]);mrho=np.stack([corr(np.broadcast_to(x,by.shape),by) for x in X],axis=1);np.savez_compressed(P/'results/marginal-fit-spread-sensitivity.npz',canonical=e.canonical.to_numpy(),binding_magnitude=by,folding_magnitude=fy,rho=mrho)
 measurement=[dict(method=name,rho_marginal_spread95=np.quantile(mrho[:,j],[.025,.975]).tolist(),scope='Marginal Normal independence sensitivity only; unavailable joint fit covariance prevents calibrated joint confidence.') for j,name in enumerate(METHODS)]
 # Use primary site mask for all predeclared alternate structures.
 sens=pd.read_csv(P/'results/structure-sensitivity-predictions.csv');srows=[]
 for model,g in sens.groupby('model'):
  g=e[['canonical','binding_magnitude']].merge(g,on='canonical',validate='one_to_one');srows.append(dict(model=int(model),gaussian_rho=float(corr(g.gaussian,g.binding_magnitude)),harmonic_rho=float(corr(g.harmonic,g.binding_magnitude)),gaussian_rank_correlation_with_primary=float(corr(g.gaussian,e.gaussian)),primary_sites_also_distal_in_own_geometry=int(g.candidate_in_own_geometry.sum()),n=n))
 pd.DataFrame(srows).to_csv(P/'results/structure-sensitivity-metrics.csv',index=False)
 signal=bool(n>=12 and rho[0]>0 and all(c['paired_site_bootstrap98_75'][0]>0 and c['delta_rho']>=.10 for c in comparisons))
 out={'status':'COMPLETED_PRESPECIFIED_RETROSPECTIVE_GRB2_EVALUATION','created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'protocol_sha256':freeze['protocol_sha256'],'prediction_sha256':receipt['prediction_sha256'],'source_Table7_sha256':sha(P/'local-only/supplementary-table-7.xlsx'),'raw_single_substitution_rows':len(df),'valid_common_variant_rows':len(good),'variant_exclusion_counts':df.loc[~df.valid,'exclusion'].value_counts().to_dict(),'common_mask_n':n,'common_mask_canonical':e.canonical.astype(int).tolist(),'structural_candidates':int(d.structural_candidate.sum()),'site_exclusion_counts':d.loc[~d.eligible,'exclusion'].value_counts().to_dict(),'min10_coverage_excluded_structural_sites':d.loc[d.structural_candidate&~d.eligible,'canonical'].astype(int).tolist(),'metrics':metrics,'gaussian_minus_controls':comparisons,'strong_within_assay_incremental_signal_gate':signal,'minimum12_sites_for_inference':n>=12,'ohm_seed_rho':{str(seed):float(corr(e['ohm_seed'+str(seed)],y)) for seed in [11,29,47]},'marginal_fit_spread_sensitivity':measurement,'structure_sensitivity':srows,'uncertainty_scope':'Paired resampling of sites conditional on observed aggregates; mutations not pseudoreplicates. Spatial dependence and one small domain limit inferential calibration and generalization. Fit spreads are not independent experimental replicates.','biology_scope':'One inferred GAB2-binding endpoint in GRB2 C-SH3; no supported negative pockets, no biological-family superiority, no original nonlinear fidelity or quantum claim.','evaluation_seconds':time.monotonic()-start,'peak_RSS_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024)}
 (P/'results/summary.json').write_text(json.dumps(out,indent=2));print(json.dumps({k:v for k,v in out.items() if k not in ['structure_sensitivity','metrics','marginal_fit_spread_sensitivity']},indent=2));print(pd.DataFrame(metrics)[['method','rho','top5_canonical','folding_partial_rho']].to_string(index=False))
if __name__=='__main__':main()
