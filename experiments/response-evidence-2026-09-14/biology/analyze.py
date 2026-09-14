"""Exploratory incremental signal audit; frozen protocol and immutable assay inputs."""
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:
    os.environ[key]='1'
from pathlib import Path
import hashlib,json,datetime,time,sys,platform
import numpy as np
import pandas as pd
import scipy
from scipy.stats import rankdata,spearmanr

ROOT=Path(__file__).resolve().parent
NBOOT=20000
PARTIAL_NAMES=['gaussian_degree_distance','harmonic_degree_distance','ohm_degree_distance','gaussian_degree_distance_harmonic','gaussian_degree_distance_folding','harmonic_degree_distance_folding','gaussian_minus_harmonic_degree_distance']
MODEL_FEATURES={'geometry':['degree','negative_distance'],'geometry_gaussian':['degree','negative_distance','gaussian'],'geometry_harmonic':['degree','negative_distance','harmonic'],'geometry_harmonic_gaussian':['degree','negative_distance','harmonic','gaussian']}
METHODS=['gaussian','harmonic','ohm','degree','negative_distance']

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def clean(o):
    if isinstance(o,dict):return {k:clean(v) for k,v in o.items()}
    if isinstance(o,list):return [clean(x) for x in o]
    if isinstance(o,np.ndarray):return clean(o.tolist())
    if isinstance(o,(np.integer,np.floating)):return clean(o.item())
    if isinstance(o,float) and not np.isfinite(o):return None
    return o

def corr(a,b):
    a=a-a.mean(axis=0);b=b-b.mean(axis=0)
    den=np.sqrt(np.sum(a*a,axis=0)*np.sum(b*b,axis=0))
    return np.divide(np.sum(a*b,axis=0),den,out=np.full(np.shape(den),np.nan),where=den>1e-10)

def partial_stats(a):
    # Columns: Gaussian, harmonic, Ohm, degree, negative distance, folding, outcome.
    r=rankdata(a,axis=0,method='average');n=len(r)
    Z=np.column_stack([np.ones(n),r[:,3:5]])
    values=r[:,[0,1,2,6]];res=values-Z@np.linalg.lstsq(Z,values,rcond=None)[0]
    base=corr(res[:,:3],res[:,[-1]])
    ZH=np.column_stack([Z,r[:,1]])
    gh=r[:,[0,6]];gh=gh-ZH@np.linalg.lstsq(ZH,gh,rcond=None)[0]
    ZF=np.column_stack([Z,r[:,5]])
    gf=r[:,[0,1,6]];gf=gf-ZF@np.linalg.lstsq(ZF,gf,rcond=None)[0]
    folding=corr(gf[:,:2],gf[:,[-1]])
    return np.array([*base,float(corr(gh[:,0],gh[:,1])),*folding,base[0]-base[1]])

def ecdf(train,query):
    s=np.sort(train)
    return (np.searchsorted(s,query,side='left')+np.searchsorted(s,query,side='right'))/(2*len(s))

def oof(e,all_positions):
    folds={int(c):i for i,part in enumerate(np.array_split(all_positions,5)) for c in part}
    e=e.copy();e['fold']=[folds[int(c)] for c in e.canonical]
    records=[];n=len(e)
    for name in MODEL_FEATURES:e[name+'_prediction']=np.nan
    e['test_outcome_train_percentile']=np.nan
    for fold in range(5):
        train=np.flatnonzero(e.fold.to_numpy()!=fold);test=np.flatnonzero(e.fold.to_numpy()==fold)
        assert len(train)>=5 and len(test)>0
        ytrain=ecdf(e.outcome.to_numpy()[train],e.outcome.to_numpy()[train]);ytest=ecdf(e.outcome.to_numpy()[train],e.outcome.to_numpy()[test]);e.loc[test,'test_outcome_train_percentile']=ytest
        foldrec={'fold':fold,'train_canonical':e.iloc[train].canonical.astype(int).tolist(),'test_canonical':e.iloc[test].canonical.astype(int).tolist(),'models':{}}
        for name,features in MODEL_FEATURES.items():
            xtrain=np.column_stack([ecdf(e[f].to_numpy()[train],e[f].to_numpy()[train]) for f in features]);xtest=np.column_stack([ecdf(e[f].to_numpy()[train],e[f].to_numpy()[test]) for f in features]);mu=xtrain.mean(axis=0);sd=xtrain.std(axis=0);sd=np.where(sd>1e-12,sd,1)
            z=(xtrain-mu)/sd;zt=(xtest-mu)/sd;intercept=ytrain.mean();coef=np.linalg.solve(z.T@z+np.eye(len(features)),z.T@(ytrain-intercept));pred=intercept+zt@coef;e.loc[test,name+'_prediction']=pred
            foldrec['models'][name]={'features':features,'alpha':1.,'train_means':mu.tolist(),'train_sd':sd.tolist(),'intercept':float(intercept),'coef':coef.tolist(),'test_mse':float(np.mean((ytest-pred)**2))}
        records.append(foldrec)
    metrics={}
    for name in MODEL_FEATURES:
        pred=e[name+'_prediction'].to_numpy();assert np.isfinite(pred).all();e[name+'_squared_error']=(e.test_outcome_train_percentile-pred)**2
        metrics[name]={'oof_spearman':float(spearmanr(pred,e.outcome).statistic),'oof_percentile_mse':float(e[name+'_squared_error'].mean()),'fold_mse':[r['models'][name]['test_mse'] for r in records]}
    return e,records,metrics

def boot_summary(a,names):
    d={}
    for j,name in enumerate(names):
        v=a[:,j];v=v[np.isfinite(v)];d[name]={'valid_draws':len(v),'undefined_draws':NBOOT-len(v),'percentile95':np.quantile(v,[.025,.975]).tolist() if len(v) else None,'percentile97_5':np.quantile(v,[.0125,.9875]).tolist() if len(v) else None}
    return d

def main(outdir):
    started=time.monotonic();outdir.mkdir(parents=True,exist_ok=True)
    freeze=json.loads((ROOT/'protocol-freeze.json').read_text());assert sha(ROOT/'PROTOCOL.md')==freeze['protocol_sha256']
    manifest=json.loads((ROOT/'input-manifest.json').read_text());out={'status':'COMPLETED_EXPLORATORY_REANALYSIS_OF_PREVIOUSLY_OPENED_ASSAYS','created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'protocol_sha256':freeze['protocol_sha256'],'script_sha256':sha(__file__),'assays':{}}
    for target in ['KRAS','GRB2']:
        path=ROOT/manifest[target]['copied_path'];assert sha(path)==manifest[target]['sha256'];d=pd.read_csv(path)
        if target=='KRAS':
            e=d[d.primary_eligible].copy().reset_index(drop=True).rename(columns={'binding_weighted_abs_ddG':'outcome','folding_weighted_abs_ddG':'folding','negative_receiver_distance':'negative_distance'});seed=2026091403;expected=[.29854157744065996,.21945765615490387,-.007028742210011813,.3481149590917138,.23027748348849267]
            assert len(e)==110 and ((e.binding_count>=10)&e.sequence_match&~e.receiver&(e.receiver_CA_distance_A>=6)).all()
            endpoint='Inverse-variance-weighted mean absolute KRAS RAF1-binding ddG, kcal/mol'
        else:
            e=d[d.eligible].copy().reset_index(drop=True).rename(columns={'binding_magnitude':'outcome','folding_magnitude':'folding'});seed=2026091404;expected=[.38051948051948054,.33506493506493507,-.23766233766233766,.480627267493178,.23506493506493506]
            assert len(e)==21 and ((e.mutation_count>=10)&~e.receiver&(e.receiver_CA_distance_A>=6)).all()
            endpoint='Equal-weight mean absolute GRB2 C-SH3 GAB2-binding ddG, kcal/mol'
        assert e.canonical.is_unique and np.isfinite(e[['outcome','folding']+METHODS]).all().all()
        observed=[float(spearmanr(e[m],e.outcome).statistic) for m in METHODS];assert np.allclose(observed,expected,atol=1e-14,rtol=0)
        positions=np.sort(d.canonical.astype(int).to_numpy());e,folds,predmetrics=oof(e,positions);e['sequence_block10']=((e.canonical-positions[0])//10).astype(int)
        a=e[['gaussian','harmonic','ohm','degree','negative_distance','folding','outcome']].to_numpy();point=partial_stats(a);blocks=[np.flatnonzero(e.sequence_block10.to_numpy()==g) for g in sorted(e.sequence_block10.unique())];rng=np.random.default_rng(seed)
        base=e.geometry_squared_error.to_numpy();gauss=e.geometry_gaussian_squared_error.to_numpy();harm=e.geometry_harmonic_squared_error.to_numpy();both=e.geometry_harmonic_gaussian_squared_error.to_numpy()
        bootout={}
        for kind in ['site','sequence_block10']:
            stats=np.full((NBOOT,len(PARTIAL_NAMES)),np.nan);losses=np.full((NBOOT,4),np.nan);lens=np.zeros(NBOOT,dtype=np.int16);inds=np.full((NBOOT,(len(e) if kind=='site' else max(len(b) for b in blocks)*len(blocks))),-1,dtype=np.int16)
            for b in range(NBOOT):
                idx=rng.integers(0,len(e),size=len(e)) if kind=='site' else np.concatenate([blocks[i] for i in rng.integers(0,len(blocks),size=len(blocks))]);lens[b]=len(idx);inds[b,:len(idx)]=idx
                stats[b]=partial_stats(a[idx]);m0=base[idx].mean();mg=gauss[idx].mean();mh=harm[idx].mean();mb=both[idx].mean();losses[b]=[m0-mg,(m0-mg)/m0,mh-mb,(mh-mb)/mh]
            np.savez_compressed(outdir/f'{target}-{kind}-bootstrap.npz',partial_names=np.array(PARTIAL_NAMES),partial=stats,loss_names=np.array(['baseline_minus_plus_gaussian_mse','relative_mse_reduction_gaussian','harmonic_minus_plus_gaussian_mse','relative_mse_reduction_after_harmonic']),loss_contrasts=losses,site_indices=inds,sample_lengths=lens,canonical=e.canonical.to_numpy())
            bootout[kind]={'partial':boot_summary(stats,PARTIAL_NAMES),'fixed_oof_loss':boot_summary(losses,['baseline_minus_plus_gaussian_mse','relative_mse_reduction_gaussian','harmonic_minus_plus_gaussian_mse','relative_mse_reduction_after_harmonic'])}
        shortlist=[]
        for method in METHODS:
            top=e.sort_values([method,'canonical'],ascending=[False,True]).head(5);shortlist.append({'method':method,'canonical':top.canonical.astype(int).tolist(),'mean_binding_magnitude_kcal_per_mol':float(top.outcome.mean()),'mean_folding_magnitude_kcal_per_mol':float(top.folding.mean())})
        e.to_csv(outdir/f'{target}-per-site.csv',index=False);(outdir/f'{target}-fold-models.json').write_text(json.dumps(clean(folds),indent=2)+'\n')
        mp=predmetrics;delta=mp['geometry']['oof_percentile_mse']-mp['geometry_gaussian']['oof_percentile_mse'];dh=mp['geometry_harmonic']['oof_percentile_mse']-mp['geometry_harmonic_gaussian']['oof_percentile_mse']
        fold_delta=np.array(mp['geometry']['fold_mse'])-np.array(mp['geometry_gaussian']['fold_mse'])
        r={'n':len(e),'endpoint':endpoint,'input_sha256':manifest[target]['sha256'],'mask_canonical':e.canonical.astype(int).tolist(),'prior_marginal_rho_reproduced':dict(zip(METHODS,observed)),'partial_spearman':dict(zip(PARTIAL_NAMES,point)),'fold_test_n':[len(f['test_canonical']) for f in folds],'sequence_block10_sizes':[len(b) for b in blocks],'oof_prediction':mp,'oof_gaussian_addition':{'mse_difference_positive_favours_gaussian':delta,'relative_mse_reduction':delta/mp['geometry']['oof_percentile_mse'],'delta_spearman':mp['geometry_gaussian']['oof_spearman']-mp['geometry']['oof_spearman'],'fold_mse_difference':fold_delta.tolist(),'folds_improved':int((fold_delta>0).sum())},'oof_gaussian_after_harmonic':{'mse_difference_positive_favours_gaussian':dh,'relative_mse_reduction':dh/mp['geometry_harmonic']['oof_percentile_mse']},'bootstrap':bootout,'shortlists':shortlist}
        out['assays'][target]=r
        print(json.dumps(clean({'target':target,'partial':r['partial_spearman'],'oof':r['oof_gaussian_addition'],'blocks':r['sequence_block10_sizes']})),flush=True)
    out['runtime_seconds']=time.monotonic()-started
    (outdir/'summary.json').write_text(json.dumps(clean(out),indent=2,allow_nan=False)+'\n')
    (outdir/'environment.json').write_text(json.dumps({'python':sys.version,'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__,'platform':platform.platform(),'thread_limits':{k:os.environ[k] for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']}},indent=2)+'\n')

if __name__=='__main__':
    main(Path(sys.argv[1]) if len(sys.argv)>1 else ROOT/'results')
