"""Post-freeze retrospective validation. No parameters are selected here."""
from pathlib import Path
import json,csv,datetime
import numpy as np
from scipy.stats import spearmanr,hypergeom
from prepare import ROOT,atoms,sha
from run_pilot import ranking

def main():
    freeze=ROOT/'prediction-freeze.json';assert freeze.exists(),'Predictions must be frozen before MOV reference inspection.'
    for rel,h in json.loads(freeze.read_text())['files'].items():assert sha(ROOT/rel)==h,rel
    model=np.load(ROOT/'model/static-model.npz');ids=model['canonical'];rec=model['receiver'];candidate=model['candidate'];idx=np.flatnonzero(candidate);out=ROOT/'evaluation';out.mkdir(exist_ok=True)
    mapping=json.loads((ROOT/'raw/6OIM-mapping.json').read_text())['6oim']['UniProt']['P01116']['mappings'];mapping=[m for m in mapping if m['chain_id']=='A'];assert len(mapping)==1
    mp=mapping[0];assert mp['unp_start']==mp['start']['author_residue_number']==1
    aa=atoms(ROOT/'raw/6OIM.pdb');lig=np.array([a['xyz'] for a in aa if a['resname']=='MOV' and a['record']=='HETATM']);assert len(lig)>0
    input_aa=atoms(ROOT/'raw/4OBE.pdb');input_name={a['auth_seq_id']:a['resname'] for a in input_aa if a['record']=='ATOM' and a['atom']=='CA'};ref_name={a['auth_seq_id']:a['resname'] for a in aa if a['record']=='ATOM' and a['atom']=='CA'}
    distance=np.full(len(ids),np.nan);mutations=[]
    for i,resid in enumerate(ids):
        heavy=np.array([a['xyz'] for a in aa if a['record']=='ATOM' and a['auth_seq_id']==resid])
        if len(heavy):distance[i]=np.sqrt(np.sum((heavy[:,None]-lig[None])**2,axis=-1)).min()
        if resid in ref_name and ref_name[resid]!=input_name[resid]:mutations.append({'canonical_residue':int(resid),'input':input_name[resid],'reference':ref_name[resid]})
    labels=distance<=5;missing=ids[~np.isfinite(distance)].tolist();eligible_labels=labels&candidate
    primary=np.load(ROOT/'results/biquadratic-d2-n33-e4.npz');sd=np.load(ROOT/'model/all-mode-harmonic-scales.npz')['harmonic_sd'];full=np.load(ROOT/'results/analytic-harmonic-d492.npz')
    methods={'Biquadratic d2 grid':primary['score'],'Harmonic d2 same grid':np.load(ROOT/'results/harmonic-d2-n33-e4.npz')['score'],'Distance-Hookean d2 same grid':np.load(ROOT/'results/distance_hookean-d2-n33-e4.npz')['score'],'Harmonic all 492 modes':full['score'],'Equilibrium d2 biquadratic':ranking(primary['equilibrium'],rec,candidate,ids)[0],'Equilibrium all-mode harmonic':ranking(full['equilibrium'],rec,candidate,ids)[0]}
    structural=np.load(ROOT/'results/structural-baselines.npz')
    methods.update({'Graph heat kernel':structural['graph_heat_kernel'],'Contact degree':structural['degree'],'Receiver inverse distance':structural['receiver_inverse_distance']})
    # An Ohm companion may be supplied later, always on the identical candidate mask.
    ohm=ROOT/'external/ohm-scores.json'
    if ohm.exists():
        oo=json.loads(ohm.read_text());methods['Ohm pinned upstream']=np.array([oo['score_by_canonical'][str(int(r))] for r in ids])
    # Predeclared matched random sets preserve degree/receiver-distance joint tertile counts.
    cuts_degree=np.quantile(model['degree'][idx],[1/3,2/3]);cuts_distance=np.quantile(model['receiver_distance'][idx],[1/3,2/3]);cell=3*np.searchsorted(cuts_degree,model['degree'],side='right')+np.searchsorted(cuts_distance,model['receiver_distance'],side='right')
    rng=np.random.default_rng(20260912);draws=[]
    for _ in range(10000):
        chosen=[]
        for c in range(9):
            count=int(np.sum(eligible_labels&(cell==c)));pool=np.flatnonzero(candidate&(cell==c))
            if count:chosen.extend(rng.choice(pool,count,replace=False).tolist())
        draws.append(chosen)
    draws=np.array(draws,dtype=int);rows=[];fullrows=[]
    for name,score in methods.items():
        order=idx[np.lexsort((ids[idx],-score[idx]))];top=order[:5];hits=int(labels[top].sum());null=score[draws].mean(axis=1);actual=float(score[eligible_labels].mean())
        rows.append({'method':name,'top5':ids[top].tolist(),'reference_heavy_atom_distances_A':distance[top].tolist(),'hits_at_5A':hits,'precision_at_5':hits/5,'eligible_known_site_mean_score':actual,'matched_random_mean_score':float(null.mean()),'matched_random_2.5pct':float(np.quantile(null,.025)),'matched_random_97.5pct':float(np.quantile(null,.975)),'matched_random_fraction_at_least_known_site':float((np.sum(null>=actual)+1)/(len(null)+1)),'descriptive_hypergeometric_tail':float(hypergeom.sf(hits-1,len(idx),int(eligible_labels.sum()),5)),'spearman_to_primary':float(spearmanr(score[idx],primary['score'][idx]).statistic)})
        for rank,i in enumerate(order,1):fullrows.append({'method':name,'rank':rank,'canonical_residue':int(ids[i]),'input_resname':input_name[int(ids[i])],'score':float(score[i]),'MOV_reference_min_heavy_atom_distance_A':float(distance[i]),'MOV_contact_5A':bool(labels[i]),'receiver_CA_distance_A':float(model['receiver_distance'][i]),'contact_degree':int(model['degree'][i])})
    with (out/'rankings.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(fullrows[0]));w.writeheader();w.writerows(fullrows)
    sensitivity=[]
    for path in sorted((ROOT/'results').glob('biquadratic-d*.npz')):
        q=np.load(path);meta=json.loads(path.with_suffix('.json').read_text());sensitivity.append({'case':path.stem,'max_abs_C_difference_from_primary':float(np.max(np.abs(q['C']-primary['C']))),'spearman_from_primary':float(spearmanr(q['score'][idx],primary['score'][idx]).statistic),'top5':ids[q['order'][:5]].tolist(),'boundary_probability':meta['boundary_probability']})
    compression=[]
    for d in [1,2,4,8,16,32,64,492]:
        q=np.load(ROOT/f'results/analytic-harmonic-d{d}.npz');block=np.ix_(idx,rec);compression.append({'d':d,'max_abs_C_difference_from_all_mode':float(np.max(np.abs(q['C']-full['C']))),'frobenius_relative_error':float(np.linalg.norm(q['C']-full['C'])/np.linalg.norm(full['C'])),'candidate_to_receiver_max_abs_C_difference':float(np.max(np.abs(q['C'][block]-full['C'][block]))),'candidate_to_receiver_relative_frobenius_error':float(np.linalg.norm(q['C'][block]-full['C'][block])/np.linalg.norm(full['C'][block])),'spearman_to_all_mode':float(spearmanr(q['score'][idx],full['score'][idx]).statistic),'top5':ids[q['order'][:5]].tolist()})
    mask=[]
    for threshold in [4,6,8]:
        can=model['receiver_distance']>=threshold;can[rec]=False;score,order=ranking(primary['C'],rec,can,ids);mask.append({'minimum_CA_distance_A':threshold,'candidate_count':int(can.sum()),'top5':ids[order[:5]].tolist(),'precision_at_5':float(labels[order[:5]].mean())})
    result={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'retrospective engineering pilot; no prospective biological validation','prediction_freeze_sha256':sha(freeze),'evaluation_addendum_sha256':sha(ROOT/'evaluation-protocol-addendum.json'),'reference':'6OIM chain A MOV','mapping':'PDBe SIFTS P01116 author1 maps canonical1; resolved atom names compared at every modeled position','reference_vs_input_mutations':mutations,'missing_reference_residues':missing,'all_reference_contact_residues_5A':ids[labels].tolist(),'eligible_reference_contact_residues':ids[eligible_labels].tolist(),'excluded_reference_contact_residues':ids[labels&~candidate].tolist(),'candidate_count':int(candidate.sum()),'receiver_count':len(rec),'methods':rows,'nonfunctional_pocket_controls':'NOT EXECUTED: no independently established nonfunctional pocket supplied; unannotated regions are not negatives','matched_random_residues':{'draws':10000,'seed':20260912,'degree_tertile_cuts':cuts_degree.tolist(),'receiver_distance_tertile_cuts_A':cuts_distance.tolist(),'interpretation':'descriptive conditional random-set distribution, not biological negative pockets or prospective uncertainty'},'grid_sensitivity':sensitivity,'harmonic_coordinate_compression':compression,'candidate_mask_sensitivity':mask,'claim_limit':'A converged d2 grid remains only a reduced-coordinate calculation. Failure against the exact all-mode harmonic reference rejects d2 as a validated protein compression, irrespective of any pocket hits.'}
    (out/'evaluation.json').write_text(json.dumps(result,indent=2)+'\n');np.savez_compressed(out/'reference-mapping.npz',canonical=ids,contact=labels,distance_A=distance,ligand_coordinates_A=lig,matched_random_nodes=draws)
    print(json.dumps({'candidates':result['candidate_count'],'reference_contacts':result['all_reference_contact_residues_5A'],'methods':[{'method':r['method'],'top5':r['top5'],'P@5':r['precision_at_5']} for r in rows],'compression':compression},indent=2))

if __name__=='__main__':main()
