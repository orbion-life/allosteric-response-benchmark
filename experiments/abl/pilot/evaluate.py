#!/usr/bin/env python3
"""Retrospective AY7 contact evaluation after the prediction freeze exists."""
from pathlib import Path
import copy,csv,hashlib,json,sys
import numpy as np
from scipy.stats import rankdata,spearmanr
from Bio.PDB import MMCIFParser,MMCIF2Dict
from Bio.Data.PDBData import protein_letters_3to1
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'vendor'))
from matched_null import plan,configurations
from run_pilot import ranking
from run import sha,dump

def main():
    out=ROOT/'evaluation';out.mkdir(exist_ok=True)
    freeze=json.loads((ROOT/'prediction-freeze.json').read_text())
    assert all(sha(ROOT/k)==v for k,v in freeze['files'].items())
    model=np.load(ROOT/'model/static-model.npz'); inp=json.loads((ROOT/'raw/ABL-input.json').read_text());rows=inp['residues'];ids=model['canonical'];candidate=model['candidate'];idx=np.flatnonzero(candidate);receiver=model['receiver']
    primary=np.load(ROOT/'results/biquadratic-d2-n33-e4.npz')
    methods={'biquadratic_d2_n33':primary['score']}
    for key,name in [('harmonic_d2_n33','harmonic-d2-n33-e4'),('Hookean_d2_n33','distance_hookean-d2-n33-e4'),('analytic_harmonic_d2','analytic-harmonic-d2'),('all_mode_harmonic','analytic-harmonic-d750'),('biquadratic_d2_n65','biquadratic-d2-n65-e4')]:methods[key]=np.load(ROOT/f'results/{name}.npz')['score']
    methods['equilibrium_d2_n33']=ranking(primary['equilibrium'],receiver,candidate,ids)[0]
    full=np.load(ROOT/'results/analytic-harmonic-d750.npz');methods['all_mode_equilibrium']=ranking(full['equilibrium'],receiver,candidate,ids)[0]
    structural=np.load(ROOT/'results/structural-baselines.npz')
    for key in ['graph_heat_kernel','degree','receiver_inverse_distance']:methods[key]=structural[key]
    # Primary source mapping uses label-sequence numbers, not an assumed author offset.
    cif=MMCIF2Dict.MMCIF2Dict(ROOT/'raw/5MO4.cif');struct=MMCIFParser(QUIET=True).get_structure('5MO4',ROOT/'raw/5MO4.cif')[0]
    fields=['auth_asym_id','auth_seq_id','pdbx_PDB_ins_code','label_seq_id','label_asym_id'];atommap={}
    for ac,an,ins,lab,lc in zip(*(cif['_atom_site.'+k] for k in fields)):
        if lab not in ('.','?'):atommap[(ac,int(an),'' if ins in ('.','?') else ins)]=(int(lab),lc)
    sifts=json.loads((ROOT/'raw/5MO4-sifts.json').read_text())['5mo4']['UniProt']['P00519']['mappings']
    seq=json.loads((ROOT/'raw/P00519.json').read_text())['sequence']['value']
    ligand=[r for r in struct['A'] if r.resname=='AY7']; assert len(ligand)==1
    ligcoords=np.array([a.coord for a in ligand[0] if a.element not in ('H','D')],float)
    assert len(ligcoords)>0
    observed={};variants=[]
    expected_heavy={'ALA':5,'ARG':11,'ASN':8,'ASP':8,'CYS':6,'GLN':9,'GLU':9,'GLY':4,'HIS':10,'ILE':8,'LEU':8,'LYS':9,'MET':8,'PHE':11,'PRO':7,'SER':6,'THR':7,'TRP':14,'TYR':12,'VAL':7}
    for r in struct['A']:
        if r.id[0]!=' ':continue
        key=('A',r.id[1],r.id[2].strip())
        if key not in atommap:continue
        lab,lc=atommap[key]
        match=[m for m in sifts if m['chain_id']=='A' and m['start']['residue_number']<=lab<=m['end']['residue_number']]
        assert len(match)<=1
        if not match:continue
        m=match[0]; pos=m['unp_start']+lab-m['start']['residue_number']
        if pos not in ids:continue
        heavy=[a for a in r if a.element not in ('H','D')]
        if not heavy:continue
        d=float(np.linalg.norm(np.array([a.coord for a in heavy],float)[:,None]-ligcoords[None],axis=2).min())
        ob=protein_letters_3to1.get(r.resname,'X');ex=seq[pos-1]
        record={'canonical':pos,'author':r.id[1],'label_seq':lab,'resname':r.resname,'expected_aa':ex,'observed_aa':ob,'distance_AY7_A':d,'state':'contact' if d<=5 else 'observed_noncontact','resolved_heavy_atom_count':len(heavy),'expected_heavy_atom_count_without_terminal_OXT':expected_heavy.get(r.resname),'has_CA':'CA' in r}
        assert pos not in observed;observed[pos]=record
        if ex!=ob:variants.append(record)
    labels=[]
    for i,row in enumerate(rows):
        pos=int(ids[i]);l=observed.get(pos,{'canonical':pos,'distance_AY7_A':None,'state':'unknown'})
        labels.append({**l,'prediction_candidate':bool(candidate[i]),'evaluation_eligible':bool(candidate[i]) and pos in observed})
    labelmap={r['canonical']:r for r in labels};positive=[r['canonical'] for r in labels if r['evaluation_eligible'] and r['state']=='contact']
    evalrows=copy.deepcopy(rows)
    for r in evalrows:r['eligible_6A']=r['eligible_6A'] and r['canonical'] in observed
    pl=plan(evalrows,set(positive))
    for unused in ['first_holm_threshold','sample_space_large_enough_for_first_holm_threshold']:pl.pop(unused,None)
    configs=list(configurations(pl,positive))
    # Same random configurations are used for every score and every paired statistic.
    mapidx={int(p):i for i,p in enumerate(ids)}; configidx=np.array([[mapidx[p] for p in c] for c in configs]);posidx=np.array([mapidx[p] for p in positive])
    def nulltest(score):
        observed_mean=float(np.mean(score[posidx]));null=score[configidx].mean(axis=1);extreme=int(np.sum(null>=observed_mean));exact=pl['M']<=100000
        return {'observed_contact_mean':observed_mean,'null_mean':float(null.mean()),'observed_minus_null_mean':observed_mean-float(null.mean()),'p':extreme/len(null) if exact else (extreme+1)/(len(null)+1),'exact':exact,'null_q025_median_q975':np.quantile(null,[.025,.5,.975]).tolist()}
    records={};percentiles={}
    for name,score in methods.items():
        order=idx[np.lexsort((ids[idx],-score[idx]))];top=[labelmap[int(p)] for p in ids[order[:5]]]
        hits=sum(t['state']=='contact' for t in top);unknown=sum(t['state']=='unknown' for t in top)
        pscore=np.full(len(ids),np.nan);pscore[idx]=rankdata(score[idx],method='average')/len(idx);percentiles[name]=pscore
        records[name]={'top5':top,'known_hits':hits,'unknown_top5':unknown,'P_at5_lower_bound':hits/5,'P_at5_upper_bound':(hits+unknown)/5,'precision_among_observed_top5':hits/(5-unknown) if unknown<5 else None,'mean_contact_percentile':float(np.mean(pscore[posidx])),'score_enrichment_null':nulltest(score),'percentile_enrichment_null':nulltest(pscore),'spearman_vs_primary_frozen_candidates':float(spearmanr(score[idx],methods['biquadratic_d2_n33'][idx]).statistic)}
    paired={name:nulltest(percentiles['biquadratic_d2_n33']-pscore) for name,pscore in percentiles.items() if name not in ['biquadratic_d2_n33','biquadratic_d2_n65']}
    # Sensitivity retains original prediction-universe cutpoints but removes unknowns BEFORE pools.
    original=plan(rows,set(positive));origgroups=[]
    from math import comb
    M=1
    for g in original['groups']:
        members=[p for p in g['candidates'] if p in observed];m=sum(p in positive for p in members);M*=comb(len(members),m);origgroups.append({**g,'candidates':members,'n':len(members),'m':m})
    alt={**original,'groups':origgroups,'M':M,'candidate_count':sum(len(g['candidates']) for g in origgroups)}
    altconfigs=list(configurations(alt,positive));aidx=np.array([[mapidx[p] for p in c] for c in altconfigs]);score=methods['biquadratic_d2_n33'];null=score[aidx].mean(axis=1);mean=float(score[posidx].mean());ext=int(np.sum(null>=mean));altresult={'cutpoints':alt['cuts'],'M':M,'draws':len(null),'p':ext/len(null) if M<=100000 else (1+ext)/(1+len(null)),'unknowns_in_pools':False}
    masks={}
    for cutoff in [4,8]:
        cand=np.array([r[f'eligible_{cutoff}A'] for r in rows]);ii=np.flatnonzero(cand);order=ii[np.lexsort((ids[ii],-methods['biquadratic_d2_n33'][ii]))]
        masks[str(cutoff)]={'candidate_count':int(cand.sum()),'top5':[labelmap[int(p)] for p in ids[order[:5]]],'interpretation':'Geometric sensitivity only; primary universe and inference remain unchanged.'}
    matrixrows=[]
    for i,r in enumerate(rows):
        matrixrows.append({'canonical':r['canonical'],'author':r['auth_seq'],'receiver':r['receiver'],'prediction_candidate':bool(candidate[i]),'evaluation_eligible':labelmap[int(ids[i])]['evaluation_eligible'],'label_state':labelmap[int(ids[i])]['state'],'reference_distance_AY7_A':labelmap[int(ids[i])]['distance_AY7_A'],**{k:float(v[i]) for k,v in methods.items()}})
    with (out/'all-residue-scores.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=matrixrows[0]);w.writeheader();w.writerows(matrixrows)
    dump(out/'reference-labels.json',{'reference':'5MO4 author A; AY7','ligand_author_number':ligand[0].id[1],'reference_sha256':sha(ROOT/'raw/5MO4.cif'),'mapping_sha256':sha(ROOT/'raw/5MO4-sifts.json'),'variants':variants,'unknown_input_positions':[r['canonical'] for r in labels if r['state']=='unknown'],'unknown_prediction_candidates':[r['canonical'] for r in labels if r['state']=='unknown' and r['prediction_candidate']],'partial_heavy_atom_residues':[r for r in observed.values() if r['expected_heavy_atom_count_without_terminal_OXT'] and r['resolved_heavy_atom_count']<r['expected_heavy_atom_count_without_terminal_OXT']],'labels':labels})
    dump(out/'evaluation.json',{'prediction_freeze_sha256':sha(ROOT/'prediction-freeze.json'),'protocol_sha256':sha(ROOT/'preanalysis-protocol.json'),'prediction_candidates':len(idx),'evaluable_candidates':pl['candidate_count'],'eligible_known_contacts':positive,'methods':records,'paired_primary_minus_control_percentile_effects':paired,'matched_null':{**pl,'draws':len(configs),'seed':20260912,'unknowns_in_pools':False,'assumption':'Conditional exchangeability of contact labels given frozen degree/distance/exposure strata; not an experimental biological null.','multiplicity':'Raw retrospective conditional p values only. No four-target or confirmatory significance claim.'},'original_prediction_cutpoint_sensitivity':altresult,'candidate_distance_sensitivity':masks,'scope':'Known-bound-site contact recovery on one engineered, ligand-conditioned domain; neither functional validation nor independent generalization. Unverified nonfunctional pockets remain unresolved.'})
    np.savez_compressed(out/'null-configurations.npz',canonical_sets=np.array(configs,dtype=int))
    print(json.dumps({'candidate_count':len(idx),'evaluable':pl['candidate_count'],'positive':positive,'unknown':[r['canonical'] for r in labels if r['state']=='unknown'],'top5':[(n,[x['canonical'] for x in r['top5']],r['known_hits']) for n,r in records.items()],'paired':paired},indent=2))
if __name__=='__main__':main()
