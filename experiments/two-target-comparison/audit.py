#!/usr/bin/env python3
"""Reproducible retrospective multiplicity audit; Python standard library only.

Read the corrected unknown-aware inputs. Preserve planned but unrun MYH7 slots
at p=1 without inventing a measured effect or a completed experiment.
"""
from pathlib import Path
import argparse,csv,hashlib,json,math
P=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def holm(pvalues):
    if not all(math.isfinite(p) and 0<=p<=1 for p in pvalues):raise ValueError('Finite p values in [0,1] are required.')
    order=sorted(range(len(pvalues)),key=lambda i:(pvalues[i],i));out=[None]*len(order);previous=0.
    for rank,i in enumerate(order):
        previous=max(previous,(len(order)-rank)*pvalues[i]);out[i]=min(1.,previous)
    return out
def select(records,field,value):
    chosen=[r for r in records if r['null_family']=='degree_distance_exposure' and r['universe']=='evaluable123' and r[field]==value]
    if len(chosen)!=1:raise ValueError(f'Expected exactly one corrected record for {value!r}; found {len(chosen)}.')
    return chosen[0]
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,default=P/'results');args=parser.parse_args();out=args.output;out.mkdir(parents=True,exist_ok=True)
    receipts=json.loads((P/'input-receipts.json').read_text())
    for r in receipts:
        if sha(P/r['input'])!=r['sha256']:raise ValueError('Input hash mismatch: '+r['input'])
    kn=json.loads((P/'inputs/kras-null-comparisons.json').read_text());kp=json.loads((P/'inputs/kras-paired-comparisons.json').read_text());a=json.loads((P/'inputs/abl-evaluation.json').read_text())
    assert a['prediction_candidates']==193 and a['evaluable_candidates']==175 and a['matched_null']['unknowns_in_pools'] is False
    k=select(kn,'method','Biquadratic d2 grid');ab=a['methods']['biquadratic_d2_n33']['score_enrichment_null']
    primary=[{'target':'KRAS','status':'retrospective_measured','raw_p':k['raw_p'],'effect_kind':'observed_minus_matched_null_mean_score','effect':k['observed_minus_null_mean'],'source_record':'degree_distance_exposure / evaluable123 / Biquadratic d2 grid'},
             {'target':'ABL','status':'retrospective_measured','raw_p':ab['p'],'effect_kind':'observed_minus_matched_null_mean_score','effect':ab['observed_minus_null_mean'],'source_record':'methods.biquadratic_d2_n33.score_enrichment_null'},
             {'target':'MYH7','status':'not_run_conservative_correction_slot','raw_p':1.,'effect_kind':'not_measured','effect':None,'source_record':None}]
    paired=[]
    for target in ['KRAS','ABL','MYH7']:
        for comparator,kname,aname in [('harmonic','Harmonic d2 same grid','harmonic_d2_n33'),('distance_Hookean','Distance-Hookean d2 same grid','Hookean_d2_n33')]:
            if target=='KRAS':
                record=select(kp,'comparison','Biquadratic minus '+kname);raw=record['raw_p'];effect=record['observed_mean'];centered=record['observed_minus_null_mean'];status='retrospective_measured';denominator=record['percentile_rank_denominator']
            elif target=='ABL':
                record=a['paired_primary_minus_control_percentile_effects'][aname];raw=record['p'];effect=record['observed_contact_mean'];centered=record['observed_minus_null_mean'];status='retrospective_measured';denominator=193
            else:raw=1.;effect=None;centered=None;status='not_run_conservative_correction_slot';denominator=None
            paired.append({'target':target,'comparator':comparator,'status':status,'raw_p':raw,'mean_primary_minus_control_percentile':effect,'mean_percentile_difference_in_percentage_points':None if effect is None else 100*effect,'effect_minus_matched_null_mean':centered,'percentile_rank_denominator':denominator})
    for family in [primary,paired]:
        for row,adjusted in zip(family,holm([r['raw_p'] for r in family])):
            row['holm_adjusted_p']=adjusted;row['below_0p05_in_this_retrospective_family']=row['status']=='retrospective_measured' and adjusted<=.05
    result={'status':'retrospective_statistical_audit','alpha':.05,'primary_family_size':3,'primary_hypothesis':'Mean primary score at known reference contacts exceeds its degree/distance/exposure-matched label null.','primary':primary,'paired_family_size':6,'paired_hypothesis':'Mean primary-minus-control percentile rank at known contacts exceeds the corresponding matched-label null.','paired':paired,'unrun_slot_policy':'MYH7 uses p=1 only for multiplicity correction; no effect or experimental result is imputed.','unknown_policy':'Input prediction ranks remain frozen; unknown reference residues are excluded from evaluation and null pools. KRAS ranks126/evaluable123; ABL ranks193/evaluable175.','interpretation':'KRAS retains conditional primary enrichment, whereas ABL fails transfer. No paired comparison passes Holm correction at0.05. This does not demonstrate nonlinear, finite-time, functional or quantum advantage.','statistical_limits':'Conditional label-exchangeability assumptions, retrospective known references and no spatial pocket-contiguity matching limit biological inference. These adjusted p values are not prospective performance estimates or a correction across every exploratory analysis.','source_receipts_sha256':sha(P/'input-receipts.json'),'code_sha256':sha(__file__)}
    (out/'audit.json').write_text(json.dumps(result,indent=2)+'\n')
    for name,family in [('primary',primary),('paired',paired)]:
        with(out/(name+'.csv')).open('w')as f:w=csv.DictWriter(f,fieldnames=list(family[0]));w.writeheader();w.writerows(family)
    print(json.dumps({'primary':[(r['target'],r['raw_p'],r['holm_adjusted_p']) for r in primary],'paired':[(r['target'],r['comparator'],r['mean_primary_minus_control_percentile'],r['holm_adjusted_p']) for r in paired]},indent=2))
if __name__=='__main__':main()
