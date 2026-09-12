#!/usr/bin/env python3
"""Conditional matched-label randomization, with an explicit finite sample space.
This is a retrospective statistical model, not a biological randomization experiment.
"""
from pathlib import Path
from itertools import combinations,product
from math import comb
import argparse,csv,hashlib,json
import numpy as np
P=Path(__file__).resolve().parent

def plan(rows,positive):
 U=[r for r in rows if r['eligible_6A']];assert set(positive)<={r['canonical']for r in U};cuts={k:sorted(set(float(v)for v in np.quantile([r[k]for r in U],[1/3,2/3])))for k in ['degree','receiver_distance_A']}
 strata={}
 for r in U:
  key=tuple(int(np.searchsorted(cuts[k],r[k],side='right'))for k in ['degree','receiver_distance_A'])+(int(r['rsa']>=.2),)
  strata.setdefault(key,[]).append(r['canonical'])
 groups=[];M=1
 for k,v in sorted(strata.items()):
  m=sum(i in positive for i in v);M*=comb(len(v),m);groups.append({'key':k,'candidates':v,'n':len(v),'m':m})
 return {'M':M,'cuts':cuts,'groups':groups,'candidate_count':len(U),'positive_count':len(positive),'smallest_exact_p':1/M,'first_holm_threshold':.05/3,'sample_space_large_enough_for_first_holm_threshold':M>=60}

def configurations(pl,positive,seed=20260912):
 if pl['M']<=100000:
  for values in product(*(combinations(g['candidates'],g['m'])for g in pl['groups'])):yield tuple(sorted(i for v in values for i in v))
 else:
  rng=np.random.Generator(np.random.PCG64(seed));seen={tuple(sorted(positive))};count=0
  while count<10000:
   x=tuple(sorted(int(i)for g in pl['groups']for i in rng.choice(g['candidates'],g['m'],replace=False)))
   if x in seen:continue
   seen.add(x);count+=1;yield x

def test(rows,positive,scores,seed=20260912):
 positive=tuple(sorted(set(positive)))
 pl=plan(rows,set(positive))
 if not positive:raise ValueError('At least one eligible positive label is required.')
 universe={r['canonical']for r in rows if r['eligible_6A']}
 if not universe<=scores.keys():raise ValueError('Scores must cover every eligible candidate.')
 if not all(np.isfinite(scores[i])for i in universe):raise ValueError('Every eligible score must be finite.')
 obs=float(np.mean([scores[i]for i in positive]));null=np.array([np.mean([scores[i]for i in c])for c in configurations(pl,positive,seed)],float)
 extreme=int(np.sum(null>=obs));exact=pl['M']<=100000
 p=extreme/len(null)if exact else(1+extreme)/(1+len(null))
 return {**pl,'observed':obs,'p':p,'exact':exact,'unique_configurations_evaluated':len(null),'seed':seed,'null_quantiles':dict(zip(['q025','median','q975'],map(float,np.quantile(null,[.025,.5,.975])))),'mean_null':float(null.mean()),'observed_minus_null_mean':obs-float(null.mean())}

def holm(pvalues):
 order=sorted(range(len(pvalues)),key=lambda i:pvalues[i]);out=[None]*len(pvalues);prev=0
 for rank,i in enumerate(order):prev=max(prev,(len(pvalues)-rank)*pvalues[i]);out[i]=min(1,prev)
 return out

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--target',choices=['KRAS','ABL'],required=True);parser.add_argument('--scores');parser.add_argument('--score-column',default='ohm_score');parser.add_argument('--out');args=parser.parse_args()
 inp=json.loads((P/f'prepared/{args.target}-input.json').read_text());lab=json.loads((P/f'evaluation/{args.target}-contacts.json').read_text());positive=lab['eligible_canonical'];pl=plan(inp['residues'],set(positive))
 if args.scores:
  with open(args.scores)as f: rows=list(csv.DictReader(f))
  positions=[int(r['canonical'])for r in rows]
  if len(positions)!=len(set(positions)):raise ValueError('Duplicate canonical positions in score CSV.')
  scores={int(r['canonical']):float(r[args.score_column])for r in rows}
  result=test(inp['residues'],positive,scores)
 else:result=pl
 result['input_sha256']=hashlib.sha256((P/f'prepared/{args.target}-input.json').read_bytes()).hexdigest();result['labels_sha256']=hashlib.sha256((P/f'evaluation/{args.target}-contacts.json').read_bytes()).hexdigest()
 text=json.dumps(result,indent=2)+'\n'
 if args.out:Path(args.out).write_text(text)
 print(text)
if __name__=='__main__':main()
