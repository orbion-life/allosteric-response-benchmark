#!/usr/bin/env python3
"""Read designated validation ligands only after input-freeze.json exists.
Writes complete contact mappings and exclusions; does not alter inputs or rank scores.
"""
from pathlib import Path
import json,hashlib
import numpy as np
from Bio.PDB import MMCIFParser
from build_validation_inputs import map_residues,heavy
P=Path(__file__).resolve().parent;O=P/'evaluation';O.mkdir(exist_ok=True)
freeze=P/'prepared/input-freeze.json';assert freeze.exists()
for target,pdb,ligand in [('KRAS','6OIM','MOV'),('ABL','5MO4','AY7')]:
 inp=json.loads((P/f'prepared/{target}-input.json').read_text());domain={r['canonical']:r for r in inp['residues']if r['auth_chain']==inp['configuration']['rank_chain']}
 model=MMCIFParser(QUIET=True).get_structure(pdb,P/f'sources/{pdb}.cif')[0];mapping=map_residues(pdb,model)
 ligand_res=[r for r in model['A']if r.resname==ligand];assert len(ligand_res)==1
 atoms=heavy(ligand_res[0]);rows=[]
 for r in model['A']:
  if r.id[0]!=' 'or 'CA'not in r:continue
  d=float(min(np.linalg.norm(a.coord-b.coord)for a in heavy(r)for b in atoms));m=mapping.get(('A',r.id[1],r.id[2].strip()))
  if d<=5:
   n=m['canonical']if m else None;x=domain.get(n)
   reason='eligible'if x and x['eligible_6A']else('not mapped'if not m else 'outside modeled domain'if not x else 'functional receiver'if x['receiver']else 'inside distal exclusion')
   rows.append({'auth_chain':'A','auth_seq':r.id[1],'resname':r.resname,'canonical':n,'distance_A':d,'eligible':bool(x and x['eligible_6A']),'exclusion_reason':reason,'reference_sequence_match':m['observed_aa']==m['expected_aa']if m else None})
 out={'target':target,'reference':pdb,'ligand':ligand,'ligand_auth_seq':ligand_res[0].id[1],'threshold_A':5.0,'input_freeze_sha256':hashlib.sha256(freeze.read_bytes()).hexdigest(),'reference_sha256':hashlib.sha256((P/f'sources/{pdb}.cif').read_bytes()).hexdigest(),'contacts':rows,'eligible_canonical':[r['canonical']for r in rows if r['eligible']],'status':'retrospective designated-ligand labels; not an independent protein cohort'}
 (O/f'{target}-contacts.json').write_text(json.dumps(out,indent=2)+'\n');print(target,'all',len(rows),'eligible',out['eligible_canonical'],'excluded',[(r['canonical'],r['exclusion_reason'])for r in rows if not r['eligible']])
for target,reason in [('MYH7','6C1H is rat myosin-Ib, not cardiac MYH7; no eligible designated ligand truth'),('MYC','No designated ligand-contact truth; consensus/docking are descriptive')]:
 (O/f'{target}-contacts.json').write_text(json.dumps({'target':target,'status':'unresolved','reason':reason,'eligible_canonical':[]},indent=2)+'\n')
