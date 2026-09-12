#!/usr/bin/env python3
"""Derive versioned residue maps, functional receivers and masks from cached primary inputs.
No response score or allosteric reference is read by this preparation command.
"""
from pathlib import Path
import csv, hashlib, json
import numpy as np
from Bio.PDB import MMCIFParser, MMCIF2Dict, PDBIO, Select, ShrakeRupley
from Bio.Data.PDBData import protein_letters_3to1
ROOT=Path(__file__).resolve().parent
SRC=ROOT/'sources'; OUT=ROOT/'prepared'; OUT.mkdir(exist_ok=True)
CONFIG={
 'KRAS': {'pdb':'4OBE','chains':{'A':('P01116',1,166)},'sequence_accession':'P01116-2','rank_chain':'A','receiver_kind':'ligand','ligand':'GDP','ligand_chain':'A'},
 'ABL': {'pdb':'1OPL','chains':{'A':('P00519',242,493)},'rank_chain':'A','receiver_kind':'ligand','ligand':'P16','ligand_chain':'A'},
 'MYH7': {'pdb':'5TBY','chains':{'A':('P12883',85,778)},'rank_chain':'A','receiver_kind':'canonical','receiver':[178,179,180,181,182,183,184,185]},
 'MYC': {'pdb':'1NKP','chains':{'A':('P01106',368,449),'B':('P61244',23,102)},'rank_chain':'A','receiver_kind':'dna','dna_chains':['F','G']}
}
# Tien et al. 2013 theoretical maxima, DOI10.1371/journal.pone.0080635/Table1.
MAXASA={'ALA':129,'ARG':274,'ASN':195,'ASP':193,'CYS':167,'GLU':223,'GLN':225,'GLY':104,'HIS':224,'ILE':197,'LEU':201,'LYS':236,'MET':224,'PHE':240,'PRO':159,'SER':155,'THR':172,'TRP':285,'TYR':263,'VAL':174}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def heavy(r):return [a for a in r if a.element not in ('H','D')]
def map_residues(pdb, model):
 d=MMCIF2Dict.MMCIF2Dict(SRC/f'{pdb}.cif')
 fields=['auth_asym_id','auth_seq_id','pdbx_PDB_ins_code','label_seq_id','label_asym_id']
 atommap={}
 for row in zip(*(d['_atom_site.'+k] for k in fields)):
  ac,an,ins,lab,lc=row
  if lab not in ('.','?'): atommap[(ac,int(an),''if ins in('.','?')else ins)]=(int(lab),lc)
 sifts=json.loads((SRC/f'{pdb}-sifts.json').read_text())[pdb.lower()]['UniProt']
 maps={}
 for c in model:
  for r in c:
   if r.id[0]!=' 'or 'CA'not in r:continue
   key=(c.id,r.id[1],r.id[2].strip())
   if key not in atommap:continue
   label,lc=atommap[key]
   for up,rec in sifts.items():
    for m in rec['mappings']:
     if m['chain_id']==c.id and m['start']['residue_number']<=label<=m['end']['residue_number']:
      pos=m['unp_start']+label-m['start']['residue_number']
      seq=json.loads((SRC/f'{up}.json').read_text())['sequence']['value']
      sequence_accession=up
      if up=='P01116':
       seq=''.join((SRC/'P01116-2.fasta').read_text().splitlines()[1:]);sequence_accession='P01116-2'
      maps[key]={'uniprot':up,'sequence_accession':sequence_accession,'canonical':pos,'label_seq':label,'label_chain':lc,'resname':r.resname,'expected_aa':seq[pos-1],'observed_aa':protein_letters_3to1.get(r.resname,'X')}
 return maps
def main():
 summary={}
 for target,cfg in CONFIG.items():
  model=MMCIFParser(QUIET=True).get_structure(cfg['pdb'],SRC/(cfg['pdb']+'.cif'))[0]
  maps=map_residues(cfg['pdb'],model);chosen=[];rows=[]
  for chain,(up,beg,end)in cfg['chains'].items():
   for r in model[chain]:
    key=(chain,r.id[1],r.id[2].strip());m=maps.get(key)
    if m and m['uniprot']==up and beg<=m['canonical']<=end:
     m['sequence_match']=m['expected_aa']==m['observed_aa'];m['variant_note']=''
     if not m['sequence_match']:
      assert (target,key,m['expected_aa'],m['observed_aa'])==('ABL',('A',382,''),'D','N'),(target,key,m)
      m['variant_note']='D363N (author382); engineered mutation verified in 1OPL CIF'
     chosen.append(r);rows.append({'index':len(rows),'auth_chain':chain,'auth_seq':r.id[1],'insertion_code':r.id[2].strip(),**m})
  if cfg['receiver_kind']=='ligand': atoms=[a for r in model[cfg['ligand_chain']]if r.resname==cfg['ligand']for a in heavy(r)]
  elif cfg['receiver_kind']=='dna': atoms=[a for c in cfg['dna_chains']for r in model[c]for a in heavy(r)]
  else:atoms=[]
  receiver=[]
  for k,(r,row)in enumerate(zip(chosen,rows)):
   if row['auth_chain']!=cfg['rank_chain']:continue
   yes=(row['canonical']in cfg.get('receiver',[]))if cfg['receiver_kind']=='canonical'else min(np.linalg.norm(a.coord-b.coord)for a in heavy(r)for b in atoms)<=4.0
   if yes:receiver.append(k)
  coords=np.array([r['CA'].coord for r in chosen],float)
  dmat=np.linalg.norm(coords[:,None,:]-coords[None,:,:],axis=2)
  contact=(dmat<=10)&(dmat>0)
  for a in range(len(rows)-1):
   if rows[a]['auth_chain']==rows[a+1]['auth_chain']and rows[a+1]['canonical']==rows[a]['canonical']+1: contact[a,a+1]=contact[a+1,a]=True
  degree=contact.sum(axis=1);mind=dmat[:,receiver].min(axis=1)
  selected={(r.parent.id,r.id)for r in chosen}
  selected_atoms={(r.parent.id,r.id,a.name,a.get_altloc())for r in chosen for a in r if a.element not in ('H','D')}
  class DomainSelect(Select):
   def accept_residue(self,r):return (r.parent.id,r.id)in selected
   def accept_atom(self,a):return int((a.parent.parent.id,a.parent.id,a.name,a.get_altloc())in selected_atoms)
  io=PDBIO();io.set_structure(model);io.save(str(OUT/f'{target}-domain.pdb'),DomainSelect())
  clean=MMCIFParser(QUIET=True) # PDB output uses selected alt A; input MMcif rows use highest occupancy, audited separately.
  from Bio.PDB import PDBParser
  sasamodel=PDBParser(QUIET=True).get_structure(target,OUT/f'{target}-domain.pdb')[0]
  ShrakeRupley(probe_radius=1.4,n_points=960).compute(sasamodel,level='R')
  for i,row in enumerate(rows):
   rr=sasamodel[row['auth_chain']][(' ',row['auth_seq'],row['insertion_code']or' ')]
   row.update(receiver=i in receiver,degree=int(degree[i]),receiver_distance_A=float(mind[i]),sasa_A2=float(rr.sasa),rsa=float(rr.sasa)/MAXASA[row['resname']])
   for cutoff in [4,6,8]:row[f'eligible_{cutoff}A']=bool(row['auth_chain']==cfg['rank_chain']and i not in receiver and mind[i]>=cutoff)
   row['surface_0p20']=row['rsa']>=.2
  out={'target':target,'configuration':cfg,'source_sha256':sha(SRC/f"{cfg['pdb']}.cif"),'sifts_sha256':sha(SRC/f"{cfg['pdb']}-sifts.json"),'N':len(rows),'receivers':[rows[i]for i in receiver],'residues':rows,'scope':'input-only geometric preparation; no response prediction or reference contacts read'}
  (OUT/f'{target}-input.json').write_text(json.dumps(out,indent=2)+'\n')
  with(OUT/f'{target}-residue-map.csv').open('w')as f:
   w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
  summary[target]={'N':len(rows),'receiver_canonical':[rows[i]['canonical']for i in receiver],'receiver_auth':[rows[i]['auth_seq']for i in receiver],'eligible6':sum(r['eligible_6A']for r in rows),'input_sha256':sha(OUT/f'{target}-input.json')}
 (ROOT/'prepared'/'input-freeze.json').write_text(json.dumps({'script_sha256':sha(Path(__file__)),'targets':summary},indent=2)+'\n')
 print(json.dumps(summary,indent=2))

if __name__=="__main__": main()
