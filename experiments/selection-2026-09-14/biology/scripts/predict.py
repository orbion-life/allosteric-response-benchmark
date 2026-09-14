"""Fixed GRB2 predictions; this script never reads assay outcomes."""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import sys,json,hashlib,time,resource,subprocess,datetime,platform
import numpy as np,pandas as pd
from scipy.linalg import eigh
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
P=Path(__file__).resolve().parents[1];sys.path.insert(0,str(P/'vendor'))
import gaussian_protein as gp
from run_pilot import harmonic_covariances

def sha(f):return hashlib.sha256(Path(f).read_bytes()).hexdigest()

def model_from(r0,receiver,canonical,out):
 start=time.monotonic();n=len(r0);dist=np.linalg.norm(r0[:,None,:]-r0[None,:,:],axis=2)
 edges=np.array([(i,j) for i in range(n) for j in range(i+1,n) if dist[i,j]<=10 or canonical[j]-canonical[i]==1],int)
 degree=np.bincount(edges.ravel(),minlength=n);K=np.zeros((3*n,3*n))
 for i,j in edges:
  v=r0[i]-r0[j];v=v/np.linalg.norm(v);h=np.outer(v,v);si=slice(3*i,3*i+3);sj=slice(3*j,3*j+3);K[si,si]+=h;K[sj,sj]+=h;K[si,sj]-=h;K[sj,si]-=h
 vals,B=eigh(K);keep=vals>1e-8
 if np.sum(~keep)!=6:raise ArithmeticError('Model does not have exactly six rigid null modes.')
 lam=vals[keep];B=B[:,keep]
 for j in range(B.shape[1]):
  if B[np.argmax(abs(B[:,j])),j]<0:B[:,j]*=-1
 rec=np.flatnonzero(receiver);recdist=dist[:,rec].min(axis=1);candidate=(recdist>=6)&~receiver
 out.mkdir(parents=True,exist_ok=True);arrays=dict(r0=r0,edges=edges,degree=degree,Hessian=K,B=B,eigenvalues=lam,receiver=rec,candidate=candidate,canonical=canonical,receiver_distance=recdist)
 np.savez_compressed(out/'static-model.npz',**arrays)
 return arrays,{'N':n,'contacts':len(edges),'internal_dimensions':len(lam),'null_modes':int(np.sum(~keep)),'tau':float(1/lam[0]),'candidate_count':int(candidate.sum()),'build_seconds':time.monotonic()-start}

def solve(r0,receiver,canonical,out,label):
 start=time.monotonic();m,meta=model_from(r0,receiver,canonical,out);tau=meta['tau'];h0,ht=harmonic_covariances(r0,m['edges'],m['B'],m['eigenvalues'],tau);sd=np.sqrt(np.diag(h0));assert np.all(sd>0)
 np.savez_compressed(out/'harmonic-scales.npz',harmonic_sd=sd);hc=(ht-h0)/np.outer(sd,sd);np.savez_compressed(out/'harmonic-response.npz',C=hc,covariance=h0,delayed_covariance=ht,harmonic_sd=sd)
 check,cost=gp.run_model(out/'static-model.npz',out/'harmonic-scales.npz',out/'gaussian',label=label)
 assert check['scientific_status']=='PASS'
 gc=np.load(out/'gaussian/response-matrices.npz')['C'][1];scores={k:np.sqrt(np.mean(C[m['receiver'],:]**2,axis=0)) for k,C in [('gaussian',gc),('harmonic',hc)]}
 receipt=dict(meta,gaussian_numerical_checks=check,gaussian_cost=cost,seconds=time.monotonic()-start,model_sha256=sha(out/'static-model.npz'));(out/'execution-receipt.json').write_text(json.dumps(receipt,indent=2));return m,scores,receipt

def main():
 start=time.monotonic();freeze=json.loads((P/'protocol-freeze.json').read_text());assert sha(P/'PROTOCOL.md')==freeze['protocol_sha256']
 z=np.load(P/'inputs/assay-domain-input.npz');assert len(z['r0'])==56 and np.array_equal(z['canonical'],np.arange(159,215));receiver=z['receiver'];canonical=z['canonical'];m,scores,receipt=solve(z['r0'],receiver,canonical,P/'results/primary-model','GRB2-1GFC')
 # Generate same heavy-atom domain for upstream Ohm, with canonical numbering.
 lines=[]
 for line in (P/'inputs/1GFC.pdb').read_text().splitlines():
  if line.startswith('ENDMDL'):break
  if line.startswith('ATOM  ') and line[21]=='A' and 3<=int(line[22:26])<=58 and line[76:78].strip() not in ['H','D']:
   canon=int(line[22:26])+156;lines.append(line[:22]+f'{canon:4d}'+line[26:])
 (P/'ohm/grb2.pdb').write_text('\n'.join(lines)+'\nTER\nEND\n');commands=[]
 def run(argv):
  st=time.monotonic();r=subprocess.run(argv,cwd=P/'ohm',capture_output=True,text=True,timeout=90);commands.append(dict(argv=argv,seconds=time.monotonic()-st,returncode=r.returncode,stdout=r.stdout,stderr=r.stderr));assert r.returncode==0,r.stderr
 run([str(P/'ohm/bin/contacts'),'grb2.pdb','-oi','grb2.ind','-oc','grb2.mat','-c','3.4'])
 ind=pd.read_csv(P/'ohm/grb2.ind',sep=r'\s+',header=None,names=['node','chain','residue']);ids=np.array([int(x.split('/')[0]) for x in ind.residue]);assert np.array_equal(ids,canonical) and np.array_equal(ind.node,np.arange(1,57))
 sites='+'.join(str(int(i)) for i in ind.node[receiver]);frame=pd.DataFrame(dict(canonical=canonical,construct_position=canonical-158,receiver=receiver,receiver_CA_distance_A=m['receiver_distance'],structural_candidate=m['candidate'],degree=m['degree'],negative_distance=-m['receiver_distance'],**scores))
 for name,seed in [('seed11',11),('seed29',29),('seed47',47),('repeat11',11)]:
  run([str(P/'ohm/bin/diffuse'),'aci','grb2.mat',sites,name+'.nodes',name+'.bonds','-n','10000','-c','0.05','-norm','-a','4.5','-seed',str(seed),'-mdist',name+'.mdist'])
  arr=np.loadtxt(P/'ohm'/f'{name}.nodes');assert arr.shape==(56,) and np.isfinite(arr).all()
  if name!='repeat11':frame['ohm_'+name]=arr
 assert (P/'ohm/seed11.nodes').read_bytes()==(P/'ohm/repeat11.nodes').read_bytes();frame['ohm']=frame[['ohm_seed11','ohm_seed29','ohm_seed47']].mean(axis=1)
 frame.to_csv(P/'results/predictions-before-label-join.csv',index=False);receipt.update(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),protocol_sha256=freeze['protocol_sha256'],prediction_sha256=sha(P/'results/predictions-before-label-join.csv'),receiver_canonical=canonical[receiver].tolist(),candidate_canonical=canonical[m['candidate']].tolist(),ohm_commands=commands,ohm_seconds=sum(c['seconds'] for c in commands),primary_whole_score_seconds=time.monotonic()-start,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024),environment={'python':sys.version,'numpy':np.__version__,'platform':platform.platform()},scope='Predictions complete before any assay join; no outcomes read by this script.')
 (P/'results/prediction-receipt.json').write_text(json.dumps(receipt,indent=2));print('PRIMARY_COMPLETE',receipt['candidate_canonical'],receipt['primary_whole_score_seconds'],flush=True)
 # All20 alternate NMR models predeclared; no best-model selection.
 d=MMCIF2Dict(str(P/'inputs/1GFD.cif'));xs=[];secondary=[]
 for model in range(1,21):
  coords=[]
  for i,(atm,ch,rnum,mnum) in enumerate(zip(d['_atom_site.label_atom_id'],d['_atom_site.auth_asym_id'],d['_atom_site.auth_seq_id'],d['_atom_site.pdbx_PDB_model_num'])):
   if atm=='CA' and ch=='A' and 3<=int(rnum)<=58 and int(mnum)==model:coords.append((int(rnum),[float(d['_atom_site.Cartn_'+a][i]) for a in 'xyz']))
  coords.sort();r0=np.array([r[1] for r in coords]);assert len(r0)==56
  out=P/'results/structure-sensitivity'/f'1GFD-model{model:02d}'
  try:
   ms,sc,rc=solve(r0,receiver,canonical,out,f'GRB2-1GFD-{model}');secondary.append({'model':model,'status':'PASS','seconds':rc['seconds'],'tau':rc['tau']});xs.append(pd.DataFrame({'model':model,'canonical':canonical,**sc,'candidate_in_own_geometry':ms['candidate']}))
  except Exception as e:secondary.append({'model':model,'status':'FAIL','error':str(e)})
 pd.concat(xs).to_csv(P/'results/structure-sensitivity-predictions.csv',index=False)
 total={'status':'COMPLETE','whole_primary_and_sensitivity_seconds':time.monotonic()-start,'primary_seconds':receipt['primary_whole_score_seconds'],'models':secondary,'peak_RSS_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024),'predictions_never_read_outcomes':True};(P/'results/whole-prediction-cost.json').write_text(json.dumps(total,indent=2));print('ALL_PREDICTIONS_COMPLETE',total['whole_primary_and_sensitivity_seconds'],flush=True)
if __name__=='__main__':main()
