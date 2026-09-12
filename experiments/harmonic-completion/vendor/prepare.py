"""Create the prespecified KRAS static contact model before any scoring."""
from pathlib import Path
import json,hashlib,datetime,itertools,math
import numpy as np
from scipy.linalg import eigh
ROOT=Path(__file__).resolve().parent

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def atoms(path, chain='A'):
    chosen={}
    for line in Path(path).read_text().splitlines():
        if line.startswith('ENDMDL'):break
        if line[:6].strip() not in ['ATOM','HETATM'] or line[21]!=chain:continue
        el=line[76:78].strip() or line[12:16].strip()[0]
        if el in ['H','D']:continue
        key=(line[:6].strip(),int(line[22:26]),line[26],line[17:20].strip(),line[12:16].strip())
        alt=line[16];occ=float(line[54:60] or 0); score=(occ,alt==' ',alt=='A')
        if key not in chosen or score>chosen[key][0]:chosen[key]=(score,{'record':key[0],'auth_seq_id':key[1],'ins_code':key[2].strip(),'resname':key[3],'atom':key[4],'xyz':[float(line[k:k+8]) for k in [30,38,46]],'altloc':alt.strip(),'occupancy':occ,'element':el,'chain':chain})
    return [v[1] for v in chosen.values()]

def polynomial(r0,edges,B,kappa=1):
    n=len(r0);d=B.shape[1];deg=np.bincount(edges.ravel(),minlength=n)
    powers=[p for order in [2,3,4] for p in itertools.product(range(order+1),repeat=d) if sum(p)==order]
    A=np.zeros((n,len(powers)))
    for i,j in edges:
        a=r0[i]-r0[j];T=B[3*i:3*i+3]-B[3*j:3*j+3];l2=np.dot(a,a)
        v=2*np.einsum('a,ab->b',a,T);Q=np.einsum('ai,aj->ij',T,T);co={}
        for u in range(d):
            p=[0]*d;p[u]=1;co[tuple(p)]=v[u]
            for w in range(u,d):
                p=[0]*d;p[u]+=1;p[w]+=1;co[tuple(p)]=Q[u,w]*(1 if u==w else 2)
        full={p:0. for p in powers}
        for p,c in co.items():
            for q,z in co.items():full[tuple(a+b for a,b in zip(p,q))]+=c*z*kappa/(8*l2)
        for k,p in enumerate(powers):A[i,k]+=full[p]/deg[i];A[j,k]+=full[p]/deg[j]
    return A,powers

def gaussian_sd(A,powers,lam,beta=1):
    def gm(p):
        if any(i%2 for i in p):return 0.
        return math.prod(math.prod(range(1,k,2))*(1/(beta*l))**(k/2) for k,l in zip(p,lam))
    means=np.array([gm(p) for p in powers]);cov=np.array([[gm(tuple(a+b for a,b in zip(p,q))) for q in powers] for p in powers])-np.outer(means,means)
    return np.sqrt(np.einsum('ip,pq,iq->i',A,cov,A)),means,cov

if __name__=='__main__':
    (ROOT/'model').mkdir(exist_ok=True)
    protocol={'status':'PREANALYSIS; no response scores or MOV contact labels computed at creation','created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'input':'4OBE','input_chain':'A','input_sha256':sha(ROOT/'raw/4OBE.pdb'),'reference':'6OIM','reference_chain':'A','reference_ligand':'MOV','scope':'resolved alpha-carbons mapping to P01116 canonical positions1–166; not the hypervariable tail','mapping':'PDBe SIFTS author-number range and actual amino-acid comparison; report all mutations/missing labels','alternate_locations':'maximum occupancy, then blank, then A; one coordinate per atom','receiver':'input GDP heavy-atom contacts at distance<=4.0Å; exclude water/Mg; GDP removed before mechanics','candidate_primary':'nonreceiver residues whose minimum C-alpha distance to receiver nodes is >=6Å; no surface filter','candidate_sensitivity_A':[4,8],'contacts':'all C-alpha pairs<=10Å plus resolved consecutive canonical backbone neighbours; no artificial gap bridging','kappa':1.0,'beta':1.0,'mu':1.0,'units':'angstrom coordinates; uncalibrated model energy/time; kappa model energy/Å², beta inverse model energy','horizon':'tau=1/(mu*lowest_positive_harmonic_eigenvalue)','normalization':'fixed quartic Ei, analytic all-mode harmonic standard deviations common to all primary model comparisons; own-coordinate harmonic scales exported separately for circuit parity','primary_case':{'d':2,'n':33,'extent_slowest_harmonic_sigma':4,'energy':'biquadratic'},'grid_cases':[[1,4,4],[1,33,4],[1,65,4],[2,17,4],[2,33,4],[2,65,4],[2,41,5],[2,49,6]],'energies':['biquadratic','harmonic','distance_hookean'],'analytic_harmonic_dimensions':[1,2,4,8,16,32,64,'all'],'primary_metric':'receiver-RMS of C; top5 unique candidate nodes, descending score, canonical index tie break','finite_time_ablation':'same-grid equilibrium covariance; also all-mode analytic equilibrium covariance','classical_baselines':['same-grid harmonic','same-grid distance-Hookean','all-mode analytic harmonic','graph heat kernel at time1/lowest_nonzero_graph_eigenvalue','receiver inverse distance','contact degree'],'validation':'only after matrices and ranks are saved: sequence-map reference MOV heavy-atom contact residues<=5Å, report all top5 distances, P@5 and all excluded pocket labels','random_comparison':'finite matched universe, independent target-known-site score comparison under predeclared matching if feasible; unannotated pockets are not negatives','negative_pockets':'unresolved until independently evidenced controls are supplied; never invent negatives','sensitivity_acceptance':'max absolute C difference<=0.001 for grid/domain and same top5; coordinate difference<=0.002 against resolved harmonic all-mode reference; report these engineering thresholds and their failures','physical_diagnostics':'probability any modeled contact strain>0.2, expected max contact strain, strong compression ratio<sqrt(0.5), noncontact CA overlap<2Å, boundary mass; no calibrated biology claims','random_seed':20260912,'runtime_limit_seconds_per_case':120,'interpretation':'retrospective engineering pilot, with benchmark/reference known from challenge; no prospective validation or demonstrated generalization','changes':'any changes require a dated amendment before rerun, retaining failed outputs'}
    target=ROOT/'preanalysis-protocol.json'
    if not target.exists():target.write_text(json.dumps(protocol,indent=2)+'\n')
    (ROOT/'preanalysis-protocol.sha256').write_text(sha(target)+'\n')
    aa=atoms(ROOT/'raw/4OBE.pdb');ca=sorted([a for a in aa if a['record']=='ATOM' and a['atom']=='CA' and 1<=a['auth_seq_id']<=166],key=lambda a:a['auth_seq_id']);r0=np.array([a['xyz'] for a in ca]);N=len(ca)
    seqmap=json.loads((ROOT/'raw/4OBE-mapping.json').read_text())['4obe']['UniProt']['P01116']['mappings'];assert any(m['chain_id']=='A' and m['unp_start']==m['start']['author_residue_number']==1 for m in seqmap)
    ids=np.array([a['auth_seq_id'] for a in ca]);resindex={int(k):i for i,k in enumerate(ids)};gdp=np.array([a['xyz'] for a in aa if a['resname']=='GDP' and a['record']=='HETATM']);assert len(gdp)>0
    receivers=[]
    for resid,i in resindex.items():
        heavy=np.array([a['xyz'] for a in aa if a['record']=='ATOM' and a['auth_seq_id']==resid]);dist=np.sqrt(np.sum((heavy[:,None]-gdp[None])**2,axis=-1)).min()
        if dist<=4:receivers.append(i)
    dist=np.sqrt(np.sum((r0[:,None]-r0[None])**2,axis=-1));edges=np.array([(i,j) for i in range(N) for j in range(i+1,N) if dist[i,j]<=10 or ids[j]-ids[i]==1],dtype=int);deg=np.bincount(edges.ravel(),minlength=N);K=np.zeros((3*N,3*N))
    for i,j in edges:
        v=r0[i]-r0[j];v/=np.linalg.norm(v);h=np.outer(v,v)
        K[3*i:3*i+3,3*i:3*i+3]+=h;K[3*j:3*j+3,3*j:3*j+3]+=h;K[3*i:3*i+3,3*j:3*j+3]-=h;K[3*j:3*j+3,3*i:3*i+3]-=h
    vals,vecs=eigh(K);keep=vals>1e-8;lam=vals[keep];B=vecs[:,keep]
    for j in range(B.shape[1]):
        if B[np.argmax(np.abs(B[:,j])),j]<0:B[:,j]*=-1
    assert np.sum(~keep)==6,('unexpected null modes',vals[:10]);assert deg.min()>0
    recdist=dist[:,receivers].min(axis=1);candidates=(recdist>=6);candidates[receivers]=False
    payload={'status':'static model; before scoring/reference-pocket extraction','protocol_sha256':sha(target),'input':'4OBE','chain':'A','uniprot':'P01116','N':N,'contacts':len(edges),'dimensions':len(lam),'rigid_null_modes':6,'coordinates_A':r0.tolist(),'edges_zero_based':edges.tolist(),'kappa':1,'beta':1,'mu':1,'tau':float(1/lam[0]),'receiver_node_indices':receivers,'receiver_canonical_residues':ids[receivers].tolist(),'primary_candidate_nodes':np.flatnonzero(candidates).tolist(),'canonical_residues':ids.tolist(),'node_map':[dict(node_index=i,canonical_residue=int(a['auth_seq_id']),**a) for i,a in enumerate(ca)],'eigenvector_convention':'orthonormal Cartesian Hessian columns sorted ascending positive eigenvalue; largest absolute element positive','input_ligand_conditioning':'4OBE coordinates determined with GDP; deleting GDP does not undo conditioning','input_resolution_A':json.loads((ROOT/'raw/4OBE-entry.json').read_text())['rcsb_entry_info']['resolution_combined'],'source_receipts':'../raw/source-receipts.json'}
    (ROOT/'model/static-model.json').write_text(json.dumps(payload,indent=2)+'\n');np.savez_compressed(ROOT/'model/static-model.npz',r0=r0,edges=edges,degree=deg,canonical=ids,B=B,eigenvalues=lam,Hessian=K,receiver=np.array(receivers),candidate=candidates,receiver_distance=recdist)
    for d in [1,2]:
        A,powers=polynomial(r0,edges,B[:,:d]);sh,pm,pc=gaussian_sd(A,powers,lam[:d]);np.savez_compressed(ROOT/f'model/coordinate-d{d}.npz',r0=r0,edges=edges,B=B[:,:d],eigenvalues=lam[:d],A=A,powers=np.array(powers),harmonic_sd=sh,monomial_mean=pm,monomial_covariance=pc,canonical=ids,receiver=np.array(receivers),candidate=candidates)
    print(json.dumps({'N':N,'edges':len(edges),'internal_modes':len(lam),'tau':payload['tau'],'receivers':payload['receiver_canonical_residues'],'candidates':int(candidates.sum()),'protocol_sha256':payload['protocol_sha256']}))
