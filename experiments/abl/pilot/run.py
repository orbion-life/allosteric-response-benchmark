#!/usr/bin/env python3
"""Portable, locked ABL pilot. Prediction does not open the bound reference.

Usage: python run.py [--output PATH]. The default output is this package.
The two original numerical modules are copied without alteration in vendor/.
"""
from pathlib import Path
import argparse, csv, datetime, hashlib, json, resource, sys, time
import numpy as np
import scipy
from scipy.linalg import eigh
from Bio import __version__ as bioversion
from Bio.PDB import MMCIFParser

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'vendor'))
from prepare import polynomial, gaussian_sd
from run_pilot import harmonic_covariances, gaussian_time_covariance, grid_case, mm, ranking

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,obj): Path(p).write_text(json.dumps(obj,indent=2)+'\n')
def now(): return datetime.datetime.now(datetime.timezone.utc).isoformat()
def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--output',type=Path,default=ROOT); args=parser.parse_args()
    out=args.output.resolve()
    for folder in ['model','results','checks']: (out/folder).mkdir(parents=True,exist_ok=True)
    start=time.monotonic()
    assert sha(ROOT/'preanalysis-protocol.json')==(ROOT/'preanalysis-protocol.sha256').read_text().strip()
    for receipt in json.loads((ROOT/'source-receipts.json').read_text()):
        assert sha(ROOT/receipt['path'])==receipt['sha256'],receipt['path']
    proto=json.loads((ROOT/'preanalysis-protocol.json').read_text())
    inp=json.loads((ROOT/'raw/ABL-input.json').read_text()); rows=inp['residues']
    structure=MMCIFParser(QUIET=True).get_structure('1OPL',ROOT/'raw/1OPL.cif')[0]
    residues=[structure[r['auth_chain']][(' ',r['auth_seq'],r['insertion_code'] or ' ')] for r in rows]
    r0=np.array([r['CA'].coord for r in residues],float); N=len(r0)
    canonical=np.array([r['canonical'] for r in rows]); assert N==252
    assert canonical.tolist()==list(range(242,494))
    assert all(r['auth_seq']==r['canonical']+19 for r in rows)
    variants=[r for r in rows if not r['sequence_match']]
    assert len(variants)==1 and (variants[0]['canonical'],variants[0]['expected_aa'],variants[0]['observed_aa'])==(363,'D','N')
    ligand=[a for r in structure['A'] if r.resname=='P16' for a in r if a.element not in ('H','D')]
    assert ligand
    receiver=np.array([i for i,r in enumerate(residues) if min(np.linalg.norm(a.coord-b.coord) for a in r if a.element not in ('H','D') for b in ligand)<=4])
    assert canonical[receiver].tolist()==[r['canonical'] for r in inp['receivers']]
    dist=np.linalg.norm(r0[:,None]-r0[None],axis=2)
    edges=np.array([(i,j) for i in range(N) for j in range(i+1,N) if dist[i,j]<=10 or canonical[j]-canonical[i]==1])
    degree=np.bincount(edges.ravel(),minlength=N); receiver_distance=dist[:,receiver].min(axis=1)
    candidate=np.array([r['eligible_6A'] for r in rows]); assert int(candidate.sum())==193
    assert np.array_equal(degree,[r['degree'] for r in rows])
    assert np.array_equal(candidate,receiver_distance>=6)
    assert np.allclose(receiver_distance,[r['receiver_distance_A'] for r in rows],atol=1e-12,rtol=0)
    Hessian=np.zeros((3*N,3*N))
    for i,j in edges:
        v=r0[i]-r0[j]; v/=np.linalg.norm(v); h=np.outer(v,v)
        Hessian[3*i:3*i+3,3*i:3*i+3]+=h; Hessian[3*j:3*j+3,3*j:3*j+3]+=h
        Hessian[3*i:3*i+3,3*j:3*j+3]-=h; Hessian[3*j:3*j+3,3*i:3*i+3]-=h
    vals,vecs=eigh(Hessian); keep=vals>1e-8; lam=vals[keep]; B=vecs[:,keep]
    assert np.sum(~keep)==6 and len(lam)==750
    for j in range(B.shape[1]):
        if B[np.argmax(np.abs(B[:,j])),j]<0: B[:,j]*=-1
    model=dict(r0=r0,edges=edges,degree=degree,canonical=canonical,B=B,eigenvalues=lam,Hessian=Hessian,receiver=receiver,candidate=candidate,receiver_distance=receiver_distance)
    np.savez_compressed(out/'model/static-model.npz',**model)
    tau=1/lam[0]
    dump(out/'model/static-model.json',{'protocol_sha256':sha(ROOT/'preanalysis-protocol.json'),'input':'1OPL','chain':'A','uniprot':'P00519','N':N,'contacts':len(edges),'dimensions':len(lam),'null_modes':int(np.sum(~keep)),'tau':float(tau),'smallest_eigenvalues':lam[:8].tolist(),'lambda3_minus_lambda2_over_lambda2':float((lam[2]-lam[1])/lam[1]),'receivers':canonical[receiver].tolist(),'candidate_count':int(candidate.sum()),'variants':variants,'coordinate_convention':'BioPython highest-occupancy CIF atom selection; original float32 coordinates promoted to float64, identical to frozen mask preparation','scope':'Ligand-conditioned isolated kinase domain; SH3/SH2, water and ligands absent from mechanics. Not full protein dynamics.'})
    timing={'preparation_seconds':time.monotonic()-start}
    A,powers=polynomial(r0,edges,B[:,:2]); powers=np.array(powers)
    ownsd,mean,moncov=gaussian_sd(A,powers,lam[:2])
    np.savez_compressed(out/'model/coordinate-d2.npz',A=A,powers=powers,harmonic_sd=ownsd,monomial_mean=mean,monomial_covariance=moncov,B=B[:,:2],eigenvalues=lam[:2])
    for d in [2,len(lam)]:
        begin=time.monotonic(); c0,ct=harmonic_covariances(r0,edges,B[:,:d],lam[:d],tau)
        np.savez_compressed(out/f'results/analytic-harmonic-d{d}.npz',covariance=c0,time_covariance=ct,R=-(c0-ct),sd=np.sqrt(np.diag(c0)))
        timing[f'analytic_d{d}_seconds']=time.monotonic()-begin
        print('analytic',d,round(timing[f'analytic_d{d}_seconds'],2),flush=True)
        if d==2:
            Wick0=mm(mm(A,moncov),A.T); Wickt=mm(mm(A,gaussian_time_covariance(powers,lam[:2],tau)),A.T)
            error=float(max(np.max(np.abs(c0-Wick0)),np.max(np.abs(ct-Wickt))))
            assert error<1e-10,error
            dump(out/'checks/independent-Wick.json',{'max_absolute_covariance_error':error,'passed':True,'checks':'Independent polynomial Gaussian/Wick moments versus contact-displacement Hermite contraction on actual two ABL modes.'})
    full=dict(np.load(out/'results/analytic-harmonic-d750.npz')); sd=full['sd']; assert np.all(sd>0)
    np.savez_compressed(out/'model/all-mode-harmonic-scales.npz',harmonic_sd=sd,dimensions=750)
    for d in [2,750]:
        path=out/f'results/analytic-harmonic-d{d}.npz'; z=dict(np.load(path));z['C']=z['R']/np.outer(sd,sd); z['equilibrium']=z['covariance']/np.outer(sd,sd)
        z['score'],z['order']=ranking(z['C'],receiver,candidate,canonical);np.savez_compressed(path,**z)
    metas=[]
    for d,n,extent in proto['grid_cases']:
        for energy in proto['energies']:
            meta,z=grid_case(model,d,n,extent,energy,sd); name=f'{energy}-d{d}-n{n}-e{extent}'
            np.savez_compressed(out/f'results/{name}.npz',**z);dump(out/f'results/{name}.json',meta);metas.append(meta)
            print(name,round(meta['seconds'],2),meta['top5_canonical'],flush=True)
    adj=np.zeros((N,N));adj[edges[:,0],edges[:,1]]=1;adj[edges[:,1],edges[:,0]]=1
    L=np.diag(adj.sum(axis=1))-adj; lv,lvectors=eigh(L);tg=1/lv[1];heat=mm(lvectors*np.exp(-tg*lv),lvectors.T)
    np.savez_compressed(out/'results/structural-baselines.npz',graph_heat_kernel=np.sqrt(np.mean(heat[:,receiver]**2,axis=1)),degree=degree.astype(float),receiver_inverse_distance=1/np.maximum(receiver_distance,1e-12),graph_time=tg,graph_heat_kernel_matrix=heat)
    primary=np.load(out/'results/biquadratic-d2-n33-e4.npz');fine=np.load(out/'results/biquadratic-d2-n65-e4.npz'); h2=np.load(out/'results/analytic-harmonic-d2.npz'); hall=np.load(out/'results/analytic-harmonic-d750.npz')
    block=np.ix_(np.flatnonzero(candidate),receiver)
    diagnostics={'grid_receiver_block_max_C_difference':float(np.max(np.abs(primary['C']-fine['C'])[block])),'grid_full_matrix_max_C_difference':float(np.max(np.abs(primary['C']-fine['C']))),'grid_top5_set_identical':set(primary['order'][:5])==set(fine['order'][:5]),'coordinate_receiver_block_max_C_difference':float(np.max(np.abs(h2['C']-hall['C'])[block])),'coordinate_full_matrix_max_C_difference':float(np.max(np.abs(h2['C']-hall['C']))),'coordinate_relative_receiver_block_Frobenius_error':float(np.linalg.norm((h2['C']-hall['C'])[block])/np.linalg.norm(hall['C'][block])),'all_grid_cases':metas}
    diagnostics['grid_tolerance_passed']=diagnostics['grid_receiver_block_max_C_difference']<=.001 and diagnostics['grid_top5_set_identical']
    diagnostics['coordinate_tolerance_passed']=diagnostics['coordinate_receiver_block_max_C_difference']<=.002
    dump(out/'results/diagnostics.json',diagnostics)
    timing.update(total_calculation_seconds=time.monotonic()-start,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024))
    assert timing['total_calculation_seconds']<=900 and timing['peak_RSS_bytes']<=8000000000,timing
    dump(out/'checks/runtime.json',{'timing':timing,'environment':{'python':sys.version,'numpy':np.__version__,'scipy':scipy.__version__,'biopython':bioversion,'platform':sys.platform},'scope':'PDB/CIF parsing, Hessian, both harmonic references, six grid cases and graph/geometric controls. Evaluation and replay are separate.','resource_limits_passed':True})
    dump(out/'prediction-freeze.json',{'frozen_utc':now(),'protocol_sha256':sha(ROOT/'preanalysis-protocol.json'),'script_sha256':sha(__file__),'no_bound_reference_opened_for_prediction':True,'files':{str(p.relative_to(out)):sha(p) for folder in ['model','results'] for p in sorted((out/folder).iterdir()) if p.is_file()}})
    print('FROZEN',sha(out/'prediction-freeze.json'),timing,flush=True)
if __name__=='__main__': main()
