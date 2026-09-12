"""Deterministic finite-grid pilot and analytic all-mode harmonic control.

Run after prepare.py. Does not inspect the bound reference or pocket labels.
The all-mode control uses the SAME quartic residue observables as the pilot.
"""
from pathlib import Path
import json, time, hashlib, datetime, itertools, functools
import numpy as np
from scipy.linalg import eigh
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import expm_multiply
from scipy.special import expit, logsumexp
from scipy.stats import spearmanr
from prepare import ROOT, sha, polynomial, gaussian_sd

def mm(a,b):
    # Explicit contractions also avoid spurious Apple Accelerate FPE warnings.
    if a.ndim==1:return np.einsum('i,ij->j',a,b,optimize=False)
    if b.ndim==1:return np.einsum('ij,j->i',a,b,optimize=False)
    return np.einsum('...ij,...jk->...ik',a,b,optimize=False)

def contracted_covariance(a,b,S,T,C):
    """Hermite orders 1–4 of (2 a.x + x.x)^2; leading axes broadcast."""
    tr=lambda m:np.trace(m,axis1=-2,axis2=-1)
    bil=lambda x,M,y:np.einsum('...i,...ij,...j->...',x,M,y,optimize=False)
    h=4*(tr(S)[...,None]*a+2*np.einsum('...ij,...j->...i',S,a))
    k=4*(tr(T)[...,None]*b+2*np.einsum('...ij,...j->...i',T,b))
    F=4*a[...,None]*a[...,None,:]+2*tr(S)[...,None,None]*np.eye(3)+4*S
    G=4*b[...,None]*b[...,None,:]+2*tr(T)[...,None,None]*np.eye(3)+4*T
    Ct=np.swapaxes(C,-1,-2);CC=mm(C,Ct)
    return (bil(h,C,k)+2*tr(mm(mm(mm(F,C),G),Ct))+32*bil(a,C,b)*tr(CC)
            +64*bil(a,mm(CC,C),b)+8*tr(CC)**2+16*tr(mm(CC,CC)))

def harmonic_covariances(r0,edges,B,lam,tau,beta=1):
    """Exact quartic-observable OU covariance, contracted in 3-D per edge."""
    n=len(r0);m=len(edges);a=r0[edges[:,0]]-r0[edges[:,1]];coef=1/(8*np.sum(a*a,axis=1))
    deg=np.bincount(edges.ravel(),minlength=n);W=np.zeros((n,m))
    for e,(i,j) in enumerate(edges):W[i,e]=1/deg[i];W[j,e]=1/deg[j]
    V=np.einsum('ik,k,jk->ij',B,1/(beta*lam),B,optimize=False).reshape(n,3,n,3).transpose(0,2,1,3)
    S=V[edges[:,0],edges[:,0]]-V[edges[:,0],edges[:,1]]-V[edges[:,1],edges[:,0]]+V[edges[:,1],edges[:,1]]
    out=[]
    for dt in [0,tau]:
        Q=V if dt==0 else np.einsum('ik,k,jk->ij',B,np.exp(-lam*dt)/(beta*lam),B,optimize=False).reshape(n,3,n,3).transpose(0,2,1,3)
        # Multiply by residue incidence immediately, retaining only N × M output.
        KW=np.zeros((m,n))
        for start in range(0,m,32):
            idx=np.arange(start,min(start+32,m));x=edges[idx,0,None];y=edges[idx,1,None];u=edges[None,:,0];v=edges[None,:,1]
            C=Q[x,u]-Q[x,v]-Q[y,u]+Q[y,v]
            K=contracted_covariance(a[idx,None],a[None],S[idx,None],S[None],C)*coef[idx,None]*coef[None]
            for node in range(n):KW[idx,node]=K[:,np.any(edges==node,axis=1)].mean(axis=1)
        out.append(np.array([KW[np.any(edges==node,axis=1)].mean(axis=0) for node in range(n)]))
    return out

def gaussian_time_covariance(powers,lam,tau):
    d=len(lam);V=np.diag(1/lam);Ct=np.diag(np.exp(-lam*tau)/lam);S=np.block([[V,Ct],[Ct,V]])
    @functools.lru_cache(None)
    def wick(ids):
        if not ids:return 1.
        if len(ids)%2:return 0.
        return sum(S[ids[0],ids[j]]*wick(ids[1:j]+ids[j+1:]) for j in range(1,len(ids)))
    def ids(p):return tuple(i for i,k in enumerate(p) for _ in range(k))
    means=np.array([wick(ids(tuple(p)+tuple([0]*d))) for p in powers])
    return np.array([[wick(ids(tuple(p)+tuple(q))) for q in powers] for p in powers])-np.outer(means,means)

def ranking(C,rec,candidate,canonical):
    score=np.sqrt(np.mean(C[:,rec]**2,axis=1));idx=np.flatnonzero(candidate)
    order=idx[np.lexsort((canonical[idx],-score[idx]))]
    return score,order

def grid_case(model,d,n,extent,energy,common_sd):
    begin=time.monotonic();r0=model['r0'];edges=model['edges'];B=model['B'][:,:d];lam=model['eigenvalues'][:d];tau=1/model['eigenvalues'][0]
    A,powers=polynomial(r0,edges,B);powers=np.array(powers);axis=np.linspace(-extent/np.sqrt(lam[0]),extent/np.sqrt(lam[0]),n);delta=axis[1]-axis[0]
    inds=np.indices((n,)*d).reshape(d,-1).T;q=axis[inds];G=len(q)
    phi=np.prod(q[:,None,:]**powers[None,:,:],axis=2);E=mm(phi,A.T)
    Ufull=np.zeros(G);Uhook=np.zeros(G);maxstrain=np.zeros(G);compression=np.zeros(G,dtype=bool)
    rr=r0[None]+mm(q,B.T).reshape(G,len(r0),3)
    for i,j in edges:
        l2=np.sum((r0[i]-r0[j])**2);r2=np.sum((rr[:,i]-rr[:,j])**2,axis=1)
        Ufull+=(r2-l2)**2/(8*l2);ratio=np.sqrt(r2/l2);Uhook+=.5*l2*(ratio-1)**2
        maxstrain=np.maximum(maxstrain,np.abs(ratio-1));compression|=ratio<np.sqrt(.5)
    U={'biquadratic':Ufull,'harmonic':.5*np.sum(lam*q*q,axis=1),'distance_hookean':Uhook}[energy]
    pi=np.exp(-U-logsumexp(-U));sqrtpi=np.sqrt(pi);center=phi-mm(pi,phi);f=sqrtpi[:,None]*center
    row=[];col=[];rate=[];diag=np.zeros(G)
    flat=np.arange(G).reshape((n,)*d)
    for dim in range(d):
        lo=[slice(None)]*d;hi=lo.copy();lo[dim]=slice(0,-1);hi[dim]=slice(1,None);x=flat[tuple(lo)].ravel();y=flat[tuple(hi)].ravel()
        xy=2/delta**2*expit(-(U[y]-U[x]));yx=2/delta**2*expit(U[y]-U[x]);h=-np.sqrt(xy*yx)
        row.extend([x,y]);col.extend([y,x]);rate.extend([h,h]);np.add.at(diag,x,xy);np.add.at(diag,y,yx)
    row.append(np.arange(G));col.append(np.arange(G));rate.append(diag)
    H=coo_matrix((np.concatenate(rate),(np.concatenate(row),np.concatenate(col))),shape=(G,G)).tocsr()
    ft=expm_multiply(-tau*H,f,traceA=-tau*diag.sum());cov=mm(f.T,f);timecov=mm(f.T,ft)
    R=-mm(mm(A,cov-timecov),A.T);C=R/np.outer(common_sd,common_sd);equilibrium=mm(mm(A,cov),A.T)/np.outer(common_sd,common_sd)
    score,order=ranking(C,model['receiver'],model['candidate'],model['canonical'])
    boundary=np.any((inds==0)|(inds==n-1),axis=1)
    # Check all non-contact pairs for gross C-alpha collision; this is only a diagnostic.
    adjacency=np.eye(len(r0),dtype=bool);adjacency[edges[:,0],edges[:,1]]=True;adjacency[edges[:,1],edges[:,0]]=True;overlap=np.zeros(G,dtype=bool)
    for i in range(len(r0)):
        js=np.flatnonzero(~adjacency[i] & (np.arange(len(r0))>i))
        if len(js):overlap|=np.any(np.sum((rr[:,i,None]-rr[:,js])**2,axis=-1)<4,axis=1)
    result={'energy':energy,'d':d,'n':n,'extent_slowest_sigma':extent,'delta_A':float(delta),'states':G,'tau_model_time':float(tau),'seconds':time.monotonic()-begin,'max_stationarity_residual':float(np.max(np.abs(H@sqrtpi))),'boundary_probability':float(pi[boundary].sum()),'prob_any_contact_strain_gt_0.2':float(pi[maxstrain>.2].sum()),'mean_max_contact_strain':float(np.sum(pi*maxstrain)),'prob_strong_compression':float(pi[compression].sum()),'prob_noncontact_CA_overlap_lt_2A':float(pi[overlap].sum()),'top5_canonical':model['canonical'][order[:5]].tolist(),'max_abs_C':float(np.max(np.abs(C))),'normalization':'analytic all-mode harmonic SD of fixed quartic Ei','prob_sum':float(pi.sum())}
    return result,dict(C=C,R=R,equilibrium=equilibrium,score=score,order=order,q=q,U=U,pi=pi,H_data=H.data,H_indices=H.indices,H_indptr=H.indptr,A=A,powers=powers,monomial_centered=center,monomial_covariance=cov,monomial_time_covariance=timecov,canonical=model['canonical'])

def main():
    out=ROOT/'results';out.mkdir(exist_ok=True);model=np.load(ROOT/'model/static-model.npz');tau=1/model['eigenvalues'][0];protocol=json.loads((ROOT/'preanalysis-protocol.json').read_text());checks=[];summary=[]
    # Independent polynomial moment check of the contracted formula on actual KRAS modes.
    for d in [1,2]:
        f=np.load(ROOT/f'model/coordinate-d{d}.npz');c0,ct=harmonic_covariances(model['r0'],model['edges'],model['B'][:,:d],model['eigenvalues'][:d],tau)
        pc0=mm(mm(f['A'],f['monomial_covariance']),f['A'].T);pct=mm(mm(f['A'],gaussian_time_covariance(f['powers'],model['eigenvalues'][:d],tau)),f['A'].T)
        error=max(np.max(np.abs(c0-pc0)),np.max(np.abs(ct-pct)));assert error<1e-10,(d,error)
        checks.append({'d':d,'max_covariance_error_against_independent_mode_Wick_moments':float(error)})
    for d in [1,2,4,8,16,32,64,len(model['eigenvalues'])]:
        start=time.monotonic();c0,ct=harmonic_covariances(model['r0'],model['edges'],model['B'][:,:d],model['eigenvalues'][:d],tau)
        np.savez_compressed(out/f'analytic-harmonic-d{d}.npz',covariance=c0,time_covariance=ct,R=-(c0-ct),sd=np.sqrt(np.diag(c0)))
        summary.append({'d':d,'seconds':time.monotonic()-start});print('analytic',d,flush=True)
    full=np.load(out/f"analytic-harmonic-d{len(model['eigenvalues'])}.npz");sd=full['sd'];np.savez_compressed(ROOT/'model/all-mode-harmonic-scales.npz',harmonic_sd=sd,dimensions=len(model['eigenvalues']))
    for d in [1,2,4,8,16,32,64,len(model['eigenvalues'])]:
        fn=out/f'analytic-harmonic-d{d}.npz';a=dict(np.load(fn));a['C']=a['R']/np.outer(sd,sd);a['equilibrium']=a['covariance']/np.outer(sd,sd);a['score'],a['order']=ranking(a['C'],model['receiver'],model['candidate'],model['canonical']);np.savez_compressed(fn,**a)
    (out/'analytic-checks.json').write_text(json.dumps({'checks':checks,'timings':summary},indent=2)+'\n')
    results=[]
    for d,n,extent in protocol['grid_cases']:
        for energy in protocol['energies']:
            name=f'{energy}-d{d}-n{n}-e{extent}';meta,arrays=grid_case(model,d,n,extent,energy,sd);np.savez_compressed(out/f'{name}.npz',**arrays);(out/f'{name}.json').write_text(json.dumps(meta,indent=2)+'\n');results.append(meta);print(name,round(meta['seconds'],2),meta['top5_canonical'],flush=True)
    # Equal-candidate graph and geometric baselines, saved before reference labels.
    N=len(model['r0']);adj=np.zeros((N,N));ed=model['edges'];adj[ed[:,0],ed[:,1]]=1;adj[ed[:,1],ed[:,0]]=1;L=np.diag(adj.sum(axis=1))-adj;vals,V=eigh(L);tgraph=1/vals[1];heat=mm(V*np.exp(-tgraph*vals),V.T)
    baselines={'graph_heat_kernel':np.sqrt(np.mean(heat[:,model['receiver']]**2,axis=1)),'degree':model['degree'].astype(float),'receiver_inverse_distance':1/np.maximum(model['receiver_distance'],1e-12)}
    np.savez_compressed(out/'structural-baselines.npz',**baselines,graph_time=tgraph,graph_heat_kernel_matrix=heat)
    (out/'run-summary.json').write_text(json.dumps({'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'protocol_sha256':sha(ROOT/'preanalysis-protocol.json'),'script_sha256':sha(__file__),'all_grid_cases':results,'no_bound_reference_labels_read':True},indent=2)+'\n')
    receipts={str(p.relative_to(ROOT)):sha(p) for p in sorted(out.iterdir()) if p.is_file()};(ROOT/'prediction-freeze.json').write_text(json.dumps({'frozen_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'files':receipts},indent=2)+'\n')
    print('PREDICTIONS FROZEN',sha(ROOT/'prediction-freeze.json'),flush=True)

if __name__=='__main__':main()
