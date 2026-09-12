"""Independent finite-intervention check and optional complete replay comparison."""
from pathlib import Path
import argparse,json
import numpy as np
from scipy.linalg import eigh
from scipy.special import expit,logsumexp
from run import static_model,grid,ROOT

def intervention_check():
    protocol=json.loads((ROOT/'preanalysis.json').read_text());records=[]
    for k in protocol['kappa']:
        r0,edges,lam,B,sd=static_model(protocol,k)
        n=5;extent=4;axis=np.linspace(-extent/np.sqrt(lam[0]),extent/np.sqrt(lam[0]),n)
        ids=np.indices((n,n,n)).reshape(3,-1).T;q=axis[ids]
        xyz=r0[None]+np.einsum('ik,jk->ij',q,B,optimize=False).reshape(len(q),3,3)
        energy=np.zeros((len(q),3))
        for i,j in edges:
            l2=np.sum((r0[i]-r0[j])**2);r2=np.sum((xyz[:,i]-xyz[:,j])**2,axis=1)
            u=k*(r2-l2)**2/(8*l2);energy[:,i]+=u/2;energy[:,j]+=u/2
        U=np.sum(energy,axis=1) # Every contact appears twice, each divided by degree2.
        pi=np.exp(-U-logsumexp(-U));spacing=axis[1]-axis[0]
        neighbours=np.sum(np.abs(ids[:,None]-ids[None]),axis=-1)==1
        times=np.array([1/lam[0]])
        exact,_=grid(r0,edges,lam,B,k,3,n,extent,'biquadratic',times,sd)
        errors=[]
        for h in [1e-4,5e-5]:
            response=np.zeros((3,3))
            for i in range(3):
                expectations=[]
                for sign in [-1,1]:
                    u=U+sign*h*energy[:,i]
                    L=2/spacing**2*expit(-(u[None]-u[:,None]))*neighbours
                    L[np.diag_indices(len(L))]=-L.sum(axis=1)
                    # Dense reversible diagonalization avoids a known local
                    # Accelerate warning inside scipy.linalg.expm's matmul.
                    # Evolve the original equilibrium under the perturbed rates;
                    # do not evaluate any covariance/linear-response identity.
                    H=-np.sqrt(L*L.T)
                    H[np.diag_indices(len(H))]=-np.diag(L)
                    vals,V=eigh(H)
                    log_pi0=-U-logsumexp(-U);log_pih=-u-logsumexp(-u)
                    initial=np.exp(log_pi0-.5*log_pih)
                    coeff=np.einsum('ji,j->i',V,initial,optimize=False)
                    evolved=np.exp(.5*log_pih)*np.einsum('ij,j->i',V,np.exp(-times[0]*vals)*coeff,optimize=False)
                    assert abs(evolved.sum()-1)<1e-9 and evolved.min()>-1e-12
                    expectations.append(np.einsum('i,ij->j',evolved,energy,optimize=False))
                response[i]=(expectations[1]-expectations[0])/(2*h)
            error=float(np.max(np.abs(response/np.outer(sd,sd)-exact['C'][0])))
            assert error<1e-7,(k,h,error);errors.append(dict(h=h,max_normalized_C_error=error))
        records.append(dict(kappa=k,states=len(q),errors=errors))
    return records

def main(reference,repeat,output):
    result={'finite_intervention_checks':intervention_check(),'meaning':'Independent perturbation of generator rates, using dense probability evolution from unperturbed equilibrium; not another use of the covariance identity.'}
    if repeat:
        counts=0;maximum=0
        for p in sorted(reference.glob('*.npz')):
            with np.load(p) as a,np.load(repeat/p.name) as b:
                assert set(a.files)==set(b.files)
                for key in a.files:
                    x,y=a[key],b[key];assert x.shape==y.shape and np.isfinite(x).all() and np.isfinite(y).all()
                    delta=float(np.max(np.abs(x-y)));maximum=max(maximum,delta)
                    if key in ['C','equilibrium','delayed']:assert delta<=1e-10,(p.name,key,delta)
                    else:assert np.allclose(x,y,atol=1e-12,rtol=1e-10),(p.name,key)
                    counts+=1
        assert len(list(reference.glob('*.npz')))==81
        result.update(replayed_arrays=counts,replayed_archives=81,max_absolute_array_difference=maximum)
    result['status']='PASS';output.parent.mkdir(exist_ok=True,parents=True);output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--reference',type=Path,default=ROOT/'results');p.add_argument('--repeat',type=Path);p.add_argument('--output',type=Path,default=ROOT/'checks/verification.json');a=p.parse_args();main(a.reference,a.repeat,a.output)
