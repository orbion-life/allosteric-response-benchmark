"""Independent equation checks plus a repeat of the primary calculation."""
import json,datetime,itertools
import numpy as np
from scipy.linalg import expm
from scipy.special import expit
from prepare import ROOT,sha,atoms
from run_pilot import contracted_covariance,grid_case,mm

def main():
    rng=np.random.default_rng(20260912);z,w=np.polynomial.hermite.hermgauss(5);inds=np.array(list(itertools.product(range(5),repeat=6)));nodes=np.sqrt(2)*z[inds];weights=np.prod(w[inds]/np.sqrt(np.pi),axis=1);checks=[]
    for index in range(6):
        T=rng.normal(size=(6,6))*.17;S=mm(T,T.T)+np.eye(6)*.02;a=rng.normal(size=3)*2;b=rng.normal(size=3)*2;v=mm(nodes,np.linalg.cholesky(S).T)
        X=2*np.einsum('ij,j->i',v[:,:3],a)+np.sum(v[:,:3]**2,axis=1);Y=2*np.einsum('ij,j->i',v[:,3:],b)+np.sum(v[:,3:]**2,axis=1)
        exact=float(contracted_covariance(a,b,S[:3,:3],S[3:,3:],S[:3,3:]));quad=float(np.sum(weights*X**2*Y**2)-np.sum(weights*X**2)*np.sum(weights*Y**2));err=abs(exact-quad);assert err<1e-10
        checks.append({'fixture':index,'Hermite_contraction':exact,'independent_6D_Gauss_Hermite_quadrature':quad,'absolute_difference':err})
    m=np.load(ROOT/'model/static-model.npz');sd=np.load(ROOT/'model/all-mode-harmonic-scales.npz')['harmonic_sd'];meta,again=grid_case(m,2,33,4,'biquadratic',sd);old=np.load(ROOT/'results/biquadratic-d2-n33-e4.npz');differences={key:float(np.max(np.abs(again[key]-old[key]))) for key in ['C','R','pi','U','score']};assert max(differences.values())<1e-12;assert np.array_equal(again['order'],old['order'])
    # Direct finite-difference perturbation of the Markov generator on the tiny fixture.
    tiny=np.load(ROOT/'results/biquadratic-d1-n4-e4.npz');U=tiny['U'];q=tiny['q'];delta=q[1,0]-q[0,0];E=mm(tiny['monomial_centered'],tiny['A'].T);p=tiny['pi'];tau=1/m['eigenvalues'][0]
    def expectation(h,i):
        u=U+h*E[:,i];L=np.zeros((4,4))
        for x in range(3):
            y=x+1;L[x,y]=2/delta**2*expit(-(u[y]-u[x]));L[y,x]=2/delta**2*expit(u[y]-u[x])
        L-=np.diag(L.sum(axis=1));return mm(mm(p,expm(tau*L)),E)
    fd=[]
    # h≥0 follows the proposed intervention; the tiny finite grid is bounded.
    for i in [3,59,156]:
        h=1e-4;derivative=(expectation(h,i)-expectation(0,i))/h;err=float(np.max(np.abs(derivative-tiny['R'][:,i])));assert err<1e-7
        fd.append({'perturbed_canonical_residue':int(m['canonical'][i]),'one_sided_h':h,'max_absolute_response_difference':err})
    result={'verified_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'primary_repetition':{'case':'biquadratic d2 n33 extent4','maximum_absolute_array_differences':differences,'identical_full_candidate_order':True},'contact_covariance_formula':checks,'direct_finite_field_response':fd,'actual_KRAS_d1_d2_Wick_checks':'results/analytic-checks.json','scope':'These checks validate the implementation and exact harmonic contraction. They do not validate a reduced protein model biologically.','source_sha256':{x:sha(ROOT/x) for x in ['prepare.py','run_pilot.py','evaluate.py','verify.py']}}
    (ROOT/'verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
