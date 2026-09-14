"""Frozen Gaussian and harmonic responses, independently checked by Wick moments."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
from pathlib import Path
from functools import lru_cache
import sys,json,time,hashlib,platform
import numpy as np
import scipy
from scipy.linalg import expm
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'source/vendor'))
from model import build_model,analytic_gaussian,energies,polynomial_values
sys.path.insert(0,str(ROOT/'source/gaussian'))
from experiment import OUKernel

def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def wick_kernel(m,Sigma,t,sd):
    # Independent coordinate-polynomial moments, rather than edge covariance contractions.
    scale=np.sqrt(m['beta']*m['eigenvalues'])
    cov=scale[:,None]*Sigma*scale[None,:]
    cross=scale[:,None]*(expm(-t*np.linalg.inv(Sigma))@Sigma)*scale[None,:]
    joint=np.block([[cov,cross],[cross.T,cov]])
    @lru_cache(None)
    def moment(p):
        n=sum(p)
        if n==0:return 1.
        if n%2:return 0.
        i=next(i for i,k in enumerate(p) if k)
        q=list(p);q[i]-=1;s=0.
        for j,k in enumerate(q):
            if k:
                z=q.copy();z[j]-=1;s+=k*joint[i,j]*moment(tuple(z))
        return s
    powers=m['powers'];a=m['coefficients']
    mono=np.array([moment(tuple(p)+(0,0,0)) for p in powers]);mu=np.einsum('ij,j->i',a,mono)
    raw=np.array([[moment(tuple(p)+tuple(q)) for q in powers] for p in powers])
    value=np.einsum('ai,ij,bj->ab',a,raw,a)-np.outer(mu,mu)
    return value/np.outer(sd,sd)

def evaluate():
    start=time.perf_counter();p=json.loads((ROOT/'protocol.json').read_text());m=build_model()
    ts=np.array(p['all_times_sorted']);physical=ts/(m['mu']*m['eigenvalues'][0]);m['times']=physical
    h=analytic_gaussian(m);saved=np.load(ROOT/'inputs/frozen-gaussian-result.npz')
    old=np.load(ROOT/'inputs/historical/controls.npz');Sigma=saved['Sigma'].copy()
    assert np.array_equal(Sigma,old['Sigma_Gaussian'])
    model={'r0':m['r0'],'edges':m['edges'],'B':m['basis'],'eigenvalues':m['eigenvalues']}
    ou=OUKernel(model,Sigma,h['sd']);g0=ou.kernel(0,False);gK=np.array([ou.kernel(t,False) for t in physical]);gC=gK-g0
    gaussian_wick=np.array([wick_kernel(m,Sigma,t,h['sd']) for t in [0,*physical]])
    native_sigma=np.diag(1/(m['beta']*m['eigenvalues']))
    harmonic_wick=np.array([wick_kernel(m,native_sigma,t,h['sd']) for t in [0,*physical]])
    indices=[p['all_times_sorted'].index(t) for t in p['original_times']]
    points=np.array([[i,j,k] for i in [-.7,0,1.3] for j in [-.7,0,1.3] for k in [-.7,0,1.3]])
    checks={'saved_scale_maxabs':float(np.max(abs(h['sd']-saved['sd']))),'gaussian_static_historical_maxabs':float(np.max(abs(g0-old['G0_Gaussian']))),'gaussian_response_historical_maxabs':float(np.max(abs(gC[indices]-old['C_Gaussian']))),'harmonic_static_historical_maxabs':float(np.max(abs(h['G0']-old['G0_harmonic']))),'harmonic_response_historical_maxabs':float(np.max(abs(h['C'][indices]-old['C_harmonic']))),'gaussian_Wick_maxabs':float(np.max(abs(gaussian_wick-np.concatenate([g0[None],gK])))),'harmonic_Wick_maxabs':float(np.max(abs(harmonic_wick-np.concatenate([h['G0'][None],h['K']])))),'independent_polynomial_energy_maxabs':float(np.max(abs(energies(m,points)[1]-polynomial_values(m,points))))}
    ok=max(checks.values())<=1e-10
    arrays={'times_over_tau':ts,'physical_times':physical,'G0_Gaussian':g0,'K_Gaussian':gK,'C_Gaussian':gC,'G0_harmonic':h['G0'],'K_harmonic':h['K'],'C_harmonic':h['C'],'Sigma_Gaussian':Sigma,'sd':h['sd'],'r0':m['r0'],'edges':m['edges'],'basis':m['basis'],'eigenvalues':m['eigenvalues'],'Gaussian_Wick_K_with_t0':gaussian_wick,'harmonic_Wick_K_with_t0':harmonic_wick}
    rec={'status':'PASS' if ok else 'FAIL','checks':checks,'independent_check_tolerance':1e-10,'model_changed_or_refitted':False,'normalization':'original exact whole-harmonic standard deviations','source_sha256':digest(__file__),'protocol_sha256':digest(ROOT/'protocol.json'),'Sigma_source_sha256':digest(ROOT/'inputs/frozen-gaussian-result.npz'),'seconds':time.perf_counter()-start,'environment':{'python':sys.version,'numpy':np.__version__,'scipy':scipy.__version__,'platform':platform.platform()}}
    return arrays,rec

if __name__=='__main__':
    z,r=evaluate();np.savez_compressed(ROOT/'results/controls.npz',**z);(ROOT/'results/controls-verification.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
    if r['status']!='PASS':raise RuntimeError('Frozen control check failed')
