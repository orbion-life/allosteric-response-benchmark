"""Independent Gaussian/harmonic control audit on the fixed nonlinear triangle.

No random sampling, refitting or model-source imports. Five Gauss-Hermite nodes
per standard-normal coordinate give 5**6 = 15,625 nodes for the joint initial
and innovation variables. Direct quartic contact energies make a product of
observables degree at most eight; this quadrature is exact for that polynomial
in exact arithmetic. Joint-time covariance uses the full saved Sigma through
its symmetric eigendecomposition, including off-diagonal correlations.

Run python3 verify_controls.py. Writes only beside this audit script.
"""
from pathlib import Path
import hashlib, itertools, json, platform, sys, time, warnings
import numpy as np
from numpy.polynomial.hermite import hermgauss
warnings.filterwarnings('error', category=RuntimeWarning)
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
SOURCE=ROOT/'nonlinear'
ORDER=5
TOLERANCE=1e-10

def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def mm(a,b):
    return np.einsum('ij,jk->ik',a,b,optimize=False)

def main():
    start=time.perf_counter()
    frozen=json.loads((SOURCE/'protocol-freeze.json').read_text())
    protocol=json.loads((SOURCE/'protocol.json').read_text())
    physical=json.loads((SOURCE/'source/vendor/protocol.json').read_text())
    assert physical['beta']==physical['mu']==physical['kappa']==1
    mismatches=[name for name,h in frozen['input_hashes'].items() if digest(SOURCE/name)!=h]
    assert not mismatches,mismatches
    assert digest(SOURCE/'protocol.json')==frozen['protocol_sha256']
    with np.load(SOURCE/'results/controls.npz') as src:
        z={k:src[k].copy() for k in src.files}
    with np.load(SOURCE/'inputs/frozen-gaussian-result.npz') as src:
        saved=src['Sigma'].copy()
    assert np.array_equal(saved,z['Sigma_Gaussian'])
    assert np.array_equal(z['times_over_tau'],protocol['all_times_sorted'])
    r0,edges,B,sd=z['r0'],z['edges'],z['basis'],z['sd']
    assert B.shape==(9,3) and edges.shape==(3,2)
    degrees=np.bincount(edges.ravel(),minlength=3)
    x,w=hermgauss(ORDER)
    indices=np.array(list(itertools.product(range(ORDER),repeat=6)))
    nodes=np.sqrt(2)*x[indices]
    weights=np.prod(w[indices]/np.sqrt(np.pi),axis=1)
    def energy(q):
        r=r0[None]+mm(q,B.T).reshape(-1,3,3)
        E=np.zeros((len(q),3))
        for i,j in edges:
            length2=np.sum((r0[i]-r0[j])**2)
            u=(np.sum((r[:,i]-r[:,j])**2,axis=1)-length2)**2/(8*length2)
            E[:,i]+=u/degrees[i]
            E[:,j]+=u/degrees[j]
        return E
    results=[];arrays={}
    for label,Sigma in [('Gaussian',z['Sigma_Gaussian']),('harmonic',np.diag(1/z['eigenvalues']))]:
        values,vectors=np.linalg.eigh(Sigma)
        assert values[0]>0
        root=vectors*np.sqrt(values)
        q0=mm(nodes[:,:3],root.T)
        e0=energy(q0);mean0=np.einsum('i,ij->j',weights,e0,optimize=False)
        f0=(e0-mean0)/sd
        G0=mm(f0.T,weights[:,None]*f0)
        allK=[];errors=[]
        for k,t in enumerate(z['physical_times']):
            decay=np.exp(-t/values)
            # q_t=exp(-t Sigma^-1)q_0+independent N(0,Sigma*(I-exp(-2t Sigma^-1))).
            qt=mm(nodes[:,:3],(root*decay).T)+mm(nodes[:,3:],(root*np.sqrt(-np.expm1(-2*t/values))).T)
            et=energy(qt);meant=np.einsum('i,ij->j',weights,et,optimize=False)
            ft=(et-meant)/sd
            K=mm(f0.T,weights[:,None]*ft)
            allK.append(K);errors.append(float(abs(K-z['K_'+label][k]).max()))
        allK=np.array(allK);C=allK-G0
        row={'model':label,'G0_max_abs_difference':float(abs(G0-z['G0_'+label]).max()),
             'K_max_abs_difference':max(errors),'C_max_abs_difference':float(abs(C-z['C_'+label]).max()),
             'per_time_K_max_abs_difference':errors,'minimum_Sigma_eigenvalue':float(values[0])}
        results.append(row)
        arrays.update({label+'_G0':G0,label+'_K':allK,label+'_C':C})
    passed=all(max(row[k] for k in ['G0_max_abs_difference','K_max_abs_difference','C_max_abs_difference'])<=TOLERANCE for row in results)
    np.savez_compressed(HERE/'independent-control-moments.npz',**arrays,times_over_tau=z['times_over_tau'])
    receipt={'status':'PASS' if passed else 'FAIL','scope':'Fixed three-coordinate Gaussian and harmonic control evaluation only; no original nonlinear fidelity or protein claim',
             'method':'Direct contact energy evaluation with six-dimensional product Gauss-Hermite quadrature',
             'quadrature_order_per_variable':ORDER,'independent_standard_normal_variables':6,'quadrature_points':len(weights),
             'maximum_polynomial_degree_of_integrand':8,'weight_sum':float(weights.sum()),
             'joint_covariance':'Full Sigma eigensystem; q0 covariance Sigma, qt covariance Sigma, cross-covariance exp(-t Sigma^-1)Sigma',
             'randomness':'None; deterministic quadrature; no seed is used','fixed_beta':1,'fixed_mu':1,'fixed_kappa':1,
             'times_over_tau':z['times_over_tau'].tolist(),'original_times_retained':all(t in z['times_over_tau'] for t in protocol['original_times']),
             'saved_covariance_unchanged':True,'sealed_input_mismatches':mismatches,'tolerance':TOLERANCE,'results':results,
             'source_hashes':{'audit_script':digest(__file__),'protocol.json':digest(SOURCE/'protocol.json'),
                             'protocol-freeze.json':digest(SOURCE/'protocol-freeze.json'),'controls.npz':digest(SOURCE/'results/controls.npz'),
                             'saved_Gaussian_result':digest(SOURCE/'inputs/frozen-gaussian-result.npz'),
                             'independent-control-moments.npz':digest(HERE/'independent-control-moments.npz')},
             'verified_protocol_input_hashes':frozen['input_hashes'],
             'environment':{'python':sys.version,'numpy':np.__version__,'platform':platform.platform()},'seconds':time.perf_counter()-start}
    (HERE/'controls-check.json').write_text(json.dumps(receipt,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:receipt[k] for k in ['status','quadrature_points','results','seconds']},indent=2))
    assert passed

if __name__=='__main__':main()
