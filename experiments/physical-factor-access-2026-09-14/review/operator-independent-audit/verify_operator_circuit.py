"""Independent raw-array and QPY audit; imports no production/operator modules."""
import os
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[name]='2'
from pathlib import Path
import itertools, math, json, hashlib, time, sys
import numpy as np
from numpy.polynomial.hermite import hermgauss
from scipy.linalg import expm
from scipy.special import ive
from qiskit import QuantumCircuit, qpy
from qiskit.quantum_info import Pauli
from qiskit_aer import AerSimulator
from threadpoolctl import threadpool_limits, threadpool_info

HERE=Path(__file__).resolve().parent
BASE=HERE.parents[1]
OP=BASE/'operator'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load(p):
    with np.load(p) as z:return {k:z[k].copy() for k in z.files}
def occupation_counts(row,d):return np.bincount(row[row>=0],minlength=d)
def circuit(p):
    with p.open('rb') as f:return qpy.load(f)[0]
def gate_counts(q):return {'qubits':q.num_qubits,'operations':q.size(),'depth':q.depth(),'gate_counts':{str(k):int(v) for k,v in q.count_ops().items()}}

def main():
    start=time.perf_counter()
    s=load(OP/'results/small-physical-model.npz');f=load(OP/'results/small-physical-factor.npz')
    summary=json.loads((OP/'results/quantum-summary.json').read_text())
    full=load(OP/'inputs/static-model.npz');rank=load(OP/'inputs/operator.npz')['Hr'].shape[0]
    d=s['Sigma'].shape[0];alphas=np.array([occupation_counts(x,d) for x in s['occupations']])
    assert d==2 and len(alphas)==14 and rank==1261
    coordinates_match=np.array_equal(s['r0'],full['r0'][:3])
    contacts_match=np.array_equal(s['edges'],full['edges'][(full['edges']<3).all(axis=1)])
    assert coordinates_match and contacts_match and s['canonical'].tolist()==[242,243,244]
    # Direct original quartic energies integrated with a 5x5 GH product rule.
    x,w=hermgauss(5);indexes=np.array(list(itertools.product(range(5),repeat=d)))
    nodes=np.sqrt(2)*x[indexes];weights=np.prod(w[indexes]/np.sqrt(np.pi),axis=1)
    root=s['covariance_vectors']*np.sqrt(s['covariance_values'])
    q=nodes@root.T;r=s['r0'][None,:,:]+(q@s['B'].T).reshape(len(q),3,3)
    degree=np.bincount(s['edges'].ravel(),minlength=3);energies=np.zeros((len(q),3))
    for i,j in s['edges']:
        rest2=np.sum((s['r0'][i]-s['r0'][j])**2)
        energy=(np.sum((r[:,i]-r[:,j])**2,axis=1)-rest2)**2/(8*rest2)
        for residue in (i,j):energies[:,residue]+=energy/degree[residue]
    centered=(energies-weights@energies)/s['sd']
    hs=np.ones((5,len(nodes),d));hs[1]=nodes
    for k in range(1,4):hs[k+1]=nodes*hs[k]-k*hs[k-1]
    basis=np.array([np.prod(np.array([hs[k,:,j]/math.sqrt(math.factorial(k)) for j,k in enumerate(alpha)]),axis=0) for alpha in alphas]).T
    h=basis.T@(weights[:,None]*centered)
    coefficient_error=float(np.max(abs(h-s['h'])))
    reconstruction_error=float(np.max(abs(basis@h-centered)))
    # Independent second-quantized coordinate drift in the complete sectors.
    G=s['covariance_vectors'].T@s['force_Gamma']@s['covariance_vectors']
    H=np.zeros((len(alphas),len(alphas)))
    for a,alpha in enumerate(alphas):
        for b,beta in enumerate(alphas):
            for j in range(d):
                for k in range(d):
                    if alpha[j] and beta[k] and np.array_equal(alpha-np.eye(d,dtype=int)[j],beta-np.eye(d,dtype=int)[k]):
                        H[a,b]+=G[j,k]*math.sqrt(alpha[j]*beta[k])
    factor_gram_error=float(np.max(abs(f['D'].T@f['D']-H)))
    covariance_generator_error=float(np.max(abs(H-np.diag(s['omega']))))
    force_factor_error=float(np.max(abs(s['force_factor'].T@s['force_factor']-s['force_Gamma'])))
    gamma=float(np.sum(f['D']*f['D']));gamma_error=abs(gamma-float(f['gamma']))
    assert coefficient_error<1e-12 and reconstruction_error<1e-11 and factor_gram_error<1e-11 and force_factor_error<1e-12
    assert covariance_generator_error<1e-10 and gamma_error<1e-12
    # Replay all physical factor columns and compare projected amplitudes.
    loader=circuit(OP/'results/physical-factor-loader.qpy');walk=circuit(OP/'results/physical-factor-walk.qpy');controlled=circuit(OP/'results/physical-factor-controlled-walk.qpy')
    component_counts={name:gate_counts(qc) for name,qc in [('loader',loader),('walk',walk),('controlled_walk',controlled)]}
    for name,counts in component_counts.items():
        for key,value in counts.items():assert summary[name][key]==value,(name,key)
    ncol=math.ceil(math.log2(f['D'].shape[1]));nrow=math.ceil(math.log2(f['D'].shape[0]));dc,dr=2**ncol,2**nrow
    sim=AerSimulator(method='statevector',device='CPU',precision='double',max_parallel_threads=2,max_parallel_experiments=1)
    projected_factor=[];projected_walk=[]
    for column in range(len(alphas)):
        init=QuantumCircuit(loader.num_qubits)
        for bit in range(ncol):
            if column&(1<<bit):init.x(bit)
        for template,output,index in [(loader,projected_factor,np.arange(dr)*dc),(walk,projected_walk,np.arange(dc))]:
            qc=init.compose(template);qc.save_statevector()
            result=sim.run(qc).result();assert result.success
            output.append(np.asarray(result.data(0)['statevector'])[index])
    projected_factor=np.array(projected_factor).T;projected_walk=np.array(projected_walk).T
    padded=np.zeros((dr,len(alphas)));padded[:len(f['D'])]=f['D']/math.sqrt(gamma)
    expected_walk=np.eye(dc)[:,:len(alphas)];expected_walk[:len(alphas)]-=2*H/gamma
    loader_error=float(np.max(abs(projected_factor-padded)));walk_error=float(np.max(abs(projected_walk-expected_walk)))
    assert loader_error<1e-12 and walk_error<1e-12
    # Independently evaluate the truncated Chebyshev series and replay its QPY.
    L=float(np.linalg.norm(h[:,0])*np.linalg.norm(h[:,1]));static=float(h[:,0]@h[:,1]);tau=float(s['tau']);Xoperator=np.eye(14)-2*H/gamma
    responses=[]
    for recorded in summary['LCU_responses']:
        ratio=recorded['time_over_tau'];t=ratio*tau;a=gamma*t/2
        c=np.array([ive(0,a)]+[2*ive(k,a) for k in range(1,33)])
        tails=np.maximum(0,1-np.cumsum(c));degree=int(np.flatnonzero(L*tails<=2.5e-7)[0])
        assert degree==recorded['chosen_degree']
        alpha=float(c[:degree+1].sum());T0=np.eye(14);P=c[0]*T0
        if degree:
            T1=Xoperator.copy();P+=c[1]*T1
            for k in range(2,degree+1):T2=2*Xoperator@T1-T0;P+=c[k]*T2;T0,T1=T1,T2
        direct=float(h[:,0]@expm(-t*H)@h[:,1]-static)
        fitted=float(np.sum(h[:,0]*np.exp(-t*s['omega'])*h[:,1])-static)
        truncated=float(h[:,0]@P@h[:,1]-static)
        qc=circuit(OP/f'results/physical-response-{ratio}.qpy');counts=gate_counts(qc)
        for key,value in counts.items():assert recorded['native_circuit'][key]==value
        row={'time_over_tau':ratio,'native_circuit':counts,'chosen_degree':degree,'norm_L':L,'static_subtraction':static,'LCU_alpha':alpha,'independent_factor_response':direct,'independent_fitted_Gaussian_response':fitted,'independent_truncated_polynomial_response':truncated,'tail_response_bound':float(L*(1-alpha)),'truncation_actual_error':abs(truncated-direct),'frozen_resource_cap_pass':counts['operations']<=200000 and counts['qubits']<=20}
        assert abs(recorded['norm_restoration_L']-L)<1e-12 and abs(recorded['static_subtraction']-static)<1e-12
        assert abs(recorded['response_reference_factor']-direct)<1e-12 and abs(recorded['response_reference_fitted_covariance']-fitted)<1e-12
        if row['frozen_resource_cap_pass']:
            qc.save_expectation_value(Pauli('X'),[0],label='independent_X')
            st=time.perf_counter();result=sim.run(qc).result();assert result.success
            X=float(np.real(result.data(0)['independent_X']));restored=L*alpha*X-static
            row.update({'ideal_X':X,'ideal_p_plus':(1+X)/2,'norm_restored_QPY_response':restored,'QPY_vs_truncated_polynomial_error':abs(restored-truncated),'QPY_vs_fitted_Gaussian_error':abs(restored-fitted),'response_replay_seconds':time.perf_counter()-st})
            assert abs(restored-truncated)<1e-11 and abs(restored-fitted)<1e-6
        else:
            row['replay']='Not simulated because it exceeds the frozen native-operation cap.'
            assert recorded['status']=='FAIL_FROZEN_COMPILED_RESOURCE_CAP'
        responses.append(row)
    old=BASE/'review/independent-math-audit'
    mapping={'operator-Hermite-comparison.json':{'script':old/'compare_operator_hermite.py','independent_helper':old/'hermite_quadrature_check.py','operator_source':OP/'source/core.py','controls':BASE/'nonlinear/results/controls.npz','snapshot_model':BASE/'nonlinear/inputs/frozen-gaussian-result.npz'},'small-Hermite-factor-check.json':{'script':old/'hermite_quadrature_check.py','controls':BASE/'nonlinear/results/controls.npz','snapshot_model':BASE/'nonlinear/inputs/frozen-gaussian-result.npz','output':old/'small-Hermite-factor.npz'}}
    receipts={}
    for filename,paths in mapping.items():
        prior=json.loads((old/filename).read_text());checks={k:sha(p)==prior['source_hashes'][k] for k,p in paths.items()}
        receipts[filename]={'status':prior['status'],'current_source_hash_matches':checks,'scope':prior['scope']}
        assert all(checks.values())
    paths=[Path(__file__),OP/'source/core.py',OP/'source/run_classical.py',OP/'source/run_quantum.py',OP/'vendor/multiplexed.py',OP/'results/quantum-summary.json',OP/'results/small-physical-model.npz',OP/'results/small-physical-factor.npz',*sorted((OP/'results').glob('*.qpy'))]
    receipt={'status':'PASS','scope':'Independent direct quartic-energy quadrature, complete14-term small ABL polynomial generator, saved component/response QPY replay and prior audit provenance. No full-rank1261 quantum oracle or nonlinear validation claim.','threads_maximum':2,'source_coordinates_and_contacts_match':{'coordinates':coordinates_match,'contacts':contacts_match,'canonical':s['canonical'].tolist()},'small_model':{'residues':3,'coordinates':2,'Hermite_states':14,'factor_rows':len(f['D']),'full_ABL_saved_response_rank':rank},'checks':{'direct_GH_coefficients_maxabs':coefficient_error,'direct_contact_reconstruction_maxabs':reconstruction_error,'physical_force_Gram_maxabs':force_factor_error,'polynomial_factor_Gram_vs_independent_ladder_generator_maxabs':factor_gram_error,'polynomial_generator_vs_covariance_rates_maxabs':covariance_generator_error,'factor_Frobenius_gamma_error':gamma_error,'QPY_all14_factor_columns_maxabs':loader_error,'QPY_all14_walk_columns_maxabs':walk_error},'component_circuits':component_counts,'responses':responses,'prior_independent_audit_receipts':receipts,'current_source_sha256':{str(path.relative_to(BASE)):sha(path) for path in paths},'software':{'python':sys.version,'numpy':np.__version__},'threadpools':threadpool_info(),'wall_seconds':time.perf_counter()-start}
    (HERE/'operator-circuit-audit.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:receipt[k] for k in ['status','checks','responses','wall_seconds']},indent=2),flush=True)
if __name__=='__main__':
    with threadpool_limits(limits=2):main()
