"""Check saved ABL scope, response arithmetic and declared template costs."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[key]='2'
from pathlib import Path
import math,json,hashlib
import numpy as np
from scipy.linalg import cholesky,solve
from scipy.special import ive
from threadpoolctl import threadpool_limits
HERE=Path(__file__).resolve().parent;BASE=HERE.parents[1];OP=BASE/'operator'
def load(p):
    with np.load(p) as z:return {k:z[k].copy() for k in z.files}
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    protocol=json.loads((OP/'protocol.json').read_text())
    for name,spec in protocol['input']['files'].items():assert sha(OP/'inputs'/name)==spec['sha256']
    saved=load(OP/'inputs/operator.npz');sigma=load(OP/'inputs/fitted-covariance.npz')['covariance'];factor=load(OP/'results/physical-factor.npz');responses=load(OP/'results/matched-responses.npz')
    classical=json.loads((OP/'results/classical-summary.json').read_text());resources=json.loads((OP/'results/hermite-access.json').read_text());cost=json.loads((OP/'results/complete-cost-accounting.json').read_text())
    N=saved['B'].shape[1];rank=len(saved['Hr']);d=len(sigma);dimensions=math.comb(d+4,4)-1;rows=1<<(dimensions-1).bit_length();columns=1<<(rank-1).bit_length()
    dimensions_check={'residues':N,'coordinates':d,'saved_reduced_rank':rank,'Hermite_dimension':dimensions,'explicit_Hermite_by_rank_factor_bytes':8*dimensions*rank,'padded_generic_table_bytes':8*rows*columns,'padded_generic_table_loader_qubits':int(math.log2(rows)+math.log2(columns)+1),'dense_saved_Hr_factor_bytes':8*rank*rank,'dense_saved_Hr_loader_qubits':2*int(math.log2(columns))+1}
    assert dimensions_check['Hermite_dimension']==resources['full_centered_Hermite_dimension']
    assert dimensions_check['explicit_Hermite_by_rank_factor_bytes']==resources['full_factor_float64_bytes']
    assert dimensions_check['padded_generic_table_loader_qubits']==resources['generic_table_loader_qubits']
    stationarity=float(np.max(abs(factor['Gamma']-solve(sigma,np.eye(d),assume_a='pos'))))
    inverse_identity=float(np.max(abs(factor['Gamma']@factor['Sigma_physical']-np.eye(d))))
    errors={'force_Gamma_vs_inverse_covariance_maxabs':stationarity,'factor_covariance_inverse_identity_maxabs':inverse_identity,'reduced_vs_full_Gaussian_maxabs':float(np.max(abs(responses['reduced_spectral']-responses['full_Gaussian']))),'matrix_free_vs_reduced_maxabs':float(np.max(abs(responses['reduced_spectral']-responses['reduced_matrix_free']))),'factor_consistent_vs_original_Gaussian_maxabs':float(np.max(abs(responses['factor_consistent_Gaussian']-responses['full_Gaussian'])))}
    assert errors['reduced_vs_full_Gaussian_maxabs']==classical['reduction_error_maxabs_vs_full_Gaussian']
    assert errors['matrix_free_vs_reduced_maxabs']==classical['matrix_free_error_maxabs']
    assert errors['factor_consistent_vs_original_Gaussian_maxabs']==classical['factor_covariance_response_change_maxabs']
    assert stationarity<1e-10 and inverse_identity<1e-10
    C=cholesky(saved['Hr'],lower=False);gamma=float(np.sum(C*C));assert np.isclose(gamma,cost['dense_reduced_fallback_factor']['gamma_Frobenius_squared'],rtol=1e-14,atol=1e-12)
    gamma_error=abs(gamma-float(np.trace(saved['Hr'])))
    Lmax=float(np.diag(saved['G0']).max());eps=cost['error_budget_planning_only']['statistical_response_estimation'];Q=N*(N+1)//2*4
    single=math.ceil(2*Lmax**2*math.log(2/.05)/eps**2);simultaneous=math.ceil(2*Lmax**2*math.log(2*Q/.05)/eps**2)
    assert single==cost['readout_planning_only']['independent_shots_single_query_95pct_Hoeffding']
    assert simultaneous==cost['readout_planning_only']['per_query_shots_95pct_simultaneous_symmetric_matrices']
    assert simultaneous*Q==cost['readout_planning_only']['total_independent_shots_95pct_simultaneous_symmetric_matrices']
    planning=[]
    for row in cost['full_ABL_dense_factor_query_plan']:
        a=gamma*row['time']/2;coeff=np.array([ive(0,a)]+[2*ive(k,a) for k in range(1,8193)])
        degree=int(np.flatnonzero(Lmax*np.maximum(0,1-np.cumsum(coeff))<=.0005)[0]);walks=(1<<degree.bit_length())-1
        assert degree==row['Chebyshev_degree_for_planning_budget'] and walks==row['binary_SELECT_controlled_walk_calls']
        assert 2*walks*columns**2==row['generic_Ry_angles_in_SELECT_lower_component_count']
        planning.append({'time_over_tau':row['time_over_tau'],'degree':degree,'controlled_walk_calls':walks})
    source_paths=[Path(__file__),OP/'protocol.json',OP/'results/complete-cost-accounting.json',OP/'results/classical-summary.json',OP/'results/matched-responses.npz',OP/'results/physical-factor.npz',OP/'results/hermite-access.json',OP/'source/run_costs.py']
    record={'status':'PASS','scope':'Direct raw-array arithmetic and the declared unoptimized template planning formulas; no independent timing benchmark or full quantum execution.','maximum_threads':2,'dimensions':dimensions_check,'array_checks':errors,'independent_dense_Cholesky_Gram_maxabs':float(np.max(abs(C.T@C-saved['Hr']))),'independent_dense_factor_gamma_minus_trace_Hr':gamma_error,'independent_shot_counts':{'single_95pct_sufficient':single,'per_entry_for_joint95pct':simultaneous,'entries':Q,'total':Q*simultaneous},'planning_degrees_verified':planning,'source_sha256':{str(p.relative_to(BASE)):sha(p) for p in source_paths}}
    (HERE/'operator-scope-and-cost-audit.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps({k:v for k,v in record.items() if k!='source_sha256'},indent=2))
if __name__=='__main__':
    with threadpool_limits(limits=2):main()
