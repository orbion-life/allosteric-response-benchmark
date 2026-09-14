"""Measured shared preprocessing and explicitly hypothetical full-query costs."""
from core import *
from scipy.linalg import cholesky
from scipy.special import ive
from threadpoolctl import threadpool_limits
import resource

def main():
    out=ROOT/'results';start=time.perf_counter();model=dict(np.load(ROOT/'inputs/static-model.npz'));saved=dict(np.load(ROOT/'inputs/operator.npz'));Sigma=np.load(ROOT/'inputs/fitted-covariance.npz')['covariance'];N=len(model['r0']);d=len(Sigma);r=len(saved['Hr'])
    st=time.perf_counter();K=np.zeros((3*N,3*N))
    for i,j in model['edges']:
        a=model['r0'][i]-model['r0'][j];block=np.outer(a,a)/np.dot(a,a)
        for u,su in [(i,1),(j,-1)]:
            for v,sv in [(i,1),(j,-1)]:K[3*u:3*u+3,3*v:3*v+3]+=su*sv*block
    assemble=time.perf_counter()-st;st=time.perf_counter();lam,V=eigh(K);eigentime=time.perf_counter()-st
    # The frozen protocol keeps original B0/eigenvalue thresholds: compare its
    # operator restriction, not an arbitrary newly rotated degenerate basis.
    basis_residual=float(np.max(abs(mm(model['B'].T,mm(K,model['B']))-np.diag(model['eigenvalues']))))
    st=time.perf_counter();obj=gp.CovarianceObjective(model['r0'],model['edges'],model['B'],model['eigenvalues']);setup=time.perf_counter()-st
    def progress(x):print('Gaussian refit',x,flush=True)
    st=time.perf_counter();Y,fit=gp.newton_cg(obj,progress=progress);fittime=time.perf_counter()-st
    Snew=Y/np.sqrt(np.outer(model['eigenvalues'],model['eigenvalues']));delta=float(np.max(abs(Snew-Sigma)))
    st=time.perf_counter();C=cholesky(saved['Hr'],lower=False);choltime=time.perf_counter()-st;cholerror=float(np.max(abs(mm(C.T,C)-saved['Hr'])));gamma=float(np.sum(C*C));rp=2**math.ceil(math.log2(r));norms=np.sqrt(np.diag(saved['G0']));Lmax=float(max(norms)**2)
    cost={'shared_classical_preprocessing':{'harmonic_cartesian_assembly_seconds':assemble,'dense_cartesian_harmonic_eigh_seconds':eigentime,'frozen_basis_harmonic_identity_maxabs':basis_residual,'Gaussian_fit_objective_setup_seconds':setup,'Gaussian_refit_seconds':fittime,'Gaussian_refit':fit,'fitted_covariance_replay_maxabs':delta,'frozen_B0_kept_for_deterministic_access':True,'covariance_fit_common_to_classical_and_quantum_Gaussian_routes':True,'original_snapshot_and_whitening_construction':'Precomputed common input; NOT retimed in this run. It requires Gaussian kernels and dense whitening. Prior campaign provenance is retained; query timings alone are not end-to-end model construction.'},
      'dense_reduced_fallback_factor':{'scope':'A measured dense Cholesky fallback, not a cheap physical force oracle. It explicitly constructs the saved reduced operator first.','seconds':choltime,'factor_float64_bytes':C.nbytes,'factor_Gram_maxabs':cholerror,'gamma_Frobenius_squared':gamma,'gamma_over_Hr_spectral_norm':gamma/float(saved['lambda_'][-1]),'factor_shape':list(C.shape),'padded_dimension':rp,'generic_loader_qubits':2*int(math.log2(rp))+1,'generic_loader_Ry_angles_per_UL':rp*rp,'estimated_native_operations_per_UL_before_optimization':2*rp*rp-3,'generic_loader_statevector_complex128_bytes':16*2**(2*int(math.log2(rp))+1),'status':'NOT_COMPILED_EXCEEDS_FROZEN_QUBIT_AND_GATE_CAP','observable_state_loading_Ry_angles_per_state':rp-1,'observable_state_loading_notes':'Generic real binary tree for one selected column of B, norm calculated classically; neither loading nor dense preprocessing is free.'},
      'error_budget_planning_only':{'reduction':5e-5,'force_and_covariance_difference':1e-8,'semigroup_approximation':5e-4,'combined_two_state_preparation_and_static_term':2.5e-4,'statistical_response_estimation':1e-3,'sum':.00180001,'target':.002,'scope':'A hypothetical engineering budget, not a passed ABL quantum error certificate. Reduction and physical covariance checks only measured at frozen4times; full coherent state/circuit implementation absent.'}}
    times=np.array([.001,.03,.1,1])*float(saved['tau']);estimate=[]
    for t in times:
        a=gamma*t/2;coeff=np.array([ive(0,a)]+[2*ive(k,a) for k in range(1,8193)]);tail=np.maximum(0,1-np.cumsum(coeff));ok=np.flatnonzero(Lmax*tail<=5e-4);degree=int(ok[0]) if len(ok) else None;selectcalls=(2**math.ceil(math.log2(degree+1))-1) if degree is not None else None
        estimate.append({'time':float(t),'time_over_tau':float(t/float(saved['tau'])),'gamma_t_over2':a,'Chebyshev_degree_for_planning_budget':degree,'binary_SELECT_controlled_walk_calls':selectcalls,'tail_response_bound':float(Lmax*tail[degree]) if degree is not None else None,'generic_UL_invocations_in_SELECT':2*selectcalls if selectcalls is not None else None,'generic_Ry_angles_in_SELECT_lower_component_count':int(2*selectcalls*rp*rp) if selectcalls is not None else None,'count_exclusions':['additional controls and reflections','observable preparations','coefficient loading','routing and fault tolerance','repeat measurements'],'scope':'Extrapolated exact generic-loader formula; no full ABL circuit compiled.'})
    Q=N*N*len(times);Qs=N*(N+1)//2*len(times);eps=.001
    shots=lambda delta:math.ceil(2*Lmax*Lmax*math.log(2/delta)/(eps*eps))
    cost['full_ABL_dense_factor_query_plan']=estimate;cost['readout_planning_only']={'max_observable_norm_product':Lmax,'response_relation':'C=L*alpha*(2p_plus-1)-G0; alpha<=1, use alpha=1 for conservative shots','independent_shots_single_query_95pct_Hoeffding':shots(.05),'full_matrix_four_times_entries':Q,'symmetry_distinct_four_times_entries':Qs,'per_query_shots_95pct_simultaneous_symmetric_matrices':shots(.05/Qs),'total_independent_shots_95pct_simultaneous_symmetric_matrices':Qs*shots(.05/Qs),'assumptions':'Exact independent Bernoulli outputs, no bias/noise, allocated absolute response statistical error. Conservative sufficient bound, not a minimum or measured QPU cost. Adaptive amplitude estimation is a separate unimplemented coherent protocol.','cached_classical_query_results_are_also_reused':True}
    cost['wall_seconds']=time.perf_counter()-start;cost['peak_RSS_bytes']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss;dump(out/'complete-cost-accounting.json',cost);print(json.dumps(cost,indent=2),flush=True)

if __name__=='__main__':
    with threadpool_limits(limits=4):main()
