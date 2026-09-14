"""Frozen ABL direct-factor diagnostic. Run from any working directory."""
from core import *
import platform, resource, datetime, importlib.metadata
from threadpoolctl import threadpool_limits, threadpool_info

def timed(fn,repeats=3):
    fn();ts=[]
    for _ in range(repeats):
        st=time.perf_counter();value=fn();ts.append(time.perf_counter()-st)
    return value,{'seconds':ts,'median_seconds':float(np.median(ts)),'warm_call_excluded':True}

def main():
    out=ROOT/'results';start=time.perf_counter();protocol=json.loads((ROOT/'protocol.json').read_text())
    for name,spec in protocol['input']['files'].items():
        assert sha(ROOT/'inputs'/name)==spec['sha256'],name
    model=dict(np.load(ROOT/'inputs/static-model.npz'));Sigma=np.load(ROOT/'inputs/fitted-covariance.npz')['covariance'];sd=np.load(ROOT/'inputs/all-mode-harmonic-scales.npz')['harmonic_sd'];saved=dict(np.load(ROOT/'inputs/operator.npz'))
    N=len(sd);d=len(Sigma);rank=saved['Hr'].shape[0];tau=float(saved['tau']);times=tau*np.array(protocol['outputs']['times_over_tau'])
    candidate=np.flatnonzero(model['candidate']);candidate=candidate[np.argsort(model['canonical'][candidate])][:4]
    receiver=model['receiver'][np.argsort(model['canonical'][model['receiver']])][:2];pairs=np.array([(i,j) for i in candidate for j in receiver])
    records={'scope':'Actual ABL fixed Gaussian response; no original-nonlinear fidelity or quantum advantage inference','N':N,'d':d,'rank':rank,'tau':tau,'times':times.tolist(),'times_over_tau':(times/tau).tolist(),'pairs_zero_based':pairs.tolist(),'pairs_canonical':model['canonical'][pairs].tolist(),'fixed_sd_sha256':sha(ROOT/'inputs/all-mode-harmonic-scales.npz')}
    print('ABL force factor',flush=True);fac=force_factor(model,Sigma);records['physical_factor']=fac['stats'];dump(out/'physical-factor.json',records)
    np.savez_compressed(out/'physical-factor.npz',blocks=fac['blocks'],Gamma=fac['Gamma'],Sigma_physical=fac['Sigma_physical'])
    print('ABL selected Hermite coefficient rows',flush=True);hermite=Hermite(model,Sigma,sd);rng=np.random.default_rng(20260914);selected=[]
    for k in range(1,5):
        rows=[tuple([0]*k),tuple([d-1]*k)]
        while len(rows)<16:
            ix=tuple(sorted(rng.integers(0,d,size=k).tolist()))
            if ix not in rows:rows.append(ix)
        selected.extend(rows)
    nodes=tau*np.array([0,.0001,.001,.01,.1,1]);Hs=[];Ps=[];Fs=[];omegas=[];row_times=[]
    for ix in selected:
        st=time.perf_counter();h=hermite.coeff(ix);psi,f,om=hermite.row(ix,saved['W'],nodes);row_times.append(time.perf_counter()-st);Hs.append(h);Ps.append(psi);Fs.append(f);omegas.append(om)
    Dfull=math.comb(d+4,4)-1;padded=2**math.ceil(math.log2(Dfull));columnpad=2**math.ceil(math.log2(rank));factorbytes=8*Dfull*rank
    norms2=np.diag(saved['G0']);samplefactor=np.asarray(Fs)
    resources={'degree_counts':{str(k):math.comb(d+k-1,k) for k in range(1,5)},'full_centered_Hermite_dimension':Dfull,'full_observable_coefficient_float64_bytes':8*Dfull*N,'full_snapshot_embedding_float64_bytes':factorbytes,'full_factor_float64_bytes':factorbytes,'dense_covariance_bytes':Sigma.nbytes,'dense_B0_bytes':model['B'].nbytes,'dense_whitening_bytes':saved['W'].nbytes,'Hermite_edge_transform_float64_bytes':hermite.A.nbytes,'selected_rows':len(selected),'selected_coefficient_and_embedding_seconds':row_times,'coefficient_setup_seconds_including_covariance_eigh':hermite.setup_seconds,'selected_rows_are_not_full_identity_verification':True,'sequential_full_row_time_extrapolation_seconds':float(np.median(row_times)*Dfull),'sequential_full_row_time_extrapolation_scope':'Empirical 64-row median including duplicate coefficient call and dense W multiply; no asymptotic lower bound','padded_factor_rows':padded,'padded_factor_columns':columnpad,'generic_table_loader_qubits':int(math.log2(padded)+math.log2(columnpad)+1),'generic_loader_Ry_angles_per_UL':int(padded-1+padded*(columnpad-1)+1),'generic_loader_table_float64_bytes':8*padded*columnpad,'actual_saved_reduced_factor_Frobenius_squared':float(np.trace(saved['Hr'])),'normalization_scope':'Trace(Hr) equals exact full sqrt(omega)Psi Frobenius squared if the embedding identity holds; full ABL embedding identity not checked here. It does not construct coherent row norms.','selected_factor_rows_squared_norm':float(np.sum(samplefactor*samplefactor)),'gate':'FAIL_GENERIC_TABLE_ACCESS','gate_reasons':['The explicit full table exceeds the 8GiB memory budget','The fixed generic loader exceeds the 20-qubit statevector budget','A compact coherent Hermite row-norm and whitening/state preparation circuit is not implemented'],'response_observable_norm_squared_minmax':[float(norms2.min()),float(norms2.max())],'full_output_scalar_count_at_four_times':N*N*len(times),'symmetric_scalar_count_at_four_times':N*(N+1)//2*len(times),'polynomial_norm_bound_from_coordinate_drift':4*float(eigh(fac['Gamma'],eigvals_only=True)[-1])}
    records['Hermite_access_resources']=resources;dump(out/'hermite-access.json',resources)
    np.savez_compressed(out/'selected-hermite-rows.npz',indices=np.array([list(ix)+[-1]*(4-len(ix)) for ix in selected]),coefficients=Hs,embedding=Ps,factor_rows=Fs,omega=omegas,nodes=nodes)
    print('ABL matched classical baselines',flush=True)
    st=time.perf_counter();lam,Q=eigh(saved['Hr']);weights=mm(Q.T,saved['B']);setup=time.perf_counter()-st;G0=saved['G0'];B=saved['B']
    def spectral_full(t):return mm(weights.T,np.exp(-t*lam)[:,None]*weights)-G0
    def spectral_pairs(t,ps):return np.array([float(np.sum(np.exp(-t*lam)*weights[:,i]*weights[:,j])-G0[i,j]) for i,j in ps])
    spec={};Cs=[]
    for t in times:
        val,full=timed(lambda:spectral_full(t));Cs.append(val);_,single=timed(lambda:spectral_pairs(t,pairs[:1]));_,few=timed(lambda:spectral_pairs(t,pairs));spec[str(t)]={'single':single,'eight':few,'full':full}
    Cs=np.array(Cs);records['cached_spectral']={'fresh_eigh_and_all_observable_projection_seconds':setup,'by_time':spec}
    def hact(x):return mm(saved['Hr'],np.asarray(x))
    Hop=LinearOperator((rank,rank),matvec=hact,rmatvec=hact,matmat=hact,dtype=float);traceH=float(np.trace(saved['Hr']));actions=[];actionC=[]
    for ti,t in enumerate(times):
        row={'time':float(t)}
        for label,js in [('single',[int(pairs[0,1])]),('eight',sorted(set(pairs[:,1].tolist()))),('full',list(range(N)))]:
            st=time.perf_counter();E=expm_multiply(-t*Hop,B[:,js],traceA=-t*traceH);dt=time.perf_counter()-st
            if label=='full':val=mm(B.T,E)-G0;actionC.append(val);err=np.max(abs(val-Cs[ti]))
            else:
                ps=pairs[:1] if label=='single' else pairs
                val=np.array([float(mm(B[:,i],E[:,js.index(j)])-G0[i,j]) for i,j in ps]);err=np.max(abs(val-spectral_pairs(t,ps)))
            row[label]={'seconds':dt,'max_abs_error_vs_cached_spectral':float(err),'right_hand_sides':len(js)}
        actions.append(row);print('matrix-free action',t/tau,flush=True)
    records['matrix_free_saved_Hr']=actions;dump(out/'classical-progress.json',records)
    st=time.perf_counter();kernel=OUKernel(model,Sigma,sd);kernel_setup=time.perf_counter()-st
    st=time.perf_counter();K0=kernel.kernel(0,False);zero_seconds=time.perf_counter()-st
    analytic=[];refs=[]
    for ti,t in enumerate(times):
        row={'time':float(t)}
        for label,ps in [('single',pairs[:1]),('eight',pairs)]:
            def calc():return np.array([pair_kernel(kernel,t,int(i),int(j))-pair_kernel(kernel,0,int(i),int(j)) for i,j in ps])
            val,timing=timed(calc);row[label]=timing;row[label]['max_abs_vs_reduced_response']=float(np.max(abs(val-spectral_pairs(t,ps))))
        st=time.perf_counter();ref=kernel.kernel(t,False)-K0;dt=time.perf_counter()-st;refs.append(ref)
        row['full']={'seconds':dt,'max_abs_reduced_response_error':float(np.max(abs(ref-Cs[ti])))};analytic.append(row);print('analytic Gaussian',t/tau,flush=True)
    refs=np.array(refs);records['analytic_Gaussian']={'covariance_eigh_and_kernel_setup_seconds':kernel_setup,'shared_zero_time_full_kernel_seconds':zero_seconds,'by_time':analytic,'normalization':'same fixed original harmonic SD as saved Hr; t=0 scalar reference recomputed in timed single/few queries'}
    print('factor-derived drift action and covariance shift',flush=True);Dop=drift_action(fac,model);traceD=float(np.trace(fac['Gamma']));drift_results=[];physCs=[]
    st=time.perf_counter();pkernel=OUKernel(model,fac['Sigma_physical'],sd);physical_setup=time.perf_counter()-st
    st=time.perf_counter();pK0=pkernel.kernel(0,False);physical_zero=time.perf_counter()-st
    # The factor matvec is tested on one deterministic physical covariance column.
    # Full kernels use the factor-derived covariance eigensystem, explicitly charged.
    for ti,t in enumerate(times):
        st=time.perf_counter();y=expm_multiply(-t*Dop,fac['Sigma_physical'][:,0],traceA=-t*traceD);dt=time.perf_counter()-st
        expected=mm(pkernel.vectors*(pkernel.values*np.exp(-t/pkernel.values)),pkernel.vectors[0])
        st=time.perf_counter();pC=pkernel.kernel(t,False)-pK0;kt=time.perf_counter()-st;physCs.append(pC)
        drift_results.append({'time':float(t),'single_covariance_column_factor_exponential_action_seconds':dt,'column_max_abs_error_vs_factor_covariance_eigen_action':float(np.max(abs(y-expected))),'physical_full_Wick_seconds':kt,'full_fixed_normalization_response_change':float(np.max(abs(pC-refs[ti])))})
    records['physical_drift_response']={'factor_consistent_covariance_kernel_setup_seconds':physical_setup,'zero_time_kernel_seconds':physical_zero,'by_time':drift_results,'qualification':'Changes both covariance and OU drift to the force-factor-consistent covariance; no claim this repairs the original nonlinear model.'}
    print('actual ABL small physical model',flush=True);sm,sS,ssd,sfit=physical_small_model(model);sf=force_factor(sm,sS);sh=Hermite(sm,sS,ssd);occ=occupations(len(sS));hh=np.array([sh.coeff(ix) for ix in occ]);omega=np.array([sum(sh.rates[list(ix)]) for ix in occ]);sk=OUKernel(sm,sS,ssd);checks=[]
    for t in np.array([0,.01,.1])*sk.tau:
        exact=mm(hh.T,np.exp(-t*omega)[:,None]*hh);reference=sk.kernel(t,False);checks.append({'time':float(t),'kernel_maxabs':float(np.max(abs(exact-reference)))})
    np.savez_compressed(out/'small-physical-model.npz',**sm,Sigma=sS,sd=ssd,h=hh,omega=omega,occupations=np.array([list(ix)+[-1]*(4-len(ix)) for ix in occ]),force_factor=sf['F'],force_Gamma=sf['Gamma'],covariance_vectors=sh.V,covariance_values=sh.values,tau=sk.tau)
    records['small_physical_model']={'r0':sm['r0'].tolist(),'edges':sm['edges'].tolist(),'canonical':sm['canonical'].tolist(),'coordinate_dimension':len(sS),'Hermite_dimension':len(occ),'fit':sfit,'force_factor':sf['stats'],'coefficient_kernel_checks':checks}
    records['reduction_error_maxabs_vs_full_Gaussian']=float(np.max(abs(Cs-refs)));records['matrix_free_error_maxabs']=float(np.max(abs(np.array(actionC)-Cs)));records['factor_covariance_response_change_maxabs']=float(np.max(abs(np.array(physCs)-refs)))
    records['protocol_sha256']=sha(ROOT/'protocol.json');records['wall_seconds']=time.perf_counter()-start;records['peak_RSS_bytes']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    np.savez_compressed(out/'matched-responses.npz',times=times,pairs=pairs,reduced_spectral=Cs,reduced_matrix_free=np.array(actionC),full_Gaussian=refs,factor_consistent_Gaussian=np.array(physCs),G0=G0,reference_G0=K0)
    records['status']='COMPLETE_DIAGNOSTIC_NO_FULL_PROTEIN_QUANTUM_ORACLE';dump(out/'classical-summary.json',records)
    env={'UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'platform':platform.platform(),'python':sys.version,'packages':{k:importlib.metadata.version(k) for k in ['numpy','scipy','qiskit','qiskit-aer','threadpoolctl']},'threadpools':threadpool_info(),'sources':{str(p.relative_to(ROOT)):sha(p) for p in list((ROOT/'source').glob('*.py'))+list((ROOT/'vendor').glob('*.py'))}}
    dump(out/'environment.json',env);print(json.dumps({k:records[k] for k in ['status','reduction_error_maxabs_vs_full_Gaussian','matrix_free_error_maxabs','factor_covariance_response_change_maxabs','wall_seconds','peak_RSS_bytes']},indent=2),flush=True)

if __name__=='__main__':
    with threadpool_limits(limits=4):main()
