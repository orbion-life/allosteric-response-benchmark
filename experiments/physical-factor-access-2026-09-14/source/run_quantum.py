"""Bounded direct physical-factor circuit; no dense unitary is materialized."""
from core import *
from qiskit import QuantumCircuit, transpile, qpy
from qiskit_aer import AerSimulator
from qiskit.quantum_info import Pauli
from scipy.special import ive
from threadpoolctl import threadpool_limits
import signal, resource
import multiplexed as mux

def occupation_counts(ix,d):
    return tuple(ix.count(j) for j in range(d))

def polynomial_factor(factor,V,occ):
    """Physical gradient lowering factor in normalized Hermite coordinates."""
    d=V.shape[0];base=[()]+occupations(d,3);base_lookup={occupation_counts(x,d):i for i,x in enumerate(base)}
    F=mm(factor,V);D=np.zeros((len(base)*len(F),len(occ)))
    for col,ix in enumerate(occ):
        counts=list(occupation_counts(ix,d))
        for j in range(d):
            if counts[j]:
                lower=counts.copy();lower[j]-=1;row=base_lookup[tuple(lower)]
                D[row*len(F):(row+1)*len(F),col]+=math.sqrt(counts[j])*F[:,j]
    return D,base

def native(qc):return transpile(qc,basis_gates=['u','cx'],optimization_level=1,seed_transpiler=20260914)
def counts(qc):return {'qubits':qc.num_qubits,'operations':qc.size(),'depth':qc.depth(),'gate_counts':{str(k):int(v) for k,v in qc.count_ops().items()},'statevector_complex128_bytes':16*2**qc.num_qubits}

def simulator(qc):
    """Bounded Aer CPU ideal statevector run with an interrupting time cap."""
    def expired(*_):raise TimeoutError('90-second simulation cap')
    old=signal.signal(signal.SIGALRM,expired);signal.setitimer(signal.ITIMER_REAL,90)
    st=time.perf_counter()
    try:
        result=AerSimulator(method='statevector',device='CPU',precision='double',max_parallel_threads=4,max_parallel_experiments=1).run(qc).result()
        if not result.success:raise ArithmeticError(str(result.status))
        return result.data(0),time.perf_counter()-st
    finally:signal.setitimer(signal.ITIMER_REAL,0);signal.signal(signal.SIGALRM,old)

def main():
    start=time.perf_counter();out=ROOT/'results';small=dict(np.load(out/'small-physical-model.npz'));d=len(small['Sigma']);occ=occupations(d);h=small['h'];D,lower=polynomial_factor(small['force_factor'],small['covariance_vectors'],occ)
    H=mm(D.T,D);gamma=float(np.sum(D*D));target=np.diag(small['omega']);delta=float(np.max(abs(H-target)));tau=float(small['tau']);norm=np.linalg.norm(h,axis=0)
    records={'scope':'First3 actual ABL residues and first2 harmonic coordinates; complete degree1–4 polynomial observables. Not a full ABL response circuit.','canonical':small['canonical'].tolist(),'coordinate_dimension':d,'Hermite_dimension':len(occ),'factor_shape':list(D.shape),'factor_nonzero_entries':int(np.count_nonzero(D)),'factor_bytes':D.nbytes,'gamma_constructed_Frobenius_squared':gamma,'factor_Gram_vs_covariance_generator_maxabs':delta,'factor_vs_covariance_generator_tau_scaled_maxabs':tau*delta,'gamma_over_generator_spectral_norm':gamma/float(eigh(H,eigvals_only=True)[-1]),'observable_norms':norm.tolist(),'classical_covariance_eigendecomposition_required':True,'physical_force_roots_dimension':3,'full_reduced_Hr_eigendecomposition_used_for_circuit':False,'precision_target_normalized_response':1e-6}
    ul,pin,pout,m=mux.loader_circuit(D,gamma,False);uln=native(ul);walk=QuantumCircuit(ul.num_qubits,name='physical_factor_walk');walk.compose(ul,inplace=True);mux.reflection(walk,m['col']+[m['flag']]);walk.compose(ul.inverse(),inplace=True);mux.reflection(walk,m['row']+[m['flag']],input_reflection=True);wn=native(walk)
    records['loader']=counts(uln);records['walk']=counts(wn);records['metadata']={k:v for k,v in m.items() if k not in ['readout']}
    if uln.num_qubits>20 or uln.size()>200000:raise ArithmeticError('Small physical loader exceeds frozen cap')
    with open(out/'physical-factor-loader.qpy','wb') as f:qpy.dump(uln,f)
    with open(out/'physical-factor-walk.qpy','wb') as f:qpy.dump(wn,f)
    projected=[];walk_projected=[];simtimes=[];dc=m['padded_columns'];dr=m['padded_rows'];good=np.arange(dr)*dc
    # Exactly14 factor columns, never a dense4096×4096 unitary.
    for col in range(len(occ)):
        initial=QuantumCircuit(uln.num_qubits)
        for j in range(m['ncol']):
            if col&(1<<j):initial.x(j)
        c=initial.compose(uln);c.save_statevector();data,dt=simulator(c);projected.append(np.asarray(data['statevector'])[good]);simtimes.append(dt)
        c=initial.compose(wn);c.save_statevector();data,dt=simulator(c);walk_projected.append(np.asarray(data['statevector'])[:dc]);simtimes.append(dt)
    actual=np.array(projected).T;actual_walk=np.array(walk_projected).T;lp=np.zeros((dr,len(occ)));lp[:len(D)]=D;wp=np.eye(dc)[:,:len(occ)];wp[:len(occ)]-=2*H/gamma
    records['compiled_component_checks']={'all_factor_columns_maxabs':float(np.max(abs(actual-lp/np.sqrt(gamma)))),'all_projected_walk_columns_maxabs':float(np.max(abs(actual_walk-wp))),'simulation_seconds':simtimes,'validation_method':'28 separate statevector-column runs, no full unitary matrix'}
    np.savez_compressed(out/'small-physical-factor.npz',D=D,H=H,gamma=gamma,projected_factor_columns=actual,projected_walk_columns=actual_walk)
    print('small physical factor compiled/ideal checked',records['compiled_component_checks'],flush=True)
    # Native controlled walk is compiled once, with its control as bit0.
    cwn=native(mux.controlled_walk(D,gamma));records['controlled_walk']=counts(cwn)
    with open(out/'physical-factor-controlled-walk.qpy','wb') as f:qpy.dump(cwn,f)
    results=[];i,j=0,1;L=float(norm[i]*norm[j]);static=float(mm(h[:,i],h[:,j]));evals,evecs=eigh(H)
    for ratio in [.01,.1]:
        t=ratio*tau;a=gamma*t/2;coeff=np.array([ive(0,a)]+[2*ive(k,a) for k in range(1,33)]);tails=np.maximum(0,1-np.cumsum(coeff));eligible=np.flatnonzero(L*tails<=2.5e-7);degree=int(eligible[0]) if len(eligible) else None
        row={'time_over_tau':ratio,'time':t,'pair_canonical':small['canonical'][[i,j]].tolist(),'norm_restoration_L':L,'static_subtraction':static,'a_gamma_t_over2':a,'degree_cap32_tail_response_bound':float(L*tails[-1]),'chosen_degree':degree,'response_reference_factor':float(mm(h[:,i],mm(evecs*np.exp(-t*evals),mm(evecs.T,h[:,j])))-static),'response_reference_fitted_covariance':float(np.sum(h[:,i]*np.exp(-t*small['omega'])*h[:,j])-static)}
        if degree is None:
            row['status']='FAIL_CHEBYSHEV_DEGREE_CAP';results.append(row);continue
        coeff=coeff[:degree+1];alpha=float(sum(coeff));nb=max(1,math.ceil(math.log2(degree+1)));cvalues=np.zeros(2**nb);cvalues[:len(coeff)]=np.sqrt(coeff/alpha)
        total=1+uln.num_qubits+nb+1;readout=0;system=list(range(1,1+uln.num_qubits));lcu=list(range(1+uln.num_qubits,1+uln.num_qubits+nb));andbit=total-1
        qc=QuantumCircuit(total,name='physical_factor_LCU_response');qc.h(readout);obs=np.zeros((2,dc));obs[0,:len(h)]=h[:,i]/norm[i];obs[1,:len(h)]=h[:,j]/norm[j]
        qc.compose(mux.real_loader(obs,[1+x for x in m['col']],[readout],total),inplace=True)
        prep=mux.real_loader(cvalues[None,:],lcu,[],total);qc.compose(prep,inplace=True)
        for bit,power in zip(lcu,[2**k for k in range(nb)]):
            qc.ccx(readout,bit,andbit)
            for _ in range(power):qc.compose(cwn,[andbit]+system,inplace=True)
            qc.ccx(readout,bit,andbit)
        qc.compose(prep.inverse(),inplace=True)
        nqc=native(qc);row['LCU_alpha']=alpha;row['truncation_response_bound']=float(L*(1-alpha));row['controlled_walk_calls']=2**nb-1;row['native_circuit']=counts(nqc);row['unnormalized_response_from_probability']='C=L*alpha*(2p_plus-1)-static; p_plus is X=+1 Hadamard probability'
        with open(out/f'physical-response-{ratio}.qpy','wb') as f:qpy.dump(nqc,f)
        if nqc.num_qubits>20 or nqc.size()>200000:
            row['status']='FAIL_FROZEN_COMPILED_RESOURCE_CAP';results.append(row);print('LCU resource gate',ratio,row['native_circuit'],flush=True);continue
        try:
            nqc.save_expectation_value(Pauli('X'),[readout],label='X');data,dt=simulator(nqc);X=float(np.real(data['X']));measured=L*alpha*X-static
            row['ideal_simulation_seconds']=dt;row['ideal_X_expectation']=X;row['ideal_plus_probability']=(1+X)/2;row['norm_restored_response']=measured;row['error_vs_factor_response']=abs(measured-row['response_reference_factor']);row['error_vs_fitted_Gaussian_response']=abs(measured-row['response_reference_fitted_covariance']);row['status']='PASS' if row['error_vs_fitted_Gaussian_response']<=1e-6 else 'FAIL_RESPONSE_ACCURACY'
        except TimeoutError as exc:row['status']='FAIL_FROZEN_SIMULATION_TIME_CAP';row['reason']=str(exc)
        results.append(row);print('LCU response',ratio,row,flush=True);dump(out/'quantum-progress.json',{**records,'LCU_responses':results})
    records['LCU_responses']=results;records['status']='COMPLETE_BOUNDED_DIRECT_FACTOR_DIAGNOSTIC';records['full_ABL_oracle_status']='NOT_IMPLEMENTED_RESOURCE_GATE_FAILED';records['hardware_execution']=False;records['shot_cost_scope']='These are ideal compiled simulations. A sampled estimate needs uncertainty analysis and coherent hardware feasibility; no advantage from simulator runtimes.';records['wall_seconds']=time.perf_counter()-start;records['peak_RSS_bytes']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    dump(out/'quantum-summary.json',records);print(json.dumps(records,indent=2),flush=True)

if __name__=='__main__':
    with threadpool_limits(limits=4):main()
