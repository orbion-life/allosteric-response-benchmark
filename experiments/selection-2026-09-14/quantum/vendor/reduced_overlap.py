"""Small actual-Galerkin-operator Chebyshev overlap circuit; no expm gate.

Dense block synthesis is deliberately explicit and charged. This is a
compatibility canary, not a sparse oracle or a quantum-advantage claim.
"""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[k]='1'
import argparse, hashlib, json, math, time, sys, resource
from pathlib import Path
import numpy as np
import scipy
from scipy.linalg import eigh, expm
from scipy.special import ive
import qiskit
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, transpile, qasm3
from qiskit.circuit.library import StatePreparation, UnitaryGate
from qiskit.quantum_info import Statevector, Operator
from qiskit.transpiler import CouplingMap


def coefficients(z, allocation=1e-8):
    a=[float(ive(0,z))]
    for k in range(1,4097):
        a.append(float(2*ive(k,z)))
        if 1-sum(a)<allocation:
            return np.asarray(a), max(0.,1-sum(a))
    raise RuntimeError('Polynomial degree ceiling exceeded')


def circuit(H,B,t,i,j,tail_budget=1e-8):
    r=len(H); n=max(1,int(math.ceil(math.log2(r)))); d=2**n
    ev,V=eigh(H)
    if min(ev)<-1e-9: raise ValueError('H is not positive semidefinite')
    alpha=max(1e-12,float(max(ev))*1.000000000001)
    A=np.eye(d); A[:r,:r]-=H/alpha
    av,Q=eigh(A)
    if min(av)<-1e-10 or max(av)>1+1e-10: raise ValueError('Invalid contraction')
    S=(Q*np.sqrt(np.maximum(0,1-av**2)))@Q.T
    U=np.block([[A,S],[S,-A]])
    R=np.diag([1.]*d+[-1.]*d)
    unitary_error=float(np.max(abs(U.T@U-np.eye(2*d))))
    norms=np.linalg.norm(B,axis=0)
    coeff,tail=coefficients(alpha*t,tail_budget/max(1,norms[i]*norms[j]))
    K=len(coeff)-1; nc=max(1,int(math.ceil(math.log2(K+1))))
    if K>512: raise RuntimeError('Canary degree exceeds 512')
    b=QuantumRegister(1,'readout'); c=QuantumRegister(nc,'coefficient'); x=QuantumRegister(n,'observable'); a=QuantumRegister(1,'block')
    qc=QuantumCircuit(b,c,x,a)
    left=np.zeros(d);right=np.zeros(d);left[:r]=B[:,i]/norms[i];right[:r]=B[:,j]/norms[j]
    amp=np.zeros(2**nc);amp[:len(coeff)]=np.sqrt(coeff/coeff.sum())
    pj=StatePreparation(right.tolist());pi=StatePreparation(left.tolist());pc=StatePreparation(amp.tolist())
    # Two controls apply R*U iff both bits are1. U is a dense ACTUAL H block,
    # never an inserted time-evolution or whole-estimator unitary.
    w=QuantumCircuit(n+1,name='block_walk')
    w.append(UnitaryGate(U,label='dense_U_A'),list(range(n+1)))
    w.z(n)
    ccw=w.to_gate().control(2)
    qc.h(b[0]);qc.append(pj.control(),list(b)+list(x));qc.append(pc.control(),list(b)+list(c))
    for bit in range(nc):
        for repeat in range(2**bit): qc.append(ccw,[b[0],c[bit]]+list(x)+list(a))
    qc.append(pc.inverse().control(),list(b)+list(c));qc.append(pi.inverse().control(),list(b)+list(x));qc.h(b[0])
    scale=float(norms[i]*norms[j]*coeff.sum());static=float(B[:,i]@B[:,j])
    # Independent three-term polynomial, rather than the quantum circuit object.
    T0=np.eye(d); poly=coeff[0]*T0
    if K:
        T1=A.copy();poly+=coeff[1]*T1
        for k in range(2,K+1):
            T2=2*A@T1-T0;poly+=coeff[k]*T2;T0,T1=T1,T2
    predicted=float(norms[i]*norms[j]*(left@poly@right)-static)
    metadata=dict(rank=r,padded_rank=d,alpha=alpha,time=float(t),z=alpha*t,degree=K,
        coefficient_qubits=nc,block_qubits=1,observable_qubits=n,readout_qubits=1,
        logical_qubits=qc.num_qubits,walk_calls=2**nc-1,unitary_error=unitary_error,
        polynomial_tail=tail,response_tail_bound=float(norms[i]*norms[j]*tail),
        scale=scale,static=static,polynomial_response=predicted,
        preprocessing='Dense eigensolves, dense square root and dense block synthesis are explicitly included; no efficient oracle assumed.')
    return qc,metadata


def expected(qc):
    readouts=[qc.find_bit(inst.qubits[0]).index for inst in qc.data if inst.operation.name=='measure']
    if len(readouts)>1: raise ValueError('Expected one readout')
    bit=readouts[0] if readouts else 0
    sv=Statevector.from_instruction(qc.remove_final_measurements(inplace=False))
    p=sv.probabilities([bit]);return float(p[0]-p[1])


def counts(qc):
    ops={str(k):int(v) for k,v in qc.count_ops().items()}
    return dict(qubits=qc.num_qubits,depth=qc.depth(),cx=ops.get('cx',0),u=ops.get('u',0),ops=ops,
                two_qubit_depth=qc.depth(lambda inst:len(inst.qubits)==2))


def run(inp,out):
    start=time.perf_counter();out.mkdir(parents=True,exist_ok=False)
    z=np.load(inp); H=z['Hr'];B=z['B'];times=z['times']
    if H.shape[0]>16: raise ValueError('Actual-operator canary rank ceiling16')
    records=[];result=[];direct=[];buildtime=0;simtime=0;classicaltime=0
    for ti,t in enumerate(times):
        at=time.perf_counter();exact=B.T@(expm(-float(t)*H)-np.eye(len(H)))@B;classicaltime+=time.perf_counter()-at
        direct.append(exact);C=np.zeros_like(exact)
        for i in range(B.shape[1]):
            for j in range(i,B.shape[1]):
                at=time.perf_counter();qc,m=circuit(H,B,float(t),i,j);buildtime+=time.perf_counter()-at
                at=time.perf_counter();chi=expected(qc);simtime+=time.perf_counter()-at
                C[i,j]=C[j,i]=m['scale']*chi-m['static']
                m.update(time_index=ti,left=i,right=j,statevector_expectation=chi,response=float(C[i,j]),
                         direct_response=float(exact[i,j]),absolute_error=float(abs(C[i,j]-exact[i,j])))
                records.append(m)
        result.append(C)
        (out/'progress.json').write_text(json.dumps(dict(completed_times=ti+1,records=records),indent=2)+'\n')
    # Compile complete largest-time off-diagonal circuit, including all loading,
    # coefficient preparation/unpreparation, controlled walk powers and readout.
    qc,m=circuit(H,B,float(times[-1]),0,min(1,B.shape[1]-1));compile_records=[]
    qc.add_register(ClassicalRegister(1,'result'));qc.measure(0,0)
    for topology in ('all_to_all','line'):
        at=time.perf_counter();kw=dict(basis_gates=['u','cx'],optimization_level=1,seed_transpiler=1729,num_processes=1)
        if topology=='line': kw.update(coupling_map=CouplingMap.from_line(qc.num_qubits,bidirectional=True),initial_layout=list(range(qc.num_qubits)),routing_method='sabre')
        compiled=transpile(qc,**kw);elapsed=time.perf_counter()-at
        if set(compiled.count_ops())-{'u','cx','measure'}: raise AssertionError('Unsynthesized gate')
        cnt=counts(compiled)
        if sum(cnt['ops'].values())>200000: raise RuntimeError('Compiled operation cap exceeded')
        at=time.perf_counter();value=expected(compiled);sec=time.perf_counter()-at
        original=expected(qc)
        compile_records.append(dict(topology=topology,compile_seconds=elapsed,simulation_seconds=sec,
             compiled_expectation=value,uncompiled_expectation=original,parity_error=abs(value-original),readout_physical_qubit=[compiled.find_bit(inst.qubits[0]).index for inst in compiled.data if inst.operation.name=='measure'][0],**cnt))
        qasm3.dump(compiled,(out/f'complete-overlap-{topology}.qasm').open('w'))
    response=np.array(result);reference=np.array(direct)
    np.savez_compressed(out/'response-matrices.npz',Hr=H,B=B,times=times,C=response,C_classical=reference)
    maxerr=float(np.max(abs(response-reference)));parity=max(x['parity_error'] for x in compile_records)
    maxscale=max(x['scale'] for x in records);mcount=len(records)
    shot_eps=.002;delta=.05
    shots=int(math.ceil(2*maxscale**2*math.log(2*mcount/delta)/shot_eps**2))
    receipt=dict(engineering_status='COMPLETE',scientific_status='PASS' if maxerr<1e-7 and parity<1e-9 else 'FAIL',
        scope='Actual finite-triangle Galerkin matrix; exact noiseless statevector compatibility only. No protein quantum run or quantum speedup.',
        input_path=str(inp),input_sha256=hashlib.sha256(inp.read_bytes()).hexdigest(),
        versions=dict(python=sys.version,numpy=np.__version__,scipy=scipy.__version__,qiskit=qiskit.__version__),
        max_response_error=maxerr,max_compiled_expectation_parity_error=parity,
        all_pair_circuit_build_seconds=buildtime,all_pair_statevector_seconds=simtime,
        classical_all_matrix_expm_seconds=classicaltime,complete_circuit_compilation=compile_records,
        overlap_records=records,measurement_plan=dict(precision_per_response=shot_eps,family_failure_probability=delta,
        distinct_unordered_overlaps=mcount,shots_per_overlap_sufficient=shots,total_shots_sufficient=shots*mcount,
        scope='Hoeffding simultaneous sufficient bound for noiseless Bernoulli sampling; excludes noise bias; not measured/minimum shots or ranking certificate.'),
        worker_wall_seconds=time.perf_counter()-start,worker_peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        whole_pipeline_exclusions='Protein equilibrium/moment integration, hardware pulses and actual shots are not executed; their costs are additional.',
        noise_status='NOT EVALUATED for this new circuit; prior grid-fixture noise results remain separate.')
    (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n'); print(json.dumps({k:receipt[k] for k in ['engineering_status','scientific_status','max_response_error','worker_wall_seconds']},indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);p.add_argument('--output',type=Path,required=True);args=p.parse_args();run(args.input,args.output)
