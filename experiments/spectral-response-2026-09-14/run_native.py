"""Build all spectral components and verify selected actual Grover circuits."""
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:
    os.environ[key]='4'
import argparse,hashlib,json,math,resource,time
from pathlib import Path
import numpy as np
from scipy.linalg import eigh,expm
from qiskit import QuantumCircuit,qpy
from qiskit.quantum_info import Pauli
from qiskit_aer import AerSimulator
from spectral import observable,attenuation,preparation,grover,counts
ROOT=Path(__file__).resolve().parent


def write(path,obj):path.write_text(json.dumps(obj,indent=2)+'\n')
def save(qc,path):
    with path.open('wb') as f:qpy.dump(qc,f)


def run(name):
    started=time.perf_counter();out=ROOT/'results'/name;out.mkdir(parents=True,exist_ok=False)
    z=np.load(ROOT/'inputs'/f'{name}.npz');H=z['Hr'] if name=='triangle' else z['H'];B=z['B'];times=z['times']
    n=(len(H)-1).bit_length();N=B.shape[1];norm=np.linalg.norm(B,axis=0);scale=float(max(norm)**2)
    tick=time.perf_counter();ev,V=eigh(H);teigh=time.perf_counter()-tick
    assert ev[0]>=0,ev[0]
    tick=time.perf_counter();VB=V.T@B;ttransform=time.perf_counter()-tick
    tick=time.perf_counter();C=np.array([(VB.T*np.expm1(-t*ev))@VB for t in times]);tspectral=time.perf_counter()-tick
    tick=time.perf_counter();K=np.array([(VB.T*np.exp(-t*ev))@VB for t in times]);tkernel=time.perf_counter()-tick
    tick=time.perf_counter();Cex=np.array([B.T@(expm(-t*H)-np.eye(len(H)))@B for t in times]);texpm=time.perf_counter()-tick
    if name=='triangle':projection=float(np.max(abs(C-z['reference_C'])))
    else:projection=json.loads((ROOT/'inputs/kras-prior-receipt.json').read_text())['reference_checks']['projection_to_exact_Gaussian']
    eps_response=.002-projection-1e-8-2e-6
    assert eps_response>0
    tick=time.perf_counter();preps=[observable(VB[:,i],n) for i in range(N)];filters=[attenuation(ev,t,n) for t in times];tcomponents=time.perf_counter()-tick
    (out/'components').mkdir()
    for i,p in enumerate(preps):save(p,out/'components'/f'prepare-{i}.qpy')
    for i,p in enumerate(filters):save(p,out/'components'/f'filter-{i}.qpy')
    np.savez_compressed(out/'model.npz',H=H,B=B,eigenvalues=ev,eigenvectors=V,transformed_B=VB,norm=norm,times=times,C=C,K=K,C_expm=Cex)
    pc=[counts(p) for p in preps];fc=[counts(p) for p in filters]
    cases=[(ti,0,1) for ti in range(3)]+[(1,0,0),(1,N-2,N-1)]
    backend=AerSimulator(method='statevector',max_parallel_threads=4,max_parallel_experiments=1,max_memory_mb=12000)
    checks=[]
    for index,(ti,i,j) in enumerate(cases):
        A=preparation(preps[i],preps[j],filters[ti]);Q,S0=grover(A);p=(1+K[ti,i,j]/(norm[i]*norm[j]))/2
        save(A,out/f'A-case{index}.qpy');save(Q,out/f'Q-case{index}.qpy')
        if index==0:save(S0,out/'zero-reflection.qpy')
        powers=[0,1,2,4] if name=='triangle' else ([0,1] if index==1 else [0])
        for k in powers:
            qc=QuantumCircuit(Q.num_qubits);qc.compose(A,range(A.num_qubits),inplace=True)
            for _ in range(k):qc.compose(Q,inplace=True)
            qc.save_expectation_value(Pauli('Z'),[0],label='z')
            tick=time.perf_counter();result=backend.run(qc,shots=1).result();assert result.success
            observed=(1+float(result.data(0)['z']))/2;expected=math.sin((2*k+1)*math.asin(math.sqrt(min(1,max(0,p)))))**2
            error=abs(observed-expected);assert error<1e-8,(name,index,k,observed,expected)
            checks.append(dict(case=index,time_index=ti,i=i,j=j,grover_power=k,active_qubits=A.num_qubits,executed_qubits=qc.num_qubits,
                               actual_A=counts(A),actual_Q=counts(Q),expected_probability=expected,observed_probability=observed,error=error,seconds=time.perf_counter()-tick))
            print(name,index,k,error,flush=True)
    result=dict(status='PASS_ACTUAL_NATIVE_SPECTRAL_AND_GROVER_CHECKS',rank=len(H),residues=N,padded_rank=2**n,
                input_sha256=hashlib.sha256((ROOT/'inputs'/f'{name}.npz').read_bytes()).hexdigest(),
                source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                active_qubits=n+2,Grover_total_qubits=2*(n+2)-2,clean_reflection_helpers=n,
                exact_classical=dict(eigh_seconds=teigh,transform_all_B_seconds=ttransform,three_full_response_matrices_seconds=tspectral,
                                     three_full_kernel_matrices_seconds=tkernel,three_expm_comparator_seconds=texpm,spectral_vs_expm_max=float(np.max(abs(C-Cex)))),
                inherited_projection_error=projection,response_sampling_budget=eps_response,amplitude_epsilon=eps_response/(2*scale),scale=scale,
                preparations=pc,filters=fc,zero_reflection=counts(S0),checks=checks,component_construction_seconds=tcomponents,
                maximum_probability_error=max(c['error'] for c in checks),wall_seconds=time.perf_counter()-started,
                peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                scope='All observable and attenuation circuits actually constructed; representative exact native A and Q executed by noiseless local statevector including clean reflection helpers. Deterministic expectations, no QPU or measured sampling run. New exact spectral representation uses classical eigendecomposition and transformed B.')
    write(out/'native-receipt.json',result);print(name,'COMPLETE',result['wall_seconds'],flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('model',choices=['triangle','kras']);run(parser.parse_args().model)
