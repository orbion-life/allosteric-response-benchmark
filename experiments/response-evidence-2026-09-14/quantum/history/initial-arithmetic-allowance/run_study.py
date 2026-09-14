"""Reproducible local compilation, ideal validation, noise sensitivity and shot budgets."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='4'
from pathlib import Path
import json,sys,time,hashlib,math,platform
import numpy as np
from scipy.linalg import eigh,expm
from scipy.stats import beta
from qiskit import QuantumCircuit,transpile,qpy
from qiskit.transpiler import CouplingMap
from qiskit.quantum_info import Statevector,Pauli
from qiskit_aer import AerSimulator
from qiskit_aer.noise import NoiseModel,depolarizing_error
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'vendor'))
from spectral import observable,attenuation,preparation,counts

def write(p,v):p.write_text(json.dumps(v,indent=2)+'\n')
def cp(k,n,a):return (0.0 if k==0 else float(beta.ppf(a/2,k,n-k+1)),1.0 if k==n else float(beta.ppf(1-a/2,k+1,n-k)))
def native_qasm(qc,readout):
    physical=[7,8,9,10,13,14]
    lines=['OPENQASM 3.0;','bit[1] b;','#pragma braket verbatim','box {']
    for ins in qc.data:
        qs=[physical[qc.find_bit(q).index] for q in ins.qubits]
        if ins.operation.name=='r':lines.append(f'  prx({float(ins.operation.params[0]):.17g}, {float(ins.operation.params[1]):.17g}) ${qs[0]};')
        elif ins.operation.name=='cz':lines.append(f'  cz ${qs[0]}, ${qs[1]};')
        else:raise ValueError(ins.operation.name)
    lines+=['}',f'b[0] = measure ${physical[readout]};']
    return '\n'.join(lines)+'\n'

def run():
    start=time.perf_counter();out=ROOT/'results';out.mkdir(exist_ok=False);(out/'circuits').mkdir()
    z=np.load(ROOT/'triangle.npz');H=z['Hr'];B=z['B'];times=z['times'];N=B.shape[1]
    tick=time.perf_counter();ev,V=eigh(H);VB=V.T@B;norm=np.linalg.norm(B,axis=0);G0=B.T@B
    Cex=np.array([B.T@(expm(-t*H)-np.eye(len(H)))@B for t in times]);classical_seconds=time.perf_counter()-tick
    scale=float(max(norm)**2);projection=float(np.max(abs(Cex-z['reference_C'])))
    prep=[observable(VB[:,i],4) for i in range(N)]
    cmap=CouplingMap([[a,b] for i,j in [(0,1),(1,2),(2,3),(1,4),(2,5)] for a,b in [(i,j),(j,i)]])
    cases=[];ideal=[];compiled=[]
    # Full nontrivial fixture uses each diagonal and off-diagonal pair; identity controls only distinct pairs.
    queries=[(f'data-t{ti}-{i}-{j}',float(t),i,j,ti,False) for ti,t in enumerate(times) for i in range(N) for j in range(i,N)]
    queries += [(f'control-t0-{i}-{j}',0.,i,j,-1,True) for i in range(N) for j in range(i+1,N)]
    for label,t,i,j,ti,iscontrol in queries:
        A=preparation(prep[i],prep[j],attenuation(ev,t,4))
        tick=time.perf_counter();c=transpile(A,basis_gates=['r','cz'],coupling_map=cmap,initial_layout=list(range(6)),optimization_level=3,seed_transpiler=1729,num_processes=1);elapsed=time.perf_counter()-tick
        readout=c.layout.final_index_layout()[0]
        # Every two-qubit operation must fit the declared local graph.
        for inst in c.data:
            if len(inst.qubits)==2:assert tuple(c.find_bit(q).index for q in inst.qubits) in cmap.get_edges()
        p=float(Statevector.from_instruction(c).probabilities([readout])[0]);expected=float((1+VB[:,i]@(np.exp(-t*ev)*VB[:,j])/(norm[i]*norm[j]))/2)
        response=(2*p-1)*norm[i]*norm[j]-G0[i,j]
        reference=0. if iscontrol else Cex[ti,i,j]
        item=dict(label=label,time=t,i=i,j=j,time_index=ti,identity_control=iscontrol,qubits=c.num_qubits,mode_rank=len(H),probability=p,expected_probability=expected,probability_error=abs(p-expected),response=float(response),reference_response=float(reference),response_error=float(abs(response-reference)),response_probability_multiplier=float(2*norm[i]*norm[j]),active_readout_local=readout,active_readout_device=[7,8,9,10,13,14][readout],original=counts(A),compiled=counts(c),compile_seconds=elapsed,final_index_layout=c.layout.final_index_layout(),qasm_file=f'circuits/{label}.qasm')
        with (out/'circuits'/f'{label}.qpy').open('wb') as f:qpy.dump(c,f)
        (out/'circuits'/f'{label}.qasm').write_text(native_qasm(c,readout))
        cases.append(item);compiled.append(c)
        print(label,'CZ',c.count_ops().get('cz'),'depth',c.depth(),'err',item['response_error'],flush=True)
    write(out/'circuit_results.json',cases)
    scenarios=[('ideal',0.,0.),('optimistic',1e-5,1e-4),('low_noise',1e-4,1e-3),('moderate_noise',1e-3,.01),('high_noise',.01,.01),('readout_only',0.,.01)]
    noise_results=[]
    for name,p2,readerror in scenarios:
        nm=NoiseModel()
        if p2:
            nm.add_all_qubit_quantum_error(depolarizing_error(p2,2),['cz']);nm.add_all_qubit_quantum_error(depolarizing_error(p2/10,1),['r'])
        backend=AerSimulator(method='density_matrix',noise_model=nm,max_parallel_threads=4,max_parallel_experiments=1)
        row=[]
        for c,item in zip(compiled,cases):
            test=c.copy();test.save_expectation_value(Pauli('Z'),[item['active_readout_local']],label='z')
            # Aer r gate supported; no further transpilation that would alter native noise placement.
            result=backend.run(test,shots=1).result();assert result.success
            raw=(1+float(result.data(0)['z']))/2;p=readerror+(1-2*readerror)*raw
            response=(2*p-1)*norm[item['i']]*norm[item['j']]-G0[item['i'],item['j']]
            row.append(dict(label=item['label'],probability=p,response=float(response),bias_to_reduced_response=float(abs(response-item['reference_response'])),identity_control=item['identity_control']))
        datum=dict(name=name,two_qubit_depolarizing_parameter=p2,one_qubit_depolarizing_parameter=p2/10,readout_flip_probability=readerror,scope='Assumed independent stochastic noise, not device calibration',maximum_data_response_bias=max(x['bias_to_reduced_response'] for x in row if not x['identity_control']),maximum_identity_control_bias=max(x['bias_to_reduced_response'] for x in row if x['identity_control']),queries=row)
        noise_results.append(datum);write(out/'noise_sensitivity.json',noise_results);print('noise',name,datum['maximum_data_response_bias'],flush=True)
    # Every saved interval has exact binomial coverage for its own stationary noisy probability.
    rng=np.random.default_rng(20260914);sampling=[];J=18
    for shots in [16384,65536]:
        for scenario in noise_results:
            data=[x for x in scenario['queries'] if not x['identity_control']];items=[]
            for item,d in zip(cases[:J],data):
                k=int(rng.binomial(shots,d['probability']));lo,hi=cp(k,shots,.05/J);mult=item['response_probability_multiplier'];offset=-norm[item['i']]*norm[item['j']]-G0[item['i'],item['j']]
                rlo=mult*lo+offset;rhi=mult*hi+offset;mid=(rlo+rhi)/2
                items.append(dict(label=item['label'],shots=shots,zeros=k,probability_interval=[lo,hi],response_interval=[float(rlo),float(rhi)],response_midpoint=float(mid),response_halfwidth=float((rhi-rlo)/2),noisy_probability_covered=bool(lo<=d['probability']<=hi),ideal_response_covered=bool(rlo<=item['reference_response']<=rhi),response_midpoint_error=float(abs(mid-item['reference_response']))))
            sampling.append(dict(scenario=scenario['name'],shots_per_query=shots,family_failure_probability=.05,queries=items,maximum_response_halfwidth=max(x['response_halfwidth'] for x in items),noisy_probabilities_covered=sum(x['noisy_probability_covered'] for x in items),ideal_responses_covered=sum(x['ideal_response_covered'] for x in items),maximum_response_midpoint_error=max(x['response_midpoint_error'] for x in items)))
    write(out/'sampling_simulations.json',sampling)
    budgets=[]
    for shots in [16384,65536]:
        # Rigorous simultaneous sufficient halfwidth for sample mean based on Hoeffding, CP used on actual counts.
        wh=2*scale*math.sqrt(math.log(2*J/.05)/(2*shots))
        ctrlshots=4*8192+3*16384;datashots=18*shots
        tasks=18*math.ceil(shots/20000)+7
        budgets.append(dict(shots_per_data_query=shots,total_data_shots=datashots,control_shots=ctrlshots,total_shots=datashots+ctrlshots,tasks=tasks,price_USD_at_verified_tariff=(datashots+ctrlshots)*.00145+tasks*.3,rigorous_statistical_response_halfwidth=wh,projection_error=projection,arithmetic_compile_allowance=1e-8,remaining_bias_allowance_for_0_05=.05-projection-1e-8-wh))
    target=.002;physical=.0005;numeric=1e-8;remaining=target-projection-physical-numeric;prob=remaining/(2*scale);n=math.ceil(math.log(2*J/.05)/(2*prob**2));ctrlshots=4*8192+3*16384;task=18*math.ceil(n/20000)+7
    finalbudget=dict(total_response_target=target,projection_error=projection,numerical_allowance=numeric,assumed_qualified_physical_bias=physical,statistical_budget=remaining,sufficient_shots_per_query=n,query_count=J,total_data_shots=J*n,total_shots_with_diagnostic_controls=J*n+ctrlshots,tasks=task,price_USD_at_verified_tariff=(J*n+ctrlshots)*.00145+task*.3,qualification='Sampling bound assumes physical bias separately established at <=0.0005; these diagnostic controls do not establish that bias for all circuit errors')
    write(out/'shot_budgets.json',dict(diagnostic=budgets,final_precision=finalbudget))
    summary=dict(status='LOCAL_NATIVE_COMPILATION_IDEAL_AND_ASSUMED_NOISE_STUDY_COMPLETE',platform='IQM Garnet via Amazon Braket, documented topology subset; not a live calibration',qubits=6,rank=9,sites=3,queries=18,controls=7,physical_positions=[7,8,9,10,13,14],native_gate_set=['PRx','CZ'],data_CZ_range=[min(x['compiled']['ops']['cz'] for x in cases[:18]),max(x['compiled']['ops']['cz'] for x in cases[:18])],data_PRx_range=[min(x['compiled']['ops']['r'] for x in cases[:18]),max(x['compiled']['ops']['r'] for x in cases[:18])],data_depth_range=[min(x['compiled']['depth'] for x in cases[:18]),max(x['compiled']['depth'] for x in cases[:18])],maximum_ideal_probability_error=max(x['probability_error'] for x in cases),maximum_ideal_response_error=max(x['response_error'] for x in cases),inherited_projection_error=projection,maximum_response_scale=scale,classical_preprocessing_and_reference_seconds=classical_seconds,wall_seconds=time.perf_counter()-start,versions=dict(python=platform.python_version(),numpy=np.__version__),files_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [ROOT/'PROTOCOL.md',ROOT/'run_study.py',ROOT/'triangle.npz',ROOT/'vendor/spectral.py',ROOT/'vendor/vendor/multiplexed.py']},noise_scenarios=[{k:v for k,v in s.items() if k!='queries'} for s in noise_results],shot_budgets=budgets,final_precision=finalbudget,no_qpu_execution=True)
    write(out/'summary.json',summary);print(json.dumps(summary,indent=2))
if __name__=='__main__':run()
