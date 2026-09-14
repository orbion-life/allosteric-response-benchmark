"""Bounded numerical MPS probes of actual full33-bit Clifford+T phase bank."""
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[key]='1'
import hashlib,json,math,time
from pathlib import Path
import numpy as np
from qiskit import QuantumCircuit,qpy
from qiskit.quantum_info import Pauli
from qiskit_aer import AerSimulator
ROOT=Path(__file__).resolve().parent

def main():
    start=time.perf_counter();b=33
    with (ROOT/'results/phase-gradient-bank-b33.qpy').open('rb') as f:bank=qpy.load(f)[0]
    with (ROOT/'results/gradient-preparation-b33.qpy').open('rb') as f:prep=qpy.load(f)[0]
    cases=[('low_and_high_bits',0x123456789,False)]
    simulator=AerSimulator(method='matrix_product_state',max_parallel_threads=1,max_parallel_experiments=1,max_parallel_shots=1,max_memory_mb=12000,
        matrix_product_state_truncation_threshold=1e-20,matrix_product_state_max_bond_dimension=256,mps_log_data=True,seed_simulator=1729)
    records=[]
    for label,word,coherent in cases:
        qc=QuantumCircuit(bank.num_qubits);target=2*b+1;helper=2*b;scratch=2*b+2
        for k in range(b):
            if word>>k&1:qc.x(k)
        if coherent:qc.h(b-1)
        qc.compose(prep,list(range(b,2*b)),inplace=True);qc.compose(bank,inplace=True)
        for axis in ['X','Y','Z']:qc.save_expectation_value(Pauli(axis),[target],label='target_'+axis)
        if coherent:qc.save_expectation_value(Pauli('X'),[b-1],label='word_coherence')
        for bit in range(b):qc.save_expectation_value(Pauli('Z'),[bit],label='word_'+str(bit))
        qc.save_expectation_value(Pauli('Z'),[helper],label='helper');qc.save_expectation_value(Pauli('Z'),[scratch],label='scratch')
        qc.compose(prep.inverse(),list(range(b,2*b)),inplace=True)
        for bit in range(b):qc.save_expectation_value(Pauli('Z'),[b+bit],label='gradient_return_'+str(bit))
        result=simulator.run(qc,shots=1).result();assert result.success;values={k:float(v) for k,v in result.data(0).items()}
        angle=4*math.pi*word/2**b;expected={'target_X':math.sin(angle),'target_Y':0.,'target_Z':math.cos(angle),'helper':1.,'scratch':1.}
        if coherent:expected['word_coherence']=-1.
        for bit in range(b):expected['word_'+str(bit)]=0. if coherent and bit==b-1 else (1. if not(word>>bit&1) else -1.)
        for bit in range(b):expected['gradient_return_'+str(bit)]=1.
        error=max(abs(values[k]-v) for k,v in expected.items());assert error<1e-8
        records.append(dict(case=label,word=word,coherent_word=coherent,maximum_expectation_error=error,observed=values,expected=expected,metadata=result.results[0].metadata))
        (ROOT/'results/low-word-MPS-probe.json').write_text(json.dumps(dict(status='RUNNING' if len(records)<1 else 'PASS',scope='Numerical full69-qubit bank probes, fourdeclaredinputs; deterministicexpectations, nothardware/shotsorfull182-qubitestimator. MPSbond256/truncation1e-20 stated; no noise.',cases=records,wall_seconds=time.perf_counter()-start),indent=2)+'\n')
    print(json.dumps({'status':'PASS','seconds':time.perf_counter()-start,'maximum_error':max(x['maximum_expectation_error'] for x in records)},indent=2))

if __name__=='__main__':main()
