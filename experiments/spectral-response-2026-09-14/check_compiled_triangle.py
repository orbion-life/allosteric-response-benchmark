"""Execute the emitted triangle Clifford+T A and one full Grover iteration."""
import json,math,time
from pathlib import Path
from qiskit import QuantumCircuit,qpy
from qiskit.quantum_info import Pauli
from qiskit_aer import AerSimulator
ROOT=Path(__file__).resolve().parent
def load(path):
    with path.open('rb') as f:return qpy.load(f)[0]

def main():
    out=ROOT/'results/triangle';A=load(out/'fault-tolerant/A-case1.qpy');Q=load(out/'fault-tolerant/Q-case1.qpy')
    native=json.loads((out/'native-receipt.json').read_text());ft=json.loads((out/'fault-tolerant/receipt.json').read_text())
    original=next(x['expected_probability'] for x in native['checks'] if x['case']==1 and x['grover_power']==0)
    backend=AerSimulator(method='statevector',max_parallel_threads=1,max_memory_mb=12000);checks=[];start=time.perf_counter()
    for k in [0,1]:
        qc=QuantumCircuit(Q.num_qubits);qc.compose(A,range(A.num_qubits),inplace=True)
        for _ in range(k):qc.compose(Q,inplace=True)
        qc.save_expectation_value(Pauli('Z'),[0],label='z');r=backend.run(qc,shots=1).result();assert r.success
        p=(1+float(r.data(0)['z']))/2
        if k==0:base=p;expected=original
        else:expected=math.sin((2*k+1)*math.asin(math.sqrt(base)))**2
        error=abs(p-expected)
        assert error<1e-8
        checks.append(dict(grover_power=k,executed_qubits=qc.num_qubits,probability=p,expected=expected,error=error,operations=len(qc.data)-1))
    receipt=dict(status='PASS_ACTUAL_CLIFFORD_T_EXECUTION',checks=checks,wall_seconds=time.perf_counter()-start,
                 scope='Noiseless numerical statevector of emitted ten-qubit circuit at k0 andk1. No high-power adaptive schedule execution or QPU result.')
    (out/'compiled-FT-check.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))

if __name__=='__main__':main()
