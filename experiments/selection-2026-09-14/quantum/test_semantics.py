"""Independent operator checks of lookup, workspace echo and binary rotations."""
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[key]='1'
import json,math,time
from pathlib import Path
import numpy as np
from qiskit.quantum_info import Operator
from qiskit import transpile
from select_swap import lookup,dirty_lookup,rotation_oracle
from synthesis import exact_clifford_t,counts,save,Synthesizer
ROOT=Path(__file__).resolve().parent

def main():
    start=time.perf_counter();records=[];(ROOT/'results/semantic-circuits').mkdir(exist_ok=True)
    for label,words,b,lam,dirty in [('clean_select',[1,2,3,0,2,1,0,3],2,1,False),('clean_swap',[1,2,3,0],2,2,False),('dirty_echo',[1,2,3,0],2,2,True)]:
        qc,meta=dirty_lookup(words,b,lam) if dirty else lookup(words,b,lam)
        native=exact_clifford_t(qc);U=Operator(native).data;n=meta['address_qubits'];err=0.
        # All basis states, including arbitrary target/dirty values; prefix ancillas start0.
        for x in range(len(words)):
            for data in range(2**(lam*b)):
                for output in range(2**b if dirty else 1):
                    i=x+(data<<n)+(output<<(qc.num_qubits-b) if dirty else 0)
                    if dirty:j=x+(data<<n)+((output^words[x])<<(qc.num_qubits-b))
                    else:
                        lanes=[(data>>(k*b))&((1<<b)-1) for k in range(lam)]
                        for k in range(lam):lanes[k]^=words[(x//lam)*lam+k]
                        for level in range((lam-1).bit_length()):
                            if x>>level&1:
                                step=2**level
                                for q in range(0,lam,2*step):lanes[q],lanes[q+step]=lanes[q+step],lanes[q]
                        j=x+(sum(a<<(k*b) for k,a in enumerate(lanes))<<n)
                    target=np.zeros(len(U),complex);target[j]=1
                    err=max(err,float(np.max(abs(U[:,i]-target))))
        assert err<1e-9
        if not dirty:
            raw=transpile(qc,basis_gates=['h','s','sdg','t','tdg','x','z','cx'],optimization_level=0,num_processes=1)
            literal=4*sum(v for k,v in qc.count_ops().items() if k in ['rccx','rccx_dg'])+7*qc.count_ops().get('cswap',0)
            assert literal==meta['lookup_T_exact_before_optimization']
            assert counts(raw)['T']<=literal
            assert counts(native)['T']<=meta['lookup_T_exact_before_optimization']
        save(native,ROOT/'results/semantic-circuits'/f'{label}.qpy')
        records.append(dict(label=label,semantic_max_error=err,metadata=meta,compiled=counts(native)))
    # A negative angle and2pi check conditional signs. Use representable low-bit values.
    theta=np.array([0.,-math.pi,2*math.pi,math.pi/2]);qc,meta,words=rotation_oracle(theta,3,2)
    U=Operator(qc).data;n=2;targetbit=qc.num_qubits-1;err=0.
    for x,angle in enumerate(theta):
        rot=np.array([[np.cos(angle/2),-np.sin(angle/2)],[np.sin(angle/2),np.cos(angle/2)]])
        for bit in range(2):
            v=np.zeros(len(U),complex);v[x]=rot[0,bit];v[x+(1<<targetbit)]=rot[1,bit]
            err=max(err,float(np.max(abs(U[:,x+(bit<<targetbit)]-v))))
    assert err<1e-9
    records.append(dict(label='signed_rotation_and_clean_uncompute',semantic_max_error=err,words=words.tolist()))
    synth=Synthesizer(1e-10,ROOT/'results/synthesis-fixture')
    for angle in [.137,-1.193,math.pi/4,2*math.pi]:synth.rz(angle)
    result=dict(status='PASS',seconds=time.perf_counter()-start,checks=records,rotation_synthesis=synth.records)
    (ROOT/'results/semantic-tests.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
