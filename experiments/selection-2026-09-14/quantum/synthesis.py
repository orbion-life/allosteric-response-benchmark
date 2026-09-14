"""Actual deterministic Clifford+T strings, operator-norm checked per rotation."""
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[key]='1'
import hashlib,json,math,sys,time
from pathlib import Path
import numpy as np
from qiskit import QuantumCircuit,transpile,qpy
from qiskit.quantum_info import Operator
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'runtime'))
from pygridsynth.gridsynth import gridsynth_gates

def counts(qc):
    ops={str(k):int(v) for k,v in qc.count_ops().items()}
    return dict(qubits=qc.num_qubits,operations=sum(ops.values()),ops=ops,T=ops.get('t',0)+ops.get('tdg',0),CX=ops.get('cx',0),
                depth=qc.depth(),T_depth=qc.depth(filter_function=lambda x:x.operation.name in ['t','tdg']))

class Synthesizer:
    def __init__(self,epsilon,path):
        self.epsilon=float(epsilon);self.path=Path(path);self.path.mkdir(parents=True,exist_ok=True);self.cache={};self.records=[]
    def rz(self,angle):
        angle=float(angle);key=angle.hex()
        if key in self.cache:return self.cache[key]
        qc=QuantumCircuit(1);t0=time.perf_counter();multiple=round(angle/(math.pi/4))
        if abs(angle-multiple*math.pi/4)<1e-14:
            # Rz(k*pi/4)=exp(-ik*pi/8)T^k. Global phase is required.
            qc.global_phase=-multiple*math.pi/8;k=multiple%8
            if k>=4:qc.z(0);k-=4
            if k>=2:qc.s(0);k-=2
            if k:qc.t(0)
            gates=None
        else:
            gates=gridsynth_gates(theta=str(angle),epsilon=str(self.epsilon),seed=1729,dps=80,dtimeout=1000,ftimeout=1000)
            for gate in reversed(gates):
                if gate=='W':qc.global_phase+=math.pi/4
                elif gate=='H':qc.h(0)
                elif gate=='T':qc.t(0)
                elif gate=='S':qc.s(0)
                elif gate=='X':qc.x(0)
                else:raise ValueError(gate)
        reference=np.diag(np.exp(np.array([-1j,1j])*angle/2));error=float(np.linalg.norm(Operator(qc).data-reference,2))
        if error>self.epsilon+2e-14:raise AssertionError(('synthesis_error',angle,error,self.epsilon))
        record=dict(angle=angle,angle_hex=key,epsilon=self.epsilon,gates=gates,error=error,seconds=time.perf_counter()-t0,**counts(qc))
        self.cache[key]=qc;self.records.append(record)
        (self.path/'rotation-records.json').write_text(json.dumps(self.records,indent=2)+'\n')
        return qc
    def ry(self,angle):
        qc=QuantumCircuit(1);qc.sdg(0);qc.h(0);qc.compose(self.rz(angle),inplace=True);qc.h(0);qc.s(0);return qc
    def compile(self,qc):
        out=QuantumCircuit(qc.num_qubits);out.global_phase=qc.global_phase
        for inst in qc.data:
            name=inst.operation.name;bits=[qc.find_bit(q).index for q in inst.qubits];p=inst.operation.params
            if name in ['h','s','sdg','t','tdg','x','y','z','cx','id']:
                if name!='id':out.append(inst.operation,bits)
            elif name=='ry':out.compose(self.ry(float(p[0])),bits,inplace=True)
            elif name=='rz':out.compose(self.rz(float(p[0])),bits,inplace=True)
            elif name=='u':
                theta,phi,lam=map(float,p);out.global_phase+=(phi+lam)/2
                out.compose(self.rz(lam),bits,inplace=True);out.compose(self.ry(theta),bits,inplace=True);out.compose(self.rz(phi),bits,inplace=True)
            else:raise ValueError('Unsupported gate '+name)
        return out

def exact_clifford_t(qc):
    out=transpile(qc,basis_gates=['h','s','sdg','t','tdg','x','z','cx'],optimization_level=1,num_processes=1,seed_transpiler=1729)
    assert set(out.count_ops())<={'h','s','sdg','t','tdg','x','z','cx'}
    return out

def save(qc,path):
    path=Path(path)
    with path.open('wb') as f:qpy.dump(qc,f)
