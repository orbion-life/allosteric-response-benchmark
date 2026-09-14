"""Exact fixed-width controlled ripple adders and Fourier phase kickback."""
import sys
from pathlib import Path
import numpy as np
from qiskit import QuantumCircuit,transpile
from qiskit.circuit.library import CDKMRippleCarryAdder,RCCXGate
sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
from synthesis import exact_clifford_t

def controlled_adder(b):
    adder=transpile(CDKMRippleCarryAdder(b,kind='fixed'),basis_gates=['x','cx','ccx'],optimization_level=0,num_processes=1)
    assert adder.num_qubits==2*b+1
    qc=QuantumCircuit(2*b+3,name='explicit_controlled_CDKM');control=2*b+1;scratch=2*b+2
    for inst in adder.data:
        bits=[adder.find_bit(q).index for q in inst.qubits]
        if inst.operation.name=='x':qc.cx(control,bits[0])
        elif inst.operation.name=='cx':qc.ccx(control,bits[0],bits[1])
        elif inst.operation.name=='ccx':
            qc.append(RCCXGate(),[control,bits[0],scratch]);qc.ccx(scratch,bits[1],bits[2]);qc.append(RCCXGate().inverse(),[control,bits[0],scratch])
        else:raise ValueError(inst.operation.name)
    return exact_clifford_t(qc),dict(adder_boolean_ops=dict(adder.count_ops()),word_bits=b,adder_helper_qubits=1,control_helper_qubits=1)

def rotation_bank(b):
    add,meta=controlled_adder(b);qc=QuantumCircuit(add.num_qubits,name='phase_gradient_Ry_bank');target=2*b+1
    # S H maps Z to Y. In temporal order the inverse basis change comes first.
    qc.sdg(target);qc.h(target)
    qc.x(target);qc.compose(add.inverse(),inplace=True);qc.x(target);qc.compose(add,inplace=True)
    qc.h(target);qc.s(target)
    return qc,add,meta

def fourier_state(b):
    return np.exp(-2j*np.pi*np.arange(2**b)/2**b)/np.sqrt(2**b)
