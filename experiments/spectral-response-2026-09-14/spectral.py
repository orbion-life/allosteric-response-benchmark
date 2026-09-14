"""Actual controlled spectral response preparation and exact Grover reflections."""
import sys
from pathlib import Path
import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit.circuit.library import UCRYGate, RCCXGate
sys.path.insert(0, str(Path(__file__).resolve().parent / 'vendor'))
from multiplexed import real_loader


def native(qc):
    return transpile(qc, basis_gates=['ry','cx','h','x','z'], optimization_level=0,
                     seed_transpiler=1729, num_processes=1)


def observable(vector, n):
    padded = np.zeros(2**n)
    padded[:len(vector)] = vector / np.linalg.norm(vector)
    return native(real_loader(padded[None, :], list(range(1, n+1)), [], n+2, 0))


def attenuation(eigenvalues, time, n):
    f = np.zeros(2**n)
    f[:len(eigenvalues)] = np.exp(-time * eigenvalues)
    if np.min(f) < 0 or np.max(f) > 1 + 1e-12:
        raise ValueError('Spectral attenuation requires an admitted positive operator')
    angles = np.r_[np.zeros(2**n), 2 * np.arccos(np.clip(f, 0, 1))]
    qc = QuantumCircuit(n+2)
    qc.append(UCRYGate(angles.tolist()), [n+1] + list(range(1, n+1)) + [0])
    return native(qc)


def preparation(pi, pj, ft):
    qc = QuantumCircuit(pi.num_qubits, name='spectral_response_A')
    qc.h(0)
    qc.compose(pj, inplace=True)
    qc.compose(ft, inplace=True)
    qc.compose(pi.inverse(), inplace=True)
    qc.h(0)
    return qc


def zero_reflection(active):
    """I-2|0><0| on active wires, with active-2 clean restored helpers."""
    total = active + max(0, active-2)
    qc = QuantumCircuit(total, name='exact_zero_reflection')
    wires = list(range(active))
    work = list(range(active,total))
    qc.x(wires)
    controls, target = wires[:-1], wires[-1]
    qc.append(RCCXGate(), [controls[0], controls[1], work[0]])
    for i in range(2, len(controls)):
        qc.append(RCCXGate(), [work[i-2], controls[i], work[i-1]])
    qc.cz(work[len(controls)-2], target)
    for i in range(len(controls)-1, 1, -1):
        qc.append(RCCXGate().inverse(), [work[i-2], controls[i], work[i-1]])
    qc.append(RCCXGate().inverse(), [controls[0], controls[1], work[0]])
    qc.x(wires)
    return transpile(qc,basis_gates=['h','s','sdg','t','tdg','x','z','cx'],optimization_level=0,num_processes=1)


def grover(A):
    reflection = zero_reflection(A.num_qubits)
    qc = QuantumCircuit(reflection.num_qubits, name='spectral_response_Q')
    qc.global_phase = np.pi
    qc.x(0);qc.z(0);qc.x(0)  # Objective is readout zero.
    qc.compose(A.inverse(), list(range(A.num_qubits)), inplace=True)
    qc.compose(reflection, inplace=True)
    qc.compose(A, list(range(A.num_qubits)), inplace=True)
    return qc, reflection


def counts(qc):
    ops = {str(k):int(v) for k,v in qc.count_ops().items()}
    return dict(qubits=qc.num_qubits,operations=sum(ops.values()),ops=ops,
                CX=ops.get('cx',0),T=ops.get('t',0)+ops.get('tdg',0),
                rotations=ops.get('ry',0),depth=qc.depth())
