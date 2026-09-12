"""One bounded real gate compilation. Invoked by compile_primitives.py."""
import argparse
import json
import math
import resource
import platform
import time
from pathlib import Path

import numpy as np
import qiskit
from qiskit import QuantumCircuit, transpile
from qiskit.circuit.library import StatePreparation
from qiskit.quantum_info import Statevector
from scipy import sparse


def counts(circuit):
    operations = dict(circuit.count_ops())
    assert set(operations) <= {'u','cx'}
    return dict(qubits=circuit.num_qubits,one_qubit_u=int(operations.get('u',0)),
                two_qubit_cx=int(operations.get('cx',0)),depth=circuit.depth(),
                two_qubit_depth=circuit.depth(lambda item:len(item.qubits)==2),
                total_native_operations=int(sum(operations.values())))


def main(kind, index, root, output):
    started = time.perf_counter()
    arrays = np.load(root/'results/analysis-arrays.npz')
    results = json.loads((root/'results/results.json').read_text())
    padded = results['loading']['padded_states']
    n = int(math.log2(padded))
    if kind in ('observable','raw_observable'):
        values = np.zeros(padded)
        if kind == 'observable':
            values[:arrays['basis'].shape[0]] = arrays['basis'][:,index]
        else:
            primary=np.load(root/'inputs/primary.npz')
            column=np.sqrt(primary['pi'])*primary['monomial_centered'][:,index]
            values[:len(column)]=column/np.linalg.norm(column)
        controls, control_state = 1, 1
    elif kind == 'coefficient':
        coeff = arrays['fixed_coefficients' if index == 25 else 'oracle_coefficients']
        assert len(coeff)-1 == index
        values = np.zeros(2**math.ceil(math.log2(len(coeff))))
        values[:len(coeff)] = np.sqrt(coeff/coeff.sum())
        controls, control_state = 1, 1
    elif kind == 'transition_row':
        G = results['model']['states']
        P = sparse.csr_matrix((arrays['P_data'],arrays['P_indices'],arrays['P_indptr']),shape=(G,G))
        values = np.zeros(padded)
        if index < G:
            values[:G] = np.sqrt(P.getrow(index).toarray().ravel())
        else:
            values[index] = 1.
        controls, control_state = n, index
    else:
        raise ValueError(kind)
    norm = np.linalg.norm(values)
    assert abs(norm-1) < 1e-12
    values /= norm
    gate = StatePreparation(values.tolist(),label=f'{kind}_{index}')
    # Exact same generic controlled StatePreparation constructor as the public
    # circuit. It is not initialize/reset, and preserves controlled phases.
    controlled = gate.control(controls,ctrl_state=control_state,annotated=False)
    circuit = QuantumCircuit(controls+int(math.log2(len(values))))
    circuit.append(controlled,list(range(circuit.num_qubits)))
    constructor_seconds = time.perf_counter()-started
    compiled = transpile(circuit,basis_gates=['u','cx'],optimization_level=1,
                         seed_transpiler=20260912,approximation_degree=1.,num_processes=1)
    compiled_counts = counts(compiled)
    record = dict(kind=kind,index=index,qiskit_version=qiskit.__version__,
                  scope='Actual isolated controlled loading primitive; all-to-all u/cx, no device topology',
                  controls=controls,control_state=control_state,prepared_state_nonzeros=int(np.count_nonzero(values)),
                  constructor_seconds=constructor_seconds,compile_total_seconds=time.perf_counter()-started,
                  counts=compiled_counts)
    if kind != 'transition_row':
        # Coherence check tests the relative phase against the inactive b branch.
        prepare = QuantumCircuit(circuit.num_qubits)
        prepare.h(0)
        prepare.compose(compiled,inplace=True)
        actual = np.asarray(Statevector.from_instruction(prepare).data)
        expected = np.zeros(len(actual),complex)
        expected[0] = 1/math.sqrt(2)
        expected[1::2] = values/math.sqrt(2)
        record['phase_sensitive_statevector_max_error'] = float(np.max(abs(actual-expected)))
        assert record['phase_sensitive_statevector_max_error'] < 1e-9
        restored = Statevector(actual).evolve(compiled.inverse())
        expected_reset = np.zeros(len(actual),complex)
        expected_reset[:2] = 1/math.sqrt(2)
        record['coherent_inverse_max_error'] = float(np.max(abs(restored.data-expected_reset)))
        assert record['coherent_inverse_max_error'] < 1e-9
        # Native inverse exists without extra workspace, with the same gate
        # counts/depth; this count uses the inverse of the compiled unitary.
        record['native_inverse_counts'] = counts(compiled.inverse())
    else:
        record['statevector_check'] = 'Not attempted: this 22-qubit row primitive exceeds this proof scope; circuit synthesis only'
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if platform.system() != 'Darwin': rss *= 1024
    record['worker_peak_RSS_bytes'] = int(rss)
    record['worker_seconds'] = time.perf_counter()-started
    record['status'] = 'completed'
    output.write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('kind',choices=['observable','raw_observable','coefficient','transition_row'])
    parser.add_argument('index',type=int)
    parser.add_argument('--root',type=Path,default=Path(__file__).parent)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    main(args.kind,args.index,args.root,args.output)
