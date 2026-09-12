"""Explicit logical gates for the Pulsar qubitized-walk overlap estimator.

No whole-estimator unitary or matrix exponential is inserted into the circuit.
V uses row-controlled amplitude preparation; the walk uses V, SWAP, V inverse,
and the zero-Y reflection. All five estimator components are controlled by b.
"""
import math
import numpy as np
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, transpile
from qiskit.circuit.library import StatePreparation, SwapGate
from qiskit.transpiler import CouplingMap

SEED = 20260912
BASIS = ["u", "cx"]


def unit_vector(values):
    v = np.asarray(values, dtype=complex)
    length = float(np.sqrt(np.sum(abs(v)**2)))
    if length <= 1e-14:
        raise ValueError("Zero observable norm: bypass its circuit and record zero contribution")
    return v / length


def preparation(values, label):
    # StatePreparation is unitary, unlike initialize/reset. Its complex amplitudes
    # include signs and phases; the controlled operation retains its global phase.
    return StatePreparation(unit_vector(values).tolist(), label=label)


def walk_gate(P):
    P = np.asarray(P, float)
    G = len(P)
    n = int(round(math.log2(G)))
    if 2**n != G or P.shape != (G, G):
        raise ValueError("Pad unused configurations as invariant states before constructing V")
    if P.min() < -1e-13 or np.max(abs(P.sum(axis=1)-1)) > 1e-12:
        raise ValueError("P must be stochastic")
    x, y = QuantumRegister(n, "x"), QuantumRegister(n, "y")
    row = QuantumCircuit(x, y, name="V")
    for state in range(G):
        gate = preparation(np.sqrt(np.maximum(P[state], 0)), f"row{state}")
        row.append(gate.control(n, ctrl_state=state, annotated=False), list(x)+list(y))
    V = row.to_gate(label="V")
    walk = QuantumCircuit(x, y, name="W")
    walk.append(V, list(x)+list(y))
    for i in range(n):
        walk.swap(x[i], y[i])
    walk.append(V.inverse(), list(x)+list(y))
    # X...controlled-Z...X implements I-2|0><0|. The global minus sign is
    # essential because W is subsequently controlled, and gives 2|0><0|-I.
    for bit in y:
        walk.x(bit)
    if n == 1:
        walk.z(y[0])
    else:
        walk.mcp(math.pi, list(y[:-1]), y[-1])
    for bit in y:
        walk.x(bit)
    walk.global_phase += math.pi
    return walk.to_gate(label="W"), V


def factorized_controlled_walk(P):
    """Exactly control R V† S V using unconditional V/V† conjugation.

    If either control is false, S and R are identities and V†V cancels.
    If both are true, the block applies precisely R V† S V. No loading gate
    is omitted: both V and V† are synthesized on every invocation.
    """
    n = int(round(math.log2(len(P))))
    controls, x, y = QuantumRegister(2, "ctrl"), QuantumRegister(n, "x"), QuantumRegister(n, "y")
    block = QuantumCircuit(controls, x, y, name="ccW_factorized")
    _, V = walk_gate(P)
    block.append(V, list(x)+list(y))
    controlled_swap = SwapGate().control(2, annotated=False)
    for bit in range(n):
        block.append(controlled_swap, list(controls)+[x[bit], y[bit]])
    block.append(V.inverse(), list(x)+list(y))
    reflection = QuantumCircuit(n, name="R_Y")
    for bit in range(n): reflection.x(bit)
    if n == 1: reflection.z(0)
    else: reflection.mcp(math.pi, list(range(n-1)), n-1)
    for bit in range(n): reflection.x(bit)
    reflection.global_phase += math.pi
    block.append(reflection.to_gate().control(2, annotated=False), list(controls)+list(y))
    return block.to_gate(label="ccW_factorized")


def overlap_circuit(P, coefficients, left, right, name="overlap", walk_implementation="factorized"):
    coeff = np.asarray(coefficients, float)
    if coeff.min() < 0 or coeff.sum() <= 0:
        raise ValueError("Coefficient preparation requires nonnegative weights with positive sum")
    K = len(coeff)-1
    a = max(1, int(math.ceil(math.log2(K+1))))
    n = int(round(math.log2(len(P))))
    b, c, x, y = (QuantumRegister(1, "b"), QuantumRegister(a, "c"),
                  QuantumRegister(n, "x"), QuantumRegister(n, "y"))
    out = ClassicalRegister(1, "out")
    circ = QuantumCircuit(b, c, x, y, out, name=name)
    amplitudes = np.zeros(2**a)
    amplitudes[:K+1] = np.sqrt(coeff / coeff.sum())
    coef = preparation(amplitudes, "COEF")
    prep_left, prep_right = preparation(left, "PREP_i"), preparation(right, "PREP_j")
    W, V = walk_gate(P)
    if walk_implementation == "factorized":
        controlled_walk = factorized_controlled_walk(P)
    elif walk_implementation == "naive":
        controlled_walk = W.control(2, annotated=False)
    else:
        raise ValueError(walk_implementation)
    components = {}

    def component(label, add):
        block = QuantumCircuit(b, c, x, y, name=label)
        add(block)
        circ.compose(block, qubits=list(b)+list(c)+list(x)+list(y), inplace=True)
        components[label] = block

    circ.h(b[0])
    component("controlled_PREP_j", lambda q: q.append(prep_right.control(1, annotated=False), list(b)+list(x)))
    component("controlled_COEF", lambda q: q.append(coef.control(1, annotated=False), list(b)+list(c)))

    def select(q):
        # On coefficient state k, binary-controlled powers apply W^k. Inactive
        # coefficient states have zero prepared amplitude and remain coherent.
        for bit in range(a):
            for _ in range(2**bit):
                q.append(controlled_walk, [b[0], c[bit]]+list(x)+list(y))

    component("controlled_SELECT", select)
    component("controlled_COEF_inverse", lambda q: q.append(coef.inverse().control(1, annotated=False), list(b)+list(c)))
    component("controlled_PREP_i_inverse", lambda q: q.append(prep_left.inverse().control(1, annotated=False), list(b)+list(x)))
    circ.h(b[0])
    circ.measure(b[0], out[0])
    metadata = {"configuration_qubits_each": n, "coefficient_qubits": a,
                "readout_qubits": 1, "arithmetic_workspace_qubits": 0,
                "walk_calls_in_binary_SELECT": 2**a-1,
                "K": K, "retained_sum": float(coeff.sum()),
                "controlled_walk_implementation": walk_implementation,
                "logical_qubits": 2*n+a+1,
                "state_preparation": "Qiskit unitary StatePreparation, including signed amplitudes",
                "walk_construction": "controlled row PREP, SWAP, inverse row PREP, zero-Y reflection"}
    return circ, components, metadata


def compile_circuit(circuit, topology="all_to_all"):
    options = dict(basis_gates=BASIS, optimization_level=1,
                   seed_transpiler=SEED, approximation_degree=1.0, num_processes=1)
    if topology == "line":
        options.update(coupling_map=CouplingMap.from_line(circuit.num_qubits, bidirectional=True),
                       initial_layout=list(range(circuit.num_qubits)), routing_method="sabre")
    elif topology != "all_to_all":
        raise ValueError(topology)
    compiled = transpile(circuit, **options)
    unexpected = set(compiled.count_ops()) - {"u", "cx", "measure", "barrier"}
    if unexpected:
        raise AssertionError(f"Circuit contains unsynthesized operations: {unexpected}")
    if compiled.num_qubits != circuit.num_qubits:
        raise AssertionError("Compiler introduced unreported workspace")
    return compiled


def resource_counts(circuit):
    counts = dict(circuit.count_ops())
    return {"qubits": circuit.num_qubits, "one_qubit_u_gates": int(counts.get("u", 0)),
            "two_qubit_cx_gates": int(counts.get("cx", 0)),
            "measurements": int(counts.get("measure", 0)),
            "depth_including_measurement": circuit.depth(),
            "two_qubit_depth": circuit.depth(lambda item: len(item.qubits) == 2),
            "operation_counts": {k: int(v) for k,v in counts.items()},
            "gate_set_scope": "abstract universal u/cx gate basis; not a named backend pulse schedule"}


def unmeasured_and_readout(circuit):
    readout = [circuit.find_bit(inst.qubits[0]).index for inst in circuit.data
               if inst.operation.name == "measure"]
    if len(readout) != 1:
        raise ValueError("Expected one final readout measurement")
    return circuit.remove_final_measurements(inplace=False), readout[0]
