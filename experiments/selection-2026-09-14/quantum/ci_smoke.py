"""Source-only lookup/phase tests. No archive, synthesis runtime or writes required."""
import json
import math
import numpy as np
from qiskit.quantum_info import Operator
from select_swap import lookup, dirty_lookup, rotation_oracle


def main():
    errors = {}
    words = [1, 2, 3, 0]
    clean, meta = lookup(words, 2, 2)
    unitary = Operator(clean).data
    error = 0.0
    for address in range(4):
        lanes = words[2 * (address // 2):2 * (address // 2) + 2]
        if address & 1:
            lanes = lanes[::-1]
        expected = np.zeros(len(unitary), complex)
        expected[address + (lanes[0] << 2) + (lanes[1] << 4)] = 1
        error = max(error, float(np.max(abs(unitary[:, address] - expected))))
    errors['clean_lookup_with_swap'] = error

    dirty, meta = dirty_lookup(words, 2, 2)
    unitary = Operator(dirty).data
    error = 0.0
    for address in range(4):
        for dirty_data in range(16):
            for output in range(4):
                start = address + (dirty_data << 2) + (output << 6)
                finish = address + (dirty_data << 2) + ((output ^ words[address]) << 6)
                expected = np.zeros(len(unitary), complex)
                expected[finish] = 1
                error = max(error, float(np.max(abs(unitary[:, start] - expected))))
    errors['dirty_echo_restores_arbitrary_workspace_and_phase'] = error

    angles = np.array([0.0, -math.pi, 2 * math.pi, math.pi / 2])
    oracle, _, _ = rotation_oracle(angles, 3, 2)
    unitary = Operator(oracle).data
    target = oracle.num_qubits - 1
    error = 0.0
    for address, angle in enumerate(angles):
        rotation = np.array([[np.cos(angle / 2), -np.sin(angle / 2)],
                             [np.sin(angle / 2), np.cos(angle / 2)]])
        for bit in range(2):
            expected = np.zeros(len(unitary), complex)
            expected[address] = rotation[0, bit]
            expected[address + (1 << target)] = rotation[1, bit]
            error = max(error, float(np.max(abs(unitary[:, address + (bit << target)] - expected))))
    errors['signed_rotation_and_clean_uncompute'] = error
    assert max(errors.values()) < 1e-9, errors
    print(json.dumps({'status': 'PASS', 'scope': 'Source-only tiny exact semantic tests; no cost or full-size execution claim.', 'errors': errors}, indent=2))


if __name__ == '__main__':
    main()
