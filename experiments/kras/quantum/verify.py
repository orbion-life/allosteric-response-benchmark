#!/usr/bin/env python3
"""Verify phase-sensitive control factoring and persisted experiment artifacts."""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np
from qiskit import QuantumCircuit, qpy
from qiskit.quantum_info import Operator
from circuit import walk_gate, factorized_controlled_walk, overlap_circuit, compile_circuit, resource_counts
from run import readout_probabilities

HERE = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--results', type=Path, default=HERE/'results')
    args = parser.parse_args()
    result_dir = args.results.resolve()
    result = json.loads((result_dir/'quantum-results.json').read_text())
    fixture = np.load(result_dir/'fixture.npz')
    for name, expected in result['source_sha256'].items():
        assert hashlib.sha256((HERE/name).read_bytes()).hexdigest() == expected
    for item in result['fixture']['input_manifest']:
        assert hashlib.sha256((HERE/'inputs'/item['file']).read_bytes()).hexdigest() == item['sha256']
    W, _ = walk_gate(fixture['P'])
    walk_matrix = Operator(W).data
    projected_walk_errors = []
    walk_power = np.eye(16, dtype=complex)
    chebyshev = [np.eye(4), fixture['M']]
    for degree in range(2, len(fixture['coefficients'])):
        chebyshev.append(2*np.einsum('ij,jk->ik', fixture['M'], chebyshev[-1], optimize=False)-chebyshev[-2])
    for degree, target in enumerate(chebyshev):
        if degree: walk_power = np.einsum('ij,jk->ik',walk_power,walk_matrix,optimize=False)
        # X is the low register and Y the high register; Y=0 gives indices 0..3.
        projected_walk_errors.append(float(np.max(abs(walk_power[:4,:4]-target))))
    assert max(projected_walk_errors)<1e-11
    factored = factorized_controlled_walk(fixture['P'])
    # Matrix comparison includes the absolute global phase; equivalence up to
    # global phase would be insufficient inside this controlled Hadamard test.
    direct_matrix = Operator(W.control(2, annotated=False)).data
    factored_matrix = Operator(factored).data
    control_error = float(np.max(abs(direct_matrix-factored_matrix)))
    assert control_error < 1e-11
    # Explicit false-control sectors must implement identity on arbitrary XY.
    inactive_errors = []
    for control in (0,1,2):
        indices = np.arange(16)*4+control
        block = factored_matrix[np.ix_(indices,indices)]
        inactive_errors.append(float(np.max(abs(block-np.eye(16)))))
    assert max(inactive_errors) < 1e-11
    # A sign change in one observable reverses the overlap. This catches a
    # physically relevant phase being discarded by amplitude preparation.
    basis = fixture['orthonormal_basis']
    signed,_,_=overlap_circuit(fixture['P'],fixture['coefficients'],-basis[:,0],basis[:,0])
    probabilities,_=readout_probabilities(compile_circuit(signed))
    sign_error=float(abs(probabilities[0]-probabilities[1]+fixture['overlap_polynomial'][0,0]))
    assert sign_error < 1e-11
    persisted=[]
    for pair in result['compiled_pairs']:
        i,j=pair['pair']
        for label, key in [('all_to_all','all_to_all'),('line','bidirectional_line')]:
            path=result_dir/'circuits'/f'basis_{i}_{j}_{label}.qpy'
            with path.open('rb') as handle: circuit=qpy.load(handle)[0]
            assert resource_counts(circuit)==pair[key]
            probabilities,_=readout_probabilities(circuit)
            error=float(abs(probabilities[0]-probabilities[1]-pair['polynomial_reference_chi']))
            assert error<1e-11
            persisted.append({'pair':[i,j],'topology':label,'absolute_overlap_error':error,
                              'qpy_sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
    assert result['ideal']['observed_ranking_intervals_separated']
    assert result['ideal']['total_executed_ideal_shots']==151887
    assert result['fixture']['retained_algebraic_rank']==2
    assert all(not row['rank_bound_separated_with_gate_mixture_bias'] for row in result['noise_cases'])
    checks={'status':'passed','exact_control_factorization_max_entry_error':control_error,
            'projected_walk_power_vs_Chebyshev_errors':projected_walk_errors,
            'inactive_control_sector_identity_errors':inactive_errors,'signed_preparation_overlap_error':sign_error,
            'saved_native_circuits_reloaded_and_simulated':persisted,
            'results_sha256':hashlib.sha256((result_dir/'quantum-results.json').read_bytes()).hexdigest(),
            'verified_source_sha256':result['source_sha256'],
            'scope':'phase-sensitive algebra and persisted circuit repeat execution; no new protein accuracy claim'}
    (result_dir/'verification.json').write_text(json.dumps(checks,indent=2)+'\n')
    print(json.dumps(checks,indent=2))


if __name__=='__main__': main()
