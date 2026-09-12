#!/usr/bin/env python3
"""Compile and simulate the complete frozen protein-derived overlap workload.

Run in the pinned environment: python run.py [--skip-noise] [--output RESULTS]
The default run includes all three independent orthonormal-basis pairs, every
logical preparation/walk/unpreparation gate, all-to-all and line synthesis,
finite-shot ideal sampling, and per-gate depolarizing density-matrix simulation.
"""
from pathlib import Path
import argparse
import hashlib
import json
import platform
import resource
import sys
import time
import warnings
import numpy as np
import scipy
import psutil
import qiskit
import qiskit_aer
from qiskit import qpy, qasm3
from qiskit_aer import AerSimulator
from qiskit_aer.noise import NoiseModel, depolarizing_error
from scipy.linalg import expm
from circuit import overlap_circuit, compile_circuit, resource_counts, unmeasured_and_readout, SEED
from model import prepare_fixture, mm

HERE = Path(__file__).resolve().parent
warnings.filterwarnings('error', category=RuntimeWarning)
np.seterr(divide='raise', invalid='raise', over='raise')


def write_json(path, data):
    Path(path).write_text(json.dumps(data, indent=2, allow_nan=False)+'\n')


def peak_rss_bytes():
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(value if sys.platform == 'darwin' else value*1024)


def scores_and_ranking(C, fixture):
    score = np.sqrt(np.mean(C[fixture['receiver'],:]**2, axis=0))
    candidates = np.flatnonzero(fixture['candidate'])
    ranking = sorted(candidates, key=lambda i:(-score[i],fixture['canonical'][i]))
    return score, [int(fixture['canonical'][i]) for i in ranking[:5]]


def reconstruct(means, fixture, summary):
    D, sh = fixture['reconstruction'], fixture['harmonic_sd']
    delayed = summary['retained_coefficient_sum']*mm(mm(D,means),D.T)
    R = -summary['beta']*(fixture['covariance']-delayed)
    return R/(summary['beta']*np.outer(sh,sh))


def readout_probabilities(compiled, method='statevector', noise_model=None):
    bare, readout = unmeasured_and_readout(compiled)
    bare.save_probabilities([readout], label='readout')
    backend = AerSimulator(method=method, noise_model=noise_model,
                           max_parallel_threads=1, fusion_enable=False)
    start = time.perf_counter()
    result = backend.run(bare, shots=1, seed_simulator=SEED).result()
    probs = np.asarray(result.data(0)['readout'], float)
    assert np.isfinite(probs).all() and abs(probs.sum()-1)<1e-9
    return probs, time.perf_counter()-start


def save_circuit(path, circuit):
    with Path(path).with_suffix('.qpy').open('wb') as handle:
        qpy.dump(circuit,handle)
    Path(path).with_suffix('.qasm').write_text(qasm3.dumps(circuit))


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--skip-noise', action='store_true')
    parser.add_argument('--walk-implementation', choices=['factorized','naive'], default='factorized')
    parser.add_argument('--output', type=Path, default=HERE/'results')
    args=parser.parse_args()
    output=args.output.resolve();output.mkdir(parents=True,exist_ok=True)
    (output/'circuits').mkdir(exist_ok=True)
    start=time.perf_counter()
    start_rss=psutil.Process().memory_info().rss
    fixture,summary=prepare_fixture(HERE/'inputs',output)
    rank=fixture['orthonormal_basis'].shape[1]
    pairs=[(i,j) for i in range(rank) for j in range(i,rank)]
    shots=int(summary['selected_sampling']['shots_per_overlap'])
    if shots>1_000_000:
        raise RuntimeError('The predeclared local shot cap is exceeded; report the derived cost instead of silently reducing precision')
    exact_means=np.zeros((rank,rank));sample_means=np.zeros((rank,rank))
    all_circuits=[];records=[]
    for serial,(i,j) in enumerate(pairs):
        pair_start=time.perf_counter()
        logical,components,logical_meta=overlap_circuit(fixture['P'],fixture['coefficients'],
                                fixture['orthonormal_basis'][:,i],fixture['orthonormal_basis'][:,j],
                                name=f'basis_{i}_{j}', walk_implementation=args.walk_implementation)
        synth_start=time.perf_counter();compiled=compile_circuit(logical)
        synthesis_seconds=time.perf_counter()-synth_start
        line_start=time.perf_counter();line=compile_circuit(logical,'line')
        line_seconds=time.perf_counter()-line_start
        counts=resource_counts(compiled);line_counts=resource_counts(line)
        save_circuit(output/'circuits'/f'basis_{i}_{j}_all_to_all',compiled)
        save_circuit(output/'circuits'/f'basis_{i}_{j}_line',line)
        with (output/'circuits'/f'basis_{i}_{j}_logical.qpy').open('wb') as handle:qpy.dump(logical,handle)
        probs, ideal_seconds=readout_probabilities(compiled)
        line_probs, line_ideal_seconds=readout_probabilities(line)
        chi=float(probs[0]-probs[1]);reference=float(fixture['overlap_polynomial'][i,j])
        if abs(chi-reference)>1e-9 or abs(chi-(line_probs[0]-line_probs[1]))>1e-9:
            raise AssertionError('Compiled circuit readout does not match the target polynomial overlap')
        sampler=AerSimulator(method='statevector',max_parallel_threads=1,fusion_enable=False)
        shot_start=time.perf_counter()
        sample_result=sampler.run(compiled,shots=shots,seed_simulator=SEED+serial).result()
        observed=sample_result.get_counts(0)
        sampling_seconds=time.perf_counter()-shot_start
        mean=(observed.get('0',0)-observed.get('1',0))/shots
        exact_means[i,j]=exact_means[j,i]=chi
        sample_means[i,j]=sample_means[j,i]=mean
        # Isolated component counts disclose loading and controlled-walk costs.
        # They need not sum to whole-circuit optimized totals across boundaries.
        component_costs={}
        comp_start=time.perf_counter()
        for label,component in components.items():
            component_costs[label]=resource_counts(compile_circuit(component))
        component_seconds=time.perf_counter()-comp_start
        record={'pair':[i,j],'logical':logical_meta,'all_to_all':counts,'bidirectional_line':line_counts,
                'component_counts':component_costs,
                'component_counts_scope':'isolated decompositions; do not sum these to obtain optimized whole-circuit totals',
                'ideal_readout_probabilities':probs.tolist(),'line_readout_probabilities':line_probs.tolist(),
                'ideal_chi':chi,'polynomial_reference_chi':reference,'absolute_chi_error':abs(chi-reference),
                'executed_ideal_shots':shots,'seed_simulator':SEED+serial,'observed_counts':observed,
                'sample_mean':float(mean),'sampling_overlap_error':abs(mean-chi),
                'synthesis_seconds_all_to_all':synthesis_seconds,'synthesis_seconds_line':line_seconds,
                'isolated_component_synthesis_seconds':component_seconds,
                'ideal_simulation_seconds_all_to_all':ideal_seconds,'ideal_simulation_seconds_line':line_ideal_seconds,
                'finite_shot_simulation_seconds':sampling_seconds,'total_pair_seconds':time.perf_counter()-pair_start}
        records.append(record);all_circuits.append(compiled)
        write_json(output/'progress.json',{'status':'running','completed_pairs':records})
        print(json.dumps({'pair':[i,j],'qubits':counts['qubits'],'cx':counts['two_qubit_cx_gates'],
                          'depth':counts['depth_including_measurement'],'chi_error':record['absolute_chi_error'],
                          'elapsed_seconds':record['total_pair_seconds']}),flush=True)
    C_ideal=reconstruct(exact_means,fixture,summary)
    C_sample=reconstruct(sample_means,fixture,summary)
    score_ideal,top_ideal=scores_and_ranking(C_ideal,fixture)
    score_sample,top_sample=scores_and_ranking(C_sample,fixture)
    score_exact,top_exact=scores_and_ranking(fixture['C_exact'],fixture)
    score_half=fixture['sampling_score_halfwidth']+fixture['score_tail_bound']
    ordered=fixture['ranking'];top=ordered[:5];others=ordered[5:]
    observed_lower=float(np.min(score_sample[top]-score_half[top]))
    observed_upper=float(np.max(score_sample[others]+score_half[others]))
    ideal_summary={'maximum_compiled_vs_polynomial_C_error':float(np.max(abs(C_ideal-fixture['C_polynomial']))),
                   'maximum_compiled_vs_exact_grid_C_error':float(np.max(abs(C_ideal-fixture['C_exact']))),
                   'maximum_sample_vs_exact_grid_C_error':float(np.max(abs(C_sample-fixture['C_exact']))),
                   'exact_top5_canonical':top_exact,'compiled_top5_canonical':top_ideal,'sample_top5_canonical':top_sample,
                   'observed_top5_lower_bound':observed_lower,'observed_outside_upper_bound':observed_upper,
                   'observed_ranking_intervals_separated':observed_lower>observed_upper,
                   'interval_scope':'95% simultaneous sampling plus analytic truncation; same fixed grid and physical slice, no biological or physical-reduction coverage',
                   'total_executed_ideal_shots':shots*len(pairs)}
    if ideal_summary['maximum_compiled_vs_polynomial_C_error']>1e-10:
        raise AssertionError('Classical reconstruction of compiled circuit failed')
    np.savez(output/'ideal-response-matrices.npz',exact_grid_C=fixture['C_exact'],compiled_C=C_ideal,
             sampled_C=C_sample,exact_means=exact_means,sampled_means=sample_means,
             exact_scores=score_exact,sampled_scores=score_sample,score_halfwidth=score_half)
    write_json(output/'ideal-results.json',{'summary':ideal_summary,'pairs':records})
    # Same-grid classical comparator includes finite-time propagation and the
    # complete N x N reconstruction; inputs/normalization are shared.
    classical_times=[]
    for repeat in range(30):
        t0=time.perf_counter()
        direct=mm(mm(fixture['orthonormal_basis'].T,expm(-summary['tau']*fixture['H'])),fixture['orthonormal_basis'])
        classical_C=reconstruct(direct/summary['retained_coefficient_sum'],fixture,summary)
        classical_times.append(time.perf_counter()-t0)
    assert np.max(abs(classical_C-fixture['C_exact']))<1e-10
    classical={'repeats':30,'median_seconds':float(np.median(classical_times)),
               'p90_seconds':float(np.quantile(classical_times,.9)),
               'times_seconds':classical_times,'scope':'small-matrix exponential plus all-residue reconstruction, same grid; shared preprocessing excluded from both solve timings'}
    noise_records=[]
    if not args.skip_noise:
        for p in (1e-4,1e-3,1e-2):
            noise=NoiseModel()
            noise.add_all_qubit_quantum_error(depolarizing_error(p,1),['u'])
            noise.add_all_qubit_quantum_error(depolarizing_error(p,2),['cx'])
            means=np.zeros((rank,rank));times=[];channel_bounds=np.zeros_like(means)
            for (i,j),compiled,record in zip(pairs,all_circuits,records):
                probs,seconds=readout_probabilities(compiled,method='density_matrix',noise_model=noise)
                means[i,j]=means[j,i]=probs[0]-probs[1]
                gate_count=record['all_to_all']['one_qubit_u_gates']+record['all_to_all']['two_qubit_cx_gates']
                channel_bounds[i,j]=channel_bounds[j,i]=2*(-math_expm1_gate(p,gate_count))
                times.append(seconds)
            for readout in (0.,.01,.03):
                noisy_means=(1-2*readout)*means
                C_noise=reconstruct(noisy_means,fixture,summary)
                noise_score,noise_top=scores_and_ranking(C_noise,fixture)
                delta=np.abs(C_noise-C_ideal)
                # Gate-mixture bound: at least (1-p)^g weight is the ideal
                # circuit; arbitrary residual branches alter a +/-1 mean by <=2.
                bounds=np.minimum(2., channel_bounds+2*readout)
                Dabs=abs(fixture['reconstruction'])
                C_bound=summary['retained_coefficient_sum']*mm(mm(Dabs,bounds),Dabs.T)/np.outer(fixture['harmonic_sd'],fixture['harmonic_sd'])
                score_bias_bound=np.sqrt(np.mean(C_bound[fixture['receiver'],:]**2,axis=0))
                noise_record={'per_gate_depolarizing_p':p,'independent_readout_flip':readout,
                       'noisy_overlap_means':noisy_means.tolist(),
                       'max_C_device_bias_vs_ideal_circuit':float(delta.max()),
                       'max_C_error_vs_exact_grid':float(np.max(abs(C_noise-fixture['C_exact']))),
                       'top5_canonical':noise_top,'observed_top5_unchanged':noise_top==top_exact,
                       'max_score_bias':float(np.max(abs(noise_score-score_ideal))),
                       'device_C_bias_budget':.001,'passes_device_C_bias_budget':bool(delta.max()<=.001),
                       'worst_case_channel_C_bias_bound':float(C_bound.max()),
                       'rank_bound_separated_with_gate_mixture_bias':bool(np.min(score_sample[top]-score_half[top]-score_bias_bound[top]) >
                                 np.max(score_sample[others]+score_half[others]+score_bias_bound[others])),
                       'density_simulation_seconds_per_pair':times,
                       'readout_treatment':'exact affine transformation of Z expectation for independent final classical bit flips; not extra circuit execution',
                       'noise_scope':'synthetic local depolarization after every synthesized u/cx; no device calibration, T1/T2 or hardware execution'}
                noise_records.append(noise_record)
                np.savez(output/f'noise-p{p:g}-readout{readout:g}.npz',C=C_noise,score=noise_score,overlap=noisy_means)
            print(json.dumps({'noise_p':p,'density_seconds':sum(times),'last_bias':noise_records[-1]['max_C_device_bias_vs_ideal_circuit']}),flush=True)
            write_json(output/'noise-results.json',noise_records)
    operation_totals={}
    for topology in ('all_to_all','bidirectional_line'):
        operation_totals[topology]={key:int(sum(r[topology][key]*shots for r in records))
                                   for key in ('one_qubit_u_gates','two_qubit_cx_gates','measurements')}
    result={'status':'complete','scope':'actual synthesized and simulated complete protein-derived circuit; not hardware or validated physical compression',
            'environment':{'python':platform.python_version(),'platform':platform.platform(),
                           'qiskit':qiskit.__version__,'qiskit_aer':qiskit_aer.__version__,
                           'numpy':np.__version__,'scipy':scipy.__version__},
            'compiler':{'basis_gates':['u','cx'],'optimization_level':1,'approximation_degree':1.0,
                        'seed':SEED,'topologies':['all-to-all','bidirectional seven-qubit line'],
                        'all_to_all_is_hardware_calibrated':False},
            'simulator':{'max_parallel_threads':1,'fusion_enable':False,'actual_finite_shots':shots*len(pairs)},
            'fixture':summary,'ideal':ideal_summary,'compiled_pairs':records,
            'full_sampling_operation_totals':operation_totals,'classical_comparator':classical,
            'noise_cases':noise_records,'noise_was_skipped':args.skip_noise,
            'total_run_seconds':time.perf_counter()-start,'process_rss_at_start_bytes':start_rss,
            'process_peak_rss_bytes':peak_rss_bytes(),
            'memory_scope':'process high-water mark including Python, Qiskit, synthesis and simulators; not an isolated circuit-state allocation',
            'unmeasured_costs':['device pulse durations','physical T1/T2','queue time','actual hardware execution','target-scale state preparation'],
            'source_sha256':{name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in ('run.py','model.py','circuit.py')}}
    write_json(output/'quantum-results.json',result)
    print(json.dumps({'status':result['status'],'total_seconds':result['total_run_seconds'],
                      'ideal':ideal_summary,'sampling_operation_totals':operation_totals,
                      'peak_rss_bytes':result['process_peak_rss_bytes']},indent=2))


def math_expm1_gate(p,count):
    # Stable (1-p)^count - 1, with saturation at -1 at large counts.
    return float(np.expm1(count*np.log1p(-p)))


if __name__=='__main__':
    main()
