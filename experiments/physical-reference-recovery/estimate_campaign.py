"""Measured-input planning estimate; never runs a grid or propagator."""
import hashlib
import argparse
import json
import math
from pathlib import Path
from chebyshev import required_degree

ROOT=Path(__file__).resolve().parent


def read(name):
    return json.loads((ROOT/name).read_text())


def run():
    inputs=['nonlinear-primary/nonlinear-h0.5/receipt.json',
            'nonlinear-primary/nonlinear-h0.25/receipt.json',
            'nonlinear-primary/watchdog.json','capacity-primary/case/receipt.json',
            'capacity-primary/watchdog.json','supplementary-minima-preflight.json',
            'harmonic-primary/watchdog.json']
    a,b=[read(x) for x in inputs[:2]]
    outer=read(inputs[2]);cap=read(inputs[3]);capouter=read(inputs[4])
    grids=read(inputs[5])['cases'][0]['nonlinear_grids']
    profiles=[a,b]
    seconds_per_nnz=max([p['action_microbenchmark']['csr_ten_three_column_actions_seconds']/10/p['actual_csr_nonzeros'] for p in profiles]+
                        [cap['matrix_action_measurements']['csr']['maximum']/cap['actual_csr_entries']])
    overhead_taylor=max(query['seconds']/(query['counts']['vector_equivalents']/3*(p['action_microbenchmark']['csr_ten_three_column_actions_seconds']/10))
                        for p in profiles for query in p['propagation'])
    overhead_chebyshev=max(query['seconds']/(query['three_column_actions']*(p['action_microbenchmark']['csr_ten_three_column_actions_seconds']/10))
                           for p in profiles for query in p['independent_chebyshev'])
    count_per_spectral=[max(p['propagation'][j]['counts']['vector_equivalents']/p['independent_chebyshev'][j]['spectral_upper_bound'] for p in profiles) for j in range(3)]
    preparation_per_state=max([sum(p['assembly_and_preprocessing'].values())/p['states'] for p in profiles]+
                               [(cap['grid_seconds']+cap['assembly_seconds'])/cap['states']])
    fixed_case_overhead=max([job['outer_worker_wall_seconds']-p['wall_seconds'] for job,p in zip(outer['jobs'],profiles)]+
                            [capouter['outer_wall_seconds']-cap['worker_timed_seconds']])
    safety_factor=2.
    cases=[]
    for grid in grids:
        # All h=1/8 boxes are aligned subsets of the measured largest lattice.
        # Removing outside hops lowers each retained absolute row sum.
        fine=grid['spacing']==.125
        spectral=cap['spectral_gershgorin_upper_bound'] if fine else grid['spectral_upper_bound_by_gershgorin_exit_bound']
        block_seconds=grid['csr_nonzeros']*seconds_per_nnz
        predicted_vectors=[math.ceil(x*spectral) for x in count_per_spectral]
        limits=[math.ceil(safety_factor*x) for x in predicted_vectors]
        taylor=sum(predicted_vectors)/3*block_seconds*overhead_taylor
        amplification=grid['variance_normalization_amplification_upper_bound']
        degrees=[required_degree(q['time']*spectral/2,min(.5,1e-9/amplification)) for q in a['propagation']]
        cheb=sum(degrees)*block_seconds*overhead_chebyshev
        prep=grid['states']*preparation_per_state+fixed_case_overhead
        cases.append(dict(faces=grid['faces'],spacing=grid['spacing'],shape=grid['shape'],states=grid['states'],
                          csr_nonzeros=grid['csr_nonzeros'],spectral_planning_bound=spectral,
                          spectral_scope='Measured largest-grid row-sum upper bound also bounds aligned reflecting subboxes.' if fine else 'Universal exit-rate/Gershgorin upper bound; no coarser-grid propagation is inferred exact.',
                          estimated_three_column_action_seconds=block_seconds,
                          predicted_taylor_vector_equivalents_per_query=predicted_vectors,
                          proposed_taylor_vector_equivalent_stop_per_query=limits,
                          chebyshev_degrees_using_polynomial_variance_bound=degrees,
                          estimated_taylor_seconds_before_margin=taylor,
                          estimated_chebyshev_seconds_before_margin=cheb,
                          estimated_preparation_and_serialization_seconds_before_margin=prep,
                          estimated_total_seconds_with_twofold_margin=safety_factor*(taylor+cheb+prep)))
    validation_reserve=60.
    estimated=sum(c['estimated_total_seconds_with_twofold_margin'] for c in cases)+validation_reserve+read(inputs[6])['wall_seconds']
    proposed_hours=math.ceil(estimated/3600)
    admission_memory=1.25*capouter['sampled_aggregate_peak_bytes']+256_000_000
    return dict(status='PROPOSED LATER ALLOCATION; COMPLETE CAMPAIGN NOT EXECUTED',
                source_inputs_sha256={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in inputs},
                formulas=dict(csr_seconds_per_nonzero=seconds_per_nnz,
                              taylor_observed_overhead_multiplier=overhead_taylor,
                              chebyshev_observed_overhead_multiplier=overhead_chebyshev,
                              count_per_spectral_bound_per_query=count_per_spectral,
                              preparation_seconds_per_state=preparation_per_state,
                              fixed_case_overhead_seconds=fixed_case_overhead,
                              engineering_safety_factor=safety_factor,
                              validation_and_receipt_reserve_seconds=validation_reserve),
                cases=cases, estimated_full_five_case_seconds_with_margin=estimated,
                proposed_allocation_wall_hours=proposed_hours,proposed_allocation_wall_seconds=proposed_hours*3600,
                proposed_aggregate_memory_cap_bytes=6_000_000_000,
                measured_capacity_aggregate_peak_bytes=capouter['sampled_aggregate_peak_bytes'],
                memory_admission_scenario_bytes=admission_memory,
                memory_admission_scenario_pass=admission_memory<=6_000_000_000,
                memory_scenario='1.25 times the measured touched-workspace aggregate peak plus256MB. This is engineering headroom, not a proven bound on unseen propagation allocations. A6GB watchdog remains mandatory.',
                uncertainty='Taylor counts are extrapolated from two smaller cost probes, not rigorously bounded by spectral radius; norm estimation, stopping behavior, cache effects and new equilibrium weights can change cost. The twofold factor and proposed action stops are declared planning safeguards, not probabilities or performance guarantees. Chebyshev degrees use a conservative polynomial variance bound and floating-point tail evaluation. Both complete propagators and normalized component comparisons are included.',
                launch_rule='Only a separately approved dated campaign protocol may launch these five cases. Use one worker, one-thread environment controls, the stated hard wall/memory caps and per-query Taylor action stops. Retain every failure. Do not launch any unlisted finer grid or geometry automatically. The existing two cost-probe K/C refinement failure remains evidence, and this executable schedule is not promised to pass nonlinear convergence.',
                physical_limit='The enlarged box follows the two-unit margin around all currently known congruent minima. It is not proof of complete whole-plane coverage or a native protein basin.')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,default=Path('recomputed-campaign-allocation.json'))
    args=parser.parse_args()
    result=run();args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['estimated_full_five_case_seconds_with_margin','proposed_allocation_wall_hours','memory_admission_scenario_bytes','memory_admission_scenario_pass']},indent=2))
