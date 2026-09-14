"""Full-output ideal IQAE likelihood study with exact complete cost composition."""
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[key]='1'
import argparse,hashlib,json,math,resource,time
from pathlib import Path
import numpy as np
from iqae import simulated_estimates
ROOT=Path(__file__).resolve().parent


def run(name):
    start=time.perf_counter();out=ROOT/'results'/name
    native=json.loads((out/'native-receipt.json').read_text());ft=json.loads((out/'fault-tolerant/receipt.json').read_text());model={k:v for k,v in np.load(out/'model.npz').items()}
    norm=model['norm'];N=len(norm);ij=np.array([(ti,i,j) for ti in range(3) for i in range(N) for j in range(i,N)],dtype=int);number=len(ij)
    factors=np.array([norm[i]*norm[j] for _,i,j in ij]);kernel=np.array([model['K'][ti,i,j] for ti,i,j in ij]);a=(1+kernel/factors)/2
    assert np.all((a>=-1e-12)&(a<=1+1e-12));a=np.clip(a,0,1)
    delta=.05/number;epsilon=native['amplitude_epsilon'];seed=2026091401 if name=='triangle' else 2026091402
    results=simulated_estimates(a,epsilon,delta,seed)
    rows={f"{r['kind']}_{r['index']}":r for r in ft['components']}
    A_T=np.array([rows[f'prepare_{i}']['T']+rows[f'prepare_{j}']['T']+rows[f'filter_{ti}']['T'] for ti,i,j in ij],dtype=np.int64)
    A_CX=np.array([native['preparations'][i]['CX']+native['preparations'][j]['CX']+native['filters'][ti]['CX'] for ti,i,j in ij],dtype=np.int64)
    A_ops=np.array([rows[f'prepare_{i}']['operations']+rows[f'prepare_{j}']['operations']+rows[f'filter_{ti}']['operations']+2 for ti,i,j in ij],dtype=np.int64)
    A_depth=np.array([native['preparations'][i]['depth']+native['preparations'][j]['depth']+native['filters'][ti]['depth']+2 for ti,i,j in ij],dtype=np.int64)
    ref=native['zero_reflection'];refT=ref['T'];refCX=ref['CX'];refops=ref['operations']+3
    total_T=sum(int(q)*int(t)+int(g)*refT for q,t,g in zip(results['A_queries'],A_T,results['Q_queries']))
    total_CX=sum(int(q)*int(t)+int(g)*refCX for q,t,g in zip(results['A_queries'],A_CX,results['Q_queries']))
    native_depth=A_depth+results['maximum_k']*(2*A_depth+ref['depth']+3)
    coherent_T=A_T+results['maximum_k']*(2*A_T+refT)
    serial_FT_ops=A_ops+results['maximum_k']*(2*A_ops+refops)
    shots=math.ceil(math.log(2*number/.05)/(2*epsilon**2))
    sampling_T=shots*sum(map(int,A_T));sampling_CX=shots*sum(map(int,A_CX))
    cover=(a>=results['a_lower'])&(a<=results['a_upper']);complete=results['status']==1
    G0=model['B'].T@model['B'];static=np.array([G0[i,j] for _,i,j in ij])
    C_est=factors*(results['a_lower']+results['a_upper']-1)-static
    C_exact=np.array([model['C'][ti,i,j] for ti,i,j in ij]);stat_radius=factors*(results['a_upper']-results['a_lower'])
    raw={k:v for k,v in results.items() if k!='rounds'}
    raw.update(queries=ij,true_ideal_a=a,normalization=factors,A_T=A_T,A_CX=A_CX,A_operations=A_ops,A_depth_upper=A_depth,
               C_estimate=C_est,C_exact=C_exact,statistical_response_halfwidth=stat_radius)
    for ri,r in enumerate(results['rounds']):
        for key,value in r.items():raw[f'round_{ri+1}_{key}']=value
    np.savez_compressed(out/'IQAE-all-outputs.npz',**raw)
    receipt=dict(status='COMPLETE_IDEAL_LIKELIHOOD_AND_COST_STUDY' if complete.all() else 'INCOMPLETE_OR_INTERVAL_CONTRADICTION',
                 model=name,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),number_outputs=number,seed=seed,
                 confidence_family=.95,per_output_failure=delta,per_round_failure='delta_output/[r(r+1)], fresh batches',
                 amplitude_epsilon=epsilon,completed=int(complete.sum()),interval_contradictions=int(np.sum(results['status']==-1)),
                 true_ideal_amplitudes_covered=int(cover.sum()),rounds=len(results['rounds']),
                 maximum_response_point_error=float(np.max(abs(C_est-C_exact))),maximum_statistical_response_halfwidth=float(np.max(stat_radius)),
                 inherited_projection_error=native['inherited_projection_error'],
                 complete_response_error_allowance=native['inherited_projection_error']+1e-8+ft['word_rounding_response_bound']+ft['synthesis_response_bound']+float(np.max(stat_radius)),
                 total_readouts=int(results['readouts'].sum()),total_A_or_inverse_calls=int(results['A_queries'].sum()),total_Grover_calls=int(results['Q_queries'].sum()),
                 maximum_Grover_power=int(results['maximum_k'].max()),maximum_A_or_inverse_calls_one_shot=int(2*results['maximum_k'].max()+1),
                 estimated_complete_T=total_T,estimated_complete_CX=total_CX,
                 maximum_coherent_T_one_shot=int(coherent_T.max()),maximum_coherent_CX_one_shot=int(np.max(A_CX+results['maximum_k']*(2*A_CX+refCX))),maximum_native_depth_upper_one_shot=int(native_depth.max()),maximum_serial_FT_operations_upper_one_shot=int(serial_FT_ops.max()),
                 matched_spectral_sampling=dict(shots_per_output=shots,total_readouts=shots*number,complete_T=sampling_T,complete_CX=sampling_CX),
                 IQAE_to_same_spectral_sampling_T_ratio=total_T/sampling_T,IQAE_to_same_spectral_sampling_CX_ratio=total_CX/sampling_CX,
                 per_logical_T_failure_sufficient_union_bound_for_1pct_campaign=.01/total_T,
                 scope='Adaptive schedules and successes are classical ideal binomial simulations from the exact spectral model probabilities, with representative actual native circuit verification. True amplitude is not supplied to scheduler. Gate totals apply actual constructive component counts to these realized ideal schedules; not a guaranteed hardware runtime bound or QPU experiment. Approximate compiled A can slightly change realized schedules. Statistical guarantee is conditional on exact repeated A and reflections, with a separately bounded fixed encoding bias. No routing, error correction or physical gate-time estimate.',
                 exact_classical=native['exact_classical'],wall_seconds=time.perf_counter()-start,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    (out/'IQAE-receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('model',choices=['triangle','kras']);run(p.parse_args().model)
