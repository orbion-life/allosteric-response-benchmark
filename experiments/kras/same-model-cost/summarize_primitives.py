"""Join measured primary/variant preparation counts without mixing their shots."""
import json
from pathlib import Path


def main():
    p=Path(__file__).parent/'results'
    r=json.loads((p/'results.json').read_text())
    original=json.loads((p/'primitive-summary.json').read_text())['records']
    variant=json.loads((p/'raw-primitive-summary.json').read_text())['records']
    coefficient={x['index']:x for x in original if x['kind']=='coefficient'}
    families=[('weighted_monomial_SVD',[x for x in original if x['kind']=='observable']),
              ('separately_normalized_monomials',variant)]
    records=[]
    for name,observables in families:
        assert len(observables)==12 and all(x['status']=='completed' for x in observables)
        oracle=next(x for x in r['exact_readout_basis_comparison'] if x['family']==name)
        for K in [25,42]:
            coeff=coefficient[K]['counts']['two_qubit_cx']
            totals=[a['counts']['two_qubit_cx']+b['counts']['two_qubit_cx']+2*coeff
                    for i,a in enumerate(observables) for b in observables[i:]]
            record=dict(family=name,degree=K,
                minimum_modular_boundary_CX_per_overlap=min(totals),
                maximum_modular_boundary_CX_per_overlap=max(totals),
                max_phase_sensitive_PREP_state_error=max(x['phase_sensitive_statevector_max_error'] for x in observables),
                max_coherent_inverse_state_error=max(x['coherent_inverse_max_error'] for x in observables),
                scope='Actual separately compiled PREP/COEF and inverse modules only; SELECT/readout/routing/global optimization excluded; not a whole-circuit count or lower bound')
            if K==42:
                shots=oracle['oracle_sufficient_total_shots']
                record.update(oracle_sufficient_shots=shots,
                    shots_times_minimum_modular_boundary_CX=shots*min(totals),
                    shots_times_maximum_modular_boundary_CX=shots*max(totals))
            records.append(record)
    result=dict(records=records,status='Matched measured primitive counts and their own basis-specific planning budgets',
                scope='No primary complete circuits or quantum shots executed')
    (p/'modular-boundary-cost.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__': main()
