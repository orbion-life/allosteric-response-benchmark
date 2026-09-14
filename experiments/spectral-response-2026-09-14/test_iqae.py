"""Independent finite-binomial confidence checks and fixed-amplitude trials."""
import json,math,time
from pathlib import Path
import numpy as np
from scipy.stats import binom
from iqae import cp_interval,simulated_estimates,choose_power
ROOT=Path(__file__).resolve().parent


def main():
    start=time.perf_counter();out=ROOT/'results';out.mkdir(exist_ok=True);checks=[]
    for delta in [1e-3,1e-6]:
        counts=np.arange(257);lo,hi=cp_interval(counts,256,delta)
        grid=np.linspace(0,1,1001);mass=binom.pmf(counts[:,None],256,grid[None,:])
        coverage=np.sum(mass*((grid[None,:]>=lo[:,None])&(grid[None,:]<=hi[:,None])),axis=0)
        assert coverage.min()>=1-delta-1e-12
        checks.append(dict(delta=delta,minimum_exact_binomial_coverage=float(coverage.min()),nominal=1-delta,grid_points=len(grid),shots=256))
    amps=np.repeat([0,.001,.01,.1,.25,.5,.75,.9,.99,.999,1.],32)
    result=simulated_estimates(amps,.001,.01,2026091403)
    covered=(amps>=result['a_lower'])&(amps<=result['a_upper']);completed=result['status']==1
    # Coverage outcomes are reported, not required to be perfect by a 99% procedure.
    assert completed.all(),np.unique(result['status'],return_counts=True)
    assert np.max(result['a_upper']-result['a_lower'])<=.002+1e-15
    payload={k:v for k,v in result.items() if k!='rounds'};payload['true_amplitudes']=amps
    np.savez_compressed(out/'IQAE-fixed-amplitude-tests.npz',**payload)
    for r in result['rounds']:
        low=r['theta_lower'];high=r['theta_upper'];m=r['multipliers'];branch=r['branches']
        assert np.all(low*m>=branch*np.pi/2-1e-12)
        assert np.all(high*m<=(branch+1)*np.pi/2+1e-12)
    receipt=dict(status='PASS_CONFIDENCE_IMPLEMENTATION_CHECKS',exact_CP_checks=checks,fixed_amplitude_trials=len(amps),
                 nominal_individual_confidence=.99,covered=int(covered.sum()),completed=int(completed.sum()),
                 maximum_halfwidth=float(np.max((result['a_upper']-result['a_lower'])/2)),
                 scope='Ideal classical binomial simulations. Finite empirical coverage is a diagnostic, not a replacement for the conditional CP union-bound proof.',seconds=time.perf_counter()-start)
    (out/'IQAE-test-receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))


if __name__=='__main__':main()
