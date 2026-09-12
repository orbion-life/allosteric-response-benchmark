"""Analyse every locked case; accuracy failures are data, not test failures."""
from pathlib import Path
import argparse,csv,hashlib,json
import numpy as np
ROOT=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def maximum(a):return float(np.max(np.abs(a)))
def run(results,output):
    protocol=json.loads((ROOT/'preanalysis.json').read_text())
    receipt=json.loads((results/'run.json').read_text())
    assert receipt['status']=='COMPLETE' and len(receipt['cases'])==72
    output.mkdir(parents=True,exist_ok=True)
    def get(k,d,n,e,kind):return dict(np.load(results/f'{kind}-k{k}-d{d}-n{n}-e{e}.npz'))
    rows=[];checks=[]
    for k in protocol['kappa']:
        hf=dict(np.load(results/f'analytic-k{k}-d3.npz'))
        for n in protocol['grid_n']:
            full=get(k,3,n,4,'biquadratic');hfg=get(k,3,n,4,'harmonic')
            for d in protocol['retained_dimensions']:
                nl=get(k,d,n,4,'biquadratic');hg=get(k,d,n,4,'harmonic')
                completed=hf['C']+nl['C']-hg['C']
                samegrid=hfg['C']+nl['C']-hg['C']
                # Decompose the correction into its static and delayed parts.
                eq=hf['equilibrium']+nl['equilibrium']-hg['equilibrium']
                delayed=hf['delayed']+nl['delayed']-hg['delayed']
                assert maximum(completed-(delayed-eq))<1e-10
                assert maximum((samegrid-full['C'])-((nl['C']-hg['C'])-(full['C']-hfg['C'])))<1e-10
                if d==3:assert maximum(samegrid-full['C'])<1e-10
                for i,tau in enumerate(protocol['times_over_tau']):
                    rows.append(dict(kappa=k,d=d,n=n,time_over_tau=tau,
                                     restriction_error=maximum(nl['C'][i]-full['C'][i]),
                                     analytic_completion_error=maximum(completed[i]-full['C'][i]),
                                     same_grid_completion_error=maximum(samegrid[i]-full['C'][i]),
                                     full_harmonic_grid_vs_analytic_error=maximum(hfg['C'][i]-hf['C'][i]),
                                     completion_static_error=maximum(eq-full['equilibrium']),
                                     completion_delayed_error=maximum(delayed[i]-full['delayed'][i])))
        full33=get(k,3,33,4,'biquadratic');full65=get(k,3,65,4,'biquadratic');full81=get(k,3,81,5,'biquadratic')
        for d in protocol['retained_dimensions']:
            a=get(k,d,33,4,'biquadratic');b=get(k,d,65,4,'biquadratic');c=get(k,d,81,5,'biquadratic')
            ha=get(k,d,33,4,'harmonic');hb=get(k,d,65,4,'harmonic');hc=get(k,d,81,5,'harmonic')
            for i,tau in enumerate(protocol['times_over_tau']):
                error=maximum(hf['C'][i]+b['C'][i]-hb['C'][i]-full65['C'][i])
                grid_delta=max(maximum(full65['C'][i]-full33['C'][i]),maximum((b['C']-hb['C'])[i]-(a['C']-ha['C'])[i]))
                domain_delta=max(maximum(full81['C'][i]-full65['C'][i]),maximum((c['C']-hc['C'])[i]-(b['C']-hb['C'])[i]))
                checks.append(dict(kappa=k,d=d,time_over_tau=tau,completion_error_n65=error,
                                   grid_max_difference=grid_delta,domain_max_difference=domain_delta,
                                   representation_pass=error<=.002,grid_pass=grid_delta<=.001,domain_pass=domain_delta<=.001,
                                   joint_pass=error<=.002 and grid_delta<=.001 and domain_delta<=.001))
    summary=dict(protocol_sha256=sha(ROOT/'preanalysis.json'),script_sha256=sha(__file__),run_receipt_sha256=sha(results/'run.json'),
                 finite_grid_rows=rows,acceptance_checks=checks,
                 overall_dimension_acceptance={str(d):all(r['joint_pass'] for r in checks if r['d']==d) for d in [1,2,3]},
                 status='Every fixed case retained; engineering acceptance is distinct from algebra/replay verification.',
                 limits=['A fully retained three-coordinate finite grid is not a certified continuum reference.',
                         'A refinement difference is an empirical diagnostic, not a rigorous continuum error bound.',
                         'Exact harmonic restoration and d=3 same-grid restoration are algebraic identities, not biological evidence.',
                         'Stiffness/time regimes are synthetic and were not selected after observing performance.'])
    (output/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    for name,data in [('all-errors.csv',rows),('acceptance.csv',checks)]:
        with (output/name).open('w',newline='') as stream:
            w=csv.DictWriter(stream,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
    print(json.dumps({'acceptance':summary['overall_dimension_acceptance'],'d2_n65':[r for r in rows if r['d']==2 and r['n']==65]},indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--results',type=Path,default=ROOT/'results');p.add_argument('--output',type=Path,default=ROOT/'analysis');a=p.parse_args();run(a.results,a.output)
