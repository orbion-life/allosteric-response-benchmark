#!/usr/bin/env python3
"""One cached-input replay with wall time/RSS; no scientific settings changed."""
from pathlib import Path
import argparse
import datetime
import hashlib
import json
import platform
import resource
import runpy
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parent


def rss():
    value=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(value if sys.platform=='darwin' else value*1024)


def child(script):
    start=time.perf_counter()
    runpy.run_path(str(ROOT/script),run_name='__main__')
    record={'script':script,'in_process_script_wall_seconds':time.perf_counter()-start,
            'peak_process_rss_bytes':rss(),'peak_method':'resource.getrusage(RUSAGE_SELF).ru_maxrss; macOS bytes or Linux KiB converted to bytes',
            'exit_status':0}
    (ROOT/f'cost-{script}.json').write_text(json.dumps(record,indent=2)+'\n')


def parent(reference, destination):
    date=datetime.datetime.now(datetime.timezone.utc).isoformat()
    phases=[]
    total_start=time.perf_counter()
    for script in ('prepare.py','run_pilot.py'):
        phase_start=time.perf_counter()
        with (ROOT/f'cost-{script}.log').open('w') as output:
            result=subprocess.run([sys.executable,str(ROOT/'measure_cost.py'),'--child',script],
                                  cwd=ROOT,stdout=output,stderr=subprocess.STDOUT)
        elapsed=time.perf_counter()-phase_start
        if result.returncode:raise RuntimeError(f'{script} failed; consult local phase log')
        phase=json.loads((ROOT/f'cost-{script}.json').read_text())
        phase['fresh_process_wall_seconds_including_startup_and_exit']=elapsed
        phases.append(phase)
        print(json.dumps(phase),flush=True)
    total_wall=time.perf_counter()-total_start
    # Verification is deliberately outside the measured workflow.
    import numpy as np
    import scipy
    verification=[]
    patterns=['model/static-model.npz','model/coordinate-d1.npz','model/coordinate-d2.npz',
              'model/all-mode-harmonic-scales.npz']
    patterns += ['results/'+p.name for p in sorted((ROOT/'results').glob('*.npz'))]
    for name in patterns:
        actual=np.load(ROOT/name);expected=np.load(reference/name)
        assert actual.files==expected.files,(name,'keys')
        fields={}
        for key in actual.files:
            a,b=actual[key],expected[key]
            assert a.shape==b.shape,(name,key,'shape')
            if a.dtype.kind in 'iubUS':
                assert np.array_equal(a,b),(name,key,'exact')
                fields[key]={'exact_equal':True,'max_absolute_difference':0.0}
            else:
                error=float(np.max(abs(a-b))) if a.size else 0.0
                assert np.allclose(a,b,atol=1e-10,rtol=1e-10,equal_nan=False),(name,key,error)
                fields[key]={'exact_equal':bool(np.array_equal(a,b)),'max_absolute_difference':error}
        verification.append({'file':name,'fields':fields,'passed':True})
    hessian=json.loads((ROOT/'cost-hessian.json').read_text())
    analytic=json.loads((ROOT/'results/analytic-checks.json').read_text())
    full=[x for x in analytic['timings'] if x['d']==492][0]
    protocol=json.loads((ROOT/'preanalysis-protocol.json').read_text())
    provenance=json.loads((ROOT/'instrumentation-provenance.json').read_text())
    if sys.platform=='darwin':
        cpu=subprocess.run(['/usr/sbin/sysctl','-n','machdep.cpu.brand_string'],capture_output=True,text=True).stdout.strip()
        memory=int(subprocess.run(['/usr/sbin/sysctl','-n','hw.memsize'],capture_output=True,text=True).stdout.strip())
    else:
        import os
        cpu=platform.processor() or platform.machine()
        memory=int(os.sysconf('SC_PAGE_SIZE')*os.sysconf('SC_PHYS_PAGES'))
    record={'status':'completed_once_and_verified','measured_utc':date,
      'workload':'cached-input prepare.py followed by run_pilot.py, each in a fresh process',
      'execution_count':{'prepare.py':1,'run_pilot.py':1},
      'machine':{'architecture':platform.machine(),'processor_model':cpu,'physical_memory_bytes':int(memory),
                 'os':platform.system(),'os_version':platform.mac_ver()[0] if sys.platform=='darwin' else platform.release()},
      'software':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__},
      'phases':phases,'total_sequential_wall_seconds':total_wall,
      'total_scope':'from first process launch through second process exit, including script imports, output serialization and phase-control overhead; cached input copying/install/verification excluded',
      'workflow_peak_process_rss_bytes':max(p['peak_process_rss_bytes'] for p in phases),
      'workflow_peak_scope':'maximum of sequential worker process high-water RSS, not a sum; includes interpreter/libraries and all arrays in each worker; controller memory excluded',
      'shared_preprocessing':{
        'prepare_phase_wall_seconds_including_startup':phases[0]['fresh_process_wall_seconds_including_startup_and_exit'],
        'prepare_phase_scope':'cached PDB parsing, sequence mapping, receiver definition, contact graph, Hessian/eigenbasis, d1/d2 observable polynomials/scales and saved model artifacts',
        'contact_hessian_eigenbasis_wall_seconds':hessian['wall_seconds'],
        'contact_hessian_eigenbasis_scope':hessian['scope'],
        'full492_harmonic_covariance_and_normalizer_wall_seconds':full['seconds'],
        'full492_harmonic_scope':'full equal-time/delayed covariance, raw response, harmonic SD extraction and save; a combined computation, not isolated SD-only time',
        'all8_harmonic_dimension_checks_wall_seconds':sum(x['seconds'] for x in analytic['timings']),
        'all8_harmonic_scope':'sum of per-dimension covariance/delayed-covariance/SD/save timers; excludes two earlier independent Wick checks and later normalization/ranking pass',
        'double_counting_warning':'These shared preparation and harmonic calculations are included in the total workflow and must not be added again to that total. They supply both classical and quantum model inputs.',
        'isolated_normalizer_only_time_measured':False},
      'included':['cached-input structural/model preparation','two independent harmonic moment checks','eight analytic harmonic dimensions','all24 finite-grid energy/coordinate/domain cases','built-in grid physical diagnostics','graph/geometric/degree baselines','model/response/ranking array serialization','prediction-freeze hashing'],
      'excluded':['package installation','copying cached inputs into fresh directory','network fetches','external Ohm baseline execution','reference-pocket evaluation and matched-resampling analysis','plots and PDF generation','independent post-run array verification','quantum compilation and simulation','hardware execution'],
      'verification':{'status':'passed','relative_tolerance':1e-10,'absolute_tolerance':1e-10,'integer_boolean_order_arrays':'exact equality required','compared_array_files':len(verification),'files':verification},
      'source_provenance':provenance,
      'protocol':{'input':protocol['input'],'reference':protocol['reference'],'grid_cases':protocol['grid_cases'],'energy_models':protocol['energies'],'analytic_dimensions':protocol['analytic_harmonic_dimensions']},
      'limitations':['one same-host observed run, not an average or performance guarantee','cached input and operating-system caches; cold disk cache not controlled','system load was not reserved','no per-subphase isolated peak RSS or pure-normalizer-only timing','small explicit timer instrumentation is included in observed preparation time'],
      'procedure_code_sha256':hashlib.sha256((ROOT/'measure_cost.py').read_bytes()).hexdigest()}
    destination.parent.mkdir(parents=True,exist_ok=True)
    payload=json.dumps(record,indent=2,allow_nan=False)+'\n'
    destination.write_text(payload)
    destination.with_name('full-pilot-cost-public.json').write_text(payload)
    print(json.dumps({'status':record['status'],'total_wall_seconds':total_wall,
                      'peak_worker_rss_bytes':record['workflow_peak_process_rss_bytes'],
                      'verified_array_files':len(verification),'shared_preprocessing':record['shared_preprocessing']},indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--child',choices=['prepare.py','run_pilot.py'])
    parser.add_argument('--reference',type=Path)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    if args.child:child(args.child)
    else:
        if (ROOT/'cost-prepare.py.json').exists() or (ROOT/'cost-run_pilot.py.json').exists():
            raise RuntimeError('Measured workflow has already been launched; this task permits one execution only')
        parent(args.reference.resolve(),args.output.resolve())
