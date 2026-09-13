"""Proposed later five-case run. Default is a non-executing admission preview."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import numpy as np
import psutil

ROOT=Path(__file__).resolve().parent


def main(protocol_path,output,execute):
    p=json.loads(protocol_path.read_text())
    for name,expected in {**p['executing_source_sha256'],**p['admission_evidence_sha256']}.items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=expected:
            raise ValueError(f'Continuation admission input changed: {name}')
    if not p['memory_admission_scenario_pass']:
        raise RuntimeError('Memory admission scenario fails.')
    if not execute:
        print(json.dumps(dict(status='PREVIEW ONLY; NO WORKERS LAUNCHED',protocol_id=p['protocol_id'],
                              cases=p['cases'],allocation=p['allocation'],warning=p['runtime_estimate_scope']),indent=2));return
    output.mkdir(parents=True,exist_ok=False)
    env=dict(os.environ)
    for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):env[key]='1'
    start=time.perf_counter();peak=0;records=[];arrays=[]
    for i,case in enumerate(p['cases']):
        case_path=output/f'case-{i+1}'
        stage_start=time.perf_counter();stop=None
        with (output/f'case-{i+1}.log').open('w') as log:
            child=subprocess.Popen([sys.executable,'-W','error',str(ROOT/'continuation_worker.py'),'--protocol',str(protocol_path),
                                    '--case-index',str(i),'--output',str(case_path)],env=env,stdout=log,stderr=subprocess.STDOUT)
            while child.poll() is None:
                try:
                    proc=psutil.Process(child.pid)
                    rss=sum(x.memory_info().rss for x in [proc,*proc.children(recursive=True)] if x.is_running())
                except (psutil.NoSuchProcess,psutil.AccessDenied):rss=0
                aggregate=rss+psutil.Process().memory_info().rss;peak=max(peak,aggregate)
                elapsed=time.perf_counter()-start;stage_elapsed=time.perf_counter()-stage_start
                if elapsed>p['allocation']['wall_seconds'] or stage_elapsed>case['stage_wall_seconds'] or aggregate>p['allocation']['aggregate_memory_bytes']:
                    stop='global_wall' if elapsed>p['allocation']['wall_seconds'] else ('stage_wall' if stage_elapsed>case['stage_wall_seconds'] else 'memory')
                    child.kill();break
                time.sleep(.05)
            code=child.wait()
        records.append(dict(case_index=i,return_code=code,stop=stop,outer_case_seconds=time.perf_counter()-stage_start))
        receipt=dict(status='RUNNING',protocol_id=p['protocol_id'],protocol_sha256=hashlib.sha256(protocol_path.read_bytes()).hexdigest(),
                     cases=records,outer_wall_seconds=time.perf_counter()-start,sampled_aggregate_peak_bytes=peak)
        (output/'run.json').write_text(json.dumps(receipt,indent=2)+'\n')
        if code or stop:
            receipt['status']='RESOURCE_STOP' if stop else 'WORKER_FAILED'
            (output/'run.json').write_text(json.dumps(receipt,indent=2)+'\n');raise SystemExit(1)
        with np.load(case_path/'response.npz') as z:
            arrays.append({key:z[key] for key in ('G0','K','C')})
    comparisons=[]
    for item in p['numerical_screen_pairs']:
        i,j=item['first'],item['second']
        errors={key:float(np.max(np.abs(arrays[j][key]-arrays[i][key]))) for key in ('G0','K','C')}
        comparisons.append(dict(**item,errors=errors,pass_all_components=all(v<=p['numerical_screen_tolerance'] for v in errors.values())))
    coarse,fine=comparisons[0],comparisons[1]
    trend={key:fine['errors'][key]<=coarse['errors'][key]+p['three_level_trend_arithmetic_allowance'] for key in ('G0','K','C')}
    receipt.update(status='COMPLETE',comparisons=comparisons,three_level_nonincreasing_error_trend=trend,
                   numerical_screens_pass=all(c['pass_all_components'] for c in comparisons if c['controlling']) and all(trend.values()),
                   whole_plane_coverage_status='NOT ESTABLISHED BY THIS FINITE-BOX CAMPAIGN',
                   outer_wall_seconds=time.perf_counter()-start,sampled_aggregate_peak_bytes=peak)
    (output/'run.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--protocol',type=Path,default=ROOT/'continuation-protocol.json')
    parser.add_argument('--output',type=Path,default=Path('later-campaign'))
    parser.add_argument('--execute',action='store_true',help='Launch the separately proposed later campaign; otherwise print admission preview only.')
    args=parser.parse_args();main(args.protocol.resolve(),args.output.resolve(),args.execute)
