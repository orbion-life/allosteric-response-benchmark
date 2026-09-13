"""Independent watchdog for the separately authorized construction-only probe."""
import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import psutil

ROOT=Path(__file__).resolve().parent


def main(output):
    output.mkdir(parents=True,exist_ok=False)
    p=json.loads((ROOT/'capacity-protocol.json').read_text())
    env=dict(os.environ)
    for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):env[k]='1'
    started=time.perf_counter();peak=0;stop=None
    with (output/'worker.log').open('w') as log:
        child=subprocess.Popen([sys.executable,'-W','error',str(ROOT/'capacity_probe.py'),'--output',str(output/'case')],env=env,stdout=log,stderr=subprocess.STDOUT)
        while child.poll() is None:
            try:
                process=psutil.Process(child.pid)
                rss=sum(x.memory_info().rss for x in [process,*process.children(recursive=True)] if x.is_running())
            except (psutil.NoSuchProcess,psutil.AccessDenied):rss=0
            aggregate=rss+psutil.Process().memory_info().rss;peak=max(peak,aggregate)
            if time.perf_counter()-started>p['allowance']['wall_seconds'] or aggregate>p['allowance']['aggregate_memory_bytes']:
                stop='time' if time.perf_counter()-started>p['allowance']['wall_seconds'] else 'memory';child.kill();break
            time.sleep(.02)
        code=child.wait()
    receipt=dict(status='COMPLETE' if code==0 and stop is None else 'FAILED',return_code=code,watchdog_stop=stop,
                 capacity_protocol_sha256=hashlib.sha256((ROOT/'capacity-protocol.json').read_bytes()).hexdigest(),
                 outer_wall_seconds=time.perf_counter()-started,sampled_aggregate_peak_bytes=peak,
                 cap_seconds=p['allowance']['wall_seconds'],cap_bytes=p['allowance']['aggregate_memory_bytes'],
                 scope='Outer fresh worker invocation, including imports and receipt serialization; dependency installation excluded.')
    (output/'watchdog.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
    if code or stop:raise SystemExit(1)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();main(args.output.resolve())
