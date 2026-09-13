"""Serial fresh-worker execution with aggregate RSS and wall-clock watchdog."""
import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time
import psutil
from model import ROOT, protocol


def main(stage, output):
    output.mkdir(parents=True, exist_ok=False)
    p = protocol()
    env = dict(os.environ)
    for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS'):
        env[key] = '1'
    if stage == 'calibration':
        jobs = [('calibration', ['calibration.py', '--output', str(output/'calibration')])]
        cap = p['execution']['calibration_wall_seconds']
        memory_cap = 512_000_000
    else:
        jobs = [(f'nonlinear-h{h:g}', ['nonlinear_profile.py', '--output', str(output/f'nonlinear-h{h:g}'), '--spacing', str(h)]) for h in (.5, .25)]
        cap = p['execution']['nonlinear_profile_wall_seconds']
        memory_cap = p['execution']['aggregate_memory_bytes']
    started = time.perf_counter()
    receipt = dict(status='RUNNING', started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                   protocol_sha256=hashlib.sha256((ROOT/'protocol.json').read_bytes()).hexdigest(),
                   platform=dict(system=platform.system(), release=platform.release(), machine=platform.machine(), python=platform.python_version()),
                   stage=stage, workers=1, BLAS_threads=1, cap_seconds=cap, memory_cap_bytes=memory_cap, jobs=[])
    peak = 0
    for name, args in jobs:
        with (output/(name+'.log')).open('w') as log:
            child = subprocess.Popen([sys.executable, '-W', 'error', *args], cwd=ROOT, env=env, stdout=log, stderr=subprocess.STDOUT)
            sampled = 0
            killed = None
            while child.poll() is None:
                try:
                    process = psutil.Process(child.pid)
                    rss = sum(x.memory_info().rss for x in [process, *process.children(recursive=True)] if x.is_running())
                except (psutil.NoSuchProcess, psutil.AccessDenied):
                    rss = 0
                aggregate = rss+psutil.Process().memory_info().rss
                sampled = max(sampled, rss)
                peak = max(peak, aggregate)
                if time.perf_counter()-started > cap or aggregate > memory_cap:
                    killed = 'time' if time.perf_counter()-started > cap else 'memory'
                    child.kill()
                    break
                time.sleep(.02)
            code = child.wait()
        receipt['jobs'].append(dict(name=name, return_code=code, sampled_worker_peak_bytes=sampled, watchdog_stop=killed))
        receipt.update(wall_seconds=time.perf_counter()-started, sampled_aggregate_peak_bytes=peak)
        (output/'watchdog.json').write_text(json.dumps(receipt, indent=2)+'\n')
        if code or killed:
            receipt['status'] = 'STOPPED' if killed else 'FAILED'
            (output/'watchdog.json').write_text(json.dumps(receipt, indent=2)+'\n')
            raise SystemExit(1)
    receipt.update(status='COMPLETE', wall_seconds=time.perf_counter()-started, sampled_aggregate_peak_bytes=peak)
    (output/'watchdog.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--stage', choices=['calibration', 'profiles'], required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    main(args.stage, args.output.resolve())
