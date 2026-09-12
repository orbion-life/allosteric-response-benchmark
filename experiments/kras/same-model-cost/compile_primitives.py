"""Compile actual primary-model loading gates under predeclared local limits.

Each worker is stopped at 20 seconds or 1 GiB process-tree RSS. Full estimator
construction is deliberately not attempted: the known statevector exceeds
512 MiB even before any compiler workspace. No timings are extrapolated.
"""
import json
import argparse
import subprocess
import sys
import time
from pathlib import Path
import psutil
from analyze import CAPS


def main(output,raw_only=False):
    root=Path(__file__).parent
    output.mkdir(parents=True,exist_ok=True)
    jobs=[('observable',k) for k in range(12)]+[('coefficient',25),('coefficient',42),
          ('transition_row',544),('transition_row',0),('transition_row',1089)]
    if raw_only:
        jobs=[('raw_observable',k) for k in range(12)]
    receipts=[]
    for kind,index in jobs:
        destination=output/f'{kind}-{index}.json'
        if destination.exists():
            raise RuntimeError('Refusing to overwrite a compilation receipt; move results/primitives before a new run')
        started=time.perf_counter()
        with (output/f'{kind}-{index}.log').open('w') as log:
            p=subprocess.Popen([sys.executable,str(root/'compile_worker.py'),kind,str(index),
                                '--output',str(destination)],stdout=log,stderr=log)
            peak=0
            reason=None
            while p.poll() is None:
                try:
                    process=psutil.Process(p.pid)
                    total=process.memory_info().rss
                    for child in process.children(recursive=True):
                        try: total+=child.memory_info().rss
                        except psutil.Error: pass
                    peak=max(peak,total)
                except psutil.Error:
                    pass
                if time.perf_counter()-started > CAPS['primitive_wall_seconds']:
                    reason='Stopped at declared 20-second worker wall limit'
                elif peak > CAPS['primitive_rss_bytes']:
                    reason='Stopped at declared 1-GiB sampled process-tree RSS limit'
                if reason:
                    try:
                        process=psutil.Process(p.pid)
                        for child in process.children(recursive=True): child.kill()
                    except psutil.Error: pass
                    p.kill()
                    break
                time.sleep(.025)
            p.wait()
        record=json.loads(destination.read_text()) if destination.exists() else dict(kind=kind,index=index)
        record.update(controller_wall_seconds=time.perf_counter()-started,
                      controller_sampled_peak_tree_RSS_bytes=peak,worker_exit_code=p.returncode)
        if reason:
            record.update(status='stopped',reason=reason)
        elif p.returncode != 0:
            record.update(status='failed',reason='Worker failed; see sanitized log')
        elif record['counts']['total_native_operations'] > CAPS['primitive_native_operations']:
            record.update(status='compiled_but_operation_cap_exceeded',
                          reason='Actual compilation exceeded the declared 200,000 native-operation acceptance cap')
        destination.write_text(json.dumps(record,indent=2)+'\n')
        # Logs are diagnostic only; strip machine-local paths before public staging.
        log_path=output/f'{kind}-{index}.log'
        text=log_path.read_text().replace(str(root),'<experiment>').replace(str(Path(sys.prefix)),'<python-environment>')
        log_path.write_text(text)
        receipts.append(record)
        print(json.dumps({k:record[k] for k in ['kind','index','status','controller_wall_seconds']}),flush=True)
    summary=dict(scope='Actual isolated loading primitives, not a compiled primary overlap circuit',
                 caps=CAPS,records=receipts,
                 inference='No full-circuit CX/depth estimate is inferred from these isolated primitives.')
    summary_name='raw-primitive-summary.json' if raw_only else 'primitive-summary.json'
    (output.parent/summary_name).write_text(json.dumps(summary,indent=2)+'\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=Path(__file__).parent/'results/primitives')
    parser.add_argument('--raw-only',action='store_true',help='Compile normalized monomial readout variant only; no transition retries')
    args=parser.parse_args()
    main(args.output,args.raw_only)
