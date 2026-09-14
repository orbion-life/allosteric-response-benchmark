"""Bounded local FP64 implementation check and first-time reference pilot."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[key]='1'
from pathlib import Path
import sys, json, time, subprocess, datetime, hashlib, platform
ROOT=Path(__file__).resolve().parent

def worker(kind):
    sys.path.insert(0,str(ROOT/'source/dynamics'))
    if kind=='small':
        from check_small import check
        return check(ROOT/'local/small','numpy')
    from core import run_reference
    return run_reference(ROOT/'local/pilot',[8,7,19],.125,'numpy',times=[.0075],independent=True)

def main():
    import psutil
    import numpy, scipy
    records=[]
    for kind in ['small','pilot']:
        start=time.perf_counter(); peak=0; status='RUNNING'
        with (ROOT/'local'/f'{kind}.log').open('w') as log:
            proc=subprocess.Popen([sys.executable,str(Path(__file__).resolve()),'--worker',kind],stdout=log,stderr=subprocess.STDOUT)
            ps=psutil.Process(proc.pid)
            while proc.poll() is None:
                try: peak=max(peak,ps.memory_info().rss)
                except psutil.NoSuchProcess: pass
                if peak>8_000_000_000 or time.perf_counter()-start>300:
                    status='STOPPED_MEMORY_OR_TIME_LIMIT';proc.terminate()
                    try:proc.wait(timeout=5)
                    except subprocess.TimeoutExpired:proc.kill();proc.wait()
                    break
                time.sleep(.1)
        status=('COMPLETE' if proc.returncode==0 else 'FAILED') if status=='RUNNING' else status
        row={'kind':kind,'status':status,'returncode':proc.returncode,'outer_seconds':time.perf_counter()-start,'observed_peak_RSS_bytes':peak}
        records.append(row);print(json.dumps(row),flush=True)
        (ROOT/'local/supervisor.json').write_text(json.dumps({'records':records,'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'environment':{'python':sys.version,'numpy':numpy.__version__,'scipy':scipy.__version__,'platform':platform.platform()},'protocol_sha256':hashlib.sha256((ROOT/'protocol.json').read_bytes()).hexdigest()},indent=2)+'\n')
        if status!='COMPLETE':raise RuntimeError('Local gate failed; retained evidence')

if __name__=='__main__':
    if len(sys.argv)>1:worker(sys.argv[2])
    else:main()
