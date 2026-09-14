"""Bound each local worker and compilation, retaining all failures."""
import argparse,datetime,json,os,subprocess,sys,time
from pathlib import Path
import psutil
ROOT=Path(__file__).resolve().parent
p=argparse.ArgumentParser();p.add_argument('mode');a=p.parse_args();dest=ROOT/'results/watchdogs';dest.mkdir(exist_ok=True)
env=os.environ.copy()
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:env[key]='1'
progress=dest/(a.mode+'-progress.json');env['PULSAR_PROGRESS']=str(progress)
start=time.monotonic();began=datetime.datetime.now(datetime.timezone.utc).isoformat();peak=0;reason=None
with (dest/(a.mode+'.log')).open('w') as log:
    proc=subprocess.Popen([sys.executable,str(ROOT/'candidate_experiment.py'),a.mode],stdout=log,stderr=subprocess.STDOUT,env=env)
    while proc.poll() is None:
        children=[]
        try:
            parent=psutil.Process(proc.pid);children=parent.children(recursive=True)
            rss=sum(x.memory_info().rss for x in [parent]+children if x.is_running());peak=max(peak,rss)
        except psutil.NoSuchProcess:rss=0
        if rss>12*1024**3:reason='12GiB worker memory cap'
        if time.monotonic()-start>3600:reason='one-hour worker cap'
        try:
            stage=json.loads(progress.read_text())
            if any(x in stage['stage'] for x in ['compile','synthesis','operator_check']) and time.monotonic()-stage['monotonic']>600:reason='600-second compilation or operator-check cap'
        except (OSError,ValueError,KeyError):pass
        if reason:
            for child in children:
                try:child.kill()
                except psutil.NoSuchProcess:pass
            proc.kill();break
        time.sleep(.2)
    code=proc.wait()
receipt=dict(mode=a.mode,started_utc=began,completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),exit_code=code,termination_reason=reason,wall_seconds=time.monotonic()-start,sampled_peak_aggregate_RSS_bytes=peak)
(dest/(a.mode+'.json')).write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2));print((dest/(a.mode+'.log')).read_text()[-3000:]);sys.exit(code)
