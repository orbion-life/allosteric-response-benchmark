"""Enforce the recorded wall-clock and aggregate process-memory ceilings."""
from pathlib import Path
import json,os,signal,subprocess,sys,time
ROOT=Path(__file__).resolve().parent
protocol=json.loads((ROOT/'preanalysis.json').read_text());limits=protocol['resources']
output=Path(sys.argv[1]).resolve();start=time.monotonic();peak=0;reason=None
p=subprocess.Popen([sys.executable,str(ROOT/'operator_study.py'),'--output',str(output)],start_new_session=True)
try:
    while p.poll() is None:
        current=0
        data=subprocess.check_output(['ps','-axo','pid=,pgid=,rss='],text=True)
        for line in data.splitlines():
            pid,group,memory=map(int,line.split())
            if group==p.pid:current+=memory*1024
        peak=max(peak,current)
        if time.monotonic()-start>limits['max_total_wall_seconds']:
            reason='wall-time limit';break
        if current>limits['max_aggregate_bytes']:
            reason='aggregate-memory limit';break
        time.sleep(.25)
except BaseException as exc:
    reason='monitoring failure: '+type(exc).__name__+': '+str(exc)
finally:
    if reason is not None and p.poll() is None:
        os.killpg(p.pid,signal.SIGTERM)
        try:p.wait(timeout=5)
        except subprocess.TimeoutExpired:os.killpg(p.pid,signal.SIGKILL)
    code=p.wait()
    receipt={'status':'COMPLETE' if code==0 and reason is None else 'STOPPED','returncode':code,'limit_reason':reason,'wall_seconds':time.monotonic()-start,'sampled_peak_aggregate_bytes':peak,'sampling_seconds':.25,'limits':limits}
    output.mkdir(parents=True,exist_ok=True);(output/'watchdog.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt),flush=True)
sys.exit(0 if code==0 and reason is None else 1)
