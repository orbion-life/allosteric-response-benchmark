#!/usr/bin/env python3
"""Enforce wall-time and sampled resident-memory bounds on a calculation."""
from pathlib import Path
import json,subprocess,sys,time
p=Path(__file__).resolve().parent
start=time.monotonic(); proc=subprocess.Popen([sys.executable,str(p/'run.py'),*sys.argv[1:]])
peak=0
while proc.poll() is None:
    try:
        rss=int(subprocess.check_output(['ps','-o','rss=','-p',str(proc.pid)],text=True).strip())*1024
        peak=max(peak,rss)
    except (ValueError,subprocess.CalledProcessError): pass
    if time.monotonic()-start>900 or peak>8000000000:
        proc.kill();proc.wait();raise RuntimeError('Preanalysis resource cap exceeded; partial outputs retained.')
    time.sleep(.5)
print(json.dumps({'guard_wall_seconds':time.monotonic()-start,'sampled_peak_RSS_bytes':peak,'returncode':proc.returncode}))
sys.exit(proc.returncode)
