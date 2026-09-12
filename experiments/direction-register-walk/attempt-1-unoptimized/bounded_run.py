#!/usr/bin/env python3
"""Stop the local canary if its frozen wall/RSS bound is exceeded."""
from pathlib import Path
import subprocess,sys,time,json
root=Path(__file__).resolve().parent
start=time.monotonic();peak=0;reason=None
with (root/'execution.log').open('w') as log:
    p=subprocess.Popen([sys.executable,'-u',str(root/'run.py')],cwd=root,stdout=log,stderr=subprocess.STDOUT)
    while p.poll() is None:
        try:
            text=subprocess.check_output(['ps','-o','rss=','-p',str(p.pid)],text=True).strip()
            rss=int(text)*1024 if text else 0;peak=max(peak,rss)
        except (subprocess.CalledProcessError,ValueError):pass
        if time.monotonic()-start>1800:reason='wall cap'
        if peak>4294967296:reason='RSS cap'
        if reason:
            p.terminate()
            try:p.wait(timeout=5)
            except subprocess.TimeoutExpired:p.kill();p.wait()
            break
        time.sleep(.5)
    code=p.wait()
receipt={'exit_code':code,'watchdog_reason':reason,'seconds':time.monotonic()-start,'polled_peak_RSS_bytes':peak,'scope':'One local canary process; compilation uses num_processes=1'}
(root/'execution-receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt));sys.exit(code or (1 if reason else 0))
