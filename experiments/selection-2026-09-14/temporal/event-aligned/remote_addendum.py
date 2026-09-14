"""Explicit second frozen GPU characterization; original outputs stay unchanged."""
from pathlib import Path
import modal,json,time
ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent
image=(modal.Image.debian_slim(python_version='3.13').pip_install('numpy==2.2.6','scipy==1.15.3','cupy-cuda12x[ctk]==14.1.1','psutil==7.0.0')
 .add_local_dir(BASE/'inputs','/work/temporal/inputs').add_local_dir(BASE/'vendor','/work/temporal/vendor')
 .add_local_file(BASE/'protocol.json','/work/temporal/protocol.json').add_local_file(BASE/'input-manifest.json','/work/temporal/input-manifest.json').add_local_file(BASE/'analyze.py','/work/temporal/analyze.py')
 .add_local_dir(ROOT,'/work/temporal/event-aligned',ignore=['results','*.log','__pycache__']))
volume=modal.Volume.from_name('pulsar-temporal-20260914')
app=modal.App('pulsar-temporal-event-aligned-20260914')
@app.function(image=image,gpu='A100-40GB',cpu=3,memory=12288,timeout=1200,startup_timeout=600,retries=0,max_containers=1,scaledown_window=2,volumes={'/results':volume},serialized=True,include_source=False)
def campaign():
 import sys
 sys.path.insert(0,'/work/temporal/event-aligned');from full_kernel_addendum import run
 start=time.perf_counter();rows=[]
 for target in ['kras','abl','myc','myh7']:
  rows.append(run(target,Path('/results/event-aligned-r2')/target));volume.commit()
 receipt={'wall_seconds':time.perf_counter()-start,'records':rows};Path('/results/event-aligned-r2/campaign.json').write_text(json.dumps(receipt,indent=2));volume.commit();return receipt
@app.local_entrypoint()
def main():
 start=time.perf_counter();r=campaign.remote();r['client_wall_seconds']=time.perf_counter()-start;(ROOT/'remote-receipt.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
