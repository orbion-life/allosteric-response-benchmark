"""Explicitly invoked private single-GPU launcher; portable numerical code is separate."""
from pathlib import Path
import modal,json,time,os
ROOT=Path(__file__).resolve().parent
image=(modal.Image.debian_slim(python_version='3.13').pip_install('numpy==2.2.6','scipy==1.15.3','cupy-cuda12x[ctk]==14.1.1','psutil==7.0.0').add_local_dir(ROOT,'/work/temporal',ignore=['results','*.log','__pycache__']))
volume=modal.Volume.from_name('pulsar-temporal-20260914',create_if_missing=True)
app=modal.App('pulsar-temporal-fixed-20260914')
@app.function(image=image,gpu='A100-40GB',cpu=3,memory=12288,timeout=2700,startup_timeout=600,retries=0,max_containers=1,scaledown_window=2,volumes={'/results':volume},serialized=True,include_source=False)
def campaign():
 import sys
 sys.path.insert(0,'/work/temporal')
 from full_kernel import run
 rows=[];start=time.perf_counter()
 for target in ['kras','abl','myc','myh7']:
  rows.append(run(target,Path('/results')/target,'cupy'));volume.commit()
 receipt={'records':rows,'wall_seconds':time.perf_counter()-start};Path('/results/campaign.json').write_text(json.dumps(receipt,indent=2));volume.commit();return receipt
@app.local_entrypoint()
def main():
 start=time.perf_counter();out=ROOT/'remote-receipt.json'
 receipt=campaign.remote();receipt['client_wall_seconds']=time.perf_counter()-start;out.write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
