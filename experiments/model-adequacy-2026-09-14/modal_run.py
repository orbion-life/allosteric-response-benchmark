"""Isolated, bounded GPU reference execution; source and evidence retained."""
from pathlib import Path
import datetime, hashlib, json, time, os, traceback
import modal
HERE=Path(__file__).resolve().parent
NAME='pulsar-final-nonlinear-20260914'
VOL_NAME='pulsar-final-nonlinear-evidence-20260914'
image=modal.Image.debian_slim(python_version='3.12').pip_install('numpy==2.2.6','scipy==1.15.3','cupy-cuda12x[ctk]==14.1.1','psutil==7.0.0')
volume=modal.Volume.from_name(VOL_NAME,create_if_missing=True)
app=modal.App(NAME)


@app.function(image=image,gpu='A100-40GB',cpu=8,memory=65536,timeout=1800,retries=0,max_containers=1,serialized=True,volumes={'/results':volume})
def work(sources,job):
    import sys,resource
    start=time.perf_counter(); root=Path('/tmp/experiment');root.mkdir(exist_ok=True)
    out=Path('/results')/job['id'];out.mkdir(exist_ok=False)
    receipt={'job':job,'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
             'source_sha256':{name:hashlib.sha256(data.encode()).hexdigest() for name,data in sources.items()},
             'scope':'Original finite reflecting model; physical admission and claims assessed from independent receipts',
             'modal_task_id':os.environ.get('MODAL_TASK_ID'),'memory_requested_mib':65536,
             'GPU_requested':'A100-40GB','CPU_requested':8,'timeout_seconds':1800,'automatic_retries':0}
    try:
        for name,data in sources.items():
            rel=Path(name)
            if rel.is_absolute() or '..' in rel.parts: raise ValueError('invalid source path')
            p=root/rel;p.parent.mkdir(parents=True,exist_ok=True);p.write_text(data)
            q=out/'sources'/rel;q.parent.mkdir(parents=True,exist_ok=True);q.write_text(data)
        sys.path.insert(0,str(root/'dynamics'))
        from core import run_reference,backend
        xp=backend('cupy')
        receipt['GPU']=xp.cuda.runtime.getDeviceProperties(0)['name'].decode()
        if job.get('small_check'):
            from check_small import check
            receipt['small_check']=check(out/'gpu-small-check','cupy')
        records=[]
        if job['mode']=='reference':
            for n,case in enumerate(job['cases']):
                records.append(run_reference(out/f'case-{n}',case['faces'],case['spacing'],'cupy',
                                             times=job.get('times',[.1,1,10]),independent=job.get('independent',True)))
                volume.commit()
                xp.get_default_memory_pool().free_all_blocks()
        elif job['mode']=='field':
            from run_field import run
            records.append(run(out/'field',job['faces'],job['spacing'],'cupy',job.get('times',[.1,1,10])))
        else: raise ValueError(job['mode'])
        receipt.update(status='COMPLETE',records=records)
    except BaseException as e:
        receipt.update(status='FAIL',error=repr(e),traceback=traceback.format_exc())
    receipt.update(seconds_inside_function=time.perf_counter()-start,
                   peak_process_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                   completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
    receipt['raw_files']=[{'path':str(p.relative_to(out)),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
                          for p in sorted(out.rglob('*')) if p.is_file()]
    (out/'provider-receipt.json').write_text(json.dumps(receipt,indent=2));volume.commit()
    # Strip NumPy scalar pickle dependencies before returning to the minimal client.
    return json.loads(json.dumps(receipt))


@app.local_entrypoint()
def main(job_file:str):
    job=json.loads(Path(job_file).read_text());dest=HERE/'remote'/job['id'];dest.mkdir(parents=True,exist_ok=False)
    sources={}
    for name in ['core.py','check_small.py','protocol.json']+(['run_field.py'] if job['mode']=='field' else []):
        sources['dynamics/'+name]=(HERE/'source/dynamics'/name).read_text()
    for name in ['model.py','protocol.json','chebyshev.py']: sources['vendor/'+name]=(HERE/'source/vendor'/name).read_text()
    sources['study-protocol.json']=(HERE/'protocol.json').read_text()
    ledger={'job':job,'volume':VOL_NAME,'status':'RUNNING','started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
    (dest/'launch.json').write_text(json.dumps(ledger,indent=2));start=time.perf_counter()
    try:
        rec=work.remote(sources,job);(dest/'provider-receipt.json').write_text(json.dumps(rec,indent=2))
        downloaded=[]
        for row in rec['raw_files']:
            p=dest/'data'/row['path'];p.parent.mkdir(parents=True,exist_ok=True)
            with p.open('wb') as f:
                for chunk in volume.read_file('/'+job['id']+'/'+row['path']): f.write(chunk)
            actual=hashlib.sha256(p.read_bytes()).hexdigest()
            downloaded.append({**row,'verified':actual==row['sha256'] and p.stat().st_size==row['bytes']})
        ledger.update(status=rec['status'],downloads=downloaded,all_downloads_verified=all(x['verified'] for x in downloaded))
    except BaseException as e:
        ledger.update(status='FAIL',error=repr(e));raise
    finally:
        ledger.update(client_seconds_including_download=time.perf_counter()-start,completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
        (dest/'launch.json').write_text(json.dumps(ledger,indent=2))
    print(json.dumps({k:v for k,v in ledger.items() if k!='downloads'},indent=2))
