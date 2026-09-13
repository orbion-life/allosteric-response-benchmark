"""One sequential numerical worker under a single recorded resource allocation."""
from pathlib import Path
import datetime,hashlib,json,os,signal,subprocess,sys,time
ROOT=Path(__file__).resolve().parent

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dump(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
protocol_path=ROOT/(sys.argv[1] if len(sys.argv)>1 else 'protocol.json')
protocol=json.loads(protocol_path.read_text())
assert sha(protocol_path)==protocol_path.with_suffix('.sha256').read_text().strip()
execution_path=ROOT/protocol.get('execution_file','execution.json')
prior_seconds=protocol.get('prior_execution_seconds',0.)
for path,expected in protocol['code_sha256'].items():assert sha(ROOT/path)==expected,path
if execution_path.exists():raise RuntimeError('An execution already exists; never silently restart the allocation.')
env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',PYTHONWARNINGS='error')
for key in protocol['single_thread_environment']:env[key]='1'
start=time.monotonic();receipt={'status':'RUNNING','started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'protocol_sha256':sha(protocol_path),'prior_execution_seconds':prior_seconds,'jobs':[],'sampled_peak_aggregate_bytes':0,'limits':protocol['limits']};dump(execution_path,receipt)
try:
 for job in protocol['jobs']:
  elapsed=prior_seconds+time.monotonic()-start
  if elapsed>=protocol['limits']['wall_seconds']:raise TimeoutError('Aggregate wall-time cap')
  label=job.get('label','pinned');exe=protocol['executables'][label];argv=[exe,'-W','error',str(ROOT/'diagnostic.py'),job['mode'],'--label',label]
  if 'case' in job:argv+=['--case',str(job['case'])]
  if job.get('shared'):argv+=['--shared']
  if job.get('repair'):argv+=['--repair']
  log=ROOT/protocol.get('log_directory','logs')/('%02d-%s.log'%(len(receipt['jobs']),job['name']));log.parent.mkdir(exist_ok=True)
  row={'name':job['name'],'argv':argv,'log':str(log.relative_to(ROOT)),'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()};begin=time.monotonic();reason=None
  with log.open('w') as stream:
   proc=subprocess.Popen(argv,cwd=ROOT,env=env,stdout=stream,stderr=subprocess.STDOUT,start_new_session=True)
   while proc.poll() is None:
    usage=0
    for line in subprocess.check_output(['ps','-axo','pid=,pgid=,rss='],text=True).splitlines():
     pid,group,rss=map(int,line.split())
     if group==proc.pid or pid==os.getpid():usage+=rss*1024
    receipt['sampled_peak_aggregate_bytes']=max(receipt['sampled_peak_aggregate_bytes'],usage)
    if usage>protocol['limits']['aggregate_bytes']:reason='aggregate memory cap'
    if prior_seconds+time.monotonic()-start>protocol['limits']['wall_seconds']:reason='aggregate time cap'
    if time.monotonic()-begin>job['max_seconds']:reason='job time cap'
    if reason:
     os.killpg(proc.pid,signal.SIGTERM)
     try:proc.wait(timeout=3)
     except subprocess.TimeoutExpired:os.killpg(proc.pid,signal.SIGKILL)
     break
    time.sleep(.25)
   code=proc.wait()
  row.update(returncode=code,seconds=time.monotonic()-begin,limit_reason=reason,log_sha256=sha(log));receipt['jobs'].append(row);receipt['wall_seconds']=time.monotonic()-start;dump(execution_path,receipt);print(json.dumps(row),flush=True)
  if code or reason:raise RuntimeError('Job stopped: '+job['name'])
 receipt['status']='COMPLETE'
except BaseException as exc:
 receipt.update(status='STOPPED',error=type(exc).__name__+': '+str(exc));raise
finally:
 receipt.update(wall_seconds=time.monotonic()-start,aggregate_execution_seconds=prior_seconds+time.monotonic()-start,finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat());dump(execution_path,receipt)
 hashes={str(p.relative_to(ROOT)):sha(p) for p in sorted(ROOT.rglob('*')) if p.is_file() and p.name!='SHA256SUMS.json'}
 dump(ROOT/'SHA256SUMS.json',hashes)
