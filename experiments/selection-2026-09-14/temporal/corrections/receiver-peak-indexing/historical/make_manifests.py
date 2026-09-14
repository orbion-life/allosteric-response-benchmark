"""Record complete local evidence and explicit compact-publication omissions."""
from pathlib import Path
import hashlib,json,datetime
ROOT=Path(__file__).resolve().parent
EXCLUDED={'FULL-LOCAL-FILE-MANIFEST.json','COMPACT-PUBLICATION-MANIFEST.json','PUBLIC-FILE-MANIFEST.json'}
def sha(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def rows(paths):return [{'file':str(p.relative_to(ROOT)),'bytes':p.stat().st_size,'sha256':sha(p)} for p in sorted(set(paths))]
paths=[]
for p in ROOT.iterdir():
 if p.is_file() and p.name not in EXCLUDED and (p.suffix in ['.py','.json','.md'] or p.name=='LICENSE'):paths.append(p)
for folder in ['vendor','inputs','results','tiny-field','figures','event-aligned']:
 for p in (ROOT/folder).rglob('*'):
  if p.is_file() and '__pycache__' not in p.parts and p.name!='progress.json' and p.suffix in ['.py','.json','.npz','.npy','.png','.pdf','.md']:paths.append(p)
full=rows(paths);full_doc={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'Complete current local scientific evidence, including reproducible derived waveforms; not a publication receipt','file_count':len(full),'total_bytes':sum(x['bytes'] for x in full),'files':full}
(ROOT/'FULL-LOCAL-FILE-MANIFEST.json').write_text(json.dumps(full_doc,indent=2)+'\n')
omitted=[x for x in full if '/results/' in '/'+x['file'] and Path(x['file']).name.endswith('.npz') and (Path(x['file']).name=='step.npz' or Path(x['file']).name.startswith('pulse-'))]
omit_paths={x['file'] for x in omitted};published=[x for x in full if x['file'] not in omit_paths]
p=ROOT/'FULL-LOCAL-FILE-MANIFEST.json';published.append({'file':p.name,'bytes':p.stat().st_size,'sha256':sha(p)})
doc={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'Reviewed compact publication selection; no files were published by this branch','public_file_count':len(published),'public_total_bytes':sum(x['bytes'] for x in published),'scope':'All immutable raw full-kernel arrays, fixed model/operator inputs, source/protocol/hash receipts, numerical/decision summaries, direct tiny-field raw evidence and verification records are retained.','omission_reason':'Derived step/pulse matrices, scores and rankings are algebraic outputs of the included raw kernels and fixed operators; exact original hashes are retained here and in the full local manifest. No raw reference kernel or original scientific failure is omitted.','regeneration_command':'python replay_waveforms.py --output NEW_REPLAY_DIRECTORY','verification_scope':'Regeneration writes a separate tree, recreates waveforms with the sealed analyzers and runs both complete saved-array audits. It does not refit or recompute the full Gaussian reference kernels.','original_saved_array_audits':['saved-verification.json','event-aligned/saved-verification.json'],'portable_regenerated_waveform_audit':'portable-replay-verification.json','omitted_derived_file_count':len(omitted),'omitted_derived_total_bytes':sum(x['bytes'] for x in omitted),'omitted_derived_files':omitted,'files':sorted(published,key=lambda x:x['file'])}
(ROOT/'COMPACT-PUBLICATION-MANIFEST.json').write_text(json.dumps(doc,indent=2)+'\n');print({k:doc[k] for k in ['public_file_count','public_total_bytes','omitted_derived_file_count','omitted_derived_total_bytes']})
