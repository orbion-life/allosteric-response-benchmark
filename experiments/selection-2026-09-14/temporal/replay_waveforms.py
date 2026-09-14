"""Regenerate omitted derived arrays in a separate tree from public raw kernels.

Inputs and raw arrays are hardlinked where possible, otherwise copied. Only the
new tree receives analysis outputs. Existing archive files remain unchanged.
"""
from pathlib import Path
import os,shutil,json,argparse,subprocess,sys,hashlib,time
ROOT=Path(__file__).resolve().parent

def materialize(src,dst):
 dst.parent.mkdir(parents=True,exist_ok=True)
 try:os.link(src,dst)
 except OSError:shutil.copy2(src,dst)

def replay(destination):
 dest=Path(destination).resolve();dest.mkdir(parents=True,exist_ok=False);start=time.perf_counter()
 # Copy source, protocols and numerical receipts; reconstruct all omitted waveforms.
 for p in ROOT.glob('*.py'):shutil.copy2(p,dest/p.name)
 for p in ROOT.glob('*.json'):
  if 'MANIFEST' not in p.name and p.name not in ['summary.json']:shutil.copy2(p,dest/p.name)
 shutil.copy2(ROOT/'input-manifest.json',dest/'input-manifest.json')
 for sub in ['vendor','inputs','tiny-field']:
  for p in (ROOT/sub).rglob('*'):
   if p.is_file() and '__pycache__' not in p.parts and p.suffix in ['.py','.json','.npz']:materialize(p,dest/p.relative_to(ROOT))
 for target in ['kras','abl','myc','myh7']:
  for name in ['full-kernels.npy','receipt.json','fixture-check.json','saved-reference-canary.json']:
   materialize(ROOT/'results'/target/name,dest/'results'/target/name)
  log=dest/f'replay-{target}.log'
  with log.open('w') as f:subprocess.run([sys.executable,'analyze.py',target],cwd=dest,check=True,stdout=f,stderr=subprocess.STDOUT)
  print('original-grid waveform rebuilt',target,flush=True)
 with (dest/'verification.log').open('w') as f:subprocess.run([sys.executable,'verify.py','--mode','saved','--output','saved-verification.json'],cwd=dest,check=True,stdout=f,stderr=subprocess.STDOUT)
 if (ROOT/'event-aligned').exists():
  ev=dest/'event-aligned';ev.mkdir(exist_ok=True)
  for p in (ROOT/'event-aligned').iterdir():
   if p.is_file() and p.suffix in ['.py','.json']:shutil.copy2(p,ev/p.name)
  for target in ['kras','abl','myc','myh7']:
   for name in ['new-full-kernels.npy','receipt.json','canary.json','fixture-check.json']:
    materialize(ROOT/'event-aligned/results'/target/name,ev/'results'/target/name)
   with (ev/f'replay-{target}.log').open('w') as f:subprocess.run([sys.executable,'analyze_addendum.py',target],cwd=ev,check=True,stdout=f,stderr=subprocess.STDOUT)
   print('event-aligned waveform rebuilt',target,flush=True)
  with (ev/'verification.log').open('w') as f:subprocess.run([sys.executable,'verify_addendum.py','--mode','saved','--output','saved-verification.json'],cwd=ev,check=True,stdout=f,stderr=subprocess.STDOUT)
 receipt={'status':'PASS','scope':'Fresh derived waveform reconstruction and offline verification from package raw kernels; not full Gaussian kernel recomputation','seconds':time.perf_counter()-start,'output':str(dest)}
 (dest/'portable-replay-receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);a=p.parse_args();replay(a.output)
