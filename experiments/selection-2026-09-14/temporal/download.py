"""Fetch completed campaign targets; verify hashes before reusing local large arrays."""
from pathlib import Path
import modal,json,hashlib
ROOT=Path(__file__).resolve().parent
v=modal.Volume.from_name('pulsar-temporal-20260914')
def sha(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main():
 out=ROOT/'results';out.mkdir(exist_ok=True);entries=[e for e in v.listdir('/',recursive=True) if e.type.name=='FILE']
 completed={e.path.lstrip('/').split('/')[0] for e in entries if e.path.endswith('/receipt.json')}
 for target in sorted(completed):
  raw=b''.join(v.read_file('/'+target+'/receipt.json'));receipt=json.loads(raw)
  for e in entries:
   rel=e.path.lstrip('/')
   if not rel.startswith(target+'/'):continue
   dest=out/rel;dest.parent.mkdir(parents=True,exist_ok=True)
   if dest.exists() and dest.stat().st_size==e.size:
    if dest.name=='full-kernels.npy' and sha(dest)==receipt['full_kernel_sha256']:continue
    if dest.name!='full-kernels.npy':continue
   with dest.open('wb') as f:
    for b in v.read_file(e.path):f.write(b)
   print(rel,dest.stat().st_size,flush=True)
  assert sha(out/target/'full-kernels.npy')==receipt['full_kernel_sha256']
if __name__=='__main__':main()
