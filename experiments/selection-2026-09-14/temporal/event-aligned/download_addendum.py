from pathlib import Path
import modal,json,hashlib
ROOT=Path(__file__).resolve().parent;v=modal.Volume.from_name('pulsar-temporal-20260914')
def sha(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main():
 entries=[e for e in v.listdir('/event-aligned-r2',recursive=True) if e.type.name=='FILE'];done={e.path.lstrip('/').split('/')[1] for e in entries if e.path.endswith('/receipt.json')}
 for target in sorted(done):
  receipt=json.loads(b''.join(v.read_file('/event-aligned-r2/'+target+'/receipt.json')))
  for e in entries:
   rel=e.path.lstrip('/')
   if not rel.startswith('event-aligned-r2/'+target+'/'):continue
   dest=ROOT/'results'/rel.removeprefix('event-aligned-r2/');dest.parent.mkdir(parents=True,exist_ok=True)
   if dest.exists() and dest.stat().st_size==e.size:
    if dest.name!='new-full-kernels.npy' or sha(dest)==receipt['full_kernel_sha256']:continue
   with dest.open('wb') as f:
    for b in v.read_file(e.path):f.write(b)
   print(rel,dest.stat().st_size,flush=True)
  assert sha(ROOT/'results'/target/'new-full-kernels.npy')==receipt['full_kernel_sha256']
if __name__=='__main__':main()
