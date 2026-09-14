"""Build the portable evidence asset locally. No network or Git operations."""
import datetime,hashlib,json,re,zipfile
from pathlib import Path
ROOT=Path(__file__).resolve().parent
PREFIX='pulsar-biological-quantum-resolution-2026-09-14/quantum/'

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def repository_file(rel):
    if rel.parts[0]=='history':return False
    if rel.suffix=='.py':return True
    if len(rel.parts)==1 and rel.suffix in {'.md','.json','.txt'}:return True
    return str(rel) in {'inputs/triangle.npz','inputs/provenance.json'}

def main():
    sources=sorted(p for p in ROOT.rglob('*.py') if not any(x in p.parts for x in ['history','delivery','ci-relocated-check','__pycache__']))
    seal=dict(recorded_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),status='FINAL_CURRENT_SOURCES_AND_RESULTS_FROZEN',
              canonical_IQAE_entrypoint='run_iqae.py',
              historical_seals='Earlier source seals remain unchanged; original sources and explicit corrections are retained in history. Final current code hashes are given here.',
              files={str(p.relative_to(ROOT)):sha(p) for p in sources})
    (ROOT/'final-source-seal.json').write_text(json.dumps(seal,indent=2)+'\n')
    omit={'delivery','__pycache__','ci-relocated-check','runtime'}
    paths=sorted(p for p in ROOT.rglob('*') if p.is_file() and not any(x in omit for x in p.relative_to(ROOT).parts)
                 and p.suffix not in {'.pyc','.pyo'} and p.name not in {'.DS_Store','publication-manifest.json','IQAE-kras-profile.txt'}
                 and p.name!='cached-copy-IQAE-all-outputs.npz')
    secrets=re.compile(rb'-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----|\bAKIA[0-9A-Z]{16}\b|\bgh[pousr]_[A-Za-z0-9]{30,}\b')
    for p in paths:
        if p.suffix in {'.py','.md','.json','.txt','.log'}:assert not secrets.search(p.read_bytes()),p
    files=[dict(path=str(p.relative_to(ROOT)),bytes=p.stat().st_size,sha256=sha(p),repository=repository_file(p.relative_to(ROOT))) for p in paths]
    manifest=dict(status='CURATED_PUBLIC_EVIDENCE',archive_prefix=PREFIX,
                  excludes=['Installed runtimes and Python caches','Redundant relocated CI directory','Internal process profiler dumps','Redundant cached raw copy; primary cached raw retained'],
                  minimum_source_CI=['ci_smoke.py','iqae.py','spectral.py','vendor/multiplexed.py','requirements-ci.txt'],
                  scope='Own scientific source, protocols, inputs, generated circuits/results and explicit correction history. No full articles, account inventories, credentials or personal memory.',files=files)
    index=ROOT/'publication-manifest.json';index.write_text(json.dumps(manifest,indent=2)+'\n');paths.append(index)
    dest=ROOT/'delivery';dest.mkdir(exist_ok=True);archive=dest/'pulsar-protein-quantum-cost-evidence.zip'
    assert not archive.exists(),'Preserve existing delivery; do not overwrite'
    with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for p in paths:z.write(p,PREFIX+str(p.relative_to(ROOT)))
    with zipfile.ZipFile(archive) as z:
        assert z.testzip() is None
        assert len(z.namelist())==len(set(z.namelist()))==len(paths)
        for item in files:assert hashlib.sha256(z.read(PREFIX+item['path'])).hexdigest()==item['sha256']
    result=dict(status='PASS',artifact=archive.name,bytes=archive.stat().st_size,sha256=sha(archive),members=len(paths),
                uncompressed_bytes=sum(p.stat().st_size for p in paths),CRC='PASS',every_manifest_SHA256='PASS',
                manifest_sha256=sha(index),source_seal_sha256=sha(ROOT/'final-source-seal.json'),
                repository_files=[f['path'] for f in files if f['repository']]+['publication-manifest.json'],
                source_only_CI='python -m pip install -r requirements-ci.txt; python -B ci_smoke.py',
                no_publication_performed=True)
    (dest/'verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='repository_files'},indent=2))

if __name__=='__main__':main()
