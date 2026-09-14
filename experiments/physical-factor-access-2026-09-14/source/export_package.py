"""Export identical scientific bytes with explicit metadata path redaction."""
from pathlib import Path
import hashlib,json,shutil,zipfile,datetime
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT.parent
TARGET=BASE/'repository/experiments/physical-factor-access-2026-09-14'
ZIP=BASE/'delivery/pulsar-physical-factor-access-2026-09-14.zip'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    TARGET.mkdir(parents=True,exist_ok=True);changed=[]
    for source in ROOT.rglob('*'):
        if not source.is_file():continue
        rel=source.relative_to(ROOT)
        if '__pycache__' in rel.parts or rel.parts[:2]==('checks','replay') or rel.name in ['MANIFEST.sha256','portable-check.log']:continue
        target=TARGET/rel;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source,target)
    protocol=TARGET/'protocol.json';p=json.loads(protocol.read_text())
    for name,spec in p['input']['files'].items():
        group='results' if name=='operator.npz' else 'inputs'
        spec['original']=f'prior-response-evidence/temporal/{group}/abl/{name}'
    protocol.write_text(json.dumps(p,indent=2)+'\n');changed.append({'path':'protocol.json','reason':'Private absolute provenance paths replaced by descriptive relative provenance labels. No scientific field changed.','original_sha256':sha(ROOT/'protocol.json'),'exported_sha256':sha(protocol)})
    original_sha=sha(ROOT/'protocol.json');(TARGET/'provenance/original-protocol.sha256').write_text(original_sha+'\n');(TARGET/'protocol.sha256').write_text(sha(protocol)+'\n')
    audit=TARGET/'review/independent-math-audit/operator-Hermite-comparison.json';j=json.loads(audit.read_text());j['source_under_test']='source/core.py';audit.write_text(json.dumps(j,indent=2)+'\n');changed.append({'path':str(audit.relative_to(TARGET)),'reason':'Private source path redacted; source_hashes and all scientific checks unchanged.','original_sha256':sha(ROOT/audit.relative_to(TARGET)),'exported_sha256':sha(audit)})
    large=TARGET/'review/independent-math-audit/operator-results-check.json'
    if large.exists():
        j=json.loads(large.read_text());j['source_hashes']={str(Path(k).relative_to(ROOT)):v for k,v in j['source_hashes'].items()};large.write_text(json.dumps(j,indent=2)+'\n')
        changed.append({'path':str(large.relative_to(TARGET)),'reason':'Private absolute source-path keys changed to package-relative keys; hash values and scientific results unchanged.','original_sha256':sha(ROOT/large.relative_to(TARGET)),'exported_sha256':sha(large)})
    code=[]
    for folder in ['source','vendor','checks','review']:
        for source in (ROOT/folder).rglob('*.py'):
            target=TARGET/source.relative_to(ROOT)
            assert sha(source)==sha(target)
            code.append({'path':str(source.relative_to(ROOT)),'sha256':sha(source)})
    receipt={'status':'Public export prepared; no publication performed by this script','created_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'original_protocol_sha256':original_sha,'public_protocol_sha256':sha(protocol),'path_redactions':changed,'unchanged_code':code,'scientific_inputs_unchanged':{str(s.relative_to(ROOT)):sha(s) for s in (ROOT/'inputs').glob('*.npz')},'archived_receipt_hash_scope':'Archived audit receipts refer to their original local layout and hashes. Redacted metadata paths are mapped here. Portable verification emits new receipts for the exported layout.','portable_wrapper':'checks/verify.py overrides only BASE/OP/HERE directory variables in unchanged independent audit source.'}
    (TARGET/'provenance/public-export.json').write_text(json.dumps(receipt,indent=2)+'\n')
    manifest=[]
    for p in sorted(TARGET.rglob('*')):
        if p.is_file() and p.name!='MANIFEST.sha256' and '__pycache__' not in p.parts and 'replay' not in p.relative_to(TARGET).parts:
            manifest.append(f'{sha(p)}  {p.relative_to(TARGET)}')
    (TARGET/'MANIFEST.sha256').write_text('\n'.join(manifest)+'\n')
    ZIP.parent.mkdir(exist_ok=True)
    with zipfile.ZipFile(ZIP,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for p in sorted(TARGET.rglob('*')):
            if p.is_file() and '__pycache__' not in p.parts and 'replay' not in p.relative_to(TARGET).parts:z.write(p,Path('physical-factor-access-2026-09-14')/p.relative_to(TARGET))
    print(json.dumps({'directory':str(TARGET),'zip':str(ZIP),'zip_bytes':ZIP.stat().st_size,'zip_sha256':sha(ZIP),'manifest_entries':len(manifest)},indent=2))

if __name__=='__main__':main()
