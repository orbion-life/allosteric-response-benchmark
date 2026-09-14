"""Build the small v0.10.1 correction overlay without duplicating raw kernels."""
from pathlib import Path
import argparse, hashlib, json, zipfile

ROOT=Path(__file__).resolve().parent
HISTORY=ROOT/'history/receiver-ranking-before-2026-09-15'

def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream,'sha256').hexdigest()

def entry(path):
    return dict(path=str(path.relative_to(ROOT)),bytes=path.stat().st_size,sha256=digest(path))

def dump(path,data):
    path.write_text(json.dumps(data,indent=2)+'\n')

def files():
    # A full extracted evidence tree contains large raw arrays and may contain
    # local replay outputs. Keep the overlay to the declared public files.
    top_names = {
        'DATA-AND-SOFTWARE-ATTRIBUTION.md', 'LICENSE', 'PACKAGING.json',
        'README.md', 'CORRECTION-2026-09-15.md', 'CORRECTION-MANIFEST.json',
        'analysis-reporting-corrections.json', 'analyze.py', 'campaign.py',
        'input-manifest.json', 'portable-replay-verification.json', 'protocol.json',
        'requirements.txt', 'source-seal.json', 'summary.json', 'verify_package.py',
        'ranking.py', 'test_receiver_ranking.py', 'replay_correction.py',
        'verify_correction.py', 'verify_selection_scope.py',
        'make_correction_package.py', 'selection-scope-verification.json',
        'correction-source-manifest.json',
    }
    selected = []
    for path in ROOT.rglob('*'):
        if not path.is_file() or '__pycache__' in path.parts or path.suffix == '.pyc':
            continue
        relative = path.relative_to(ROOT)
        if len(relative.parts) == 1:
            keep = path.name in top_names or path.name.startswith('PUBLIC-MANIFEST-')
        elif relative.parts[0] == 'results':
            keep = len(relative.parts) == 3 and path.name in {'analysis.json', 'receipt.json', 'construction.json', 'canary.json'}
        elif relative.parts[0] in {'history', 'vendor', 'provenance', 'correction-verification'}:
            keep = path.suffix in {'.py', '.json', '.md', '.txt'}
        else:
            keep = relative.parts[0] == 'fixtures' and path.name in {f'{t}-tau.npz' for t in ['kras', 'abl', 'myc', 'myh7']}
        if keep:
            selected.append(path)
    return sorted(selected)

def build(output):
    original_inputs=[]
    correction=json.loads((ROOT/'correction-verification/summary.json').read_text())
    for row in correction['targets']:
        target=row['target'];hashes=row['source_hashes']
        names={'metadata':f'inputs/{target}/ranking-metadata.npz','operator':f'results/{target}/operator.npz','fresh_kernels':f'results/{target}/fresh-kernels.npy','old_kernels':f'reference/results/{target}/full-kernels.npy','event_kernels':f'reference/event-aligned/results/{target}/new-full-kernels.npy'}
        original_inputs.extend(dict(target=target,path=names[key],sha256=value,status='unchanged original saved input; supplied by v0.9.0 archives') for key,value in hashes.items())
    provenance=dict(date='2026-09-15',release='v0.10.1',status='CORRECTED_FROM_ORIGINAL_SAVED_ARRAYS',sources=[entry(p) for p in files() if p.suffix=='.py' or p.name in ['protocol.json','input-manifest.json','requirements.txt']],original_inputs=original_inputs,historical_provenance=[entry(p) for p in sorted(HISTORY.rglob('*')) if p.is_file()],scientific_status='175/176 informative fresh shortlists preserved; all 288 fresh matrix errors pass; original KRAS mass-orthogonality failure and all historical bound-portability failures unchanged.')
    dump(ROOT/'correction-source-manifest.json',provenance)
    common=[p for p in files() if not (p.parent==ROOT and (p.name.startswith('PUBLIC-MANIFEST-') or p.name=='CORRECTION-MANIFEST.json'))]
    for historical_manifest in sorted(HISTORY.glob('PUBLIC-MANIFEST-*.json')):
        original=json.loads(historical_manifest.read_text())
        indexed={e['path']:e for e in original['files']}
        # Existing raw entries retain their exact original size and digest.
        # Changed compact outputs and added provenance are supplied by the overlay.
        indexed.update({str(p.relative_to(ROOT)):entry(p) for p in common})
        result=dict(release='v0.10.1',study='temporal',target=original['target'],part=original['part'],requires='Original v0.9.0 target raw archive(s), then the v0.10.1 receiver correction overlay',historical_manifest=str(historical_manifest.relative_to(ROOT)),files=[indexed[k] for k in sorted(indexed)])
        dump(ROOT/historical_manifest.name,result)
    members=[p for p in files() if p.name!='CORRECTION-MANIFEST.json']
    manifest=dict(release='v0.10.1',study='temporal-receiver-correction',date='2026-09-15',scope='Portable source, corrected summaries, independent verification, actual-data regression fixtures and preserved historical provenance. Original large arrays are reused without modification.',original_input_manifest='correction-source-manifest.json',files=[entry(p) for p in members])
    dump(ROOT/'CORRECTION-MANIFEST.json',manifest)
    output.parent.mkdir(parents=True,exist_ok=True)
    if output.exists():raise FileExistsError(output)
    with zipfile.ZipFile(output,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as archive:
        for path in members+[ROOT/'CORRECTION-MANIFEST.json']:
            archive.write(path,'temporal/'+str(path.relative_to(ROOT)))
    print(json.dumps(dict(archive=output.name,bytes=output.stat().st_size,sha256=digest(output),files=len(members)+1),indent=2))

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();build(args.output)
