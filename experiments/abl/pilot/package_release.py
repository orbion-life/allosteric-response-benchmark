#!/usr/bin/env python3
"""Create a portable scientific-only package without duplicate replay arrays."""
from pathlib import Path
import json,shutil,zipfile
from run import ROOT,sha,dump
out=ROOT/'public-package';out.mkdir(exist_ok=True)
files=[p for p in ROOT.iterdir() if p.is_file() and p.suffix in ['.py','.md','.txt','.json','.sha256']]
files.append(ROOT/'LICENSE')
for folder in ['raw','vendor','model','results','evaluation','figures','checks']:
    files.extend(p for p in (ROOT/folder).glob('*') if p.is_file() and p.suffix!='.pyc')
for src in sorted(set(files)):
    relative=src.relative_to(ROOT);target=out/relative;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src,target)
manifest={str(p.relative_to(out)):sha(p) for p in sorted(out.rglob('*')) if p.is_file() and p.name!='SHA256SUMS.json'}
dump(out/'SHA256SUMS.json',manifest)
with zipfile.ZipFile(ROOT/'project-pulsar-abl-pilot.zip','w',zipfile.ZIP_DEFLATED) as z:
    for p in sorted(out.rglob('*')):
        if p.is_file():z.write(p,Path('project-pulsar-abl-pilot')/p.relative_to(out))
print(json.dumps({'files':len(manifest),'zip_sha256':sha(ROOT/'project-pulsar-abl-pilot.zip'),'zip_bytes':(ROOT/'project-pulsar-abl-pilot.zip').stat().st_size}))
