"""Verify one extracted release evidence manifest without scientific dependencies."""
from pathlib import Path
import argparse,json,hashlib
a=argparse.ArgumentParser();a.add_argument('manifest');args=a.parse_args();root=Path(__file__).resolve().parent
m=json.loads((root/args.manifest).read_text())
for entry in m['files']:
 p=root/entry['path']
 with p.open('rb') as f:actual=hashlib.file_digest(f,'sha256').hexdigest()
 assert actual==entry['sha256'],entry['path']
 assert p.stat().st_size==entry['bytes'],entry['path']
print(f"Verified {len(m['files'])} evidence files")
