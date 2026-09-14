"""Curate this branch, stage only its new Git paths, and build a lossless ZIP.

This script performs no git operation and no network operation. Run only after
science source/results are frozen. The external ZIP receipt avoids self-hashes.
"""
import argparse
import datetime
import hashlib
import json
import re
import shutil
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent
PREFIX = 'pulsar-selection-experiments-2026-09-14/quantum/'
OMIT = {'runtime', '__pycache__', 'qa-ci-source-only'}


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def repository_file(rel):
    if rel.parts[0] == 'history':
        return False
    if rel.suffix == '.py':
        return True
    if len(rel.parts) == 1 and rel.suffix in {'.json', '.md', '.txt'}:
        return True
    if str(rel) in {'inputs/triangle.npz', 'inputs/provenance.json',
                    'results/semantic-tests.json',
                    'phase-gradient-followup/results/receipt.json',
                    'phase-gradient-followup/results/semantic-tests.json',
                    'phase-gradient-followup/results/component-probe-summary.json'}:
        return True
    return (rel.parts[0] == 'phase-gradient-followup' and len(rel.parts) == 2
            and rel.suffix in {'.json', '.md', '.txt'})


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--repository', type=Path, required=True)
    parser.add_argument('--asset', type=Path, required=True)
    args = parser.parse_args()
    repo = args.repository.resolve()
    asset = args.asset.resolve()
    assert repo.is_dir()
    paths = sorted(p for p in ROOT.rglob('*') if p.is_file()
                   and not any(x in OMIT for x in p.relative_to(ROOT).parts)
                   and p.suffix not in {'.pyc', '.pyo'}
                   and p.name not in {'.DS_Store', 'publication-manifest.json',
                                      'publication-verification.json'})
    # Preserve source seals and copied-input provenance without editing bytes.
    checks = []
    for seal in [ROOT / 'source-seal.json', ROOT / 'phase-gradient-followup/source-seal.json']:
        for name, expected in json.loads(seal.read_text())['files'].items():
            path = seal.parent / name
            assert digest(path) == expected, path
        checks.append(str(seal.relative_to(ROOT)))
    for entry in json.loads((ROOT / 'inputs/provenance.json').read_text())['files']:
        assert digest(ROOT / entry['path']) == entry['sha256'], entry['path']
    # A narrow actual-secret pattern check; only paths would be reported on failure.
    secret = re.compile(rb'-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----|\bAKIA[0-9A-Z]{16}\b|\bgh[pousr]_[A-Za-z0-9]{30,}\b')
    for path in paths:
        if path.suffix in {'.py', '.json', '.md', '.txt', '.log', '.yaml', '.yml'}:
            assert secret.search(path.read_bytes()) is None, f'Review secret pattern in {path.relative_to(ROOT)}'
    entries = [dict(path=str(p.relative_to(ROOT)), bytes=p.stat().st_size,
                    sha256=digest(p), repository=repository_file(p.relative_to(ROOT))) for p in paths]
    manifest = dict(status='CURATED_PUBLIC_BRANCH', archive_prefix=PREFIX,
                    excludes=['Installed runtime and dependencies', 'Python caches',
                              'Redundant relocated CI copy', 'Mutable publication receipt'],
                    scope='Own source, generated circuits/results, copied own prior inputs and experiment provenance. No full articles, personal memory, challenge PDFs or cloud account inventory.',
                    minimum_CI_files=['ci_smoke.py', 'select_swap.py', 'requirements-ci.txt'],
                    repository_scope='Current source, protocols, small triangle fixture and compact decision receipts. Full archived verification requires the release overlay.',
                    source_seals_checked=checks, input_provenance_checked=True,
                    files=entries)
    index = ROOT / 'publication-manifest.json'
    index.write_text(json.dumps(manifest, indent=2) + '\n')
    paths.append(index)
    destination = repo / 'experiments/selection-2026-09-14/quantum'
    destination.mkdir(parents=True, exist_ok=True)
    copied = []
    for path in paths:
        rel = path.relative_to(ROOT)
        if path == index or repository_file(rel):
            target = destination / rel
            target.parent.mkdir(parents=True, exist_ok=True)
            if target.exists():
                assert target.read_bytes() == path.read_bytes(), f'Existing staged file differs: {target}'
            else:
                shutil.copy2(path, target)
            copied.append(str(rel))
    workflow = repo / '.github/workflows/selection-quantum.yml'
    workflow_text = '''name: Selection quantum source semantics
on:
  push:
    paths:
      - 'experiments/selection-2026-09-14/quantum/**'
      - '.github/workflows/selection-quantum.yml'
  pull_request:
    paths:
      - 'experiments/selection-2026-09-14/quantum/**'
      - '.github/workflows/selection-quantum.yml'
  workflow_dispatch:
permissions:
  contents: read
jobs:
  source-semantics:
    runs-on: ubuntu-latest
    timeout-minutes: 10
    env:
      OMP_NUM_THREADS: '1'
      OPENBLAS_NUM_THREADS: '1'
      MKL_NUM_THREADS: '1'
    defaults:
      run:
        working-directory: experiments/selection-2026-09-14/quantum
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-python@v5
        with:
          python-version: '3.13'
          cache: pip
          cache-dependency-path: experiments/selection-2026-09-14/quantum/requirements-ci.txt
      - run: python -m pip install -r requirements-ci.txt
      - name: Exact lookup, dirty-workspace and signed-phase checks
        run: python ci_smoke.py
'''
    workflow.parent.mkdir(parents=True, exist_ok=True)
    if workflow.exists():
        assert workflow.read_text() == workflow_text, 'Existing workflow differs; no overwrite'
    else:
        workflow.write_text(workflow_text)
    assert not asset.exists(), 'Do not overwrite an evidence archive'
    asset.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(asset, 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=6, allowZip64=True) as archive:
        for path in paths:
            archive.write(path, PREFIX + str(path.relative_to(ROOT)))
    with zipfile.ZipFile(asset) as archive:
        bad = archive.testzip()
        assert bad is None, bad
        assert len(set(archive.namelist())) == len(paths)
        for entry in entries:
            content = archive.read(PREFIX + entry['path'])
            assert hashlib.sha256(content).hexdigest() == entry['sha256'], entry['path']
    receipt = dict(status='PASS', completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                   asset=str(asset), asset_bytes=asset.stat().st_size, asset_sha256=digest(asset),
                   archive_members=len(paths), uncompressed_bytes=sum(p.stat().st_size for p in paths),
                   CRC_test='PASS', all_manifest_member_SHA256='PASS',
                   manifest_sha256=digest(index), repository=str(destination),
                   repository_file_count=len(copied), repository_bytes=sum((destination / p).stat().st_size for p in copied),
                   repository_files=copied, workflow=str(workflow), workflow_sha256=digest(workflow),
                   no_git_commit_or_push=True, no_network_operation=True)
    (ROOT / 'publication-verification.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({k: v for k, v in receipt.items() if k != 'repository_files'}, indent=2))


if __name__ == '__main__':
    main()
