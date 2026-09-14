"""Build the pinned GPL-3.0 Ohm executables locally; no global installation.

Requires GCC 15 with C++23 support, make, patch and zlib development files.
The pristine archive and two-line portability patch are included with notices.
"""
from pathlib import Path
import argparse, subprocess, tarfile, shutil, hashlib, json, time
P = Path(__file__).resolve().parents[1]
ap = argparse.ArgumentParser()
ap.add_argument('--compiler', default='g++-15')
args = ap.parse_args()
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
archive = P / 'source/ohm-upstream-462a3b1.tar.gz'
patch = P / 'source/ohm-gcc15-portability.patch'
assert sha(archive) == '395dc1c3e561e21a65a8a2929331b3eb719ae3d66dba64dfe4c879589da140d5'
assert sha(patch) == 'd558b81d4870a6ad037958a36b1cc87cc06c50b606f787c0a696f9eb8f9f3b57'
work = P / 'replay-ohm/source'
assert not work.exists(), 'Use a fresh replay copy; preserve previous build attempts.'
work.mkdir(parents=True)
start = time.monotonic()
with tarfile.open(archive) as t:
    t.extractall(work, filter='data')
commands = []
def call(argv):
    result = subprocess.run(argv, cwd=work, text=True, capture_output=True)
    commands.append({'argv': argv, 'returncode': result.returncode,
                     'stdout': result.stdout, 'stderr': result.stderr})
    (P / 'replay-ohm/build-receipt.json').write_text(json.dumps({
        'status': 'BUILDING' if result.returncode == 0 else 'FAILED',
        'commands': commands}, indent=2))
    result.check_returncode()
call(['patch', '-p1', '-i', str(patch)])
# Relative prefix avoids upstream Makefile whitespace handling of absolute paths.
(work / '.prefix').write_text('install\n')
call(['make', '-j1', 'contacts', 'diffuse', 'CC=' + args.compiler,
      'FLAGS=-std=c++23 -Wno-everything -lstdc++ -Wwrite-strings -lm -Isrc -MMD -DNDEBUG -O3'])
out = P / 'ohm/bin'
out.mkdir(parents=True, exist_ok=True)
for name in ['contacts', 'diffuse']:
    shutil.copy2(work / 'install/bin' / name, out / name)
receipt = {'status': 'PASS', 'seconds': time.monotonic() - start,
           'upstream_commit': '462a3b1318e24ebe2061137fe88af719ef89c0ae',
           'archive_sha256': sha(archive), 'patch_sha256': sha(patch),
           'compiler': subprocess.check_output([args.compiler, '--version'], text=True),
           'binary_sha256': {n: sha(out / n) for n in ['contacts', 'diffuse']},
           'commands': commands,
           'scope': 'Local rebuild. Cross-platform random-stream or bitwise equality is not assumed.'}
(P / 'replay-ohm/build-receipt.json').write_text(json.dumps(receipt, indent=2))
print(json.dumps({k: v for k, v in receipt.items() if k != 'commands'}, indent=2))
