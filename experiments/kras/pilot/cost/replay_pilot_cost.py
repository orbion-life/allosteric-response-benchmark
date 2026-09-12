#!/usr/bin/env python3
"""Replay the frozen KRAS pilot cost procedure without modifying its sources.

Usage: python replay_pilot_cost.py --pilot PATH --work NEW_DIRECTORY
Requires sibling pilot-cost-worker.py and a complete cached-input pilot package.
"""
from pathlib import Path
import argparse
import hashlib
import json
import shutil
import subprocess
import sys

EXPECTED_SOURCE_SHA256={
    'prepare.py':'3238b9a09773c21b51b427bedfad2b05fdeb53c6722681c6c170d40944428a8b',
    'run_pilot.py':'ec5ca6691fb0568d2ca971bf61b3a9b631ed63d7b6692f503a3ba705e343e053',
}
EXPECTED_INSTRUMENTED_PREPARE_SHA256='717ec8e8175e3bce66f409f4388795cae12bf29f198e6c93a783a4dfd4a9b3ae'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def instrument_prepare(original):
    """Reconstruct the exact measurement-only change used in the recorded run."""
    old_import='import json,hashlib,datetime,itertools,math'
    start='    dist=np.sqrt(np.sum((r0[:,None]-r0[None])**2,axis=-1));edges='
    end='    recdist=dist[:,receivers].min(axis=1);'
    assert original.count(old_import)==original.count(start)==original.count(end)==1
    changed=original.replace(old_import,old_import+',time')
    changed=changed.replace(start,'    _cost_hessian_start=time.perf_counter()\n'+start)
    changed=changed.replace(end,'    _cost_hessian_seconds=time.perf_counter()-_cost_hessian_start\n'+end)
    changed+='\nif __name__=="__main__":\n    (ROOT/"cost-hessian.json").write_text(json.dumps({"wall_seconds":_cost_hessian_seconds,"scope":"contact graph, Hessian assembly, eigendecomposition, rigid-mode selection and deterministic eigenvector signs; part of shared preparation"},indent=2)+"\\n")\n'
    assert hashlib.sha256(changed.encode()).hexdigest()==EXPECTED_INSTRUMENTED_PREPARE_SHA256
    return changed


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pilot',type=Path,required=True,help='Frozen pilot directory containing scripts, raw/, model/, results/ and protocol')
    parser.add_argument('--work',type=Path,required=True,help='New separate directory; must not already exist')
    parser.add_argument('--output',type=Path,help='New timing JSON; default is WORK/full-pilot-cost.json')
    args=parser.parse_args()
    pilot,work=args.pilot.resolve(),args.work.resolve()
    if work.exists():parser.error('The work directory must not exist; use a fresh location.')
    if pilot==work or pilot in work.parents:parser.error('The work directory must be outside the supplied pilot directory.')
    if sys.platform not in ('darwin','linux'):parser.error('RSS accounting is implemented for macOS and Linux only.')
    worker=Path(__file__).resolve().with_name('pilot-cost-worker.py')
    if not worker.is_file():parser.error('The sibling pilot-cost-worker.py is required.')
    required=['prepare.py','run_pilot.py','preanalysis-protocol.json','raw','model','results']
    if any(not (pilot/name).exists() for name in required):parser.error('The supplied pilot is incomplete.')
    for name,expected in EXPECTED_SOURCE_SHA256.items():
        if digest(pilot/name)!=expected:parser.error(f'{name} does not match the frozen measured version; no workload was launched.')
    # Inputs are copied, not linked, so timing instrumentation and generated
    # artifacts cannot alter the source pilot. Model/results are not copied:
    # they are regenerated and compared against read-only originals afterward.
    work.mkdir(parents=True)
    for name in ('prepare.py','run_pilot.py','preanalysis-protocol.json'):
        shutil.copy2(pilot/name,work/name)
    shutil.copytree(pilot/'raw',work/'raw')
    changed=instrument_prepare((work/'prepare.py').read_text())
    (work/'prepare.py').write_text(changed)
    shutil.copy2(worker,work/'measure_cost.py')
    provenance={
        'original_source_sha256':EXPECTED_SOURCE_SHA256,
        'executed_source_sha256':{name:digest(work/name) for name in EXPECTED_SOURCE_SHA256},
        'instrumentation':'Only prepare.py imports time, adds one monotonic timer around existing contact/Hessian/eigenbasis statements, and writes a separate cost-hessian.json after original work. No model parameters or scientific operations changed. run_pilot.py is unchanged.',
        'protocol_sha256':digest(work/'preanalysis-protocol.json'),
        'input_sha256':{str(path.relative_to(work)):digest(path) for path in sorted((work/'raw').iterdir()) if path.is_file()},
        'launcher_sha256':digest(__file__),
        'worker_sha256':digest(worker),
    }
    (work/'instrumentation-provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    destination=args.output.resolve() if args.output else work/'full-pilot-cost.json'
    if destination.exists() or destination.with_name('full-pilot-cost-public.json').exists():
        parser.error('Output records already exist; no workload was launched.')
    # Installation, cached-input copying and source checks finish before the
    # worker's measured interval. Each scientific phase then runs once in its
    # own fresh interpreter. No host name or supplied absolute path is stored.
    result=subprocess.run([sys.executable,str(work/'measure_cost.py'),'--reference',str(pilot),
                           '--output',str(destination)],cwd=work)
    if result.returncode:raise SystemExit(result.returncode)
    for name,expected in EXPECTED_SOURCE_SHA256.items():
        assert digest(pilot/name)==expected,'Original pilot source changed unexpectedly'


if __name__=='__main__':main()
