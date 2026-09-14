"""Portable path-only wrappers around the unchanged independent audit code."""
from pathlib import Path
import argparse, importlib.util, hashlib, sys
from threadpoolctl import threadpool_limits

ROOT=Path(__file__).resolve().parents[1]
def run(name,folder='operator-independent-audit'):
    path=ROOT/'review'/folder/name
    sys.path.insert(0,str(path.parent))
    spec=importlib.util.spec_from_file_location(path.stem,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    # The archived auditors originally sat beside an `operator` directory.
    # Only location variables change; every mathematical check remains intact.
    module.BASE=ROOT;module.OP=ROOT;module.HERE=ROOT/'checks/replay';module.HERE.mkdir(exist_ok=True)
    with threadpool_limits(limits=2):module.main()

def main():
    p=argparse.ArgumentParser();p.add_argument('--quantum',action='store_true');p.add_argument('--abl',action='store_true');args=p.parse_args()
    manifest=ROOT/'MANIFEST.sha256'
    if manifest.exists():
        for line in manifest.read_text().splitlines():
            digest,name=line.split('  ',1)
            assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
        print('Export manifest verified')
    run('verify_scope_and_costs.py')
    if args.quantum:run('verify_operator_circuit.py')
    if args.abl:run('verify_operator_results.py','independent-math-audit')

if __name__=='__main__':main()
