from pathlib import Path
import json,sys
ROOT=Path(__file__).resolve().parent
import openqasm3
from openqasm3.parser import parse
rows=[]
for f in sorted((ROOT/'results/circuits').glob('*.qasm')):
    program=parse(f.read_text());rows.append(dict(file=f.name,parsed=True,statements=len(program.statements)))
r=dict(parser='openqasm3',version=openqasm3.__version__,parsed=len(rows),total=len(rows),live_device_acceptance=False,scope='OpenQASM3 grammar parse only; provider gate validation remains separate',files=rows)
(ROOT/'results/qasm-parser-validation.json').write_text(json.dumps(r,indent=2)+'\n')
print(json.dumps({k:v for k,v in r.items() if k!='files'},indent=2))
