"""Audit saved results via independent native-gate QASM replay and statistical reconstruction."""
from pathlib import Path
import re,json,hashlib,math
import numpy as np
from scipy.stats import beta
ROOT=Path(__file__).resolve().parent
out=ROOT/'results';cases=json.loads((out/'circuit_results.json').read_text());summary=json.loads((out/'summary.json').read_text())
physical=[7,8,9,10,13,14];pos={p:i for i,p in enumerate(physical)};edges={frozenset(e) for e in [(7,8),(8,9),(9,10),(8,13),(9,14)]}
checks=[]
for case in cases:
    s=(out/case['qasm_file']).read_text();state=np.zeros(64,complex);state[0]=1;nr=nc=0
    for line in s.splitlines():
        m=re.fullmatch(r'  prx\(([^,]+), ([^)]+)\) \$(\d+);',line)
        if m:
            theta,phi=float(m[1]),float(m[2]);q=pos[int(m[3])];nr+=1
            c=np.cos(theta/2);v=np.sin(theta/2);U=np.array([[c,-1j*np.exp(-1j*phi)*v],[-1j*np.exp(1j*phi)*v,c]])
            for a in range(64):
                if a&(1<<q):continue
                b=a|(1<<q);state[[a,b]]=U@state[[a,b]]
            continue
        m=re.fullmatch(r'  cz \$(\d+), \$(\d+);',line)
        if m:
            p,q=int(m[1]),int(m[2]);assert frozenset((p,q)) in edges;nc+=1
            for a in range(64):
                if (a&(1<<pos[p])) and (a&(1<<pos[q])):state[a]*=-1
    q=case['active_readout_local'];p=float(sum(abs(state[a])**2 for a in range(64) if not a&(1<<q)))
    assert abs(p-case['probability'])<2e-13
    assert nr==case['compiled']['ops']['r'] and nc==case['compiled']['ops']['cz']
    checks.append(dict(label=case['label'],qasm_probability=p,qpy_probability=case['probability'],difference=abs(p-case['probability'])))
for batch in json.loads((out/'sampling_simulations.json').read_text()):
    for r in batch['queries']:
        k,n=r['zeros'],r['shots'];a=.05/18
        lo=0. if k==0 else float(beta.ppf(a/2,k,n-k+1));hi=1. if k==n else float(beta.ppf(1-a/2,k+1,n-k))
        assert max(abs(lo-r['probability_interval'][0]),abs(hi-r['probability_interval'][1]))<1e-14
for b in summary['shot_budgets']:
    assert b['total_shots']==18*b['shots_per_data_query']+4*8192+3*16384
    assert b['tasks']==18*math.ceil(b['shots_per_data_query']/20000)+7
    assert abs(b['price_USD_at_verified_tariff']-(.00145*b['total_shots']+.3*b['tasks']))<1e-8
assert summary['maximum_ideal_response_error']<1e-6
assert len(list((out/'circuits').glob('*.qasm')))==25
result=dict(status='PASS',native_qasm_independent_replays=len(checks),maximum_native_qasm_probability_difference=max(c['difference'] for c in checks),confidence_intervals_reconstructed=216,input_sha256=hashlib.sha256((ROOT/'triangle.npz').read_bytes()).hexdigest(),noise_channels_are_assumptions=True,hardware_execution=False,checks=checks)
(out/'verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='checks'},indent=2))
