#!/usr/bin/env python3
"""Independently check the saved grid topology and replay persisted native circuits."""
from pathlib import Path
import json,hashlib,sys
import numpy as np
from qiskit import qpy
from qiskit.quantum_info import Statevector
ROOT=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    out=ROOT/'results';r=json.loads((out/'result.json').read_text());f=np.load(out/'fixture.npz')
    assert r['status']=='passed'
    assert sha(ROOT/'protocol.json')==r['protocol_sha256']
    for name,h in r['source_sha256'].items():assert sha(ROOT/name)==h
    assert sha(out/'fixture.npz')==r['fixture_sha256']
    P=f['P']
    for a in range(16):
        for b in range(16):
            permitted=(a==b) or (a<9 and b<9 and abs(a//3-b//3)+abs(a%3-b%3)==1)
            if not permitted:assert P[a,b]==0,(a,b,'non-neighbour transition')
    assert P[2,3]==0 and P[5,6]==0,'row-wrap transitions are forbidden'
    assert np.array_equal(P[9:,9:],np.eye(7))
    assert np.max(abs(P.sum(1)-1))<1e-10
    expected=float(f['left']@f['polynomial']@f['right']);records=[]
    for row in r['compiled']:
        path=out/'circuits'/f"{row['method']}_{row['topology']}.qpy"
        assert sha(path)==row['qpy_sha256']
        with path.open('rb') as h:circuit=qpy.load(h)[0]
        assert circuit.num_qubits==(11 if row['method']=='two_address' else 10)
        assert set(circuit.count_ops())<=set(['u','cx','measure','barrier'])
        assert circuit.count_ops().get('measure')==1
        probability=Statevector.from_instruction(circuit.remove_final_measurements(inplace=False)).probabilities([0])
        actual=float(probability[0]-probability[1]);target=(-1 if row['method']=='direction_negative' else 1)*expected
        error=abs(actual-target);assert error<1e-10,(path.name,error)
        records.append({'path':str(path.relative_to(ROOT)),'absolute_error':error,'qpy_sha256':sha(path)})
    result={'status':'passed','protocol_sha256':r['protocol_sha256'],'results_sha256':sha(out/'result.json'),'independent_topology_and_row_wrap_check':True,'saved_circuits_replayed':records}
    (out/'independent-replay.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
