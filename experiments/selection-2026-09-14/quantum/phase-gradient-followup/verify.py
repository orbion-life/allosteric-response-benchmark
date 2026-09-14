"""Independent gradient-preparation overlap, saved components and total-cost verification."""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[k]='1'
import hashlib,json,math,time
from pathlib import Path
import mpmath as mp
import numpy as np
from qiskit import qpy
from scipy.special import ive
ROOT=Path(__file__).resolve().parent;PARENT=ROOT.parent

def load(p):
    with Path(p).open('rb') as f:return qpy.load(f)[0]
def count(p):
    q=load(p);o=q.count_ops();return dict(T=o.get('t',0)+o.get('tdg',0),CX=o.get('cx',0),operations=sum(o.values()),qubits=q.num_qubits)

def main():
    start=time.perf_counter();mp.mp.dps=70;r=json.loads((ROOT/'results/receipt.json').read_text());b=r['precision_bits'];q=load(ROOT/'results/gradient-preparation-b33.qpy')
    # Follow actual emitted circuit gates, independent of generator's gate-string reader.
    v=[[mp.mpc(1),mp.mpc(0)] for _ in range(b)];w=mp.exp(mp.j*mp.pi/4);i=mp.j;sqrt2=mp.sqrt(2)
    for inst in q.data:
        bit=q.find_bit(inst.qubits[0]).index;a,c=v[bit];g=inst.operation.name
        if g=='h':a,c=(a+c)/sqrt2,(a-c)/sqrt2
        elif g=='t':c*=w
        elif g=='tdg':c/=w
        elif g=='s':c*=i
        elif g=='sdg':c/=i
        elif g=='x':a,c=c,a
        elif g=='z':c=-c
        else:raise ValueError(g)
        v[bit]=[a,c]
    overlap=mp.mpc(1)
    for bit,(a,c) in enumerate(v):overlap*=(a+mp.exp(2*mp.j*mp.pi*2**bit/2**b)*c)/sqrt2
    error=mp.sqrt(max(mp.mpf(0),2-2*abs(overlap)));difference=abs(error-mp.mpf(r['gradient_state_phase_aligned_distance_70_digits']));assert difference<mp.mpf('1e-50')
    wc=count(ROOT/'results/controlled-walk.qpy');gc=count(ROOT/'results/gradient-preparation-b33.qpy');pc=[count(ROOT/'results'/f'prepare-{k}.qpy') for k in range(3)]
    for key,val in wc.items():assert val==r['cost']['controlled_walk'][key]
    for key,val in gc.items():assert val==r['gradient_preparation'][key]
    z=np.load(ROOT/'results/mathematical-results.npz');B=z['B'];scale=np.outer(np.linalg.norm(B,axis=0),np.linalg.norm(B,axis=0));C=[];expectedT=0.
    for ti,t in enumerate(z['times']):
        k=np.arange(r['polynomial_degrees'][ti]+1);c=ive(k,r['alpha']*t);c[1:]*=2;C.append(scale*np.einsum('k,kij->ij',c,z['terms'][:len(c)]).real-B.T@B)
        ek=float(np.dot(k,c/c.sum()))
        for a in range(3):
            for bidx in range(a,3):expectedT+=r['sampling']['shots_per_entry']*(gc['T']+pc[a]['T']+pc[bidx]['T']+ek*wc['T'])
    assert np.max(abs(np.array(C)-z['ideal_gradient_C']))<1e-12
    assert abs(expectedT/r['cost']['expected_total_T']-1)<1e-12
    for record in r['dense_fairness']:
        path=ROOT/'results'/(record['name']+'-constant-angle-bank.qpy');cc=count(path)
        for key in cc:assert cc[key]==record['compiled_constant_bank'][key]
        choices=json.loads((ROOT/'results'/(record['name']+'-rotation-choice.json')).read_text())
        assert len(choices)==record['angles_compared'] and all(x['direct_T']<=x['implemented_constant_adder_T'] for x in choices)
    baseline=min(x['expected_total_T_retained'] for x in r['dense_fairness']);ratio=expectedT/baseline
    result=dict(status='PASS',decision='PASS_SCOPED_FOLLOWUP_HALF_T_GATE',gradient_state_distance_from_emitted_gates_70_digits=str(error),gradient_distance_reconstruction_error=str(difference),
        full_expected_T_reconstructed=expectedT,candidate_to_strongest_evaluated_dense_ratio=ratio,half_gate=ratio<=.5,
        response_error_upper=r['prepared_gradient_Gaussian_response_error_upper'],full_statevector_executed=False,maximum_order_complete_circuit_executed=False,
        wall_seconds=time.perf_counter()-start,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (ROOT/'independent-verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
