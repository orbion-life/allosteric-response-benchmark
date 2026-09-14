"""Constructive complete binary-rotation T accounting, with bounded exports."""
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[key]='1'
import argparse,hashlib,json,math,resource,sys,time
from pathlib import Path
import numpy as np
from qiskit import QuantumCircuit,qpy
from qiskit.quantum_info import Operator
ROOT=Path(__file__).resolve().parent
if os.environ.get('PULSAR_SYNTH_RUNTIME'):
    sys.path.insert(0,os.environ['PULSAR_SYNTH_RUNTIME'])
sys.path.insert(0,str(ROOT/'vendor'))
from prior_synthesis import Synthesizer,counts as detailed_counts
from spectral import grover


def load(path):
    with path.open('rb') as f:return qpy.load(f)[0]
def save(qc,path):
    with path.open('wb') as f:qpy.dump(qc,f)
def word(angle,b):return int(round((float(angle)%(4*math.pi))*2**b/(4*math.pi)))%(2**b)


def run(name):
    tick=time.perf_counter();out=ROOT/'results'/name/'fault-tolerant';out.mkdir(exist_ok=False)
    native=json.loads((out.parent/'native-receipt.json').read_text());scale=native['scale'];N=native['residues']
    nrot=2*max(x['rotations'] for x in native['preparations'])+max(x['rotations'] for x in native['filters'])
    b=math.ceil(math.log2(4*scale*nrot*math.pi/1e-6));epsilon=1e-6/(4*scale*nrot*b)
    synth=Synthesizer(epsilon,out/'synthesis');bank=[synth.rz(4*math.pi*2**j/2**b) for j in range(b)]
    bc=[detailed_counts(p) for p in bank];mats=[Operator(p).data for p in bank]
    for i,p in enumerate(bank):save(p,out/f'Rz-bit-{i}.qpy')
    def replace(qc,emit=False):
        result=QuantumCircuit(qc.num_qubits) if emit else None
        if emit:result.global_phase=qc.global_phase
        totals={key:0 for key in ['T','CX','operations']};words=[]
        for inst in qc.data:
            gate=inst.operation.name;bits=[qc.find_bit(q).index for q in inst.qubits]
            if gate=='ry':
                value=word(inst.operation.params[0],b);words.append(value)
                totals['operations']+=4
                if emit:result.sdg(bits[0]);result.h(bits[0])
                for k in range(b):
                    if value>>k&1:
                        for key in totals:totals[key]+=bc[k][key]
                        if emit:result.compose(bank[k],[bits[0]],inplace=True)
                if emit:result.h(bits[0]);result.s(bits[0])
            else:
                totals['operations']+=1;totals['CX']+=int(gate=='cx');totals['T']+=int(gate in ['t','tdg'])
                if emit:result.append(inst.operation,bits)
        return totals,np.array(words,dtype=np.uint64),result
    rows=[];raw={}
    for kind,count in [('prepare',N),('filter',3)]:
        for i in range(count):
            qc=load(out.parent/'components'/f'{kind}-{i}.qpy');cost,words,_=replace(qc)
            rows.append(dict(kind=kind,index=i,**cost,rotation_count=len(words)))
            raw[f'{kind}_{i}']=words
    np.savez_compressed(out/'all-rotation-words.npz',**raw)
    A=load(out.parent/'A-case1.qpy');cost,words,_=replace(A);exports=[]
    if cost['operations']<=2000000:
        _,_,compiled=replace(A,True);assert detailed_counts(compiled)['T']==cost['T'];save(compiled,out/'A-case1.qpy')
        Q,S0=grover(compiled);qc=detailed_counts(Q)
        if qc['operations']<=2000000:save(Q,out/'Q-case1.qpy');qstatus='EXPORTED'
        else:qstatus='OPERATION_CAP_NO_EXPORT'
        exports.append(dict(name='A-case1',status='EXPORTED',**detailed_counts(compiled)))
        exports.append(dict(name='Q-case1',status=qstatus,**qc))
    else:exports.append(dict(name='A-case1',status='OPERATION_CAP_NO_CONSTRUCTION',**cost))
    result=dict(status='COMPLETE_CONSTRUCTIVE_BINARY_SYNTHESIS_ACCOUNTING',source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                precision_bits=b,per_bit_synthesis_epsilon=epsilon,maximum_A_rotations=nrot,
                word_rounding_response_bound=4*scale*nrot*math.pi/2**b,synthesis_response_bound=4*scale*nrot*b*epsilon,
                bank_counts=bc,components=rows,exports=exports,
                scope='Actual synthesized fixed binary Rz gate strings and exact bit-word composition counts for every protein preparation/filter. Counts are constructive before optimization, not minimal T estimates. Full A and Q emitted only below operation cap. The Grover inverse reuses exact emitted inverse sequence, so fixed approximate-A bias is bounded once, rather than treated as independent amplification error.',
                wall_seconds=time.perf_counter()-tick,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ['bank_counts','components']},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('model',choices=['triangle','kras']);run(p.parse_args().model)
