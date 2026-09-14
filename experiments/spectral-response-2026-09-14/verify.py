"""Independent saved-data confidence, circuit-word and complete cost audit."""
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[key]='1'
import hashlib,json,math,time
from pathlib import Path
import numpy as np
from scipy.stats import beta
from qiskit import qpy
from qiskit.quantum_info import Operator
ROOT=Path(__file__).resolve().parent

def load(path):
    with path.open('rb') as f:return qpy.load(f)[0]

def main():
    start=time.perf_counter();records=[]
    for name in ['triangle','kras']:
        folder=ROOT/'results'/name;native=json.loads((folder/'native-receipt.json').read_text());ft=json.loads((folder/'fault-tolerant/receipt.json').read_text());receipt=json.loads((folder/'IQAE-receipt.json').read_text())
        with np.load(folder/'IQAE-all-outputs.npz') as z:data={k:z[k] for k in z.files}
        with np.load(folder/'fault-tolerant/all-rotation-words.npz') as z:words={k:z[k] for k in z.files}
        b=ft['precision_bits'];bitcounts=[];maxsynth=0.
        for bit in range(b):
            qc=load(folder/'fault-tolerant'/f'Rz-bit-{bit}.qpy');ops=qc.count_ops();bitcounts.append(ops.get('t',0)+ops.get('tdg',0))
            angle=4*math.pi*2**bit/2**b;ideal=np.diag(np.exp(np.array([-1j,1j])*angle/2));error=float(np.linalg.norm(Operator(qc).data-ideal,2))
            assert error<=ft['per_bit_synthesis_epsilon']+2e-14;maxsynth=max(maxsynth,error)
        actualT={}
        for component in ft['components']:
            label=f"{component['kind']}_{component['index']}";qc=load(folder/'components'/f"{component['kind']}-{component['index']}.qpy")
            actualwords=np.array([int(round((float(i.operation.params[0])%(4*math.pi))*2**b/(4*math.pi)))%2**b for i in qc.data if i.operation.name=='ry'],dtype=np.uint64)
            assert np.array_equal(actualwords,words[label])
            t=sum(int(bitcounts[k])*int(np.count_nonzero((actualwords>>k)&1)) for k in range(b))
            assert t==component['T'];actualT[label]=t
        indices=data['queries'];number=len(indices);lo=np.zeros(number);hi=np.full(number,np.pi/2)
        aq=np.zeros(number,dtype=np.int64);gq=aq.copy();readouts=aq.copy();maxk=aq.copy();interval_error=0.
        for r in range(1,receipt['rounds']+1):
            get=lambda key:data[f'round_{r}_{key}'];ix=get('indices');m=get('multipliers');branch=get('branches');n=get('shots');x=get('successes')
            d=.05/number/(r*(r+1));assert float(get('delta'))==d
            assert np.all(m%2==1)
            assert np.all(m*lo[ix]>=branch*np.pi/2-1e-11)
            assert np.all(m*hi[ix]<=(branch+1)*np.pi/2+1e-11)
            pl=np.where(x==0,0,beta.ppf(d/2,x,n-x+1));pu=np.where(x==n,1,beta.ppf(1-d/2,x+1,n-x))
            al=np.arcsin(np.sqrt(pl));au=np.arcsin(np.sqrt(pu));even=branch%2==0
            l=np.where(even,branch*np.pi/2+al,(branch+1)*np.pi/2-au)/m
            u=np.where(even,branch*np.pi/2+au,(branch+1)*np.pi/2-al)/m
            lo[ix]=np.maximum(lo[ix],l);hi[ix]=np.minimum(hi[ix],u)
            interval_error=max(interval_error,float(np.max(abs(lo[ix]-get('theta_lower')))),float(np.max(abs(hi[ix]-get('theta_upper')))))
            aq[ix]+=n*m;gq[ix]+=n*((m-1)//2);readouts[ix]+=n;maxk[ix]=np.maximum(maxk[ix],(m-1)//2)
        for key,v in [('A_queries',aq),('Q_queries',gq),('readouts',readouts),('maximum_k',maxk)]:assert np.array_equal(v,data[key])
        assert np.max(abs(np.sin(lo)**2-data['a_lower']))<1e-14
        assert np.max(abs(np.sin(hi)**2-data['a_upper']))<1e-14
        AT=np.array([actualT[f'prepare_{i}']+actualT[f'prepare_{j}']+actualT[f'filter_{t}'] for t,i,j in indices],dtype=np.int64)
        assert np.array_equal(AT,data['A_T'])
        CX=np.array([native['preparations'][i]['CX']+native['preparations'][j]['CX']+native['filters'][t]['CX'] for t,i,j in indices],dtype=np.int64)
        ref=load(folder/'zero-reflection.qpy').count_ops();rt=ref.get('t',0)+ref.get('tdg',0);rc=ref.get('cx',0)
        totalT=sum(int(a)*int(t)+int(g)*rt for a,t,g in zip(aq,AT,gq));totalCX=sum(int(a)*int(c)+int(g)*rc for a,c,g in zip(aq,CX,gq))
        assert totalT==receipt['estimated_complete_T'];assert totalCX==receipt['estimated_complete_CX']
        shots=math.ceil(math.log(2*number/.05)/(2*receipt['amplitude_epsilon']**2))
        assert shots==receipt['matched_spectral_sampling']['shots_per_output']
        assert shots*sum(map(int,AT))==receipt['matched_spectral_sampling']['complete_T']
        assert interval_error<1e-14
        records.append(dict(model=name,status='PASS',outputs=number,reconstructed_total_T=totalT,reconstructed_total_CX=totalCX,
                            all_raw_schedules_counts_and_intervals_reconstructed=True,interval_maximum_difference=interval_error,
                            all_component_angle_words_and_T_counts_verified=True,maximum_emitted_bit_operator_error=maxsynth,
                            full_raw_sha256=hashlib.sha256((folder/'IQAE-all-outputs.npz').read_bytes()).hexdigest()))
    result=dict(status='PASS_INDEPENDENT_SAVED_EVIDENCE_CHECKS',records=records,seconds=time.perf_counter()-start,
                scope='Independent reconstruction of saved current canonical replay. The original pre-cache raw KRAS archive was overwritten; this verifier does not claim a byte comparison to it.')
    (ROOT/'independent-verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
