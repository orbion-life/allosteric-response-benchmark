"""Actual protein tables and synthesized rotation banks; other counts are construction estimates."""
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[key]='1'
import datetime,hashlib,json,math,resource,time
from pathlib import Path
import numpy as np
from qiskit import QuantumCircuit
from select_swap import tables_for_values,quantize,choose_lambda
from synthesis import Synthesizer,counts,save
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'results/proteins'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()

def table_cost(theta,b,bank):
    words=quantize(theta,b);N=len(words);lam=choose_lambda(N,b);M=N//lam
    pop=int(np.bitwise_count(words).sum());rccx=4*max(M-2,0);fredkin=b*(lam-1)
    t=rccx*4+7*fredkin;lookup_ops=9*rccx+17*fredkin+2*max(M-1,0)+pop
    lookup_cx=3*rccx+8*fredkin+(pop if M>1 else 0)
    return words,dict(N=N,lambda_value=lam,word_bits=b,word_ones=pop,
        clean_workspace=lam*b+max(0,(M-1).bit_length()-1),dirty_echo_workspace=lam*b,
        literal_lookup_T=t,literal_lookup_operations=lookup_ops,literal_lookup_CX=lookup_cx,
        forward_rotation_uncompute_T_upper=2*t+bank['T'],
        forward_rotation_uncompute_CX_upper=2*lookup_cx+bank['CX'],
        forward_rotation_uncompute_operations_upper=2*lookup_ops+bank['operations'],
        dirty_echo_lookup_T_upper=2*t+14*b*(lam-1))

def main():
    OUT.mkdir(exist_ok=True);protocol=json.loads((ROOT/'protein-protocol.json').read_text());assert sha(Path(__file__))==protocol['source_sha256']
    start=time.perf_counter();results=[]
    for target in ['kras','abl','myc','myh7']:
        folder=OUT/target;folder.mkdir(exist_ok=True)
        arrays=ROOT/'inputs'/(target+'-measured-arrays.npz');rp=ROOT/'inputs'/(target+'-receipt.json');z=np.load(arrays,allow_pickle=False);prior=json.loads(rp.read_text())
        L,B=z['L'],z['B'];rank=len(L);Nres=B.shape[1];n=(rank-1).bit_length();d=2**n;norm=np.linalg.norm(B,axis=0);scale=float(max(norm)**2)
        query=next(x['queries'] for x in prior['polynomial_degrees'] if x['encoding']=='factor');K=4096 if any('degree' not in q for q in query) else max(q['degree'] for q in query)
        paths=4*n*K+2*n;b=math.ceil(math.log2(2*math.pi*scale*paths/1e-6));epsilon=1e-6/(2*scale*2*b*paths)
        synth=Synthesizer(epsilon,folder/'synthesis');bank=QuantumCircuit(b+1)
        for bit in range(b):
            a=synth.ry(2*math.pi*2**bit/2**b);bank.compose(a,[b],inplace=True);bank.cx(bit,b);bank.compose(a.inverse(),[b],inplace=True);bank.cx(bit,b)
        bc=counts(bank);save(bank,folder/'rotation-bank.qpy')
        pL=np.zeros((d,d));pL[:rank,:rank]=L;rho=np.linalg.norm(pL,axis=1)
        pin=tables_for_values((rho/np.linalg.norm(rho))[None,:],list(range(n+1,2*n+1)),[],0)
        pout=tables_for_values(pL,list(range(1,n+1)),list(range(n+1,2*n+1)),0)
        tabledata={};tables=[]
        for k,table in enumerate(pin+pout):
            words,cost=table_cost(table['theta'],b,bc);tabledata['words_'+str(k)]=words;tabledata['angles_'+str(k)]=table['theta'];tables.append(cost)
        np.savez_compressed(folder/'factor-angle-tables.npz',**tabledata)
        prepdata={};prepcost=[];pB=np.zeros((d,Nres));pB[:rank]=B/norm
        for i in range(Nres):
            totalT=0;totalCX=0;totalops=0;space=0
            for j,table in enumerate(tables_for_values(pB[:,i][None,:],list(range(1,n+1)),[],0)):
                words,c=table_cost(table['theta'],b,bc);prepdata[f'words_{i}_{j}']=words
                totalT+=c['forward_rotation_uncompute_T_upper'];totalCX+=c['forward_rotation_uncompute_CX_upper'];totalops+=c['forward_rotation_uncompute_operations_upper'];space=max(space,c['clean_workspace'])
            prepcost.append(dict(index=i,T_upper=totalT,CX_upper=totalCX,operations_upper=totalops,clean_workspace=space))
        np.savez_compressed(folder/'observable-angle-words.npz',**prepdata)
        # The exact controlled reflections use paired RCCX chains and one CZ each.
        walkT=2*sum(x['forward_rotation_uncompute_T_upper'] for x in tables)+16*n
        walkCX=2*sum(x['forward_rotation_uncompute_CX_upper'] for x in tables)+12*n+2
        walkops=2*sum(x['forward_rotation_uncompute_operations_upper'] for x in tables)+40*n+11
        projection=prior['reference_checks']['projection_to_exact_Gaussian'];eps=.002-projection-2e-6-2e-8
        entries=3*Nres*(Nres+1)//2;shots=math.ceil(2*scale**2*math.log(2*entries/.05)/eps**2)
        queries=[];all_expected=0;complete=True
        for q in query:
            k=q.get('degree');ek=q.get('expected_order');points=Nres*(Nres+1)//2
            expected=None if ek is None else shots*(points*ek*walkT+(Nres+1)*sum(p['T_upper'] for p in prepcost))
            queries.append(dict(time_index=q['time_index'],degree=k,status=q.get('status','ADMITTED_DEGREE'),
                expected_order=ek,expected_complete_full_matrix_T_upper=expected,
                maximum_order_pair01_native_operations_upper=None if k is None else k*walkops+prepcost[0]['operations_upper']+prepcost[1]['operations_upper']+3))
            if expected is None:complete=False
            else:all_expected+=expected
        result=dict(target=target,status='ACTUAL_TABLES_AND_ROTATION_BANK_SYNTHESIS_OTHER_COMPONENTS_UNCOMPILED',source_input_sha256=sha(arrays),source_receipt_sha256=sha(rp),
            rank=rank,residues=Nres,padded_rank=d,normalization_Gamma=prior['factor_Gamma'],normalization_penalty=prior['normalization_ratio'],precision_bits=b,
            precision_maximum_order=K,precision_scope='Up to the4096cap only where a degree was not admitted; a precision setting does not admit the missing polynomial.',
            maximum_word_rounding_response_bound=2*scale*paths*math.pi/2**b,per_Rz_epsilon=epsilon,
            logical_system_qubits=2*n+2,clean_workspace_qubits=max(x['clean_workspace'] for x in tables),dirty_qubits_used=0,
            dirty_echo_note='An alternative verified echo would restorelambda*bdirtyqubits but adds lookup/swap passes andbclean output qubits. These are calculated costs, not a measured protein circuit or free workspace.',
            rotation_bank=bc,lookup_tables=tables,preparations=prepcost,
            controlled_walk_T_upper=walkT,controlled_walk_CX_upper=walkCX,controlled_walk_native_operations_upper=walkops,
            queries=queries,full_matrix_shots_sufficient=entries*shots,shots_per_entry=shots,sampling_epsilon=eps,
            expected_all_three_time_T_upper=all_expected if complete else None,
            count_scope='Literal chosen RCCX/Fredkin/bit-rotation construction before compiler cancellations; actual data bits and allprepare/unprepare/lookup/uncompute/reflections included. Complete single-circuit serial depth is bounded above by its operation count. No native protein lookup, full walk, routing or noise simulation; not a minimal bound or an admitted execution advantage.',
            factor_tables_sha256=sha(folder/'factor-angle-tables.npz'),observable_words_sha256=sha(folder/'observable-angle-words.npz'))
        (folder/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');results.append(result)
        print(target,b,walkT,entries*shots,flush=True)
    summary=dict(status='COMPLETE_SCOPED_RESOURCE_ACCOUNTING',source_sha256=sha(Path(__file__)),protocol_sha256=sha(ROOT/'protein-protocol.json'),
                 wall_seconds=time.perf_counter()-start,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,results=results)
    (OUT/'receipt.json').write_text(json.dumps(summary,indent=2)+'\n')

if __name__=='__main__':main()
