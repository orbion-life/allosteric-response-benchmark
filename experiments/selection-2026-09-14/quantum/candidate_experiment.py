"""Complete fixed SELECT–SWAP component synthesis and same-model estimator costs."""
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[key]='1'
import argparse,datetime,hashlib,json,math,resource,sys,time
from pathlib import Path
import numpy as np
from scipy.linalg import cholesky,eigh
from qiskit import QuantumCircuit,qpy
from qiskit.circuit.library import RCCXGate
from qiskit.quantum_info import Operator,Statevector
from select_swap import lookup,quantize,tables_for_values,choose_lambda
from synthesis import Synthesizer,exact_clifford_t,counts,save
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'vendor'))
from reduced_overlap import coefficients
CAP=1_000_000

def dump(p,x):Path(p).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def mark(stage):dump(Path(os.environ.get('PULSAR_PROGRESS',str(ROOT/'progress.json'))),dict(stage=stage,monotonic=time.monotonic(),utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))

def data():
    z=np.load(ROOT/'inputs/triangle.npz',allow_pickle=False);H=z['Hr'];B=z['B'];ts=z['times']
    w,V=eigh(H);QB=V.T@B;exact=np.array([(QB.T*np.expm1(-t*w))@QB for t in ts]);norm=np.linalg.norm(B,axis=0)
    return z,H,B,ts,w,exact,norm

def sampling(norm,projection):
    scale=float(max(norm)**2);eps=.002-projection-1e-8-1e-6-1e-6-1e-8
    return dict(sampling_epsilon=eps,shots_per_entry=math.ceil(2*scale**2*math.log(36/.05)/eps**2),entries=18,
                scope='Sufficient ideal sampling plan; no actual shots or noise guarantee',family_failure_probability=.05)

def summarize_costs(walk,preps,degrees,poly,shots):
    wc=counts(walk);pc=[counts(p) for p in preps];rows=[]
    for ti,c in enumerate(poly):
        ek=float(np.dot(np.arange(len(c)),c/c.sum()))
        for i in range(3):
            for j in range(i,3):
                rows.append(dict(time_index=ti,i=i,j=j,expected_order=ek,shots=shots,
                    expected_T_per_shot=pc[i]['T']+pc[j]['T']+ek*wc['T'],
                    expected_CX_per_shot=pc[i]['CX']+pc[j]['CX']+ek*wc['CX'],
                    serial_T_depth_upper_per_shot=pc[i]['T_depth']+pc[j]['T_depth']+ek*wc['T_depth']))
    return dict(controlled_walk=wc,preparations=pc,queries=rows,
        expected_total_T=sum(x['shots']*x['expected_T_per_shot'] for x in rows),
        expected_total_CX=sum(x['shots']*x['expected_CX_per_shot'] for x in rows),
        all_query_readouts=18*shots,maximum_order=max(degrees),
        maximum_order_pair01_cached_operations=pc[0]['operations']+pc[1]['operations']+max(degrees)*wc['operations']+3,
        cost_scope='Actual synthesized component counts, exact cached composition and coefficient-weighted repetition. Complete high-order circuits may exceed cap. Two Hadamards/readout charged; T=0. No routing, error correction, distillation or hardware runtime estimate.')

def basic_complete(walk,preps,out):
    info=[]
    for k in [0,1]:
        count=counts(preps[0])['operations']+counts(preps[1])['operations']+k*counts(walk)['operations']+3
        if count>CAP:info.append(dict(order=k,status='OPERATION_CAP',cached_operations=count));continue
        qc=QuantumCircuit(walk.num_qubits,1);qc.h(0);qc.compose(preps[1],inplace=True)
        for _ in range(k):qc.compose(walk,inplace=True)
        qc.compose(preps[0].inverse(),inplace=True);qc.h(0);qc.measure(0,0)
        save(qc,out/f'complete-k{k}.qpy');info.append(dict(order=k,status='ACTUAL_CACHED_CIRCUIT_EXPORTED',**counts(qc)))
    return info

def reflection(systembits,total,work,input_reflection=False):
    qc=QuantumCircuit(total);qc.x(systembits[1:]) # readout is a positive control
    controls=systembits[:-1];target=systembits[-1]
    if len(controls)==1:qc.cz(controls[0],target)
    else:
        qc.append(RCCXGate(),[controls[0],controls[1],work[0]])
        for k in range(2,len(controls)):qc.append(RCCXGate(),[work[k-2],controls[k],work[k-1]])
        qc.cz(work[len(controls)-2],target)
        for k in range(len(controls)-1,1,-1):qc.append(RCCXGate().inverse(),[work[k-2],controls[k],work[k-1]])
        qc.append(RCCXGate().inverse(),[controls[0],controls[1],work[0]])
    qc.x(systembits[1:])
    if input_reflection:qc.z(0)
    return exact_clifford_t(qc)

def apply_table(U,table,bank,b,nbits):
    """Exact effective clean-workspace action of verified lookup/rotation/uncompute."""
    ctrl=[x-1 for x in table['controls'] if x!=0];target=table['target']-1
    theta=table['theta'][len(table['theta'])//2:];words=quantize(theta,b)
    base=np.arange(2**nbits);lo=base[(base>>target&1)==0];hi=lo+(1<<target)
    addresses=sum((((lo>>q)&1)<<i for i,q in enumerate(ctrl)),np.zeros(len(lo),dtype=int));matrices=np.empty((len(words),2,2),complex)
    X=np.array([[0,1],[1,0]],complex)
    for j,word in enumerate(words):
        mat=np.eye(2,dtype=complex)
        for bit in range(b):
            if int(word)>>bit&1:
                A=bank[bit];mat=(X@A.conj().T@X@A)@mat
        matrices[j]=mat
    a=U[lo].copy();v=U[hi].copy();g=matrices[addresses]
    U[lo]=g[:,0,0,None]*a+g[:,0,1,None]*v;U[hi]=g[:,1,0,None]*a+g[:,1,1,None]*v
    return U

def candidate():
    start=time.perf_counter();out=ROOT/'results/select_swap';out.mkdir(exist_ok=False)
    z,H,B,times,w,exact,norm=data();L=cholesky(H,lower=False);gamma=float(np.linalg.norm(L)**2);alpha=gamma/2
    poly=[coefficients(float(t*alpha),1e-8/max(1,float(max(norm)**2)))[0] for t in times];degrees=[len(c)-1 for c in poly];K=max(degrees)
    n=4;d=16;system=10;scale=float(max(norm)**2);paths=4*n*K+2*n
    b=math.ceil(math.log2(2*math.pi*scale*paths/1e-6));roundbound=2*scale*paths*math.pi/2**b
    epsrot=1e-6/(2*scale*2*b*paths)
    synth=Synthesizer(epsrot,out/'synthesis');mark('candidate:binary_rotation_synthesis')
    bank=[synth.ry(2*math.pi*2**bit/2**b) for bit in range(b)]
    bankmat=[Operator(a).data for a in bank]
    rotbank=QuantumCircuit(b+1)
    for bit,A in enumerate(bank):
        rotbank.compose(A,[b],inplace=True);rotbank.cx(bit,b);rotbank.compose(A.inverse(),[b],inplace=True);rotbank.cx(bit,b)
    save(rotbank,out/'bit-controlled-rotation-bank.qpy')
    padded=np.zeros((d,d));padded[:len(L),:len(L)]=L;rho=np.linalg.norm(padded,axis=1)
    pin=tables_for_values((rho/np.linalg.norm(rho))[None,:],list(range(5,9)),[],0)
    pout=tables_for_values(padded,list(range(1,5)),list(range(5,9)),0)
    vec=np.zeros((d,3));vec[:len(B)]=B/norm
    observable=[tables_for_values(vec[:,i][None,:],list(range(1,5)),[],0) for i in range(3)]
    all_tables=pin+pout+sum(observable,[])
    workspace=max(choose_lambda(len(t['theta']),b)*b+max(0,(len(t['theta'])-1).bit_length()-(choose_lambda(len(t['theta']),b)-1).bit_length()-1) for t in all_tables)
    workspace=max(workspace,4);total=system+workspace;table_receipts=[]
    def compile_tables(tables,label):
        result=QuantumCircuit(total)
        for k,table in enumerate(tables):
            mark('candidate:compile_lookup:'+label+str(k));words=quantize(table['theta'],b);oracle,meta=lookup(words,b)
            native=exact_clifford_t(oracle);mapping=table['controls']+list(range(system,system+oracle.num_qubits-len(table['controls'])))
            result.compose(native,mapping,inplace=True)
            data_bits=list(range(system,system+b))
            result.compose(rotbank,data_bits+[table['target']],inplace=True)
            result.compose(native.inverse(),mapping,inplace=True)
            save(native,out/f'lookup-{label}-{k}.qpy')
            np.savez_compressed(out/f'lookup-{label}-{k}-data.npz',theta=table['theta'],words=words)
            table_receipts.append(dict(label=label,index=k,controls=table['controls'],target=table['target'],metadata=meta,actual_native=counts(native),
                                       forward_and_uncompute_T=2*counts(native)['T'],rotation_bank_T=counts(rotbank)['T']))
        return result
    inp=compile_tables(pin,'row');otp=compile_tables(pout,'column')
    ul=inp.copy();ul.compose(otp.inverse(),inplace=True)
    # Retain a flag wire. Its attenuation is identicallyzero for Gamma=||L||F².
    rout=reflection([0]+list(range(1,5))+[9],total,list(range(system,total)))
    rin=reflection([0]+list(range(5,9))+[9],total,list(range(system,total)),True)
    walk=ul.copy();walk.compose(rout,inplace=True);walk.compose(ul.inverse(),inplace=True);walk.compose(rin,inplace=True)
    preps=[compile_tables(t,'observable'+str(i)) for i,t in enumerate(observable)]
    mark('candidate:effective_clean_workspace_matrix')
    Uin=np.eye(512,dtype=complex);Uout=np.eye(512,dtype=complex)
    for table in pin:apply_table(Uin,table,bankmat,b,9)
    for table in pout:apply_table(Uout,table,bankmat,b,9)
    UL=Uout.conj().T@Uin
    p_in=np.arange(16);p_out=np.arange(16)*16
    routdiag=np.ones(512);routdiag[p_out]=-1;rindiag=-np.ones(512);rindiag[p_in]=1
    W=rindiag[:,None]*(UL.conj().T@(routdiag[:,None]*UL))
    states=np.zeros((512,3),complex)
    for i,tables in enumerate(observable):
        P=np.eye(16,dtype=complex)
        for t in tables:apply_table(P,t,bankmat,b,4)
        states[:16,i]=P[:,0]
    state=states.copy();terms=[]
    for k in range(K+1):terms.append(states.conj().T@state);state=W@state
    C=np.array([np.outer(norm,norm)*np.einsum('k,kij->ij',c,np.array(terms)[:len(c)]).real-B.T@B for c in poly])
    err=float(np.max(abs(C-z['reference_C'])));projection=float(np.max(abs(exact-z['reference_C'])))
    shotplan=sampling(norm,projection);cost=summarize_costs(walk,preps,degrees,poly,shotplan['shots_per_entry'])
    save(walk,out/'controlled-walk.qpy') if cost['controlled_walk']['operations']<=CAP else None
    for i,p in enumerate(preps):save(p,out/f'prepare-{i}.qpy')
    complete=basic_complete(walk,preps,out)
    np.savez_compressed(out/'mathematical-results.npz',H=H,L=L,B=B,times=times,UL=UL,W=W,states=states,terms=np.array(terms),C=C,exact=exact,reference=z['reference_C'])
    receipt=dict(status='COMPONENT_SYNTHESIS_AND_SEMANTIC_COMPOSITION_COMPLETE',scope='Actual Clifford+T lookup, bank, reflection and preparation circuits. Full 138-plus-qubit statevector not executed; full model response follows verified clean-workspace semantic composition and measured single-qubit synthesis matrices.',
        input_sha256=sha(ROOT/'inputs/triangle.npz'),source_sha256=sha(__file__),protocol_sha256=sha(ROOT/'protocol.json'),
        rank=len(H),residues=3,Gamma=gamma,alpha=alpha,polynomial_degrees=degrees,precision_bits=b,
        per_Rz_synthesis_epsilon=epsrot,maximum_word_rounding_response_bound=roundbound,synthesis_response_bound=2*scale*2*b*paths*epsrot,
        full_Gaussian_response_error=err,response_pass=err<=.002,projected_encoding_to_H_error=float(np.max(abs(UL[np.ix_(p_out,p_in)]-padded/math.sqrt(gamma)))),
        logical_system_qubits=system,clean_workspace_qubits=workspace,total_logical_qubits=total,dirty_qubits_used=0,
        tables=table_receipts,rotation_bank=counts(rotbank),reflections=[counts(rout),counts(rin)],sampling=shotplan,cost=cost,
        complete_circuits=complete,maximum_order_complete_circuit_status='OPERATION_CAP' if cost['maximum_order_pair01_cached_operations']>CAP else 'NOT_EXPORTED',
        wall_seconds=time.perf_counter()-start,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    dump(out/'receipt.json',receipt);print(json.dumps({k:receipt[k] for k in ['status','precision_bits','full_Gaussian_response_error','clean_workspace_qubits','wall_seconds']},indent=2))

def dense(name):
    start=time.perf_counter();out=ROOT/'results'/name;out.mkdir(exist_ok=False);z,H,B,times,w,exact,norm=data()
    base=ROOT/'baseline'/name;prior=json.loads((base/'receipt.json').read_text());alpha=prior['alpha'];K=prior['maximum_order']
    def load(p):
        with p.open('rb') as f:return qpy.load(f)[0]
    original=load(base/'controlled-walk.qpy');preoriginal=[load(base/f'prepare-{i}.qpy') for i in range(3)]
    scale=float(max(norm)**2);rotation_bound=3*(K*original.count_ops().get('u',0)+2*max(p.count_ops().get('u',0) for p in preoriginal));eps=1e-6/(2*scale*rotation_bound)
    synth=Synthesizer(eps,out/'synthesis');mark('dense:compile:'+name)
    walk=synth.compile(original);preps=[synth.compile(p) for p in preoriginal]
    # Operator evaluation uses the emitted Clifford+T circuit, with its global phase retained.
    mark('dense:operator_check:'+name);W=Operator(walk).data;expected=Operator(original).data
    error=float(np.linalg.norm(W-expected,2));allowed=3*original.count_ops().get('u',0)*eps+1e-10
    assert error<=allowed
    states=[];zeros=np.zeros(64,complex);zeros[0]=1
    # Complete overlaps are simulated on sixqubits for all18outputs using emitted components.
    H0=QuantumCircuit(6);H0.h(0);hmat=Operator(H0).data
    Ps=[Operator(p).data for p in preps]
    term=np.zeros((K+1,3,3));Zdiag=np.where(np.arange(64)&1,-1.,1.)
    for j in range(3):
        state=Ps[j]@hmat@zeros
        for k in range(K+1):
            for i in range(3):
                final=hmat@Ps[i].conj().T@state;term[k,i,j]=float(np.vdot(final,Zdiag*final).real)
            state=W@state
    poly=[coefficients(float(t*alpha),1e-8/max(1,float(max(norm)**2)))[0] for t in times]
    C=np.array([np.outer(norm,norm)*np.einsum('k,kij->ij',c,term[:len(c)])-B.T@B for c in poly]);projection=float(np.max(abs(exact-z['reference_C'])))
    shotplan=sampling(norm,projection);cost=summarize_costs(walk,preps,[len(c)-1 for c in poly],poly,shotplan['shots_per_entry'])
    if cost['controlled_walk']['operations']<=CAP:save(walk,out/'controlled-walk.qpy')
    for i,p in enumerate(preps):save(p,out/f'prepare-{i}.qpy')
    complete=basic_complete(walk,preps,out);err=float(np.max(abs(C-z['reference_C'])))
    np.savez_compressed(out/'mathematical-results.npz',H=H,B=B,times=times,W=W,terms=term,C=C,exact=exact,reference=z['reference_C'])
    receipt=dict(status='ACTUAL_COMPONENT_CLIFFORD_T_AND_SIX_QUBIT_OPERATOR_REPLAY_COMPLETE',source_sha256=sha(__file__),input_sha256=sha(ROOT/'inputs/triangle.npz'),
        alpha=alpha,polynomial_degrees=[len(c)-1 for c in poly],per_Rz_epsilon=eps,maximum_Rz_count=rotation_bound,
        component_operator_error=error,component_operator_error_bound=allowed,full_Gaussian_response_error=err,response_pass=err<=.002,
        sampling=shotplan,cost=cost,complete_circuits=complete,maximum_order_complete_circuit_status='OPERATION_CAP' if cost['maximum_order_pair01_cached_operations']>CAP else 'NOT_EXPORTED',
        wall_seconds=time.perf_counter()-start,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    dump(out/'receipt.json',receipt);print(json.dumps({k:receipt[k] for k in ['status','full_Gaussian_response_error','component_operator_error','wall_seconds']},indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['select_swap','dense_norm_optimal','dense_factor_normalization']);a=p.parse_args()
    candidate() if a.mode=='select_swap' else dense(a.mode)
