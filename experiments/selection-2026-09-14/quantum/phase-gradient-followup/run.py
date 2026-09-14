"""Single declared phase-gradient follow-up using the primary exact lookup circuits."""
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[key]='1'
import hashlib,json,math,resource,sys,time
from pathlib import Path
import mpmath as mp
import numpy as np
from qiskit import QuantumCircuit,qpy
from qiskit.quantum_info import Operator
from gradient import rotation_bank
ROOT=Path(__file__).resolve().parent;PARENT=ROOT.parent;sys.path.insert(0,str(PARENT))
from synthesis import Synthesizer,counts,save
from candidate_experiment import data,apply_table,reflection,summarize_costs,sampling
from select_swap import tables_for_values,quantize
from scipy.linalg import cholesky
from reduced_overlap import coefficients

def dump(path,value):Path(path).write_text(json.dumps(value,indent=2)+'\n')
def load(path):
    with Path(path).open('rb') as f:return qpy.load(f)[0]

def mp_unitary(record):
    a=mp.mpc(1);b=mp.mpc(0);c=mp.mpc(0);d=mp.mpc(1);w=mp.exp(mp.j*mp.pi/4);r=mp.sqrt(2)
    if record['gates'] is None:
        k=round(record['angle']/(math.pi/4));return mp.exp(-mp.j*k*mp.pi/8),b,c,mp.exp(mp.j*k*mp.pi/8)
    for g in record['gates']:
        if g=='H':a,b,c,d=(a+b)/r,(a-b)/r,(c+d)/r,(c-d)/r
        elif g=='T':b*=w;d*=w
        elif g=='S':b*=mp.j;d*=mp.j
        elif g=='X':a,b=b,a;c,d=d,c
        elif g=='W':a*=w;b*=w;c*=w;d*=w
        else:raise ValueError(g)
    return a,b,c,d

def main():
    start=time.perf_counter();out=ROOT/'results';out.mkdir(exist_ok=True)
    z,H,B,times,ev,exact,norm=data();scale=float(max(norm)**2);b=33
    primary=json.loads((PARENT/'results/select_swap/receipt.json').read_text());assert b==primary['precision_bits']
    bank,adder,addermeta=rotation_bank(b);bc=counts(bank);save(bank,out/'phase-gradient-bank-b33.qpy');save(adder,out/'controlled-adder-b33.qpy')
    synth=Synthesizer(1e-6/(2*scale*b),out/'gradient-preparation-synthesis');gprep=QuantumCircuit(b)
    for bit in range(b):gprep.h(bit);gprep.compose(synth.rz(-2*math.pi*2**bit/2**b),[bit],inplace=True)
    save(gprep,out/'gradient-preparation-b33.qpy');gc=counts(gprep)
    # Exact gate strings evaluated at high precision give a product-state overlap.
    mp.mp.dps=70;overlap=mp.mpc(1)
    for bit,record in enumerate(synth.records):
        a,bb,c,d=mp_unitary(record);p0=(a+bb)/mp.sqrt(2);p1=(c+d)/mp.sqrt(2)
        phase=mp.exp(-2*mp.j*mp.pi*2**bit/2**b);overlap*=(p0+mp.conj(phase)*p1)/mp.sqrt(2)
    state_error=mp.sqrt(max(mp.mpf(0),2-2*abs(overlap)));response_gradient_bound=2*scale*float(state_error)
    assert response_gradient_bound<=1e-6
    system=10;workspace=primary['clean_workspace_qubits'];gradient=list(range(system+workspace,system+workspace+b));helpbit=gradient[-1]+1;scratch=helpbit+1;total=scratch+1
    def compile_tables(label):
        result=QuantumCircuit(total)
        tables=[x for x in primary['tables'] if x['label']==label]
        for table in tables:
            oracle=load(PARENT/'results/select_swap'/f"lookup-{label}-{table['index']}.qpy")
            mapping=table['controls']+list(range(system,system+oracle.num_qubits-len(table['controls'])))
            result.compose(oracle,mapping,inplace=True)
            result.compose(bank,list(range(system,system+b))+gradient+[helpbit,table['target'],scratch],inplace=True)
            result.compose(oracle.inverse(),mapping,inplace=True)
        return result
    pin=compile_tables('row');pout=compile_tables('column');ul=pin.copy();ul.compose(pout.inverse(),inplace=True)
    rout=reflection([0]+list(range(1,5))+[9],total,list(range(system,system+workspace)))
    rin=reflection([0]+list(range(5,9))+[9],total,list(range(system,system+workspace)),True)
    walk=ul.copy();walk.compose(rout,inplace=True);walk.compose(ul.inverse(),inplace=True);walk.compose(rin,inplace=True)
    preps=[compile_tables('observable'+str(i)) for i in range(3)]
    save(walk,out/'controlled-walk.qpy')
    for i,p in enumerate(preps):save(p,out/f'prepare-{i}.qpy')
    L=cholesky(H,lower=False);gamma=float(np.linalg.norm(L)**2);alpha=gamma/2;poly=[coefficients(float(t*alpha),1e-8/max(1,scale))[0] for t in times];degrees=[len(c)-1 for c in poly]
    # Under an idealgradient, exact adders implement the alreadyquantized Ry tables.
    padded=np.zeros((16,16));padded[:9,:9]=L;rho=np.linalg.norm(padded,axis=1)
    row=tables_for_values((rho/np.linalg.norm(rho))[None,:],list(range(5,9)),[],0);column=tables_for_values(padded,list(range(1,5)),list(range(5,9)),0)
    idealbank=[]
    for bit in range(b):
        a=2*math.pi*2**bit/2**b;idealbank.append(np.array([[math.cos(a/2),-math.sin(a/2)],[math.sin(a/2),math.cos(a/2)]],complex))
    UI=np.eye(512,dtype=complex);UO=np.eye(512,dtype=complex)
    for t in row:apply_table(UI,t,idealbank,b,9)
    for t in column:apply_table(UO,t,idealbank,b,9)
    UL=UO.conj().T@UI;rd=np.ones(512);rd[np.arange(16)*16]=-1;ri=-np.ones(512);ri[:16]=1;W=ri[:,None]*(UL.conj().T@(rd[:,None]*UL))
    vectors=np.zeros((16,3));vectors[:9]=B/norm;states=np.zeros((512,3),complex)
    for i in range(3):
        P=np.eye(16,dtype=complex)
        for t in tables_for_values(vectors[:,i][None,:],list(range(1,5)),[],0):apply_table(P,t,idealbank,b,4)
        states[:16,i]=P[:,0]
    evolved=states.copy();terms=[]
    for k in range(max(degrees)+1):terms.append(states.conj().T@evolved);evolved=W@evolved
    C=np.array([np.outer(norm,norm)*np.einsum('k,kij->ij',c,np.array(terms)[:len(c)]).real-B.T@B for c in poly]);idealerror=float(np.max(abs(C-z['reference_C'])))
    shotplan=sampling(norm,float(np.max(abs(exact-z['reference_C']))));cost=summarize_costs(walk,preps,degrees,poly,shotplan['shots_per_entry'])
    # Charge fresh gradientpreparation once for everycompleteestimator shot.
    cost['expected_total_T']+=cost['all_query_readouts']*gc['T'];cost['expected_total_CX']+=cost['all_query_readouts']*gc['CX']
    cost['maximum_order_pair01_cached_operations']+=gc['operations']
    for query in cost['queries']:
        query['expected_T_per_shot']+=gc['T'];query['expected_CX_per_shot']+=gc['CX'];query['serial_T_depth_upper_per_shot']+=gc['T_depth']
    complete=[]
    for k in [0,1]:
        qc=QuantumCircuit(total,1);qc.compose(gprep,gradient,inplace=True);qc.h(0);qc.compose(preps[1],inplace=True)
        for _ in range(k):qc.compose(walk,inplace=True)
        qc.compose(preps[0].inverse(),inplace=True);qc.h(0);qc.measure(0,0)
        if sum(qc.count_ops().values())<=1000000:save(qc,out/f'complete-k{k}.qpy');status='EXPORTED_ACTUAL_CLIFFORD_T_CIRCUIT'
        else:status='OPERATION_CAP'
        complete.append(dict(order=k,status=status,**counts(qc)))
    # Offer the same implemented constant-word adder to each dense Rz at its requiredbit precision.
    fairness=[]
    for name in ['dense_norm_optimal','dense_factor_normalization']:
        dreceipt=json.loads((PARENT/'results'/name/'receipt.json').read_text());R=dreceipt['maximum_Rz_count'];db=math.ceil(math.log2(2*math.pi*scale*R/1e-6));dbank,_,_=rotation_bank(db);dc=counts(dbank);save(dbank,out/(name+'-constant-angle-bank.qpy'))
        records=json.loads((PARENT/'results'/name/'synthesis/rotation-records.json').read_text());decisions=[]
        for record in records:
            word=int(quantize([record['angle']],db)[0]);option=0 if word==0 else dc['T']
            decisions.append(dict(angle_hex=record['angle_hex'],word=word,direct_T=record['T'],implemented_constant_adder_T=option,selected='direct_Ross_Selinger' if record['T']<=option else 'phase_gradient'))
        assert all(x['selected']=='direct_Ross_Selinger' for x in decisions)
        dump(out/(name+'-rotation-choice.json'),decisions)
        fairness.append(dict(name=name,precision_bits=db,compiled_constant_bank=dc,angles_compared=len(decisions),maximum_direct_Rz_T=max(x['direct_T'] for x in decisions),phase_gradient_selected=0,
            expected_total_T_retained=dreceipt['cost']['expected_total_T'],scope='Same implemented controlledCDKMbank with classicalwordXinitialization/uninitialization haszeroextraT forlookup. EverydirectRz is cheaper. Gradientinitialization wouldaddnonnegativecost. Specializedconstantadders, jointrotation synthesis andgloballyoptimaldensecircuits remainunassessed.'))
    dense=min(x['expected_total_T_retained'] for x in fairness);ratio=cost['expected_total_T']/dense
    np.savez_compressed(out/'mathematical-results.npz',H=H,L=L,B=B,times=times,W=W,states=states,terms=np.array(terms),ideal_gradient_C=C,exact=exact,reference=z['reference_C'])
    receipt=dict(status='FOLLOWUP_COMPONENTS_COMPILED_AND_COMPOSITION_BOUND_VERIFIED',source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),input_sha256=hashlib.sha256((PARENT/'inputs/triangle.npz').read_bytes()).hexdigest(),
        Gamma=gamma,alpha=alpha,polynomial_degrees=degrees,precision_bits=b,gradient_bank=bc,controlled_adder=counts(adder),adder_metadata=addermeta,gradient_preparation=gc,
        gradient_state_phase_aligned_distance_70_digits=str(state_error),gradient_preparation_response_error_bound=response_gradient_bound,
        gradient_initialization_policy='Fresh oncepercompleteestimator shot; reusedthrough allwalkpowers andloaders. No unpreparationneededbecausegradientisdiscarded afterreadout. Initialstate normerror givesuniform wholecircuitexpectationbound despitepossible entanglement. No reuseacrossshots.',
        ideal_gradient_Gaussian_response_error=idealerror,prepared_gradient_Gaussian_response_error_upper=idealerror+response_gradient_bound,
        response_pass=idealerror+response_gradient_bound<=.002,system_qubits=system,lookup_clean_workspace=workspace,gradient_qubits=b,adder_helper_qubits=2,total_logical_qubits=total,dirty_qubits=0,
        sampling=shotplan,cost=cost,complete_circuits=complete,full_maximum_order_circuit_status='OPERATION_CAP' if cost['maximum_order_pair01_cached_operations']>1000000 else 'NOT_EXPORTED',
        dense_fairness=fairness,candidate_to_strongest_evaluated_dense_T_ratio=ratio,half_T_cost_gate_pass=ratio<=.5,
        scope='Measured exactClifford+T fullsizecomponents andexported lowordercompletecircuits; tinyfullunitarysemantictests. Largeworkspacefullstatevector andmaximumordercompletecircuit notexecuted. Modelresponse isideal-gradient semanticcomposition plusmeasuredpreparationstateerror bound. No hardware, faultcorrection or classicaladvantageclaim.',
        wall_seconds=time.perf_counter()-start,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    dump(out/'receipt.json',receipt);print(json.dumps(receipt,indent=2))

if __name__=='__main__':main()
