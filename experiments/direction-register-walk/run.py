#!/usr/bin/env python3
"""Bounded, lookup-loaded direction-register canary; no protein-scale claim."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[key]='1'
from pathlib import Path
import sys,json,hashlib,time,resource,platform
import numpy as np
from scipy.linalg import expm
from scipy.special import expit,ive
import scipy,qiskit
from qiskit import QuantumCircuit,QuantumRegister,ClassicalRegister,qpy
from qiskit.quantum_info import Operator,Statevector
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'vendor'))
import current_circuit as old
from contraction_backend import install_backend
BACKEND_PARITY=install_backend()
TOL=1e-10
BEGIN=time.monotonic()

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,x):Path(p).write_text(json.dumps(x,indent=2)+'\n')
def maxabs(x):return float(np.max(np.abs(x)))
def mm(a,b):return np.einsum('ij,jk->ik',a,b,optimize=False)
def bounded():
    rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if sys.platform!='darwin':rss*=1024
    if time.monotonic()-BEGIN>1800 or rss>4*1024**3:raise RuntimeError('Frozen wall/memory cap reached')
    return int(rss)

def fixture():
    coords=np.array([(i,j) for i in (-1.,0.,1.) for j in (-1.,0.,1.)])
    x,y=coords.T;U=.5*(x*x+.7*y*y)+.1*x**4+.07*y**4+.05*x*x*y*y+.1*x
    P=np.eye(16);directions=np.zeros((16,8));directions[:,0]=1;edges=[]
    for a in range(9):
        i,j=divmod(a,3);P[a]=0;directions[a]=0
        for label,(di,dj) in enumerate(((1,0),(-1,0),(0,1),(0,-1)),1):
            ii,jj=i+di,j+dj
            if 0<=ii<3 and 0<=jj<3:
                b=3*ii+jj;P[a,b]=.25*expit(-(U[b]-U[a]));directions[a,label]=P[a,b]
                if label in (1,3):edges.append((a,label,b,label+1))
        P[a,a]=1-P[a].sum();directions[a,0]=P[a,a]
    pi=np.exp(-U);pi/=pi.sum();D=np.sqrt(P*P.T);H=8*(np.eye(9)-D[:9,:9])
    f=x+.25*x*y+.1*y*y;g=y*y+.3*x*y+.2*x
    vec=[]
    for z in (f,g):
        v=np.zeros(16);v[:9]=np.sqrt(pi)*(z-pi@z);v/=np.linalg.norm(v);vec.append(v)
    coeff=ive(np.arange(4),.4);coeff[1:]*=2
    poly=coeff[0]*np.eye(16)+coeff[1]*D;prev=np.eye(16);now=D.copy();cheb=[prev,now]
    for k in range(2,4):
        nxt=2*mm(D,now)-prev;poly+=coeff[k]*nxt;cheb.append(nxt);prev,now=now,nxt
    exact=np.eye(16);exact[:9,:9]=expm(-.05*H)
    return dict(coords=coords,U=U,P=P,pi=pi,directions=directions,D=D,H=H,left=vec[0],right=vec[1],coefficients=coeff,polynomial=poly/coeff.sum(),exponential=exact,chebyshev=np.array(cheb),edges=np.array(edges))

def adjacent_transposition(q,a,b,n):
    changed=a^b;assert changed and changed&(changed-1)==0
    target=changed.bit_length()-1;controls=[i for i in range(n) if i!=target]
    zeros=[i for i in controls if not (a>>i)&1]
    for i in zeros:q.x(i)
    q.mcx(controls,target)
    for i in zeros:q.x(i)

def transpose_basis(q,a,b,n):
    path=[a];current=a
    for i in range(n):
        if (a^b)>>i&1:current^=1<<i;path.append(current)
    steps=list(zip(path[:-1],path[1:]))
    for u,v in steps+steps[-2::-1]:adjacent_transposition(q,u,v,n)

def direction_gates(f):
    nx,nd=4,3
    x,d=QuantumRegister(nx,'x'),QuantumRegister(nd,'dir')
    prep=QuantumCircuit(x,d,name='A_dir')
    for state in range(9):
        gate=old.preparation(np.sqrt(f['directions'][state]),f'direction_row_{state}')
        prep.append(gate.control(nx,ctrl_state=state,annotated=False),list(x)+list(d))
    A=prep.to_gate()
    shift=QuantumCircuit(nx+nd,name='S_dir');permutation=np.arange(128)
    for a,da,b,db in f['edges']:
        u=int(a)+(int(da)<<nx);v=int(b)+(int(db)<<nx)
        transpose_basis(shift,u,v,nx+nd);permutation[u],permutation[v]=v,u
    S=shift.to_gate()
    reflection=QuantumCircuit(nd,name='R_dir')
    for i in range(nd):reflection.x(i)
    reflection.mcp(np.pi,list(range(nd-1)),nd-1)
    for i in range(nd):reflection.x(i)
    reflection.global_phase+=np.pi;R=reflection.to_gate()
    walk=QuantumCircuit(nx+nd,name='W_dir');walk.append(A,range(7));walk.append(S,range(7));walk.append(A.inverse(),range(7));walk.append(R,range(4,7));W=walk.to_gate()
    block=QuantumCircuit(2+nx+nd,name='ccW_dir')
    block.append(A,range(2,9));block.append(S.control(2,annotated=False),range(9));block.append(A.inverse(),range(2,9));block.append(R.control(2,annotated=False),[0,1,6,7,8])
    return A,S,R,W,block.to_gate(),permutation

def overlap(f,gates,negative=False):
    nx,nd,a=4,3,2;b=QuantumRegister(1,'b');c=QuantumRegister(a,'c');x=QuantumRegister(nx,'x');d=QuantumRegister(nd,'dir');out=ClassicalRegister(1,'out')
    circ=QuantumCircuit(b,c,x,d,out,name='direction_overlap_negative' if negative else 'direction_overlap')
    coef=old.preparation(np.sqrt(f['coefficients']/sum(f['coefficients'])),'COEF')
    left=old.preparation((-1 if negative else 1)*f['left'],'PREP_left');right=old.preparation(f['right'],'PREP_right')
    circ.h(b[0]);circ.append(right.control(1,annotated=False),list(b)+list(x));circ.append(coef.control(1,annotated=False),list(b)+list(c))
    for bit in range(a):
        for _ in range(2**bit):circ.append(gates[4],[b[0],c[bit]]+list(x)+list(d))
    circ.append(coef.inverse().control(1,annotated=False),list(b)+list(c));circ.append(left.inverse().control(1,annotated=False),list(b)+list(x));circ.h(b[0]);circ.measure(b[0],out[0])
    return circ

def expectation(circ):
    unitary= circ.remove_final_measurements(inplace=False)
    state=Statevector.from_instruction(unitary);p=state.probabilities([0]);return float(p[0]-p[1])

def main():
    out=ROOT/'results';out.mkdir(exist_ok=True);(out/'circuits').mkdir(exist_ok=True)
    dump(out/'contraction-backend-parity.json',BACKEND_PARITY);print('backend parity passed',flush=True)
    assert sha(ROOT/'protocol.json')==(ROOT/'protocol.sha256').read_text().strip()
    protocol=json.loads((ROOT/'protocol.json').read_text());assert sha(ROOT/'vendor/current_circuit.py')==protocol['historical_source']['sha256']
    f=fixture();np.savez_compressed(out/'fixture.npz',**f);gates=direction_gates(f);A,S,R,W,ccW,perm=gates
    checks={'stochastic_error':maxabs(f['P'].sum(1)-1),'minimum_probability':float(f['P'].min()),'detailed_balance_error':maxabs(f['pi'][:,None]*f['P'][:9,:9]-f['P'][:9,:9].T*f['pi'][None,:]),'permutation_bijective':len(set(perm.tolist()))==128,'permutation_involutive':bool(np.array_equal(perm[perm],np.arange(128)))}
    assert checks['stochastic_error']<TOL and checks['minimum_probability']>=0 and checks['detailed_balance_error']<TOL
    assert checks['permutation_bijective'] and checks['permutation_involutive']
    Sm=Operator(S).data;expectedS=np.zeros((128,128));expectedS[perm,np.arange(128)]=1
    checks['shift_full_domain_error']=maxabs(Sm-expectedS);assert checks['shift_full_domain_error']<TOL
    Am=Operator(A).data;checks['direction_preparation_error']=max(maxabs(Am[:,x]-np.eye(16)[x][None,:].repeat(8,axis=0).ravel()*np.repeat(np.sqrt(f['directions'][x]),16)) for x in range(16))
    assert checks['direction_preparation_error']<TOL
    # Every invalid label is fixed in the explicitly compiled permutation.
    invalid=[]
    for x in range(16):
        for d in range(8):
            if x>=9 or d==0 or d>=5 or f['directions'][x,d]==0:invalid.append(x+(d<<4))
    checks['fixed_boundary_padding_reserved_count']=len(invalid);assert np.array_equal(perm[invalid],invalid)
    Wm=Operator(W).data;checks['walk_unitarity_error']=maxabs(mm(Wm.conj().T,Wm)-np.eye(128))
    power=np.eye(128,dtype=complex);errors=[]
    for k in range(4):
        errors.append(maxabs(power[:16,:16]-f['chebyshev'][k]));power=mm(power,Wm)
    checks['projected_chebyshev_errors']=errors;assert max(errors)<TOL
    cc=Operator(ccW).data;target=np.eye(512,dtype=complex);active=np.arange(128)*4+3;target[np.ix_(active,active)]=Wm
    checks['factorized_control_absolute_phase_error']=maxabs(cc-target);assert checks['factorized_control_absolute_phase_error']<TOL
    checks['inactive_control_identity_errors']=[maxabs(cc[np.ix_(np.arange(128)*4+c,np.arange(128)*4+c)]-np.eye(128)) for c in range(3)]
    checks['controlled_inverse_error']=maxabs(mm(Operator(ccW.inverse()).data,cc)-np.eye(512));assert checks['controlled_inverse_error']<TOL
    ref=float(f['left']@f['polynomial']@f['right']);expref=float(f['left']@f['exponential']@f['right'])
    checks['polynomial_offdiagonal_chi']=ref;checks['dense_exponential_offdiagonal']=expref;checks['polynomial_tail_mass']=float(1-f['coefficients'].sum());checks['unnormalized_polynomial_vs_exponential_error']=abs(f['coefficients'].sum()*ref-expref)
    assert checks['unnormalized_polynomial_vs_exponential_error']<=checks['polynomial_tail_mass']+TOL
    dump(out/'algebra.json',checks);print('algebra',json.dumps(checks),flush=True);bounded()
    oldcirc,_,oldmeta=old.overlap_circuit(f['P'],f['coefficients'],f['left'],f['right'])
    logical={'direction':overlap(f,gates),'two_address':oldcirc,'direction_negative':overlap(f,gates,True)}
    compiled=[]
    for name,circ in logical.items():
        with (out/'circuits'/f'{name}_logical.qpy').open('wb') as h:qpy.dump(circ,h)
        for topology in (['all_to_all'] if name=='direction_negative' else ['all_to_all','line']):
            bounded();begin=time.monotonic();native=old.compile_circuit(circ,topology);seconds=time.monotonic()-begin
            count=old.resource_counts(native)
            if sum(count['operation_counts'].values())>1000000:raise RuntimeError('Frozen native-operation cap reached')
            path=out/'circuits'/f'{name}_{topology}.qpy'
            with path.open('wb') as h:qpy.dump(native,h)
            chi=expectation(native);wanted=-ref if name=='direction_negative' else ref;error=abs(chi-wanted);assert error<TOL,(name,topology,error)
            with path.open('rb') as h:loaded=qpy.load(h)[0]
            repeated=expectation(loaded);assert abs(repeated-wanted)<TOL
            item={'method':name,'topology':topology,'compile_seconds':seconds,'counts':count,'chi':chi,'expected_chi':wanted,'absolute_error':error,'reloaded_chi':repeated,'qpy_sha256':sha(path),'peak_RSS_bytes':bounded()}
            compiled.append(item);dump(out/'compiled-results.json',compiled);print('compiled',json.dumps(item),flush=True)
    final={'status':'passed','protocol_sha256':sha(ROOT/'protocol.json'),'source_sha256':{'run.py':sha(ROOT/'run.py'),'vendor/current_circuit.py':sha(ROOT/'vendor/current_circuit.py'),'contraction_backend.py':sha(ROOT/'contraction_backend.py')},'fixture_sha256':sha(out/'fixture.npz'),'algebra':checks,'compiled':compiled,'versions':{'python':platform.python_version(),'qiskit':qiskit.__version__,'numpy':np.__version__,'scipy':scipy.__version__},'total_seconds':time.monotonic()-BEGIN,'peak_RSS_bytes':bounded(),'scope':'Synthetic3x3 finite model with lookup loading. No primary protein cost, arithmetic loading, hardware, shot reduction or quantum advantage claim.'}
    dump(out/'result.json',final);print('finished',json.dumps({'status':'passed','seconds':final['total_seconds'],'peak_RSS_bytes':final['peak_RSS_bytes']}),flush=True)

if __name__=='__main__':
    try:main()
    except Exception as exc:
        (ROOT/'results').mkdir(exist_ok=True);dump(ROOT/'results/failure.json',{'status':'failed','type':type(exc).__name__,'message':str(exc),'elapsed_seconds':time.monotonic()-BEGIN});raise
