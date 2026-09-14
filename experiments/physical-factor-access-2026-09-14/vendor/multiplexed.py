"""One fixed real binary-tree loader with directly multiplexed controls."""
import math
import numpy as np
from qiskit import QuantumCircuit
from qiskit.circuit.library import UCRYGate
from qiskit.quantum_info import Operator

def width(x):
    n=max(1,int(math.ceil(math.log2(x))));return n,2**n

def real_loader(values, targets, context, total_qubits, readout=None):
    """Context bits index rows; higher target bits index binary-tree prefixes."""
    values=np.asarray(values)
    if np.iscomplexobj(values) and np.max(abs(values.imag))>1e-14:raise ValueError('Only real factors admitted in this follow-up')
    values=values.real
    if values.shape!=(2**len(context),2**len(targets)):raise ValueError('Padded loader dimensions')
    qc=QuantumCircuit(total_qubits,name='real_Ry_loader')
    n=len(targets)
    for k in range(n-1,-1,-1):
        higher=targets[k+1:];nprefix=2**len(higher);angles=[]
        for row in range(len(values)):
            for prefix in range(nprefix):
                start=prefix*2**(k+1);part=values[row,start:start+2**(k+1)]
                if k:
                    low=float(np.linalg.norm(part[:2**k]));high=float(np.linalg.norm(part[2**k:]))
                else:low=float(part[0]);high=float(part[1])
                angles.append(2*math.atan2(high,low))
        controls=higher+context
        if readout is not None:
            # Readout is the highest angle-index bit. Its zero branch is identity.
            angles=[0.]*len(angles)+angles;controls=controls+[readout]
        qc.append(UCRYGate(angles),[targets[k]]+controls)
    return qc

def reflection(qc,bits,readout=None,input_reflection=False):
    qc.x(bits)
    controls=bits[:-1] if readout is None else [readout]+bits[:-1]
    if controls:qc.mcp(math.pi,controls,bits[-1])
    else:qc.z(bits[-1])
    qc.x(bits)
    if input_reflection:
        if readout is None:qc.global_phase+=math.pi
        else:qc.z(readout)

def loader_circuit(L,gamma,controlled=False):
    nr,dr=width(L.shape[0]);nc,dc=width(L.shape[1]);offset=int(controlled)
    lp=np.zeros((dr,dc));lp[:L.shape[0],:L.shape[1]]=np.asarray(L).real
    if np.iscomplexobj(L) and np.max(abs(np.asarray(L).imag))>1e-14:raise ValueError('Complex factor outside follow-up')
    rho=np.linalg.norm(lp,axis=1);F=float(np.linalg.norm(rho));eta=F/math.sqrt(gamma)
    if not 0<eta<=1+1e-13:raise ValueError('Invalid attenuation')
    eta=min(eta,1.);col=list(range(offset,offset+nc));row=list(range(offset+nc,offset+nc+nr));flag=offset+nc+nr
    total=nc+nr+1+offset;readout=0 if controlled else None
    pin=real_loader((rho/F)[None,:],row,[],total,readout)
    pout=real_loader(lp,col,row,total,readout)
    ul=QuantumCircuit(total,name='controlled_mux_UL' if controlled else 'mux_UL')
    ul.compose(pin,inplace=True);ul.compose(pout.inverse(),inplace=True)
    angle=2*math.acos(eta)
    if controlled:ul.cry(angle,readout,flag)
    else:ul.ry(angle,flag)
    return ul,pin,pout,dict(nrow=nr,ncol=nc,padded_rows=dr,padded_columns=dc,F=F,eta=eta,attenuation_angle=angle,col=col,row=row,flag=flag,readout=readout)

def controlled_walk(L,gamma):
    ul,pin,pout,m=loader_circuit(L,gamma,True)
    qc=QuantumCircuit(ul.num_qubits,name='direct_controlled_mux_W');qc.compose(ul,inplace=True)
    reflection(qc,m['col']+[m['flag']],0)
    qc.compose(ul.inverse(),inplace=True)
    reflection(qc,m['row']+[m['flag']],0,True)
    return qc

def factor_encoding(L,gamma):
    ul,pin,pout,m=loader_circuit(L,gamma,False)
    UL=Operator(ul).data;dc=m['padded_columns'];dr=m['padded_rows']
    good_in=np.arange(dc);good_out=np.arange(dr)*dc
    rout=np.ones(len(UL));rout[good_out]=-1;rin=-np.ones(len(UL));rin[good_in]=1
    V=UL.conj().T@(rout[:,None]*UL);W=rin[:,None]*V
    qc=QuantumCircuit(ul.num_qubits,name='mux_W');qc.compose(ul,inplace=True)
    reflection(qc,m['col']+[m['flag']]);qc.compose(ul.inverse(),inplace=True);reflection(qc,m['row']+[m['flag']],input_reflection=True)
    lp=np.zeros((dr,dc));lp[:L.shape[0],:L.shape[1]]=np.asarray(L).real
    m['checks']={'factor_overlap':float(np.max(abs(UL[np.ix_(good_out,good_in)]-lp/math.sqrt(gamma)))),
                 'UL_unitarity':float(np.max(abs(UL.conj().T@UL-np.eye(len(UL))))),
                 'V_hermitian':float(np.max(abs(V-V.conj().T))),
                 'V_involution':float(np.max(abs(V@V-np.eye(len(V))))),
                 'projected_A':float(np.max(abs(V[np.ix_(good_in,good_in)]-(np.eye(dc)-2*lp.T@lp/gamma)))),
                 'W_unitarity':float(np.max(abs(W.conj().T@W-np.eye(len(W)))))}
    m['good_input_indices']=good_in.tolist();m['good_output_indices']=good_out.tolist()
    m['loader_scope']='Separate fixed real-Ry multiplexed development revision; row/control register widths nrow/ncol are qubit counts.'
    return qc,W,UL,m
