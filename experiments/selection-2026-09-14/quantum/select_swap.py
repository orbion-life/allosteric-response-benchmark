"""Implemented clean SELECT–SWAP and exact dirty echo, with restored prefix phases."""
import math
import numpy as np
from qiskit import QuantumCircuit
from qiskit.circuit.library import RCCXGate

def choose_lambda(N,b):
    candidates=[2**s for s in range(N.bit_length())]
    return min(candidates,key=lambda lam:(16*max(N//lam-2,0)+7*b*(lam-1),lam))

def select_network(words,b,lam):
    N=len(words);n=(N-1).bit_length();s=(lam-1).bit_length();m=n-s
    assert N==2**n and lam==2**s
    anc=max(0,m-1);qc=QuantumCircuit(n+lam*b+anc,name='prefix_SELECT')
    hi=list(range(s,n));work=list(range(n+lam*b,qc.num_qubits))
    def write(prefix,control):
        for lane in range(lam):
            value=int(words[prefix*lam+lane])
            for bit in range(b):
                if value>>bit&1:
                    if control is None:qc.x(n+lane*b+bit)
                    else:qc.cx(control,n+lane*b+bit)
    def visit(depth,prefix,control):
        if depth==m:write(prefix,control);return
        bit=hi[m-1-depth]
        for val in [0,1]:
            if val==0:qc.x(bit)
            if control is None:visit(depth+1,prefix*2+val,bit)
            else:
                target=work[depth-1]
                qc.append(RCCXGate(),[control,bit,target])
                visit(depth+1,prefix*2+val,target)
                qc.append(RCCXGate().inverse(),[control,bit,target])
            if val==0:qc.x(bit)
    if m==0:write(0,None)
    else:visit(0,0,None)
    return qc

def swap_network(n,b,lam,anc):
    qc=QuantumCircuit(n+lam*b+anc,name='exact_SWAP')
    for level in range((lam-1).bit_length()):
        step=2**level
        for start in range(0,lam,2*step):
            for bit in range(b):qc.cswap(level,n+start*b+bit,n+(start+step)*b+bit)
    return qc

def lookup(words,b,lam=None):
    if lam is None:lam=choose_lambda(len(words),b)
    select=select_network(words,b,lam);n=(len(words)-1).bit_length()
    swaps=swap_network(n,b,lam,select.num_qubits-n-lam*b)
    qc=select.copy();qc.name='clean_SELECT_SWAP';qc.compose(swaps,inplace=True)
    return qc,dict(N=len(words),address_qubits=n,word_bits=b,lambda_value=lam,
                   clean_data_qubits=lam*b,prefix_clean_qubits=select.num_qubits-n-lam*b,
                   dirty_qubits=0,lookup_T_exact_before_optimization=16*max(len(words)//lam-2,0)+7*b*(lam-1),
                   select_Rccx=4*max(len(words)//lam-2,0),swap_Fredkin=b*(lam-1))

def dirty_lookup(words,b,lam):
    clean,meta=lookup(words,b,lam);n=meta['address_qubits']
    select=select_network(words,b,lam);swaps=swap_network(n,b,lam,meta['prefix_clean_qubits'])
    qc=QuantumCircuit(clean.num_qubits+b,name='dirty_echo_SELECT_SWAP');mapping=list(range(clean.num_qubits));out=list(range(clean.num_qubits,qc.num_qubits))
    qc.compose(clean,mapping,inplace=True)
    for k in range(b):qc.cx(n+k,out[k])
    qc.compose(clean.inverse(),mapping,inplace=True)
    qc.compose(swaps,mapping,inplace=True)
    for k in range(b):qc.cx(n+k,out[k])
    qc.compose(swaps.inverse(),mapping,inplace=True)
    return qc,dict(meta,clean_data_qubits=b,dirty_qubits=lam*b)

def quantize(theta,b):
    return np.rint(np.mod(np.asarray(theta,float),4*np.pi)*(2**b/(4*np.pi))).astype(np.uint64) % np.uint64(2**b)

def rotation_oracle(theta,b,lam=None):
    words=quantize(theta,b);oracle,meta=lookup(words,b,lam);n=meta['address_qubits']
    qc=QuantumCircuit(oracle.num_qubits+1,name='angle_lookup_Ry_uncompute');target=oracle.num_qubits
    qc.compose(oracle,range(oracle.num_qubits),inplace=True)
    for bit in range(b):qc.cry(4*math.pi*2**bit/2**b,n+bit,target)
    qc.compose(oracle.inverse(),range(oracle.num_qubits),inplace=True)
    return qc,meta,words

def tables_for_values(values,targets,context,readout=None):
    values=np.asarray(values,float);n=len(targets);tables=[]
    for k in range(n-1,-1,-1):
        higher=targets[k+1:];part=values.reshape(len(values),-1,2,2**k)
        low=np.linalg.norm(part[:,:,0,:],axis=2) if k else part[:,:,0,0]
        high=np.linalg.norm(part[:,:,1,:],axis=2) if k else part[:,:,1,0]
        theta=(2*np.arctan2(high,low)).ravel();controls=higher+context
        if readout is not None:theta=np.concatenate([np.zeros_like(theta),theta]);controls=controls+[readout]
        tables.append(dict(theta=theta,controls=controls,target=targets[k]))
    return tables
