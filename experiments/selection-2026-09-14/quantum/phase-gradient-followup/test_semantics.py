"""Exhaustive small exact-adder and signed phase-gradient rotation checks."""
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[key]='1'
import json,sys,time
from pathlib import Path
import numpy as np
from qiskit.quantum_info import Operator
from gradient import rotation_bank,fourier_state
from synthesis import counts,save
ROOT=Path(__file__).resolve().parent

def main():
    start=time.perf_counter();out=ROOT/'results';out.mkdir(exist_ok=True);records=[]
    for b in [1,2,3]:
        bank,adder,meta=rotation_bank(b);A=Operator(adder).data;U=Operator(bank).data;D=2**b;control=2*b+1
        erradd=0.;errbank=0.
        for word in range(D):
            for x in range(D):
                for c in range(2):
                    i=word+(x<<b)+(c<<control);j=word+(((x+c*word)%D)<<b)+(c<<control)
                    desired=np.zeros(len(A),complex);desired[j]=1;erradd=max(erradd,float(np.max(abs(A[:,i]-desired))))
            f=fourier_state(b);angle=4*np.pi*word/D;R=np.array([[np.cos(angle/2),-np.sin(angle/2)],[np.sin(angle/2),np.cos(angle/2)]])
            for c in range(2):
                initial=np.zeros(len(U),complex);expected=np.zeros(len(U),complex)
                for x in range(D):
                    initial[word+(x<<b)+(c<<control)]=f[x]
                    for result in range(2):expected[word+(x<<b)+(result<<control)]=f[x]*R[result,c]
                errbank=max(errbank,float(np.max(abs(U@initial-expected))))
        assert max(erradd,errbank)<1e-9
        save(bank,out/f'phase-bank-b{b}.qpy');records.append(dict(bits=b,adder_error=erradd,signed_Ry_and_gradient_return_error=errbank,metadata=meta,bank=counts(bank),adder=counts(adder)))
    r=dict(status='PASS',wall_seconds=time.perf_counter()-start,cases=records)
    (out/'semantic-tests.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))

if __name__=='__main__':main()
