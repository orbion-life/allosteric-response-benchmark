"""Portable source-only confidence and actual tiny circuit checks; no writes."""
import json
import numpy as np
from scipy.stats import binom
from qiskit import QuantumCircuit
from qiskit.quantum_info import Statevector
from iqae import cp_interval,simulated_estimates
from spectral import observable,attenuation,preparation,grover

def main():
    successes=np.arange(33);lo,hi=cp_interval(successes,32,.01);p=np.linspace(0,1,101)
    coverage=np.sum(binom.pmf(successes[:,None],32,p[None,:])*((lo[:,None]<=p)&(p<=hi[:,None])),axis=0)
    assert coverage.min()>=.99-1e-12
    amplitudes=np.array([.1,.25,.5,.9]);r=simulated_estimates(amplitudes,.01,.01,2026091411)
    assert np.all(r['status']==1);assert np.max(r['a_upper']-r['a_lower'])<=.02
    ev=np.array([.3,1.]);x=np.array([1.,-2.]);y=np.array([2.,1.]);t=.7
    A=preparation(observable(x,1),observable(y,1),attenuation(ev,t,1));Q,_=grover(A)
    expected=(1+np.dot(x*np.exp(-t*ev),y)/np.linalg.norm(x)/np.linalg.norm(y))/2
    errors=[]
    for k in [0,1,2]:
        qc=QuantumCircuit(Q.num_qubits);qc.compose(A,range(A.num_qubits),inplace=True)
        for _ in range(k):qc.compose(Q,inplace=True)
        v=Statevector.from_instruction(qc);actual=float(v.probabilities([0])[0])
        target=np.sin((2*k+1)*np.arcsin(np.sqrt(expected)))**2;errors.append(abs(actual-target))
    assert max(errors)<1e-10
    print(json.dumps(dict(status='PASS',exact_CP_minimum_coverage=float(coverage.min()),
                         tiny_spectral_Grover_maximum_error=max(errors),
                         scope='Source-only confidence and actual four-qubit circuit semantics. No archived protein or complete resource-claim replay.'),indent=2))

if __name__=='__main__':main()
