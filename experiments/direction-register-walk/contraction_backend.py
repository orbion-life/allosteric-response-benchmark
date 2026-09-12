"""Verification-only local override of a Qiskit tensor contraction.

The index map follows Qiskit2.5.2 Operator._einsum_matmul; only NumPy's
contraction optimizer is enabled. The installed Qiskit files are unchanged.
The parity check runs before installing this process-local classmethod.
"""
import inspect,hashlib,warnings
import numpy as np
from qiskit.quantum_info import Operator
from qiskit.exceptions import QiskitError

def optimized(cls,tensor,mat,indices,shift=0,right_mul=False):
    if mat.ndim%2:raise QiskitError('Contracted matrix must have an even number of indices.')
    rank=tensor.ndim;it=list(range(rank))
    for j,index in enumerate(indices):it[index+shift]=rank+j
    contracted=list(reversed(range(rank,rank+len(indices))))
    free=[index+shift for index in reversed(indices)]
    im=contracted+free if right_mul else free+contracted
    result=np.einsum(tensor,it,mat,im,optimize=True)
    if not np.isfinite(result).all():raise FloatingPointError('Nonfinite optimized verification contraction')
    return result

def install_backend():
    warnings.simplefilter('error',RuntimeWarning)
    original=Operator._einsum_matmul
    source_sha=hashlib.sha256(inspect.getsource(original).encode()).hexdigest()
    rng=np.random.default_rng(20260912);results=[]
    for shape,indices,shift in [((2,)*4,[0],0),((2,)*4,[1],2),((2,)*6,[0,2],0),((2,)*6,[0,2],3)]:
        tensor=rng.normal(size=shape)+1j*rng.normal(size=shape)
        matrix=rng.normal(size=(2,)*(2*len(indices)))+1j*rng.normal(size=(2,)*(2*len(indices)))
        for right in (False,True):
            reference=original(tensor,matrix,indices,shift,right)
            actual=optimized(Operator,tensor,matrix,indices,shift,right)
            assert np.isfinite(reference).all() and np.isfinite(actual).all()
            error=float(np.max(abs(reference-actual)));assert error<1e-12,error
            results.append({'tensor_shape':list(shape),'indices':indices,'shift':shift,'right_mul':right,'max_absolute_error':error})
    Operator._einsum_matmul=classmethod(optimized)
    return {'status':'passed','legacy_qiskit_method_source_sha256':source_sha,'fixture_count':len(results),'parity':results,'change':'Identical contraction indices, optimize=True, finite checks; process-local method override only','RuntimeWarning_policy':'error'}
