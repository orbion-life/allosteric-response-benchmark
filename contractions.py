"""Inspectable alternative to authored one-/two-dimensional dense @ operations.

Dense contractions use NumPy einsum with optimize=False, avoiding the dense
matmul route. Sparse operations retain SciPy dot. This does not establish the
cause of warnings from any other route, suppress warnings, or alter source files.
"""
import numpy as np
from scipy import sparse


def contract(a, b):
    """Match @ for the vector/matrix cases used by this synthetic calculation."""
    if sparse.issparse(a):
        return a.dot(b)
    if sparse.issparse(b):
        return b.T.dot(np.asarray(a).T).T
    a, b = np.asarray(a), np.asarray(b)
    forms = {(1, 1): 'i,i->', (1, 2): 'i,ij->j',
             (2, 1): 'ij,j->i', (2, 2): 'ij,jk->ik'}
    if (a.ndim, b.ndim) not in forms:
        raise ValueError('Only one-/two-dimensional contractions are supported')
    return np.einsum(forms[(a.ndim, b.ndim)], a, b, optimize=False)
