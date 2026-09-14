"""Receiver scores with an explicit distinction between masks and index arrays."""
import numpy as np


def validated_indices(selector, size, *, kind):
    """Validate a declared encoding without coercing one encoding into another."""
    values = np.asarray(selector)
    if values.ndim != 1:
        raise ValueError('A selector must have one dimension')
    if kind == 'mask':
        if values.dtype != np.dtype(bool):
            raise TypeError('A mask must have boolean dtype')
        if len(values) != size:
            raise ValueError('A mask must match the indexed axis length')
        indices = np.flatnonzero(values)
    elif kind == 'indices':
        if not np.issubdtype(values.dtype, np.integer):
            raise TypeError('Indices must have integer dtype; a boolean mask is not an index array')
        if np.any(values < 0) or np.any(values >= size):
            raise ValueError('An index is outside the indexed axis')
        if len(np.unique(values)) != len(values):
            raise ValueError('Duplicate indices would change receiver weighting')
        indices = values.astype(np.intp, copy=False)
    else:
        raise ValueError('Declare kind as mask or indices')
    if len(indices) == 0:
        raise ValueError('A selector must include at least one entry')
    return indices


def receiver_rms(response, receiver):
    """RMS across the explicitly indexed receiver columns, for every source row."""
    response = np.asarray(response)
    if response.ndim != 2:
        raise ValueError('A response must have two dimensions')
    indices = validated_indices(receiver, response.shape[1], kind='indices')
    return np.sqrt(np.mean(response[:, indices] ** 2, axis=1))
