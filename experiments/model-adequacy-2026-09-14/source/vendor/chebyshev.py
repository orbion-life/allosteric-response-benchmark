"""Independent imaginary-time action with explicit uniform tail estimate."""
import math
import time
import numpy as np
from scipy.special import ive


def log_tail_bound(z, degree):
    """2 exp[z(cosh(s)-1)-(degree+1)s], s=asinh((degree+1)/z)."""
    if z == 0:
        return -math.inf
    n = degree+1
    # hypot(z,n)-z is evaluated without its near-equal subtraction.
    return math.log(2.)+n*n/(math.hypot(z, n)+z)-n*math.asinh(n/z)


def required_degree(z, tolerance):
    if not (z >= 0 and 0 < tolerance < 1):
        raise ValueError('Invalid degree criterion.')
    target = math.log(tolerance)
    low, high = -1, 1
    while log_tail_bound(z, high) > target:
        high *= 2
    while high-low > 1:
        mid = (low+high)//2
        if log_tail_bound(z, mid) <= target:
            high = mid
        else:
            low = mid
    return high


def action(H, F, t, upper, covariance_tolerance=1e-9):
    started = time.perf_counter()
    norms = np.linalg.norm(F, axis=0)
    amplification = float(np.max(norms)**2)
    if upper <= 0 or t < 0:
        raise ValueError('Invalid spectral interval or time.')
    z = t*upper/2
    operator_tolerance = min(.5, covariance_tolerance/max(amplification, np.finfo(float).tiny))
    degree = required_degree(z, operator_tolerance)
    coeff = ive(np.arange(degree+1), z)
    if not np.isfinite(coeff).all():
        raise FloatingPointError('Invalid scaled Bessel coefficients.')
    old = F.copy()
    result = coeff[0]*old
    calls = 0
    if degree:
        current = F-(2/upper)*(H @ F)
        calls += 1
        result += 2*coeff[1]*current
        for k in range(2, degree+1):
            new = 2*current-(4/upper)*(H @ current)-old
            calls += 1
            result += 2*coeff[k]*new
            old, current = current, new
    bound = math.exp(log_tail_bound(z, degree))
    return result, dict(degree=degree, z=z, spectral_upper_bound=upper,
                        uniform_operator_tail_bound=bound,
                        maximum_normalized_covariance_tail_bound=bound*amplification,
                        observable_norms=norms.tolist(), covariance_amplification=amplification,
                        three_column_actions=calls, vector_equivalent_products=calls*F.shape[1],
                        seconds=time.perf_counter()-started,
                        scope='Analytic truncation inequality evaluated in floating point; excludes floating-point recurrence and coefficient errors. Not an interval-arithmetic certificate.')
