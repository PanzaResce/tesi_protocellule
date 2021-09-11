import numpy as np
from scipy import integrate
from scipy.integrate._ivp.common import norm

"""
These function are monkey-pathced to the scipy library, 
both replace the corresponding function without the '_mod' at the end.
This python file needs to be imported like this:
    'from protocellule import scipy_patch'

Both patch aim to add the delta time (dt) as a parameter to the integration function passed to the solve_ivp method
The dt is needed at each step to adjust the species' concentration according to the volume's change
"""


def rk_step_mod(fun, t, y, f, h, A, B, C, K):
    """
    Line changed: K[s] = fun(t + c * h, y + dy) --> K[s] = fun([t + c * h, c*h], y + dy)
    """
    K[0] = f
    for s, (a, c) in enumerate(zip(A[1:], C[1:]), start=1):
        dy = np.dot(K[:s].T, a[:s]) * h
        K[s] = fun([t + c * h, c*h], y + dy)

    y_new = y + h * np.dot(K[:-1].T, B)
    f_new = fun([t + h, h], y_new)

    K[-1] = f_new

    return y_new, f_new


# PATCH
integrate._ivp.rk.rk_step = rk_step_mod


def select_initial_step_mod(fun, t0, y0, f0, direction, order, rtol, atol):
    """
    Line changed: f1 = fun(t0 + h0 * direction, y1) --> f1 = fun([t0 + h0 * direction, h0 * direction], y1)
    """
    if y0.size == 0:
        return np.inf

    scale = atol + np.abs(y0) * rtol
    d0 = norm(y0 / scale)
    d1 = norm(f0 / scale)
    if d0 < 1e-5 or d1 < 1e-5:
        h0 = 1e-6
    else:
        h0 = 0.01 * d0 / d1

    y1 = y0 + h0 * direction * f0
    f1 = fun([t0 + h0 * direction, h0 * direction], y1)
    d2 = norm((f1 - f0) / scale) / h0

    if d1 <= 1e-15 and d2 <= 1e-15:
        h1 = max(1e-6, h0 * 1e-3)
    else:
        h1 = (0.01 / max(d1, d2)) ** (1 / (order + 1))

    return min(100 * h0, h1)


# PATCH
integrate._ivp.common.select_initial_step = select_initial_step_mod
