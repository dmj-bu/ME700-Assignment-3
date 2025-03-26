import numpy as np


def residual_fcn(u):
    """Simple nonlinear system residual"""
    r = np.zeros(2)
    r[0] = u[0]**2 + u[1] - 5
    r[1] = u[0] + u[1]**2 - 3
    return r


def jacobian_fcn(u):
    """Jacobian of the residual"""
    J = np.zeros((2, 2))
    J[0, 0] = 2*u[0]
    J[0, 1] = 1
    J[1, 0] = 1
    J[1, 1] = 2*u[1]
    return J
