import numpy as np
import numpy.typing as npt

# TODO: obsolete?
def lagrange_function(x:npt.NDArray|float, x_i:float, d:float, M:int):
    """Compute the 1D Lagrange function on `x`.

    Args:
        x: array of points at which to compute the Lagrange function.
        x_i: a grid point : ±1/2 dx, ±3/2 dx, ±5/2 dx, ...
        d: grid spacing
        M: total number of grid points (positive and negative axis) = 2N
    Caveat:
        The function returns nan when `x == x_i` because it evaluates 0/0. Checking for this corner case implies giving
        up numpy array functions. As we assume that the user is not interested in finding the value at x_i because it
        is known to be 1. The user must avoid this. However, for `x = x_i + 1e-9` the result is very close to 1,
    """
    one_over_2N = 1/M
    A_i = (np.pi/d)*(x - x_i)
    result = one_over_2N*np.sin(A_i)/np.sin(one_over_2N*A_i)
    return result
