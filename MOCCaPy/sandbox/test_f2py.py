import numpy as np

import fill_f90 as f90
import mocca.util as util

def test_fill():
    n = 10
    a = np.zeros(n, dtype=float)
    assert (a == 0.0).all()
    a0 = 2.
    f90.fill(a, a0)
    assert (a == a0).all()

def test_zero_size_array():
    a = np.array([], dtype=float)
    a0 = 2.
    f90.fill(a, a0)
