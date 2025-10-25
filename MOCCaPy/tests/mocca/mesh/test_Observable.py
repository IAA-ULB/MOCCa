import pytest
import numpy as np

from mocca.mesh import LagrangeMesh
from mocca.mesh.observable import Observable


def test_Observable_ctor_scalar():
    n_gridpoints = 4*4*4
    for symmetry in [1,-1]:
        o = np.ones((n_gridpoints,), dtype=float)
        O = Observable(o, symmetry=symmetry)
        assert O.data.shape == (n_gridpoints,1)
        assert O.shape == (1,)
        assert O.symmetry[0] == symmetry


def test_Observable_ctor_vector():
    n_gridpoints = 4*4*4
    for symmetry in [np.array([1,1,1]), np.array([-1,-1,-1])]:
        o = np.ones((n_gridpoints,3), dtype=float)
        O = Observable(o, symmetry=symmetry)
        assert O.data.shape == (n_gridpoints,3)
        assert O.shape == (3,)
        assert (O.symmetry == symmetry).all()


def test_Observable_ctor_tensor():
    n_gridpoints = 4*4*4
    symmetries = [np.ones((3,3), dtype=int, order='F'),]
    s = np.ones((3,3), dtype=int, order='F')
    for i in range(3):
        for j in range(3):
            s[i,j] =  (-1)**(i+j)
    symmetries.append(s)
    s = np.ones((3,3), dtype=int, order='F')
    for i in range(3):
        for j in range(3):
            s[i,j] =  (-1)**(i+j+1)
    symmetries.append(s)
    s = np.ones((3,3), dtype=int, order='F')
    for i in range(3):
        for j in range(3):
            s[i,j] =  (-1) if i!=j else 1
    symmetries.append(s)

    for symmetry in symmetries:
        o = np.ones((n_gridpoints,3,3), dtype=float)
        O = Observable(o, symmetry=symmetry)
        assert O.data.shape == (n_gridpoints,9)
        assert O.shape == (3,3)
        assert (O.symmetry == symmetry.reshape((9,))).all()

