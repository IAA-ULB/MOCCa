import pytest
import numpy as np

from mocca.mesh import LagrangeMesh
from mocca.mesh.observable import Observable

from src_heph.heph_functional import rreplace


def test_Observable_ctor():
    mesh = LagrangeMesh(dim=1, M=4, d=1., reduced=False)
    Q = Observable(mesh, n_components=1)

    data = np.zeros((mesh.n_gridpoints(), 2), order='F')

def test_derive1_1D(debug=False):
    for reduced in [
        False,
        # True,
    ]:
        mesh = LagrangeMesh(dim=1, M=4, d=1., reduced=reduced)
        Q = Observable(mesh, data = mesh.basis_function(ijk=0, r=mesh.gx))
        Q.derive1()
        print("test_derive1_1D finished")

def test_derive1_2D(debug=False):
    for reduced in [
        False,
        # True,
    ]:
        mesh = LagrangeMesh(dim=2, M=4, d=1., reduced=reduced)
        Q = Observable(mesh, data=mesh.basis_function(ijk=(0,0), r=mesh.grid))
        Q.derive1()
        print("test_derive1_2D finished")