import pytest
import numpy as np

from mocca.mesh import LagrangeMesh
from mocca.mesh.observable import Observable


def test_Observable_ctor():
    mesh = LagrangeMesh(dim=1, M=4, d=1., reduced=False)
    Q = Observable(mesh, n_components=1)

    data = np.zeros((mesh.n_gridpoints(), 2), order='F')

