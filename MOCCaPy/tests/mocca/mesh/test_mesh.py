import pytest

from mocca.mesh import Cartesian3D

def test_Cartesian3D_constructor():
    mesh = Cartesian3D(n=5, d=.5)
    assert mesh.n == (5, 5, 5)

    mesh = Cartesian3D(n=(3,4,5), d=.5)
    assert mesh.n == (3, 4, 5)

    with pytest.raises(AssertionError):
        mesh = Cartesian3D(n=(3,4), d=.5)

    mesh = Cartesian3D(n=5, d=5)
    assert mesh.d == (5., 5., 5.)

    mesh = Cartesian3D(n=5, d=.5)
    assert mesh.d == (.5, .5, .5)

    mesh = Cartesian3D(n=5, d=(.6,.8,1))
    assert mesh.d == (.6, .8, 1.0)
