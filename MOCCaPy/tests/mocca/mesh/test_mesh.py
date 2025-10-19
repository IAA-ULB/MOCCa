import pytest

from mocca.mesh import LagrangeMesh

def test_LagrangeMesh_constructor():
    # valid 1D mesh constructions
    mesh = LagrangeMesh(dim=1, n=5, d=.5,)
    assert mesh.dim == 1
    assert mesh.n == (5,)
    assert mesh.d == (.5,)

    mesh = LagrangeMesh(n=(5,), d=.5,)
    assert mesh.dim == 1
    assert mesh.n == (5,)
    assert mesh.d == (.5,)

    # valid 2D mesh constructions
    mesh = LagrangeMesh(dim=2, n=5, d=.5,)
    assert mesh.dim == 2
    assert mesh.n == (5,5)
    assert mesh.d == (.5,.5)

    mesh = LagrangeMesh(dim=2, n=5, d=(.5,.6))
    assert mesh.dim == 2
    assert mesh.n == (5,5)
    assert mesh.d == (.5,.6)

    mesh = LagrangeMesh(n=(5,6), d=.5)
    assert mesh.dim == 2
    assert mesh.n == (5,6)
    assert mesh.d == (.5,.5)

    mesh = LagrangeMesh(n=(5,6), d=(.5,.6))
    assert mesh.dim == 2
    assert mesh.n == (5,6)
    assert mesh.d == (.5,.6)

    # valid 3D mesh constructions
    mesh = LagrangeMesh(dim=3, n=5, d=.5)
    assert mesh.dim == 3
    assert mesh.n == (5,5,5)
    assert mesh.d == (.5,.5,.5)

    mesh = LagrangeMesh(dim=3, n=(4,5,6), d=.5,)
    assert mesh.dim == 3
    assert mesh.n == (4,5,6)
    assert mesh.d == (.5,.5,.5)

    mesh = LagrangeMesh(dim=3, n=5, d=(.4,.5,.6))
    assert mesh.dim == 3
    assert mesh.n == (5,5,5)
    assert mesh.d == (.4,.5,.6)

    mesh = LagrangeMesh(n=(4,5,6), d=(.4,.5,.6))
    assert mesh.dim == 3
    assert mesh.n == (4,5,6)
    assert mesh.d == (.4,.5,.6)

    mesh = LagrangeMesh(n=(4,5,6), d=.5)
    assert mesh.dim == 3
    assert mesh.n == (4,5,6)
    assert mesh.d == (.5,.5,.5)

    # Invalid 1D mesh constructions
    with pytest.raises(AssertionError):
        mesh = LagrangeMesh(n=5, d=.5)
        # although this could be interpreted as a 1D mesh, the user might expect a 3D mesh.
        # The user must either provide dim=1 or n=tuple(5,)

    with pytest.raises(AssertionError):
        mesh = LagrangeMesh(dim=1, n=5, d=(.5,.5))

    with pytest.raises(AssertionError):
        mesh = LagrangeMesh(dim=1, n=(5,6), d=.5)

    # Invalid 2D mesh constructions
    with pytest.raises(AssertionError):
        mesh = LagrangeMesh(n=(4,5), d=(.4,.5,.6))

    with pytest.raises(AssertionError):
        mesh = LagrangeMesh(dim=2, n=(4,5,6), d=(.4,.5))

    with pytest.raises(AssertionError):
        mesh = LagrangeMesh(dim=2, n=5, d=(.4,.5,.6))

    # Invalid 3D mesh constructions
    with pytest.raises(AssertionError):
        mesh = LagrangeMesh(n=(4,5), d=(.4,.5,.6))

    with pytest.raises(AssertionError):
        mesh = LagrangeMesh(dim=3, n=5, d=(.5,.6))

    # Invalid 4D mesh constructions
    with pytest.raises(AssertionError):
        mesh = LagrangeMesh(dim=4, n=5, d=.5)

    with pytest.raises(AssertionError):
        mesh = LagrangeMesh(n=(4,5,6,7), d=.5)