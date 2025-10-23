import pytest
import numpy as np

from mocca.mesh import LagrangeMesh

def test_LagrangeMesh_ctor_n_d():
    """Test valid and invalid n and d parameters."""

    # valid 1D mesh constructions
    mesh = LagrangeMesh(dim=1, n=6, d=.5,)
    assert mesh.dim == 1
    assert mesh.n == (6,)
    assert mesh.d == (.5,)

    mesh = LagrangeMesh(n=(6,), d=.5,)
    assert mesh.dim == 1
    assert mesh.n == (6,)
    assert mesh.d == (.5,)

    # valid 2D mesh constructions
    mesh = LagrangeMesh(dim=2, n=6, d=.5,)
    assert mesh.dim == 2
    assert mesh.n == (6,6)
    assert mesh.d == (.5,.5)

    mesh = LagrangeMesh(dim=2, n=6, d=(.5,.6))
    assert mesh.dim == 2
    assert mesh.n == (6,6)
    assert mesh.d == (.5,.6)

    mesh = LagrangeMesh(n=(4,6), d=.5)
    assert mesh.dim == 2
    assert mesh.n == (4,6)
    assert mesh.d == (.5,.5)

    mesh = LagrangeMesh(n=(4,6), d=(.5,.6))
    assert mesh.dim == 2
    assert mesh.n == (4,6)
    assert mesh.d == (.5,.6)

    # valid 3D mesh constructions
    mesh = LagrangeMesh(dim=3, n=4, d=.5)
    assert mesh.dim == 3
    assert mesh.n == (4,4,4)
    assert mesh.d == (.5,.5,.5)

    mesh = LagrangeMesh(dim=3, n=(2,4,6), d=.5,)
    assert mesh.dim == 3
    assert mesh.n == (2,4,6)
    assert mesh.d == (.5,.5,.5)

    mesh = LagrangeMesh(dim=3, n=4, d=(.4,.5,.6))
    assert mesh.dim == 3
    assert mesh.n == (4,4,4)
    assert mesh.d == (.4,.5,.6)

    mesh = LagrangeMesh(n=(2,4,6), d=(.4,.5,.6))
    assert mesh.dim == 3
    assert mesh.n == (2,4,6)
    assert mesh.d == (.4,.5,.6)

    mesh = LagrangeMesh(n=(2,4,6), d=.5)
    assert mesh.dim == 3
    assert mesh.n == (2,4,6)
    assert mesh.d == (.5,.5,.5)

    # Invalid 1D mesh constructions
    with pytest.raises(AssertionError):
        mesh = LagrangeMesh(n=4, d=.5)
        # although this could be interpreted as a 1D mesh, the user might expect a 3D mesh.
        # The user must either provide dim=1 or n=tuple(5,)

    with pytest.raises(AssertionError):
        mesh = LagrangeMesh(dim=1, n=4, d=(.5,.5))

    with pytest.raises(AssertionError):
        mesh = LagrangeMesh(dim=1, n=(4,6), d=.5)

    # Invalid 2D mesh constructions
    with pytest.raises(AssertionError):
        mesh = LagrangeMesh(n=(4,6), d=(.4,.5,.6))

    with pytest.raises(AssertionError):
        mesh = LagrangeMesh(dim=2, n=(2,4,6), d=(.4,.5))

    with pytest.raises(AssertionError):
        mesh = LagrangeMesh(dim=2, n=4, d=(.4,.5,.6))

    # Invalid 3D mesh constructions
    with pytest.raises(AssertionError):
        mesh = LagrangeMesh(n=(4,6), d=(.4,.5,.6))

    with pytest.raises(AssertionError):
        mesh = LagrangeMesh(dim=3, n=4, d=(.5,.6))

    # Invalid 4D mesh constructions
    with pytest.raises(AssertionError):
        mesh = LagrangeMesh(dim=4, n=4, d=.5)

    with pytest.raises(AssertionError):
        mesh = LagrangeMesh(n=(2,4,6,8), d=.5)

    # invalid n
    with pytest.raises(AssertionError):
        mesh = LagrangeMesh(dim=1, n=3, d=.5) # odd n

    with pytest.raises(AssertionError):
        mesh = LagrangeMesh(dim=1, n=0, d=.5) # n == 0

    with pytest.raises(AssertionError):
        mesh = LagrangeMesh(dim=1, n=2, d=-.5) # negative d

    with pytest.raises(AssertionError):
        mesh = LagrangeMesh(dim=2, n=2, d=(-.5,.5)) # negative d

    with pytest.raises(AssertionError):
        mesh = LagrangeMesh(dim=1, n=2, d=0)  # negative d

    with pytest.raises(AssertionError):
        mesh = LagrangeMesh(dim=2, n=2, d=(0,.5)) # negative d


def test_LagrangeMesh_ctor_shift():
    """Test valid and invalid reduce parameter."""

    mesh = LagrangeMesh(dim=3, n=4, d=.5)
    assert mesh.shift == (.0, .0, .0)

    mesh = LagrangeMesh(dim=3, n=4, d=.5, shift=.1, reduce=False)
    assert mesh.shift == (.1, .1, .1)

    mesh = LagrangeMesh(dim=3, n=4, d=.5, shift=(.1,.2,.3), reduce=False)
    assert mesh.shift == (.1, .2, .3)


def test_LagrangeMesh_ctor_reduce():
    """Test valid and invalid reduce parameter."""

    mesh = LagrangeMesh(dim=3, n=4, d=.5, reduce=True)
    assert mesh.reduce == (True, True, True )

    mesh = LagrangeMesh(dim=3, n=4, d=.5, reduce=False)
    assert mesh.reduce == (False, False, False)

    mesh = LagrangeMesh(dim=3, n=4, d=.5, reduce=(True,True,False))
    assert mesh.reduce == (True, True, False)

def test_LagrangeMesh_ctor_bc():
    """Test valid and invalid bc parameter."""
    mesh = LagrangeMesh(dim=3, n=4, d=.5)
    assert mesh.antiperiodic
    assert not mesh.periodic

    mesh = LagrangeMesh(dim=3, n=4, d=.5, bc='antiperiodic')
    assert mesh.antiperiodic
    assert not mesh.periodic

    mesh = LagrangeMesh(dim=3, n=4, d=.5, bc='periodic')
    assert not mesh.antiperiodic
    assert mesh.periodic


def test_LagrangeMesh_ctor_grid_1D():
    """Test grid points for 1D grid."""
    mesh = LagrangeMesh(dim=1, n=6, d=1.)
    assert np.all(mesh.gridx == [0.5,1.5,2.5])

    mesh = LagrangeMesh(dim=1, n=6, d=.5)
    assert np.all(mesh.gridx == [.25,.75,1.25])

    mesh = LagrangeMesh(dim=1, n=6, d=1., reduce=False)
    assert np.all(mesh.gridx == [-2.5, -1.5,-0.5, 0.5,1.5,2.5])

    mesh = LagrangeMesh(dim=1, n=6, d=1., reduce=False, shift=0.5)
    assert np.all(mesh.gridx == [-3., -2., -1., 0., 1., 2.])

    mesh = LagrangeMesh(dim=1, n=6, d=.5, reduce=False)
    assert np.all(mesh.gridx == [-1.25, -.75, -.25, .25, .75, 1.25])

    mesh = LagrangeMesh(dim=1, n=6, d=.5, reduce=False, shift=0.5)
    assert np.all(mesh.gridx == [-1.75, -1.25, -.75, -.25, .25, .75])

def test_LagrangeMesh_ctor_grid_2D():
    """Test grid points for 2D grid."""
    mesh = LagrangeMesh(dim=2, n=6, d=1.)
    for j in range(3):
        assert np.all(mesh.gridx[:,j] == [0.5,1.5,2.5])
    for i in range(3):
        assert np.all(mesh.gridy[i,:] == [0.5,1.5,2.5])

def test_LagrangeMesh_ctor_grid_3D():
    """Test grid points for 2D grid."""
    mesh = LagrangeMesh(dim=3, n=4, d=1.)
    for k in range(2):
        for j in range(2):
            assert np.all(mesh.gridx[:,j,k] == [0.5,1.5])
    for k in range(2):
        for i in range(2):
            assert np.all(mesh.gridy[i,:,k] == [0.5,1.5])
    for j in range(2):
        for i in range(2):
            assert np.all(mesh.gridz[i,j,:] == [0.5,1.5])

def test_LagrangeMesh_ctor_dv():
    """Test dv computation."""
    mesh = LagrangeMesh(dim=3, n=6, d=.5)
    assert mesh.dv == .5**3 * 2**3

    mesh = LagrangeMesh(dim=3, n=6, d=.5, reduce=(True,True,False))
    assert mesh.dv == .5 ** 3 * 2 ** 2

    mesh = LagrangeMesh(dim=3, n=6, d=.5, reduce=(True, False, False))
    assert mesh.dv == .5 ** 3 * 2

    mesh = LagrangeMesh(dim=3, n=6, d=.5, reduce=False)
    assert mesh.dv == .5 ** 3


def test_LagrangeMesh_ctor_box_width():
    """Test dv computation."""
    mesh = LagrangeMesh(dim=3, n=6, d=.5)
    assert mesh.box_width == tuple(3*[6*.5])


def test_LagrangeMesh_plane_wave():
    mesh = LagrangeMesh(dim=3, n=6, d=.5, reduce=False)

    bw = mesh.box_width
    oneoversqrtbw = np.sqrt(1 / np.prod(mesh.box_width))
    twopij = 2j * np.pi
    k = np.array([0.5,1.5,2.5])
    k /= bw
    rng = np.random.default_rng()
    r = (rng.random((5,3))*2 - 1)*bw[0]
    expected = oneoversqrtbw*np.exp(twopij*r@k)

    pwc = mesh.plane_wave(k, r)
    assert np.all(pwc == expected)


def test_LagrangeMesh_reshape():
    """Test reshaping of the mesh."""
    mesh = LagrangeMesh(dim=3, n=4, d=.5)
    assert mesh.shape == (2,2,2)
    assert mesh.flat_shape == (8,)

    assert mesh.gridx.shape == (2,2,2)
    mesh.flatten()
    assert mesh.gridx.shape == (8,)
    mesh.flatten() # already flat
    assert mesh.gridx.shape == (8,)
    mesh.unflatten()
    assert mesh.gridx.shape == (2, 2, 2)
    mesh.unflatten() # already unflattened
    assert mesh.gridx.shape == (2, 2, 2)

    a = np.ones((2,2,2,5,5))
    a_flat = mesh.flatten(a)
    assert a_flat.shape == (8,5,5)
    a = mesh.unflatten(a_flat)
    assert a.shape == (2,2,2,5,5)

from numba import vectorize, float64

@vectorize([float64(float64, float64, float64)])
def norm(x, y, z):
    return np.sqrt(x*x + y*y + z*z)

def test_LagrangeMesh_apply():
    mesh = LagrangeMesh(dim=3, n=10, d=.5)
    r_xyz = mesh.apply(norm)
    for i in range(mesh.flat_shape[0]):
        assert r_xyz[i] == np.sqrt(mesh.gridx[i]**2 + mesh.gridy[i]**2 + mesh.gridz[i]**2)

def test_LagrangeMesh_integrate():
    mesh = LagrangeMesh(dim=3, n=10, d=.5)
    # integrate a constant function.
    c = np.ones_like(mesh.gridx)
    c = mesh.flatten(c)
    print(f"{c=}")
    integral_of_c = mesh.integrate(c)
    assert integral_of_c == len(c) * mesh.dv