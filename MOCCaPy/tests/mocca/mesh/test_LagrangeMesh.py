import inspect
import matplotlib.pyplot as plt
from   matplotlib import cm
import numpy as np
from   pathlib import Path
import sys
import pytest

from mocca.mesh import LagrangeMesh
from mocca.mesh.lagrange import create_mesh
from mocca.mesh.mesh_quantity import MeshQuantity
from mocca.util.timer import Timer

path2MOCCaPy = Path(__file__).parent
while not path2MOCCaPy.name == 'MOCCaPy':
    path2MOCCaPy = path2MOCCaPy.parent
    print(path2MOCCaPy)
sys.path.insert(0, str(path2MOCCaPy))

from tests.util import started_finished, started, finished

this_file = Path(__file__).resolve()
project_folder = this_file
while project_folder.name != 'tantalus_full':
    project_folder = project_folder.parent
test_folder = project_folder/f"MOCCaPy/tests/mocca/mesh"
png_folder = test_folder/"png"/this_file.stem
png_folder.mkdir(exist_ok=True)

# remove all .png files
for png in png_folder.glob('*.png'):
    png.unlink()


@started_finished
def test_LagrangeMesh_ctor_n_d():
    """Test valid and invalid n and d parameters."""

    # valid 1D mesh constructions
    mesh = LagrangeMesh(dim=1, M=6, d=.5,)
    assert mesh.dim == 1
    assert mesh.M == (6,)
    assert mesh.d == (.5,)

    mesh = LagrangeMesh(M=(6,), d=.5,)
    assert mesh.dim == 1
    assert mesh.M == (6,)
    assert mesh.d == (.5,)

    # valid 2D mesh constructions
    mesh = LagrangeMesh(dim=2, M=6, d=.5,)
    assert mesh.dim == 2
    assert (mesh.M == np.array((6,6))).all()
    assert (mesh.d == np.array([.5,.5])).all()

    mesh = LagrangeMesh(dim=2, M=6, d=(.5,.6))
    assert mesh.dim == 2
    assert mesh.M == (6,6)
    assert mesh.d == (.5,.6)

    mesh = LagrangeMesh(M=(4,6), d=.5)
    assert mesh.dim == 2
    assert mesh.M == (4,6)
    assert mesh.d == (.5,.5)

    mesh = LagrangeMesh(M=(4,6), d=(.5,.6))
    assert mesh.dim == 2
    assert mesh.M == (4,6)
    assert mesh.d == (.5,.6)

    # valid 3D mesh constructions
    mesh = LagrangeMesh(dim=3, M=4, d=.5)
    assert mesh.dim == 3
    assert mesh.M == (4,4,4)
    assert mesh.d == (.5,.5,.5)

    mesh = LagrangeMesh(dim=3, M=(2,4,6), d=.5,)
    assert mesh.dim == 3
    assert mesh.M == (2,4,6)
    assert mesh.d == (.5,.5,.5)

    mesh = LagrangeMesh(dim=3, M=4, d=(.4,.5,.6))
    assert mesh.dim == 3
    assert mesh.M == (4,4,4)
    assert mesh.d == (.4,.5,.6)

    mesh = LagrangeMesh(M=(2,4,6), d=(.4,.5,.6))
    assert mesh.dim == 3
    assert mesh.M == (2,4,6)
    assert mesh.d == (.4,.5,.6)

    mesh = LagrangeMesh(M=(2,4,6), d=.5)
    assert mesh.dim == 3
    assert mesh.M == (2,4,6)
    assert mesh.d == (.5,.5,.5)

    mesh = LagrangeMesh(M=4, d=.5) # 3D mesh is default
    assert mesh.dim == 3
    assert mesh.M == (4,4,4)
    assert mesh.d == (.5,.5,.5)


    # Invalid 1D mesh constructions
    with pytest.raises(ValueError):
        mesh = LagrangeMesh(dim=1, M=4, d=(.5,.5))

    with pytest.raises(ValueError):
        mesh = LagrangeMesh(dim=1, M=(4,6), d=.5)

    # Invalid 2D mesh constructions
    with pytest.raises(ValueError):
        mesh = LagrangeMesh(M=(4,6), d=(.4,.5,.6))

    with pytest.raises(ValueError):
        mesh = LagrangeMesh(dim=2, M=(2,4,6), d=(.4,.5))

    with pytest.raises(ValueError):
        mesh = LagrangeMesh(dim=2, M=4, d=(.4,.5,.6))

    # Invalid 3D mesh constructions
    with pytest.raises(ValueError):
        mesh = LagrangeMesh(M=(4,6), d=(.4,.5,.6))

    with pytest.raises(ValueError):
        mesh = LagrangeMesh(dim=3, M=4, d=(.5,.6))

    # Invalid 4D mesh constructions
    with pytest.raises(ValueError):
        mesh = LagrangeMesh(dim=4, M=4, d=.5)

    with pytest.raises(ValueError):
        mesh = LagrangeMesh(M=(2,4,6,8), d=.5)

    # invalid N
    with pytest.raises(ValueError):
        mesh = LagrangeMesh(dim=1, M=3, d=.5) # odd n

    with pytest.raises(ValueError):
        mesh = LagrangeMesh(dim=1, M=0, d=.5) # n == 0

    with pytest.raises(ValueError):
        mesh = LagrangeMesh(dim=1, M=2, d=-.5) # negative d

    with pytest.raises(ValueError):
        mesh = LagrangeMesh(dim=2, M=2, d=(-.5,.5)) # negative d

    with pytest.raises(ValueError):
        mesh = LagrangeMesh(dim=1, M=2, d=0)  # negative d

    with pytest.raises(ValueError):
        mesh = LagrangeMesh(dim=2, M=2, d=(0,.5)) # negative d


@started_finished
def test_LagrangeMesh_ctor_shift():
    """Test valid and invalid reduced parameter."""

    mesh = LagrangeMesh(dim=3, M=4, d=.5)
    assert mesh.shift == (.0, .0, .0)

    mesh = LagrangeMesh(dim=3, M=4, d=.5, shift=.1, reduced=False)
    assert mesh.shift == (.1, .1, .1)

    mesh = LagrangeMesh(dim=3, M=4, d=.5, shift=(.1,.2,.3), reduced=False)
    assert mesh.shift == (.1, .2, .3)


@started_finished
def test_LagrangeMesh_ctor_reduced():
    """Test valid and invalid reduced parameter."""

    mesh = LagrangeMesh(dim=3, M=4, d=.5, reduced=True)
    assert mesh.reduced == (True, True, True )

    mesh = LagrangeMesh(dim=3, M=4, d=.5, reduced=False)
    assert mesh.reduced == (False, False, False)

    mesh = LagrangeMesh(dim=3, M=4, d=.5, reduced=(True,True,False))
    assert mesh.reduced == (True, True, False)

@started_finished
def test_LagrangeMesh_ctor_bc():
    """Test valid and invalid bc parameter."""
    mesh = LagrangeMesh(dim=3, M=4, d=.5)
    assert mesh.antiperiodic
    assert not mesh.periodic

    mesh = LagrangeMesh(dim=3, M=4, d=.5, bc='antiperiodic')
    assert mesh.antiperiodic
    assert not mesh.periodic

    mesh = LagrangeMesh(dim=3, M=4, d=.5, bc='periodic')
    assert not mesh.antiperiodic
    assert mesh.periodic


@started_finished
def test_LagrangeMesh_ctor_grid_1D():
    """Test grid points for 1D grid."""
    mesh = LagrangeMesh(dim=1, M=6, d=1.)
    assert np.all(mesh.gridx == [0.5,1.5,2.5])

    mesh = LagrangeMesh(dim=1, M=6, d=.5)
    assert np.all(mesh.gridx == [.25,.75,1.25])

    mesh = LagrangeMesh(dim=1, M=6, d=1., reduced=False)
    assert np.all(mesh.gridx == [-2.5, -1.5,-0.5, 0.5,1.5,2.5])

    mesh = LagrangeMesh(dim=1, M=6, d=1., reduced=False, shift=0.5)
    assert np.all(mesh.gridx == [-3., -2., -1., 0., 1., 2.])

    mesh = LagrangeMesh(dim=1, M=6, d=.5, reduced=False)
    assert np.all(mesh.gridx == [-1.25, -.75, -.25, .25, .75, 1.25])

    mesh = LagrangeMesh(dim=1, M=6, d=.5, reduced=False, shift=0.5)
    assert np.all(mesh.gridx == [-1.75, -1.25, -.75, -.25, .25, .75])


@started_finished
def test_LagrangeMesh_ctor_grid_2D():
    """Test grid points for 2D grid."""
    mesh = LagrangeMesh(dim=2, M=6, d=1.)
    for j in range(3):
        assert np.all(mesh.gridx[:,j] == [0.5,1.5,2.5])
    for i in range(3):
        assert np.all(mesh.gridy[i,:] == [0.5,1.5,2.5])

\
@started_finished
def test_LagrangeMesh_ctor_grid_3D():
    """Test grid points for 2D grid."""
    mesh = LagrangeMesh(dim=3, M=4, d=1.)
    for k in range(2):
        for j in range(2):
            assert np.all(mesh.gridx[:,j,k] == [0.5,1.5])
    for k in range(2):
        for i in range(2):
            assert np.all(mesh.gridy[i,:,k] == [0.5,1.5])
    for j in range(2):
        for i in range(2):
            assert np.all(mesh.gridz[i,j,:] == [0.5,1.5])

@started_finished
def test_LagrangeMesh_ctor_dv():
    """Test dv computation."""
    mesh = LagrangeMesh(dim=3, M=6, d=.5)
    assert mesh.dv == .5**3 * 2**3

    mesh = LagrangeMesh(dim=3, M=6, d=.5, reduced=(True,True,False))
    assert mesh.dv == .5 ** 3 * 2 ** 2

    mesh = LagrangeMesh(dim=3, M=6, d=.5, reduced=(True, False, False))
    assert mesh.dv == .5 ** 3 * 2

    mesh = LagrangeMesh(dim=3, M=6, d=.5, reduced=False)
    assert mesh.dv == .5 ** 3


@started_finished
def test_LagrangeMesh_ctor_box_width():
    """Test dv computation."""
    mesh = LagrangeMesh(dim=3, M=6, d=.5)
    assert (mesh.box_width == tuple(3*[6*.5])).all()


@started_finished
def test_LagrangeMesh_plane_wave_1D():
    d = .5
    n = 6
    mesh = LagrangeMesh(dim=3, M=6, d=d, reduced=False)
    k = 1.5
    L = n*d
    r = mesh.gridx[:,0,0]
    pw = mesh.plane_wave(L=L, k=k, r=r)

    arg = (2*np.pi * k / L) * r
    real_part = np.sqrt(1/L) * np.cos(arg)
    imag_part = np.sqrt(1/L) * np.sin(arg)
    assert (pw[:,0] == real_part).all()
    assert (pw[:,1] == imag_part).all()


@started_finished
def test_LagrangeMesh_reshape():
    """Test reshaping of the mesh."""
    mesh = LagrangeMesh(dim=3, M=4, d=.5)
    assert mesh.mesh_shape == (2,2,2)
    assert mesh.linear_size == 8

    assert mesh.gridx.shape == (2,2,2)
    mesh.gridx = mesh.cast2linear(mesh.gridx)
    assert mesh.gridx.shape == (8,1)
    mesh.gridx = mesh.cast2linear(mesh.gridx) # already linearized
    assert mesh.gridx.shape == (8,1)
    mesh.gridx = mesh.cast2grid(mesh.gridx)
    assert mesh.gridx.shape == (2,2,2,1)
    mesh.gridx = mesh.cast2grid(mesh.gridx)   # already to mesh
    assert mesh.gridx.shape == (2,2,2,1)

    a = np.ones((2,2,2,5))
    a_linear = mesh.cast2linear(a)
    assert a_linear.shape == (8,5)
    a = mesh.cast2grid(a_linear)
    assert a.shape == (2,2,2,5)


def test_LagrangeMesh_basis_function_1D(no_plot):
    """"""
    started(inspect.stack()[0][3])
    # 1D mesh
    N = 6
    d = 1
    mesh = LagrangeMesh(M=(N,), d=d, reduced=False)
    # r = mesh.gx
    r = np.linspace(-N*d/2, N*d/2, num=241)

    if no_plot:
        pass
    else:
        for i in range(N):
            pw = mesh.basis_function(i,r)
            fig, ax = plt.subplots()
            plt.title(f'{i=}, x_i={mesh.g1D[0][i]}')
            ax.plot([-N/2,N/2],[0,0],color='k')
            ax.plot([0,0],[-1,1],color='k')
            ax.plot(mesh.g1D[0],  N*[0], 'rs')
            ax.plot(mesh.g1D[0][i], [0], 'gs')
            ax.plot(r, pw[:,0], label='real')
            ax.plot(r, pw[:,1], label='imag')
            plt.legend()
            fig.savefig(png_folder/f"1D_basis_function_{i}.png")
            plt.close(fig)

    finished(inspect.stack()[0][3])

def test_LagrangeMesh_basis_function_2D(no_plot):
    started(inspect.stack()[0][3])
    # 2D mesh
    N = 6
    d = 1
    mesh = LagrangeMesh(dim=2, M=N, d=d, reduced=False)
    # r = mesh.gx
    nr1 = 201
    y = x = np.linspace(-N*d/2, N*d/2, num=nr1)
    r = np.empty((nr1**2,2), dtype=float, order='F')
    for i in range(nr1):
        r[i*nr1:(i+1)*nr1, 0] = x
        r[i*nr1:(i+1)*nr1, 1] = y[i]

    x = r[:,0].reshape((nr1,nr1), order='F')
    y = r[:,1].reshape((nr1,nr1), order='F')
    for i in range(N):
        for j in range(N):
            ij = (i,j)
            pw = mesh.basis_function(ij,r)
            pw_real = pw[:,0].reshape((nr1,nr1), order='F')
            pw_imag = pw[:,1].reshape((nr1,nr1), order='F')

            if no_plot:
                pass
            else:
                fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
                plt.title(f'{i=}, x_i={mesh.g1D[0][i]}, y_j={mesh.g1D[1][j]} real')
                ax.plot_surface(x, y, pw_real, cmap=cm.coolwarm, linewidth=0, antialiased=False)
                fig.savefig(png_folder/f"2D_basis_function_{i}_real.png")
                # plt.show()
                plt.close(fig)

                fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
                plt.title(f'{i=}, x_i={mesh.g1D[0][i]}, y_j={mesh.g1D[1][j]} imag')
                ax.plot_surface(x, y, pw_imag, cmap=cm.coolwarm, linewidth=0, antialiased=False)
                fig.savefig(png_folder/f"2D_basis_function_{i}_imag.png")
                # plt.show()
                plt.close(fig)

    finished(inspect.stack()[0][3])


from numba import vectorize, float64

@vectorize([float64(float64, float64, float64)])
def norm(x, y, z):
    return np.sqrt(x*x + y*y + z*z)


@started_finished
def test_LagrangeMesh_apply():
    mesh = LagrangeMesh(dim=3, M=10, d=.5)
    r_xyz = mesh.apply(norm)
    for i in range(mesh.linear_size):
        assert r_xyz[i] == np.sqrt(mesh.gridx[i]**2 + mesh.gridy[i]**2 + mesh.gridz[i]**2)


@started_finished
def test_LagrangeMesh_integrate():
    mesh = LagrangeMesh(dim=3, M=10, d=.5)
    # integrate a constant function.

    q = np.ones((mesh.linear_size,))
    with pytest.raises(UserWarning):
        Q = MeshQuantity(mesh, data=q)
        integral_of_Q = Q.integrate()
        assert integral_of_Q == mesh.linear_size * mesh.dv


def test_LagrangeMesh_interpolate1D(no_plot, debug=False):
    started(inspect.stack()[0][3])

    d = 1.
    N = 3
    for reduced in [
        False,
        True,
    ]:
        str_reduced = "(reduced)" if reduced else ""
        mesh  = LagrangeMesh(dim=1, M=2*N, d=d, reduced=reduced)
        r = np.empty((mesh.linear_size,), dtype=float, order='F')
        r = mesh.g1D[0]
        rip = np.linspace(-N*d, N*d, num=61)
        for ip in range(len(rip)):
            if (reduced and rip[ip] < 0 and -rip[ip] in mesh.g1D[0]) or (rip[ip] in mesh.g1D[0]):
                rip[ip] += 1e-12 # avoid 0/0 in the Lagrange Functions

        for ii in range(2*N):
            if not mesh.reduced[0]:
                i = ii
                print(f"{reduced=} bf={i}/{2*N}")
            else:
                i = ii//2
                if ii % 2 == 0:
                    print(f"{reduced=} bf={i}/{2 * N}")
                else:
                    print(f"{reduced=} bf=-{i}/{2 * N}")

            lcpw_rgp = mesh.basis_function(i,r  )
            lcpw_rip = mesh.basis_function(i,rip)

            Q = MeshQuantity( mesh, data=lcpw_rgp, symmetry=[[1],[-1]]) # real/imag component is symmetric/skew-symmetric

            Qrip = mesh.interpolate(Q, rip)

            if no_plot:
                pass
            else:
                fig, ax = plt.subplots()
                plt.title(f"test_LagrangeMesh_interpolate1D real part $x_{i}$")
                plt.plot(rip,lcpw_rip[:,0], 'o', label=f'plane wave $x_{i}$ real')
                plt.plot(rip,    Qrip[:,0],      label=f'interpolated $x_{i}$ real')
                plt.legend()
                plt.savefig(png_folder/f"test_LagrangeMesh_interpolate1D_{i}_real{str_reduced}")
                if debug:
                    plt.show()
                plt.close(fig)

                fig, ax = plt.subplots()
                plt.title(f"test_LagrangeMesh_interpolate1D imag part $x_{i}$")
                plt.plot(rip,lcpw_rip[:,1], 'o', label=f'plane wave $x_{i}$ imag')
                plt.plot(rip,    Qrip[:,1],      label=f'interpolated $x_{i}$ imag')
                plt.legend()
                plt.savefig(png_folder/f"test_LagrangeMesh_interpolate1D_{i}_imag{str_reduced}")
                if debug:
                    plt.show()
                plt.close(fig)

            assert ((lcpw_rip - Qrip) < 1e-12).all()

    for lc in range(2):
        coeff = np.ones((2*N,),dtype=float) if lc == 0 else (
                np.random.rand(2*N)
        )
        s = "sum of bf" if lc == 0 else "random lc of bf"
        if not mesh.reduced[0]:
            lcpw_rgp = coeff[0] * mesh.basis_function(0, r  )
            lcpw_rip = coeff[0] * mesh.basis_function(0, rip)
        else:
            lcpw_rgp = coeff[0] * mesh.basis_function(0, r  ) + \
                       coeff[1] * mesh.basis_function(0, r  , ijk_sign=-1)
            lcpw_rip = coeff[0] * mesh.basis_function(0, rip) + \
                       coeff[1] * mesh.basis_function(0, rip, ijk_sign=-1)

        if not mesh.reduced[0]:
            for i in range(1,2*N):
                lcpw_rgp += coeff[i] * mesh.basis_function(i,r)
                lcpw_rip += coeff[i] * mesh.basis_function(i,rip)
        else:
            for i in range(1,N):
                lcpw_rgp += coeff[2*i  ] * mesh.basis_function(i,r  ) + \
                            coeff[2*i+1] * mesh.basis_function(i,r  , ijk_sign=-1)
                lcpw_rip += coeff[2*i  ] * mesh.basis_function(i,rip) + \
                            coeff[2*i+1] * mesh.basis_function(i,rip, ijk_sign=-1)

        Q = MeshQuantity(mesh, data=lcpw_rgp, symmetry=[[1],[-1]])
        Qrip = mesh.interpolate(Q, rip)

        if no_plot:
            pass
        else:
            fig, ax = plt.subplots()
            plt.title(f"{s}: real part")
            plt.plot(rip,lcpw_rip[:,0], 'o', label=f'plane wave $x_{i}$ real')
            plt.plot(rip,    Qrip[:,0],      label=f'interpolated $x_{i}$ real')
            plt.legend()
            plt.savefig(png_folder/f"test_LagrangeMesh_interpolate1D {s} real part")
            if debug:
                plt.show()
            plt.close(fig)

            fig, ax = plt.subplots()
            plt.title(f"{s}: imag part")
            ax.set_ylim([-1, 1])
            plt.plot(rip,lcpw_rip[:,1], 'o', label=f'plane wave $x_{i}$ real')
            plt.plot(rip,    Qrip[:,1],      label=f'interpolated $x_{i}$ real')
            plt.legend()
            plt.savefig(png_folder/f"test_LagrangeMesh_interpolate1D {s} imag part")
            if debug:
                plt.show()
            plt.close(fig)

        assert ((lcpw_rip - Qrip) < 1e-12).all()

    finished(inspect.stack()[0][3])


def test_LagrangeMesh_interpolate2D_basisfunction(no_plot, debug=False):
    """"""
    started(inspect.stack()[0][3])
    d = 1.
    N = 2
    nip = 61
    r = np.linspace(-2,2,num=nip)
    xy = create_mesh([r, r])
    xy += 1e-9
    x = xy[:,0].reshape((nip,nip), order='F')
    y = xy[:,1].reshape((nip,nip), order='F')

    for reduced in [
                   (False, False),
                   # (False, True), This case does not work for basis functions (see lagrange.md)
                   # (True, False), This case does not work for basis functions (see lagrange.md)
                   # (True, True ), This case does not work for basis functions (see lagrange.md)
                   ]:
        mesh = LagrangeMesh(dim=2, M=2*N, d=d, reduced=reduced)

        for ibfx in range(2*N):
            for ibfy in range(2*N):
                # values of the basis function at the grid points (with a very small offset)
                bf_xy = mesh.basis_function((ibfx,ibfy), xy) # the interpolated quantity
                bf_real = bf_xy[:, 0].reshape((nip, nip), order='F')
                bf_imag = bf_xy[:, 1].reshape((nip, nip), order='F')

                if no_plot:
                    pass
                else:
                    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
                    plt.title(f'BF imag x_ij=({mesh.g1D[0][ibfx]}, {mesh.g1D[1][ibfy]})')
                    ax.plot_surface(x, y, bf_real, cmap=cm.coolwarm, linewidth=0, antialiased=False)
                    fig.savefig(png_folder/"basis_function_({ibfx},{ibfy}).png")
                    if debug:
                        plt.show()
                    plt.close(fig)

                    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
                    plt.title(f'BF imag x_ij=({mesh.g1D[0][ibfx]}, {mesh.g1D[1][ibfy]})')
                    ax.plot_surface(x, y, bf_imag, cmap=cm.coolwarm, linewidth=0, antialiased=False)
                    fig.savefig(png_folder/"BF_real({ibfx},{ibfy}).png")
                    if debug:
                        plt.show()
                    plt.close(fig)

                # compute the values of the basis function at the grid points
                bf_xy_g  = mesh.basis_function((ibfx,ibfy), mesh.grid)
                # wrap the interpolated quantity in a MeshQuantity
                Q = MeshQuantity(mesh, data=bf_xy_g, symmetry=[1,-1])

                # interpolate Q at xy
                Qxy = mesh.interpolate(Q, xy)
                Qxy_real = Qxy[:, 0].reshape((nip, nip), order='F')
                Qxy_imag = Qxy[:, 1].reshape((nip, nip), order='F')


                if no_plot:
                    pass
                else:
                    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
                    plt.title(f'BF real x_ij={mesh.g1D[0][ibfx]}, y_j={mesh.g1D[1][ibfy]} interpolated')
                    ax.plot_surface(x, y, Qxy_real, cmap=cm.coolwarm, linewidth=0, antialiased=False)
                    fig.savefig(png_folder/"BF_real({ibfx},{ibfy})_interpolated.png")
                    if debug:
                        plt.show()
                    plt.close(fig)

                    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
                    plt.title(f'BF imag x_ij={mesh.g1D[0][ibfx]}, y_j={mesh.g1D[1][ibfy]} interpolated')
                    ax.plot_surface(x, y, Qxy_imag, cmap=cm.coolwarm, linewidth=0, antialiased=False)
                    fig.savefig(png_folder/"BF_real({ibfx},{ibfy})_interpolated.png")
                    if debug:
                        plt.show()
                    plt.close(fig)

                print(f"{reduced=} bf = ({ibfx}, {ibfy})")
                for i in range(len(Qxy_real)):
                    assert Qxy[i,0] == pytest.approx(bf_xy[i,0]), f"{i=} {Qxy[i,0]} != {bf_xy[i,0]}, |diff|={abs(Qxy[i,0]-bf_xy[i,0])}"
                    assert Qxy[i,1] == pytest.approx(bf_xy[i,1]), f"{i=} {Qxy[i,1]} != {bf_xy[i,1]}, |diff|={abs(Qxy[i,1]-bf_xy[i,1])}"
    finished(inspect.stack()[0][3])


def gauss(x, sigma):
    if x.shape[1] == 3:
        r2 = x[:,0]**2 + x[:,1]**2 + x[:,2]**2
    elif x.shape[1] == 2:
        r2 = x[:,0]**2 + x[:,1]**2
    else:
        r2 = x**2
    return np.exp(-0.5*r2/sigma**2)/(sigma*np.sqrt(2*np.pi))


def test_LagrangeMesh_interpolate2D_bell(no_plot, debug=False):
    """"""
    started(inspect.stack()[0][3])
    d = 1.
    N = 40
    nip = 61
    r = np.linspace(-2,2,num=nip)
    xy = create_mesh([r, r])
    xy += 1e-9
    x = xy[:,0].reshape((nip,nip), order='F')
    y = xy[:,1].reshape((nip,nip), order='F')

    sigma = 1.0
    bell = gauss(xy, sigma=sigma)
    bell_plot = bell.reshape((nip,nip), order='F')

    if no_plot:
        pass
    else:
        fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
        plt.title(f'bell')
        ax.plot_surface(x, y, bell_plot, cmap=cm.coolwarm, linewidth=0, antialiased=False)
        fig.savefig(png_folder/"bell.png")
        if debug:
            plt.show()
        plt.close(fig)

    error = {}
    for reduced in [
        (False, False),
        (False, True),
        (True, False),
        (True, True )
    ]:
        mesh = LagrangeMesh(dim=2, M=2*N, d=d, reduced=reduced)
        bell_g = gauss(mesh.grid, sigma=sigma)

        # wrap the interpolated quantity in a MeshQuantity
        Q = MeshQuantity(mesh, data=bell_g, symmetry=1)
        # interpolate Q at xy
        Qxy = mesh.interpolate(Q, xy)
        Qxy_plot = Qxy.reshape((nip, nip), order='F')

        if no_plot:
            pass
        else:
            fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
            plt.title(f'bell interpolated ({reduced=})')
            ax.plot_surface(x, y, Qxy_plot, cmap=cm.coolwarm, linewidth=0, antialiased=False)
            fig.savefig(png_folder/"bell_interpolated_({reduced=}).png")
            if debug:
                plt.show()
            plt.close(fig)

        err = np.abs(Qxy_plot-bell_plot)
        error[reduced] = err

        if no_plot:
            pass
        else:
            fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
            plt.title(f'bell interpolated error ({reduced=})')
            ax.plot_surface(x, y, err, cmap=cm.coolwarm, linewidth=0, antialiased=False)
            fig.savefig(png_folder/"bell_interpolated_error_({reduced=}).png")
            if debug:
                plt.show()
            plt.close(fig)

        # the interpolation is not particularly accurate, but all four interpolations give very comparable errors
    prev = None
    for k,v in error.items():
        if prev is not None:
            for i in range(v.shape[0]):
                assert v[i,0] == pytest.approx(prev[i,0])
        prev = v

    finished(inspect.stack()[0][3])


def test_LagrangeMesh_interpolate2D_bellx(no_plot, debug=False):
    """"""
    started(inspect.stack()[0][3])

    d = 1.
    N = 40
    nip = 61
    r = np.linspace(-2,2,num=nip)
    xy = create_mesh([r, r])
    xy += 1e-9
    x = xy[:,0].reshape((nip,nip), order='F')
    y = xy[:,1].reshape((nip,nip), order='F')

    sigma = 1.0
    bellx = gauss(xy, sigma=sigma) * xy[:,0]
    bellx_plot = bellx.reshape((nip,nip), order='F')

    if no_plot:
        pass
    else:
        fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
        plt.title(f'bellx')
        ax.plot_surface(x, y, bellx_plot, cmap=cm.coolwarm, linewidth=0, antialiased=False)
        fig.savefig(png_folder/"bellx.png")
        if debug:
            plt.show()
        plt.close(fig)

    error = {}
    for reduced in [
        (False, False),
        (False, True),
        (True, False),
        (True, True )
    ]:
        mesh = LagrangeMesh(dim=2, M=2*N, d=d, reduced=reduced)
        bellx_g = gauss(mesh.grid, sigma=sigma) * mesh.grid[:,0]

        # wrap the interpolated quantity in a MeshQuantity
        Q = MeshQuantity(mesh, data=bellx_g, symmetry=[[-1,1]])
        # interpolate Q at xy
        Qxy = mesh.interpolate(Q, xy)
        Qxy_plot = Qxy.reshape((nip, nip), order='F')

        if no_plot:
            pass
        else:
            fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
            plt.title(f'bellx interpolated ({reduced=})')
            ax.plot_surface(x, y, Qxy_plot, cmap=cm.coolwarm, linewidth=0, antialiased=False)
            fig.savefig(png_folder/"bellx_interpolated_({reduced=}).png")
            if debug:
                plt.show()
            plt.close(fig)

        err = np.abs(Qxy_plot-bellx_plot)
        error[reduced] = err

        if no_plot:
            pass
        else:
            fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
            plt.title(f'bellx interpolated error ({reduced=})')
            ax.plot_surface(x, y, err, cmap=cm.coolwarm, linewidth=0, antialiased=False)
            fig.savefig(png_folder/"bellx_interpolated_error_({reduced=}).png")
            if debug:
                plt.show()
            plt.close(fig)

        # the interpolation is not particularly accurate, but all four interpolations give very comparable errors
    prev = None
    for k,v in error.items():
        if prev is not None:
            for i in range(v.shape[0]):
                assert v[i,0] == pytest.approx(prev[i,0])
        prev = v

    finished(inspect.stack()[0][3])


@started_finished
def test_LagrangeMesh_interpolate3D_bell(debug=False):
    """"""
    d = 1.
    N = 10
    nip = 21
    r = np.linspace(-2,2,num=nip)
    xyz = create_mesh([r, r, r])
    xyz += 1e-9

    sigma = 1.0
    bell = gauss(xyz, sigma=sigma)
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    error = {}
    for reduced in [
        (False, False, False),
        (False, False, True),
        (False, True, False),
        (False, True, True ),
        (True, False, False),
        (True, False, True),
        (True, True, False),
        (True, True, True),
    ]:
        print(f"{reduced=}")
        mesh = LagrangeMesh(dim=3, M=2*N, d=d, reduced=reduced)
        bell_g = gauss(mesh.grid, sigma=sigma)

        # wrap the interpolated quantity in a MeshQuantity
        Q = MeshQuantity(mesh, data=bell_g, symmetry=1)
        # interpolate Q at xy

        with Timer('old') as t:
            Qxyz = mesh.interpolate(Q, xyz, algo='old')
        print(t)

        with Timer('new') as t:
            Qxyz = mesh.interpolate(Q, xyz, algo='new')
        print(t)
        err = np.abs(Qxyz-bell)
        error[reduced] = err
        # the interpolation is not particularly accurate, but all four interpolations give very comparable errors
        prev = None
        for k,v in error.items():
            if prev is not None:
                for i in range(v.shape[0]):
                    assert v[i,0] == pytest.approx(prev[i,0])
        prev = v


@started_finished
def test_LagrangeMesh_interpolate3D_bellx(debug=False):
    """"""
    d = 1.
    N = 10
    nip = 21
    r = np.linspace(-2,2,num=nip)
    xyz = create_mesh([r, r, r])
    xyz += 1e-9

    sigma = 1.0
    bellx = gauss(xyz, sigma=sigma) * xyz[:,0]

    error = {}
    for reduced in [
        (False, False, False),
        (False, False, True),
        (False, True, False),
        (False, True, True ),
        (True, False, False),
        (True, False, True),
        (True, True, False),
        (True, True, True),
    ]:
        print(f"{reduced=}")
        mesh = LagrangeMesh(dim=3, M=2*N, d=d, reduced=reduced)
        bellx_g = gauss(mesh.grid, sigma=sigma) * mesh.grid[:,0]

        # wrap the interpolated quantity in a MeshQuantity
        Q = MeshQuantity(mesh, data=bellx_g, symmetry=[[-1,1,1]])
        # interpolate Q at xy
        Qxyz = mesh.interpolate(Q, xyz)

        err = np.abs(Qxyz-bellx)
        error[reduced] = err
        # the interpolation is not particularly accurate, but all four interpolations give very comparable errors
    prev = None
    for k,v in error.items():
        if prev is not None:
            for i in range(v.shape[0]):
                assert v[i,0] == pytest.approx(prev[i,0])
        prev = v


@started_finished
def test_create_mesh():
    # Because create_mesh(g1D[0],gy,gz) internally calls create_mesh(g1D[0],gy), we do not need a separate test for the 2D case
    # the 1D case is trivial

    gx = np.array([1,2], dtype=float)
    gy = 10*gx
    gz = 10*gy
    nx = len(gx)
    ny = len(gy)
    nz = len(gz)
    gxyz = create_mesh([gx, gy, gz])
    for i in range(nx*ny*nz):
        assert gxyz[i,0] == gx[i   %2]
        assert gxyz[i,1] == gy[i//2%2]
        assert gxyz[i,2] == gz[i//4%2]


@started_finished
def test_LagrangeMesh_lagrange_function():
    for reduced in [False, True]:

        print(f"1D case, {reduced=}")
        mesh = LagrangeMesh(dim=1, M=6, d=1, reduced=False)
        for ix,xi in enumerate(mesh.grid):
            fgridx = mesh.lagrange_function(mesh.gridx, (ix,))
            for (f,x) in zip(fgridx,mesh.gridx):
                if x == xi:
                    assert np.isnan(f), f"{ix=} {xi=} {f=} expected nan"# corner case due to 0/0, which by l'Hopitals rule should be 1
                else:
                    assert f == pytest.approx(.0, abs=1e-14), f"{ix=} {xi=} {f=} expected 0,0"

        print(f"2D case, {reduced=}")
        mesh = LagrangeMesh(dim=2, M=6, d=1, reduced=False)
        ixy = 0
        for iy in range(mesh.M[1]):
            for ix in range(mesh.M[0]):
                fgridxy = mesh.lagrange_function(mesh.grid + 1e-12,(ix,iy))
                for ig, (f,x) in enumerate(zip(fgridxy,mesh.grid)):
                    if ixy == ig:
                        assert f == pytest.approx(1., abs=1e-14), f"{ig=} {ixy=} {f=} expected 1.0"
                    else:
                        assert f == pytest.approx(.0, abs=1e-10), f"{ig=} {ixy=} {f=} expected 0.0"
                        # print("ok")
                ixy += 1

        print(f"3D case, {reduced=}")
        mesh = LagrangeMesh(dim=3, M=6, d=1, reduced=False)
        ixyz = 0
        for iz in range(mesh.M[2]):
            for iy in range(mesh.M[1]):
                for ix in range(mesh.M[0]):
                    fgridxyz = mesh.lagrange_function(mesh.grid + 1e-12, (ix,iy,iz))
                    for ig, (f,xyz) in enumerate(zip(fgridxyz,mesh.grid)):
                        if ixyz == ig:
                            assert f == pytest.approx(1., abs=1e-10), f"{ig=} {ixyz} {f=} expected 1.0"
                        else:
                            assert f == pytest.approx(.0, abs=1e-10), f"{ig=} {ixyz} {f=} expected 0.0"
                    ixyz += 1


def test_LagrangeMesh_lagrange_function_plot2D(no_plot):
    started(inspect.stack()[0][3])

    M = 4
    d = 1
    mesh = LagrangeMesh(dim=2, M=M, d=d, reduced=False)
    nip = 61
    r = np.linspace(-2,2,num=nip)
    xy = create_mesh([r, r])
    x = xy[:,0].reshape((nip,nip), order='F')
    y = xy[:,1].reshape((nip,nip), order='F')
    for i in range(M):
        for j in range(M):
            ij = (i,j)
            lf_ij = mesh.lagrange_function(xy, ij)
            lf_ij = lf_ij.reshape((nip,nip), order='F')

            if no_plot:
                pass
            else:
                fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
                plt.title(f'{i=}, x_i={mesh.g1D[0][i]}, y_j={mesh.g1D[1][j]}')
                ax.plot_surface(x, y, lf_ij, label='real', cmap=cm.coolwarm, linewidth=0, antialiased=False)
                plt.legend()
                fig.savefig(png_folder/"2D_lagrange_function_({i},{j}).png")
                # plt.show()
                plt.close(fig)

    finished(inspect.stack()[0][3])
