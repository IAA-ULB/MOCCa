import pytest
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib import cm

from mocca.mesh import LagrangeMesh
from mocca.mesh.lagrange import create_mesh
from mocca.mesh.observable import Observable

project_folder = Path(__file__)
while project_folder.name != 'tantalus_full':
    project_folder = project_folder.parent
test_folder = project_folder/f"MOCCaPy/tests/mocca/mesh"
# remove all .png files
for png in test_folder.glob('*.png'):
    png.unlink()

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
    assert mesh.M == (6,6)
    assert (mesh.d == [.5,.5]).all()

    mesh = LagrangeMesh(dim=2, M=6, d=(.5,.6))
    assert mesh.dim == 2
    assert mesh.M == (6,6)
    assert (mesh.d == (.5,.6)).all()

    mesh = LagrangeMesh(M=(4,6), d=.5)
    assert mesh.dim == 2
    assert mesh.M == (4,6)
    assert (mesh.d == (.5,.5)).all()

    mesh = LagrangeMesh(M=(4,6), d=(.5,.6))
    assert mesh.dim == 2
    assert mesh.M == (4,6)
    assert (mesh.d == (.5,.6)).all()

    # valid 3D mesh constructions
    mesh = LagrangeMesh(dim=3, M=4, d=.5)
    assert mesh.dim == 3
    assert mesh.M == (4,4,4)
    assert (mesh.d == (.5,.5,.5)).all()

    mesh = LagrangeMesh(dim=3, M=(2,4,6), d=.5,)
    assert mesh.dim == 3
    assert mesh.M == (2,4,6)
    assert (mesh.d == (.5,.5,.5)).all()

    mesh = LagrangeMesh(dim=3, M=4, d=(.4,.5,.6))
    assert mesh.dim == 3
    assert mesh.M == (4,4,4)
    assert (mesh.d == (.4,.5,.6)).all()

    mesh = LagrangeMesh(M=(2,4,6), d=(.4,.5,.6))
    assert mesh.dim == 3
    assert mesh.M == (2,4,6)
    assert (mesh.d == (.4,.5,.6)).all()

    mesh = LagrangeMesh(M=(2,4,6), d=.5)
    assert mesh.dim == 3
    assert mesh.M == (2,4,6)
    assert (mesh.d == (.5,.5,.5)).all()

    # Invalid 1D mesh constructions
    with pytest.raises(AssertionError):
        mesh = LagrangeMesh(M=4, d=.5)
        # although this could be interpreted as a 1D mesh, the user might expect a 3D mesh.
        # The user must either provide dim=1 or M=tuple(5,)

    with pytest.raises(AssertionError):
        mesh = LagrangeMesh(dim=1, M=4, d=(.5,.5))

    with pytest.raises(AssertionError):
        mesh = LagrangeMesh(dim=1, M=(4,6), d=.5)

    # Invalid 2D mesh constructions
    with pytest.raises(AssertionError):
        mesh = LagrangeMesh(M=(4,6), d=(.4,.5,.6))

    with pytest.raises(AssertionError):
        mesh = LagrangeMesh(dim=2, M=(2,4,6), d=(.4,.5))

    with pytest.raises(AssertionError):
        mesh = LagrangeMesh(dim=2, M=4, d=(.4,.5,.6))

    # Invalid 3D mesh constructions
    with pytest.raises(AssertionError):
        mesh = LagrangeMesh(M=(4,6), d=(.4,.5,.6))

    with pytest.raises(AssertionError):
        mesh = LagrangeMesh(dim=3, M=4, d=(.5,.6))

    # Invalid 4D mesh constructions
    with pytest.raises(AssertionError):
        mesh = LagrangeMesh(dim=4, M=4, d=.5)

    with pytest.raises(AssertionError):
        mesh = LagrangeMesh(M=(2,4,6,8), d=.5)

    # invalid n
    with pytest.raises(AssertionError):
        mesh = LagrangeMesh(dim=1, M=3, d=.5) # odd n

    with pytest.raises(AssertionError):
        mesh = LagrangeMesh(dim=1, M=0, d=.5) # n == 0

    with pytest.raises(AssertionError):
        mesh = LagrangeMesh(dim=1, M=2, d=-.5) # negative d

    with pytest.raises(AssertionError):
        mesh = LagrangeMesh(dim=2, M=2, d=(-.5,.5)) # negative d

    with pytest.raises(AssertionError):
        mesh = LagrangeMesh(dim=1, M=2, d=0)  # negative d

    with pytest.raises(AssertionError):
        mesh = LagrangeMesh(dim=2, M=2, d=(0,.5)) # negative d


def test_LagrangeMesh_ctor_shift():
    """Test valid and invalid reduced parameter."""

    mesh = LagrangeMesh(dim=3, M=4, d=.5)
    assert mesh.shift == (.0, .0, .0)

    mesh = LagrangeMesh(dim=3, M=4, d=.5, shift=.1, reduced=False)
    assert mesh.shift == (.1, .1, .1)

    mesh = LagrangeMesh(dim=3, M=4, d=.5, shift=(.1,.2,.3), reduced=False)
    assert mesh.shift == (.1, .2, .3)


def test_LagrangeMesh_ctor_reduced():
    """Test valid and invalid reduced parameter."""

    mesh = LagrangeMesh(dim=3, M=4, d=.5, reduced=True)
    assert mesh.reduced == (True, True, True )

    mesh = LagrangeMesh(dim=3, M=4, d=.5, reduced=False)
    assert mesh.reduced == (False, False, False)

    mesh = LagrangeMesh(dim=3, M=4, d=.5, reduced=(True,True,False))
    assert mesh.reduced == (True, True, False)

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

def test_LagrangeMesh_ctor_grid_2D():
    """Test grid points for 2D grid."""
    mesh = LagrangeMesh(dim=2, M=6, d=1.)
    for j in range(3):
        assert np.all(mesh.gridx[:,j] == [0.5,1.5,2.5])
    for i in range(3):
        assert np.all(mesh.gridy[i,:] == [0.5,1.5,2.5])

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


def test_LagrangeMesh_ctor_box_width():
    """Test dv computation."""
    mesh = LagrangeMesh(dim=3, M=6, d=.5)
    assert (mesh.box_width == tuple(3*[6*.5])).all()


def test_LagrangeMesh_plane_wave_1D():
    d = .5
    n = 6
    mesh = LagrangeMesh(dim=3, M=6, d=d, reduced=False)
    k = 1.5
    L = n*d
    r = mesh.gridx[:,0,0]
    pw = mesh.plane_wave_1D(L=L, k=k, r=r)
    
    arg = (2*np.pi * k / L) * r
    real_part = np.sqrt(1/L) * np.cos(arg)
    imag_part = np.sqrt(1/L) * np.sin(arg)
    assert (pw[:,0] == real_part).all()
    assert (pw[:,1] == imag_part).all()


def test_LagrangeMesh_reshape():
    """Test reshaping of the mesh."""
    mesh = LagrangeMesh(dim=3, M=4, d=.5)
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


def test_LagrangeMesh_basis_function_1D():
    """"""
    # 1D mesh
    N = 6
    d = 1
    mesh = LagrangeMesh(M=(N,), d=d, reduced=False)
    # r = mesh.gx
    r = np.linspace(-N*d/2, N*d/2, num=241)

    for i in range(N):
        pw = mesh.basis_function(i,r)
        fig, ax = plt.subplots()
        plt.title(f'{i=}, x_i={mesh.gx[i]}')
        ax.plot([-N/2,N/2],[0,0],color='k')
        ax.plot([0,0],[-1,1],color='k')
        ax.plot(mesh.gx,  N*[0], 'rs')
        ax.plot(mesh.gx[i], [0], 'gs')
        ax.plot(r, pw[:,0], label='real')
        ax.plot(r, pw[:,1], label='imag')
        plt.legend()
        fig.savefig(test_folder/f"1D_basis_function_{i}.png")
        fig.clear()
        plt.close(fig)

def test_LagrangeMesh_basis_function_2D():
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

            fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
            plt.title(f'{i=}, x_i={mesh.gx[i]}, y_j={mesh.gy[j]} real')
            ax.plot_surface(x, y, pw_real, cmap=cm.coolwarm, linewidth=0, antialiased=False)
            plt.legend()
            fig.savefig(test_folder/f"2D_basis_function_{i}_real.png")
            # plt.show()
            plt.close(fig)

            fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
            plt.title(f'{i=}, x_i={mesh.gx[i]}, y_j={mesh.gy[j]} imag')
            ax.plot_surface(x, y, pw_imag, cmap=cm.coolwarm, linewidth=0, antialiased=False)
            plt.legend()
            fig.savefig(test_folder/f"2D_basis_function_{i}_imag.png")
            # plt.show()
            plt.close(fig)


from numba import vectorize, float64

@vectorize([float64(float64, float64, float64)])
def norm(x, y, z):
    return np.sqrt(x*x + y*y + z*z)


def test_LagrangeMesh_apply():
    mesh = LagrangeMesh(dim=3, M=10, d=.5)
    r_xyz = mesh.apply(norm)
    for i in range(mesh.flat_shape[0]):
        assert r_xyz[i] == np.sqrt(mesh.gridx[i]**2 + mesh.gridy[i]**2 + mesh.gridz[i]**2)


def test_LagrangeMesh_integrate():
    mesh = LagrangeMesh(dim=3, M=10, d=.5)
    # integrate a constant function.

    q = np.ones((mesh.n_gridpoints(),))
    Q= Observable(mesh,data=q)
    integral_of_Q = Q.integrate()
    assert integral_of_Q == mesh.n_gridpoints() * mesh.dv


def test_LagrangeMesh_interpolate1D(debug=False):
    d = 1.
    N = 3
    for reduced in [
        False,
        True,
    ]:
        str_reduced = "(reduced)" if reduced else ""
        mesh  = LagrangeMesh(dim=1, M=2*N, d=d, reduced=reduced)
        r = np.empty((mesh.n_gridpoints(),), dtype=float, order='F')
        r = mesh.gx
        rip = np.linspace(-N*d, N*d, num=61)
        for ip in range(len(rip)):
            if (reduced and rip[ip] < 0 and -rip[ip] in mesh.gx) or (rip[ip] in mesh.gx):
                rip[ip] += 1e-12 # avoid 0/0 in the Lagrange Functions

        for ii in range(2*N):
            if not mesh.reduced[0]:
                i = ii
                negative_axis = False
                print(f"{reduced=} bf={i}/{2*N}")
            else:
                i = ii//2
                if ii % 2 == 0:
                    negative_axis = False
                    print(f"{reduced=} bf={i}/{2 * N}")
                else:
                    negative_axis = True
                    print(f"{reduced=} bf=-{i}/{2 * N}")


            lcpw_rgp = mesh.basis_function(i,r  ,negative_axis=negative_axis)
            lcpw_rip = mesh.basis_function(i,rip,negative_axis=negative_axis)

            Q = Observable( mesh, data=lcpw_rgp, symmetry=[1,-1]) # real/imag component is symmetric/skew-symmetric
            # Q = Observable( lcpw_rgp, symmetry=1)

            Qrip = mesh.interpolate(Q, rip)

            fig, ax = plt.subplots()
            plt.title(f"test_LagrangeMesh_interpolate1D real part $x_{i}$")
            plt.plot(rip,lcpw_rip[:,0], 'o', label=f'plane wave $x_{i}$ real')
            plt.plot(rip,    Qrip[:,0],      label=f'interpolated $x_{i}$ real')
            plt.legend()
            plt.savefig(test_folder/f"test_LagrangeMesh_interpolate1D_{i}_real{str_reduced}")
            if debug:
                plt.show()
            plt.close(fig)

            fig, ax = plt.subplots()
            plt.title(f"test_LagrangeMesh_interpolate1D imag part $x_{i}$")
            plt.plot(rip,lcpw_rip[:,1], 'o', label=f'plane wave $x_{i}$ imag')
            plt.plot(rip,    Qrip[:,1],      label=f'interpolated $x_{i}$ imag')
            plt.legend()
            plt.savefig(test_folder/f"test_LagrangeMesh_interpolate1D_{i}_imag{str_reduced}")
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
                       coeff[1] * mesh.basis_function(0, r  , negative_axis=True)
            lcpw_rip = coeff[0] * mesh.basis_function(0, rip) + \
                       coeff[1] * mesh.basis_function(0, rip, negative_axis=True)

        if not mesh.reduced[0]:
            for i in range(1,2*N):
                lcpw_rgp += coeff[i] * mesh.basis_function(i,r)
                lcpw_rip += coeff[i] * mesh.basis_function(i,rip)
        else:
            for i in range(1,N):
                lcpw_rgp += coeff[2*i  ] * mesh.basis_function(i,r  ) + \
                            coeff[2*i+1] * mesh.basis_function(i,r  , negative_axis=True)
                lcpw_rip += coeff[2*i  ] * mesh.basis_function(i,rip) + \
                            coeff[2*i+1] * mesh.basis_function(i,rip, negative_axis=True)

        Q = Observable(mesh, data=lcpw_rgp, symmetry=[1,-1])
        Qrip = mesh.interpolate(Q, rip)

        fig, ax = plt.subplots()
        plt.title(f"{s}: real part")
        plt.plot(rip,lcpw_rip[:,0], 'o', label=f'plane wave $x_{i}$ real')
        plt.plot(rip,    Qrip[:,0],      label=f'interpolated $x_{i}$ real')
        plt.legend()
        plt.savefig(test_folder/f"test_LagrangeMesh_interpolate1D {s} real part")
        if debug:
            plt.show()
        plt.close(fig)

        fig, ax = plt.subplots()
        plt.title(f"{s}: imag part")
        ax.set_ylim([-1, 1])
        plt.plot(rip,lcpw_rip[:,1], 'o', label=f'plane wave $x_{i}$ real')
        plt.plot(rip,    Qrip[:,1],      label=f'interpolated $x_{i}$ real')
        plt.legend()
        plt.savefig(test_folder/f"test_LagrangeMesh_interpolate1D {s} imag part")
        if debug:
            plt.show()
        plt.close(fig)

        assert ((lcpw_rip - Qrip) < 1e-12).all()


def test_LagrangeMesh_interpolate2D_basisfunction(debug=False):
    """"""
    d = 1.
    N = 2
    nip = 61
    r = np.linspace(-2,2,num=nip)
    xy = create_mesh(r, r)
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

                fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
                plt.title(f'BF imag x_ij=({mesh.gx[ibfx]}, {mesh.gy[ibfy]})')
                ax.plot_surface(x, y, bf_real, cmap=cm.coolwarm, linewidth=0, antialiased=False)
                fig.savefig(test_folder / f"basis_function_({ibfx},{ibfy}).png")
                if debug:
                    plt.show()
                plt.close(fig)

                fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
                plt.title(f'BF imag x_ij=({mesh.gx[ibfx]}, {mesh.gy[ibfy]})')
                ax.plot_surface(x, y, bf_imag, cmap=cm.coolwarm, linewidth=0, antialiased=False)
                fig.savefig(test_folder / f"BF_real({ibfx},{ibfy}).png")
                if debug:
                    plt.show()
                plt.close(fig)

                # compute the values of the basis function at the grid points
                bf_xy_g = mesh.basis_function((ibfx,ibfy), mesh.grid)
                # wrap the interpolated quantity in an observable
                Q = Observable(mesh, data=bf_xy_g, symmetry=[1,-1])

                # interpolate Q at xy
                Qxy = mesh.interpolate(Q, xy)
                Qxy_real = Qxy[:, 0].reshape((nip, nip), order='F')
                Qxy_imag = Qxy[:, 1].reshape((nip, nip), order='F')

                fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
                plt.title(f'BF real x_ij={mesh.gx[ibfx]}, y_j={mesh.gy[ibfy]} interpolated')
                ax.plot_surface(x, y, Qxy_real, cmap=cm.coolwarm, linewidth=0, antialiased=False)
                fig.savefig(test_folder / f"BF_real({ibfx},{ibfy})_interpolated.png")
                if debug:
                    plt.show()
                plt.close(fig)

                fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
                plt.title(f'BF imag x_ij={mesh.gx[ibfx]}, y_j={mesh.gy[ibfy]} interpolated')
                ax.plot_surface(x, y, Qxy_imag, cmap=cm.coolwarm, linewidth=0, antialiased=False)
                fig.savefig(test_folder / f"BF_real({ibfx},{ibfy})_interpolated.png")
                if debug:
                    plt.show()
                plt.close(fig)

                print(f"{reduced=} bf = ({ibfx}, {ibfy})")
                for i in range(len(Qxy_real)):
                    assert Qxy[i,0] == pytest.approx(bf_xy[i,0]), f"{i=} {Qxy[i,0]} != {bf_xy[i,0]}, |diff|={abs(Qxy[i,0]-bf_xy[i,0])}"
                    assert Qxy[i,1] == pytest.approx(bf_xy[i,1]), f"{i=} {Qxy[i,1]} != {bf_xy[i,1]}, |diff|={abs(Qxy[i,1]-bf_xy[i,1])}"

def gauss(x, sigma):
    if x.shape[1] == 3:
        r2 = x[:,0]**2 + x[:,1]**2 + x[:,2]**2
    elif x.shape[1] == 2:
        r2 = x[:,0]**2 + x[:,1]**2
    else:
        r2 = x**2
    return np.exp(-0.5*r2/sigma**2)/(sigma*np.sqrt(2*np.pi))

def test_LagrangeMesh_interpolate2D_bell(debug=False):
    """"""
    d = 1.
    N = 40
    nip = 61
    r = np.linspace(-2,2,num=nip)
    xy = create_mesh(r, r)
    xy += 1e-9
    x = xy[:,0].reshape((nip,nip), order='F')
    y = xy[:,1].reshape((nip,nip), order='F')

    sigma = 1.0
    bell = gauss(xy, sigma=sigma)
    bell_plot = bell.reshape((nip,nip), order='F')
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    plt.title(f'bell')
    ax.plot_surface(x, y, bell_plot, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    fig.savefig(test_folder / "bell.png")
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

        # wrap the interpolated quantity in an observable
        Q = Observable(mesh, data=bell_g, symmetry=1)
        # interpolate Q at xy
        Qxy = mesh.interpolate(Q, xy)
        Qxy_plot = Qxy.reshape((nip, nip), order='F')

        fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
        plt.title(f'bell interpolated ({reduced=})')
        ax.plot_surface(x, y, Qxy_plot, cmap=cm.coolwarm, linewidth=0, antialiased=False)
        fig.savefig(test_folder / f"bell_interpolated_({reduced=}).png")
        if debug:
            plt.show()
        plt.close(fig)

        err = np.abs(Qxy_plot-bell_plot)
        error[reduced] = err

        fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
        plt.title(f'bell interpolated error ({reduced=})')
        ax.plot_surface(x, y, err, cmap=cm.coolwarm, linewidth=0, antialiased=False)
        fig.savefig(test_folder / f"bell_interpolated_error_({reduced=}).png")
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

def test_LagrangeMesh_interpolate2D_bellx(debug=False):
    """"""
    d = 1.
    N = 40
    nip = 61
    r = np.linspace(-2,2,num=nip)
    xy = create_mesh(r, r)
    xy += 1e-9
    x = xy[:,0].reshape((nip,nip), order='F')
    y = xy[:,1].reshape((nip,nip), order='F')

    sigma = 1.0
    bellx = gauss(xy, sigma=sigma) * xy[:,0]
    bellx_plot = bellx.reshape((nip,nip), order='F')
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    plt.title(f'bellx')
    ax.plot_surface(x, y, bellx_plot, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    fig.savefig(test_folder / "bellx.png")
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

        # wrap the interpolated quantity in an observable
        Q = Observable(mesh, data=bellx_g, symmetry=[[-1],[1]])
        # interpolate Q at xy
        Qxy = mesh.interpolate(Q, xy)
        Qxy_plot = Qxy.reshape((nip, nip), order='F')

        fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
        plt.title(f'bellx interpolated ({reduced=})')
        ax.plot_surface(x, y, Qxy_plot, cmap=cm.coolwarm, linewidth=0, antialiased=False)
        fig.savefig(test_folder / f"bellx_interpolated_({reduced=}).png")
        if debug:
            plt.show()
        plt.close(fig)

        err = np.abs(Qxy_plot-bellx_plot)
        error[reduced] = err

        fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
        plt.title(f'bellx interpolated error ({reduced=})')
        ax.plot_surface(x, y, err, cmap=cm.coolwarm, linewidth=0, antialiased=False)
        fig.savefig(test_folder / f"bellx_interpolated_error_({reduced=}).png")
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


def test_LagrangeMesh_interpolate3D_bell(debug=False):
    """"""
    d = 1.
    N = 10
    nip = 21
    r = np.linspace(-2,2,num=nip)
    xyz = create_mesh(r, r, r)
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

        # wrap the interpolated quantity in an observable
        Q = Observable(mesh, data=bell_g, symmetry=1)
        # interpolate Q at xy
        Qxyz = mesh.interpolate(Q, xyz)

        err = np.abs(Qxyz-bell)
        error[reduced] = err
        # the interpolation is not particularly accurate, but all four interpolations give very comparable errors
    prev = None
    for k,v in error.items():
        if prev is not None:
            for i in range(v.shape[0]):
                assert v[i,0] == pytest.approx(prev[i,0])
        prev = v

def test_LagrangeMesh_interpolate3D_bellx(debug=False):
    """"""
    d = 1.
    N = 10
    nip = 21
    r = np.linspace(-2,2,num=nip)
    xyz = create_mesh(r, r, r)
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

        # wrap the interpolated quantity in an observable
        Q = Observable(mesh, data=bellx_g, symmetry=[[-1],[1],[1]])
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


def test_create_mesh():
    # Because create_mesh(gx,gy,gz) internally calls create_mesh(gx,gy), we do not need a separate test for the 2D case
    # the 1D case is trivial

    gx = np.array([1,2], dtype=float)
    gy = 10*gx
    gz = 10*gy
    nx = len(gx)
    ny = len(gy)
    nz = len(gz)
    gxyz = create_mesh(gx, gy, gz)
    for i in range(nx*ny*nz):
        assert gxyz[i,0] == gx[i   %2]
        assert gxyz[i,1] == gy[i//2%2]
        assert gxyz[i,2] == gz[i//4%2]

def test_LagrangeMesh_lagrange_function():
    for reduced in [False, True]:

        print(f"1D case, {reduced=}")
        mesh = LagrangeMesh(dim=1, M=6, d=1, reduced=False)
        for ix,xi in enumerate(mesh.grid):
            fgridx = mesh.lagrange_function(mesh.gridx, ix)
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

def test_LagrangeMesh_lagrange_function_plot2D():
    M = 4
    d = 1
    mesh = LagrangeMesh(dim=2, M=M, d=d, reduced=False)
    nip = 61
    r = np.linspace(-2,2,num=nip)
    xy = create_mesh(r, r)
    x = xy[:,0].reshape((nip,nip), order='F')
    y = xy[:,1].reshape((nip,nip), order='F')
    for i in range(M):
        for j in range(M):
            ij = (i,j)
            lf_ij = mesh.lagrange_function(xy, ij)
            lf_ij = lf_ij.reshape((nip,nip), order='F')
            fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
            plt.title(f'{i=}, x_i={mesh.gx[i]}, y_j={mesh.gy[j]}')
            ax.plot_surface(x, y, lf_ij, label='real', cmap=cm.coolwarm, linewidth=0, antialiased=False)
            plt.legend()
            fig.savefig(test_folder / f"2D_lagrange_function_({i},{j}).png")
            # plt.show()
            plt.close(fig)
