import numpy as np
import scipy.sparse as sparse
import sympy
import pytest
from tabulate import tabulate
from pathlib import Path
import sys

from MOCCaPy.mocca.mean_field.solvers.generalized_poisson import \
    GeneralizedPoissonSolverMF, \
    GeneralizedPoissonSolver, \
    Laplacian1D
from MOCCaPy.scripts.mocca.util import title_line
# dia_entry_set, dia_entry_add

from mocca.mesh import LagrangeMesh, MeshQuantity
from mocca.util.timer import Timer

output = open(Path(__file__).parent / "test_generalized_poisson.py.txt", mode="w")

path2MOCCaPy = Path(__file__).parent
while not path2MOCCaPy.name == 'MOCCaPy':
    path2MOCCaPy = path2MOCCaPy.parent
    # print(path2MOCCaPy)
sys.path.insert(0, str(path2MOCCaPy))

from tests.util import started_finished, started, finished

def gauss(dim, sigma):
    """A Gaussian function $g$ centered around the origin in `dim` dimensions with standard deviation `sigma`,
    and its Laplacian $f = \Delta g$

    Args:
         dim: number of spatial dimensions.
         sigma: standard deviation.
    Returns:
        ue_lambda, fe_lambda, symmetry
        ue_lambda : function of position in dim-D space which is the analytical solution of the Poisson equation.
            Used for validating the numerical solution, and for the Dirichlet boundary conditions.
        fe_lambda : function of position in dim-D space which is the analytical solution of the Poisson equation
        symmetry: symmetry components of ue and fe:
    """
    sigma2 = sigma**2

    x,y,z = sympy.symbols('x y z', real=True)

    if dim == 1:
        ue = sympy.exp(-0.5 * (x ** 2) / sigma2)  # a gaussian curve
        fe = ue.diff(x, 2)
        ue_lambda = sympy.lambdify([x], ue.simplify(), "numpy")
        fe_lambda = sympy.lambdify([x], fe.simplify(), "numpy")
    elif dim == 2:
        ue = sympy.exp(-0.5*(x**2 + y**2)/sigma2)  # a gaussian curve
        fe = ue.diff(x, 2) + ue.diff(y, 2)
        ue_lambda = sympy.lambdify([x,y], ue.simplify(), "numpy")
        fe_lambda = sympy.lambdify([x,y], fe.simplify(), "numpy")
    elif dim == 3:
        ue = sympy.exp(-0.5 * (x ** 2 + y ** 2 + z ** 2) / sigma2)  # a gaussian curve
        fe = ue.diff(x, 2) + ue.diff(y, 2) + ue.diff(z, 2)
        ue_lambda = sympy.lambdify([x, y, z], ue.simplify(), "numpy")
        fe_lambda = sympy.lambdify([x, y, z], fe.simplify(), "numpy")
    else:
        raise ValueError(f"dim must be 1, 2, or 3, got {dim}.")

    symmetry = 1

    return ue_lambda, fe_lambda, symmetry


def x_gauss(dim, sigma):
    """A Gauss function $g$ multiplied by $x$, and therefor skew-symmetric in the x-direction,
    adn symmetric in the other directions.
    (and its Laplacian $f = \Delta x g$)

    Args:
         as for `gauss`
    Returns:
         as for `gauss`
    """
    sigma2 = sigma**2

    x,y,z = sympy.symbols('x y z', real=True)

    if dim == 1:
        ue = x * sympy.exp(-0.5 * (x ** 2) / sigma2)  # a gaussian curve
        fe = ue.diff(x, 2)
        ue_lambda = sympy.lambdify([x], ue.simplify(), "numpy")
        fe_lambda = sympy.lambdify([x], fe.simplify(), "numpy")
        symmetry = -1

    elif dim == 2:
        ue = x * sympy.exp(-0.5*(x**2 + y**2)/sigma2)  # a gaussian curve
        fe = ue.diff(x, 2) + ue.diff(y, 2)
        ue_lambda = sympy.lambdify([x,y], ue.simplify(), "numpy")
        fe_lambda = sympy.lambdify([x,y], fe.simplify(), "numpy")
        symmetry = (-1,1)

    elif dim == 3:
        ue = x * sympy.exp(-0.5 * (x ** 2 + y ** 2 + z ** 2) / sigma2)  # a gaussian curve
        fe = ue.diff(x, 2) + ue.diff(y, 2) + ue.diff(z, 2)
        ue_lambda = sympy.lambdify([x, y, z], ue.simplify(), "numpy")
        fe_lambda = sympy.lambdify([x, y, z], fe.simplify(), "numpy")
        symmetry = (-1, 1, 1)

    else:
        raise ValueError(f"dim must be 1, 2, or 3, got {dim}.")

    return ue_lambda, fe_lambda, symmetry


def y_gauss(dim, sigma):
    """A Gauss function $g$ multiplied by $y$, and therefor skew-symmetric in the y-direction,
    adn symmetric in the other directions.
    (and its Laplacian $f = \Delta x g$)

    Args:
         as for `gauss`
    Returns:
         as for `gauss`
    """
    sigma2 = sigma**2

    x,y,z = sympy.symbols('x y z', real=True)

    if dim == 1:
        raise ValueError(f"y_gauss(x,y[,z]) requires dim to be 2, or 3, got {dim}.")

    elif dim == 2:
        ue = y * sympy.exp(-0.5*(x**2 + y**2)/sigma2)  # a gaussian curve
        fe = ue.diff(x, 2) + ue.diff(y, 2)
        ue_lambda = sympy.lambdify([x,y], ue.simplify(), "numpy")
        fe_lambda = sympy.lambdify([x,y], fe.simplify(), "numpy")
        symmetry = (1,-1)

    elif dim == 3:
        ue = y * sympy.exp(-0.5 * (x ** 2 + y ** 2 + z ** 2) / sigma2)  # a gaussian curve
        fe = ue.diff(x, 2) + ue.diff(y, 2) + ue.diff(z, 2)
        ue_lambda = sympy.lambdify([x, y, z], ue.simplify(), "numpy")
        fe_lambda = sympy.lambdify([x, y, z], fe.simplify(), "numpy")
        symmetry = (1, -1, 1)

    else:
        raise ValueError(f"dim must be 1, 2, or 3, got {dim}.")

    return ue_lambda, fe_lambda, symmetry


def z_gauss(dim, sigma):
    """A Gauss function $g$ multiplied by $y$, and therefor skew-symmetric in the z-direction,
    adn symmetric in the other directions.
    (and its Laplacian $f = \Delta x g$)

    Args:
         as for `gauss`
    Returns:
         as for `gauss`
    """
    sigma2 = sigma**2

    x,y,z = sympy.symbols('x y z', real=True)

    if dim !=3 :
        raise ValueError(f"z_gauss(x,y,z) requires dim to be 3, got {dim}.")

    ue = z * sympy.exp(-0.5 * (x ** 2 + y ** 2 + z ** 2) / sigma2)  # a gaussian curve
    fe = ue.diff(x, 2) + ue.diff(y, 2) + ue.diff(z, 2)
    ue_lambda = sympy.lambdify([x, y, z], ue.simplify(), "numpy")
    fe_lambda = sympy.lambdify([x, y, z], fe.simplify(), "numpy")
    symmetry = (1, 1, -1)

    return ue_lambda, fe_lambda, symmetry


def xy_gauss(dim, sigma):
    """A Gauss function $g$ multiplied by $xy$, and therefor skew-symmetric in the x- and y-direction,
    adn symmetric in the z-direction.
    (and its Laplacian $f = \Delta x g$)

    Args:
         as for `gauss`
    Returns:
         as for `gauss`
    """
    sigma2 = sigma**2

    x,y,z = sympy.symbols('x y z', real=True)

    if dim == 1:
        raise ValueError(f"xy_gauss(x,y[,z]) requires dim to be 2, or 3, got {dim}.")

    elif dim == 2:
        ue = x * y * sympy.exp(-0.5*(x**2 + y**2)/sigma2)  # a gaussian curve
        fe = ue.diff(x, 2) + ue.diff(y, 2)
        ue_lambda = sympy.lambdify([x,y], ue.simplify(), "numpy")
        fe_lambda = sympy.lambdify([x,y], fe.simplify(), "numpy")
        symmetry = (-1,-1)

    elif dim == 3:
        ue = x * y * sympy.exp(-0.5 * (x ** 2 + y ** 2 + z ** 2) / sigma2)  # a gaussian curve
        fe = ue.diff(x, 2) + ue.diff(y, 2) + ue.diff(z, 2)
        ue_lambda = sympy.lambdify([x, y, z], ue.simplify(), "numpy")
        fe_lambda = sympy.lambdify([x, y, z], fe.simplify(), "numpy")
        symmetry = (-1, -1, 1)

    else:
        raise ValueError(f"dim must be 1, 2, or 3, got {dim}.")

    return ue_lambda, fe_lambda, symmetry


def xz_gauss(dim, sigma):
    """A Gauss function $g$ multiplied by $xz$, and therefor skew-symmetric in the x- and z-direction,
    adn symmetric in the y-direction.
    (and its Laplacian $f = \Delta x g$)

    Args:
         as for `gauss`
    Returns:
         as for `gauss`
    """
    sigma2 = sigma**2

    x,y,z = sympy.symbols('x y z', real=True)

    if dim != 3:
        raise ValueError(f"xz_gauss(x,y,z) requires dim to be 2, or 3, got {dim}.")

    elif dim == 3:
        ue = x * z * sympy.exp(-0.5 * (x ** 2 + y ** 2 + z ** 2) / sigma2)  # a gaussian curve
        fe = ue.diff(x, 2) + ue.diff(y, 2) + ue.diff(z, 2)
        ue_lambda = sympy.lambdify([x, y, z], ue.simplify(), "numpy")
        fe_lambda = sympy.lambdify([x, y, z], fe.simplify(), "numpy")
        symmetry = (-1, 1, -1)

    else:
        raise ValueError(f"dim must be 1, 2, or 3, got {dim}.")

    return ue_lambda, fe_lambda, symmetry


def yz_gauss(dim, sigma):
    """A Gauss function $g$ multiplied by $yz$, and therefor skew-symmetric in the y- and z-direction,
    and symmetric in the x-direction.
    (and its Laplacian $f = \Delta x g$)

    Args:
         as for `gauss`
    Returns:
         as for `gauss`
    """
    sigma2 = sigma**2

    x,y,z = sympy.symbols('x y z', real=True)

    if dim != 3:
        raise ValueError(f"xz_gauss(x,y,z) requires dim to be 2, or 3, got {dim}.")

    ue = y * z * sympy.exp(-0.5 * (x ** 2 + y ** 2 + z ** 2) / sigma2)  # a gaussian curve
    fe = ue.diff(x, 2) + ue.diff(y, 2) + ue.diff(z, 2)
    ue_lambda = sympy.lambdify([x, y, z], ue.simplify(), "numpy")
    fe_lambda = sympy.lambdify([x, y, z], fe.simplify(), "numpy")
    symmetry = (1, -1, -1)

    return ue_lambda, fe_lambda, symmetry


def xyz_gauss(dim, sigma):
    """A Gauss function $g$ multiplied by $xyz$, and therefor skew-symmetric in all three directions.
    (and its Laplacian $f = \Delta x g$)

    Args:
         as for `gauss`
    Returns:
         as for `gauss`
    """
    sigma2 = sigma**2

    x,y,z = sympy.symbols('x y z', real=True)

    if dim != 3:
        raise ValueError(f"xz_gauss(x,y,z) requires dim to be 2, or 3, got {dim}.")

    ue = x * y * z * sympy.exp(-0.5 * (x ** 2 + y ** 2 + z ** 2) / sigma2)  # a gaussian curve
    fe = ue.diff(x, 2) + ue.diff(y, 2) + ue.diff(z, 2)
    ue_lambda = sympy.lambdify([x, y, z], ue.simplify(), "numpy")
    fe_lambda = sympy.lambdify([x, y, z], fe.simplify(), "numpy")
    symmetry = (-1, -1, -1)

    return ue_lambda, fe_lambda, symmetry


def run_ND(
        analytical_solution, sigma,
        dim, M, h, reduced,
        bc_scale=None, stencil=3, method='direct',
        verbosity = 1,
    ):
    # Warming up Timer instance
    tmr = Timer()
    tmr.start()
    elapsed = tmr.stop()

    ue_lambda, fe_lambda, symmetry = analytical_solution(dim=dim, sigma=sigma)

    mesh = LagrangeMesh(dim=dim, M=M, d=h, reduced=reduced, collect_boundary_points=stencil)

    s = title_line(char='-', width=120)
    s += f"solution: u(r)={analytical_solution.__name__}, {sigma=}, {symmetry=}\n"
    s += f"mesh    : {dim=}, {M=}, {h=}, {reduced=}\n"
    s += f"solver  : {bc_scale=}, {stencil=}, {method=}\n"
    s += f"linear system : {mesh.linear_size}x{mesh.linear_size}"
    print(s, file=output)
    if verbosity:
        print(s)

    ue = MeshQuantity(mesh, name='u', n_components=1, symmetry=symmetry)
    ue.data[:,0] = ue_lambda(mesh.grid[:,0], mesh.grid[:,1], mesh.grid[:,2]) if (dim == 3) else \
                   ue_lambda(mesh.grid[:,0], mesh.grid[:,1]) if (dim == 2) else \
                   ue_lambda(mesh.grid[:])

    f = MeshQuantity(mesh, name='f', n_components=1, symmetry=symmetry)
    f.data[:, 0] = fe_lambda(mesh.grid[:,0], mesh.grid[:,1], mesh.grid[:,2]) if (dim == 3) else \
                   fe_lambda(mesh.grid[:,0], mesh.grid[:,1]) if (dim == 2) else \
                   fe_lambda(mesh.grid[:])
    tmr.start()
    pSolver = GeneralizedPoissonSolver(mesh, stencil=stencil, bc_scale=bc_scale)
    pSolver.assemble(f, ue_lambda, symmetry=symmetry)
    cput_asmbl = tmr.stop()


    if verbosity >= 2:
        print(pSolver)
        # print(pSolver, file=output)

    tmr.start()
    u = pSolver.solve(method=method)
    cput_solve = tmr.stop()

    diff = np.abs(ue.data[:, 0] - u)
    rmse = np.sqrt(np.sum(np.square(diff)) / M)
    mean_diff = float(np.mean(diff))
    max_diff = float(np.max(diff))
    max_rel = float(np.max(diff/u))

    s = f"{bc_scale=} {M=} : {rmse=} {mean_diff=} {max_diff=} {cput_asmbl:.5f}s {cput_solve:.5f}s"
    print(s, file=output)
    if verbosity:
        print(s)

    return rmse, mean_diff, max_diff, max_rel, cput_asmbl, cput_solve, u


@pytest.mark.skip(reason="probably obsolete")
def test_1D():
    dim = 1

    s = title_line(text='test_1D', char='-', width=120)
    print(s)
    print(s, file=output)

    tbl = []

    sigma = 2
    for analytical_solution in [gauss, x_gauss]:
        for stencil in [
            3,
            5,
        ]:
            for reduced in [
                False,
                True,
            ]:
                for bc_scale in [
                    None,
                    1e20,
                ]:
                    rmse0, mean_diff0, max_diff0 = 1e9, 1e9, 1e9
                    M = 4
                    h = 1.6
                    for iter in range(3):
                        M *= 2
                        h *= .5
                        verbosity = 2 if (iter >=0) else 1

                        rmse, mean_diff, max_diff, max_rel, cput_asmbl, cput_solve, u = run_ND(
                            analytical_solution=analytical_solution, sigma=sigma,
                            dim=dim, M=M, h=h, reduced=reduced,
                            bc_scale=bc_scale, stencil=stencil,
                            verbosity=verbosity,
                        )

                        tbl.append([analytical_solution.__name__, stencil, reduced, bc_scale, M, rmse, mean_diff, max_diff, max_rel, cput_asmbl, cput_solve])
                        assert rmse < rmse0
                        assert mean_diff < mean_diff0
                        assert max_diff < max_diff0

                        rmse0, mean_diff0, max_diff0 = rmse, mean_diff, max_diff

    s = "\n" + title_line(text='SUMMARY', char='-', width=120, above=True, below=True)
    s += tabulate(tbl
        , headers=["u(r)", "stencil", "reduced", "bc_scale", "M", "RMSE", "Mean diff", "Max diff", "Max rel", "cput_asmlb", 'cput_solve']
        , tablefmt="simple"
    )
    print(s)
    print(s, file=output)


@pytest.mark.skip(reason="probably obsolete")
def test_2D():
    dim = 2

    s = title_line(text='test_2D', char='-', width=120)
    print(s)
    print(s, file=output)

    tbl = []

    sigma = 2
    for analytical_solution in [gauss, x_gauss, y_gauss, xy_gauss]:
        for stencil in [
            3,
            5,
        ]:
            for reduced in [
                False,
                True,
                (True, False),
                (False, True),
            ]:
                for bc_scale in [None, 1e20]:
                    for h in [
                        1.6,
                        (1.6, 1.61),
                    ]:
                        rmse0, mean_diff0, max_diff0 = 1e9, 1e9, 1e9
                        M = 4
                        for iter in range(3):
                            M *= 2
                            if isinstance(h, tuple):
                                h = tuple([hi*0.5 for hi in h])
                            else:
                                h *= .5

                            verbosity = 2 if (iter == 0) else 1

                            rmse, mean_diff, max_diff, max_rel, cput_asmbl, cput_solve, u = run_ND(
                                analytical_solution=analytical_solution, sigma=sigma,
                                dim=dim, M=M, h=h, reduced=reduced,
                                bc_scale=bc_scale, stencil=stencil,
                                verbosity=verbosity,
                            )

                            tbl.append([analytical_solution.__name__, stencil, reduced, bc_scale, M, rmse, mean_diff, max_diff, max_rel, cput_asmbl, cput_solve])
                            assert rmse < rmse0
                            assert mean_diff < mean_diff0
                            assert max_diff < max_diff0

                            rmse0, mean_diff0, max_diff0 = rmse, mean_diff, max_diff
                        pass

    s = "\n" + title_line(text='SUMMARY', char='-', width=120, above=True, below=True)
    s += tabulate(tbl
        , headers=["u(r)", "stencil", "reduced", "bc_scale", "M", "RMSE", "Mean diff", "Max diff", "Max rel", "cput_asmlb", 'cput_solve']
        , tablefmt="simple"
    )
    print(s)
    print(s, file=output)


# @pytest.mark.slow
@pytest.mark.skip(reason="probably obsolete")
def test_3D():
    dim = 3

    s = title_line(text='test_3D', char='-', width=120)
    print(s)
    print(s, file=output)

    tbl = []

    sigma = 2
    for stencil in [
        3,
        5,
    ]:
        for analytical_solution in [
            xy_gauss,
            gauss,
            x_gauss,
            y_gauss,
            z_gauss,
            xz_gauss,
            yz_gauss,
            xyz_gauss,
        ]:
            for reduced in [
                False,
                True,
                (True, False, False),
                (False, True, False),
                (False, False, True),
                (True, True, False),
                (True, False, True),
                (False, True, True),
            ]:
                for bc_scale in [
                    1e20,
                    None,
                ]:
                    for h in [
                        4.,
                        (4, 4.01, 4.01),
                    ]:
                        rmse0, mean_diff0, max_diff0 = 1e9, 1e9, 1e9
                        M = 4
                        w = M*h if isinstance(h, float) else [M*hi for hi in h]

                        for iter in range(7):
                            # Doubling M increases the linear system's size way too fast in 3D.
                            M += 4
                            h = w/M if isinstance(w, float) else tuple([wi/M for wi in w])
                            verbosity = 2 if (iter == 0) else 1

                            rmse, mean_diff, max_diff, cput_asmbl, cput_solve, u = run_ND(
                                analytical_solution=analytical_solution, sigma=sigma,
                                dim=dim, M=M, h=h, reduced=reduced,
                                bc_scale=bc_scale, stencil=stencil,
                                verbosity=verbosity,
                            )

                            tbl.append([analytical_solution.__name__, stencil, reduced, 0 if (bc_scale is None) else bc_scale, h, M, rmse, mean_diff, max_diff, cput_asmbl, cput_solve])
                            assert rmse < rmse0
                            assert mean_diff < mean_diff0
                            try:
                                assert max_diff < max_diff0
                            except AssertionError:
                                print(f"Assertion failed: {max_diff} < {max_diff0}")
                                print(f"Assertion failed: {max_diff} < {max_diff0}", file=output)

                            rmse0, mean_diff0, max_diff0 = rmse, mean_diff, max_diff

                        pass

    s = "\n" + title_line(text='SUMMARY', char='-', width=120, above=True, below=True)
    s += tabulate(tbl
                  , headers=["stencil", "u(r)", "reduced", "bc_scale", "h", "M", "RMSE", "Mean diff", "Max diff", "cput_asmlb", 'cput_solve']
                  , tablefmt="simple"
                  )
    print(s)
    print(s, file=output)


# TODO: dia_entry_set and dia_entry_add are no longer used. discard or improve.
# def test_set_dia_entry():
#     D = sparse.eye(4, k=1, format='dia')
#     D44 = D.toarray()
#     assert D44[0,1] == 1
#     dia_entry_set(D, [(0, 1, 2)])
#     D44 = D.toarray()
#     assert D44[0,1] == 2
#     with pytest.raises(ValueError):
#         dia_entry_set(D, [(0, 2, 2)])
#
#     D = sparse.eye(4, k=2, format='dia')
#     D44 = D.toarray()
#     assert D44[0, 2] == 1
#     dia_entry_set(D, [(0, 2, 2)])
#     D44 = D.toarray()
#     assert D44[0, 2] == 2
#
#     D = sparse.eye(4, k=-1, format='dia')
#     D44 = D.toarray()
#     assert D44[1, 0] == 1
#     dia_entry_set(D, [(1, 0, 2),
#                       (2, 1, 3)
#                      ])
#     D44 = D.toarray()
#     assert D44[1,0] == 2
#     assert D44[2,1] == 3
#
#     dia_entry_add(D, [
#         (1, 0, 2),
#         (2, 1, 3),
#     ])
#     D44 = D.toarray()
#     assert D44[1,0] == 4
#     assert D44[2,1] == 6
#
#     D = sparse.eye(4, k=-2, format='dia')
#     D44 = D.toarray()
#     assert D44[2, 0] == 1
#     dia_entry_set(D, [
#         (2, 0, 2),
#     ])
#     D44 = D.toarray()
#     assert D44[2, 0] == 2


@pytest.mark.skip(reason="probably obsolete")
def test_collect_boundary_points_1D():
    dim = 1
    M = 10
    d = 1
    boundary_width = 1
    for reduced in [True, False]:
        mesh = LagrangeMesh(
            dim=dim, M=M, d=d, reduced=reduced,
            highest_derivative_order=0,
            boundary_width = boundary_width
        )
        if reduced:
            assert len(mesh.bp_l) == 1
            assert mesh.bp_l[0] == 6
            assert mesh.bp_xyz[0] == M/2 - d/2 + boundary_width
        else:
            assert len(mesh.bp_l) == 2
            assert mesh.bp_l[0] == 0
            assert mesh.bp_xyz[0] == -(M/2 + d/2)
            assert mesh.bp_l[1] == M+1
            assert mesh.bp_xyz[1] ==  (M/2 + d/2)
        pass


@pytest.mark.skip(reason="probably obsolete")
def test_collect_boundary_points_2D():
    dim = 2
    M = (4,6)
    d = 1
    xb = (M[0]+d)/2
    yb = (M[1]+d)/2
    x = [-xb, xb]
    y = [-yb, yb]
    boundary_width = 1

    for reduced in [
        (True, True),
        (False, False),
        (True, False),
        (False, True),
    ]:
        print(f"{reduced=}")
        mesh = LagrangeMesh(
            dim=dim, M=M, d=d, reduced=reduced,
            highest_derivative_order=0,
            boundary_width = boundary_width
        )
        nbp = mesh.bp_l.size
        if reduced == (True, True):
            assert nbp == 4
        elif reduced == (False, False):
            assert nbp == 16
        else:
            assert nbp == 8

        for i in range(nbp):
            assert mesh.bp_xyz[i,0] in x or mesh.bp_xyz[i,1] in y


@pytest.mark.skip(reason="probably obsolete")
def test_collect_boundary_points_3D():
    dim = 3
    M = (4,6,8)
    d = 1
    xb = (M[0]-d)/2
    yb = (M[1]-d)/2
    zb = (M[2]-d)/2
    x = [-xb, xb]
    y = [-yb, yb]
    z = [-zb, zb]

    for reduced in [
        (True, True, True),
        (False, False, False),
        (True, False, False),
        (False, True, False),
        (False, False, True),
        (True, True, False),
        (True, False, True),
        (False, True, True),
    ]:
        print(f"{reduced=}")
        mesh = LagrangeMesh(dim=dim, M=M, d=d, reduced=reduced)
        mesh.collect_boundary_points(stencil=3)
        nbp = mesh.bp_l.size
        nbp_expected = 18 # for all reduced
        for r in reduced:
            if not r:
                nbp_expected *= 2
        assert nbp == nbp_expected

        for i in range(nbp):
            assert mesh.bp_xyz[i,0] in x or \
                   mesh.bp_xyz[i,1] in y or \
                   mesh.bp_xyz[i,2] in z


@pytest.mark.skip(reason="probably obsolete")
def test_Laplacia1D_non_uniform():
    n = 8
    h = 2
    L1 = Laplacian1D(8,    ).toarray() # h=None, corresponds to h=1 (for L at least)
    L2 = Laplacian1D(8, h=h).toarray()
    invh2 = 1/h**2
    for i in range(n):
        for j in range(n):
            assert L1[i,j]*invh2 == pytest.approx(L2[i,j])


def run_MF_ND(
        analytical_solution, sigma,
        dim, M, h, reduced,
        stencil=3, method='cg',
        verbosity = 1,
    ):
    tmr = Timer()

    ue_lambda, fe_lambda, symmetry = analytical_solution(dim=dim, sigma=sigma)

    mesh = LagrangeMesh(dim=dim, M=M, d=h, reduced=reduced)

    s = title_line(char='-', width=120)
    s += f"solution: u(r)={analytical_solution.__name__}, {sigma=}, {symmetry=}\n"
    s += f"mesh    : {dim=}, {M=}, {h=}, {reduced=}\n"
    s += f"solver  : {stencil=}, {method=}\n"
    s += f"linear system : {mesh.linear_size}x{mesh.linear_size}"
    print(s, file=output)
    if verbosity:
        print(s)

    ue = ue_lambda(mesh.grid[:,0], mesh.grid[:,1], mesh.grid[:,2]) if (dim == 3) else \
         ue_lambda(mesh.grid[:,0], mesh.grid[:,1]) if (dim == 2) else \
         ue_lambda(mesh.grid[:])

    f  = fe_lambda(mesh.grid[:,0], mesh.grid[:,1], mesh.grid[:,2]) if (dim == 3) else \
         fe_lambda(mesh.grid[:,0], mesh.grid[:,1]) if (dim == 2) else \
         fe_lambda(mesh.grid[:])

    tmr.start()
    pSolver = GeneralizedPoissonSolverMF(mesh, stencil=stencil)
    # pSolver.ue = ue # just convenient during debugging
    pSolver.assemble_rhs(f, ue_lambda, symmetry=symmetry)
    cput_asmbl = tmr.stop()

    # if verbosity >= 2:
    #     print(pSolver)
    #     # print(pSolver, file=output)

    tmr.start()
    u,ok = pSolver.solve(method=method, rtol=1e-10)
    assert ok == 0
    cput_solve = tmr.stop()

    diff = np.abs(ue - u)
    rmse = np.sqrt(np.sum(np.square(diff)) / M)
    mean_diff = float(np.mean(diff))
    max_diff = float(np.max(diff))
    max_rel = float(np.max(diff/u))

    s = f"{M=} : {rmse=} {mean_diff=} {max_diff=} {cput_asmbl:.5f}s {cput_solve:.5f}s"
    print(s, file=output)
    if verbosity:
        print(s)

    return rmse, mean_diff, max_diff, max_rel, cput_asmbl, cput_solve, u


def test_MF_ND():
    """test GeneralizedPoissonSolverLO"""
    dim = 3
    if dim == 1:
        reduced_cases = [
            False,
            True,
        ]
        h_cases = [
            1.6,
        ]
        analytical_solution_cases = [
            gauss,
            x_gauss,
        ]
        n_iter = 4
    elif (dim == 2):
        reduced_cases = [
            True,
            False,
            (True, False),
            (False, True),
        ]
        h_cases = [
            1.6,
            (1.6, 1.601)
        ]
        analytical_solution_cases = [
            gauss,
            x_gauss,
            y_gauss,
            xy_gauss,
        ]
        n_iter = 4
    else: # dim ==3
        reduced_cases = [
            False,
            True,
            (True, False, False),
            (False, True, False),
            (False, False, True),
            (True, True, False),
            (True, False, True),
            (False, True, True),
        ]
        h_cases = [
            4.,
            (4, 4.01, 4.01),
        ]
        analytical_solution_cases = [
            gauss,
            x_gauss,
            y_gauss,
            z_gauss,
            xy_gauss,
            xz_gauss,
            yz_gauss,
            xyz_gauss,
        ]
        n_iter = 4


    tbl = []
    method = "gmres"
    sigma = 2

    s = title_line(text=f'test_MF_{dim}D', char='-', width=120)
    print(s)
    print(s, file=output)
    for analytical_solution in [
        gauss,
        # x_gauss,
    ]:
        for stencil in [
            3,
            5,
        ]:
            for reduced in reduced_cases:
                for h in h_cases:
                    if dim <= 2:
                        M = 6
                    else:
                        M = 4
                        w = M * h if isinstance(h, float) else [M * hi for hi in h]

                    rmse0, mean_diff0, max_diff0 = 1e9, 1e9, 1e9
                    for iter in range(n_iter):
                        if dim <= 2:
                            M *= 2
                            h = h*0.5 if isinstance(h, float) else tuple([hi * 0.5 for hi in h])
                        else:
                            M += 4
                            h = w / M if isinstance(w, float) else tuple([wi / M for wi in w])
                        verbosity = 2 if (iter >=0) else 1

                        rmse, mean_diff, max_diff, max_rel, cput_asmbl, cput_solve, u = run_MF_ND(
                            analytical_solution=analytical_solution, sigma=sigma,
                            dim=dim, M=M, h=h, reduced=reduced,
                            stencil=stencil, method=method,
                            verbosity=verbosity,
                        )

                        tbl.append([analytical_solution.__name__, stencil, reduced, M, rmse, mean_diff, max_diff, max_rel, cput_asmbl, cput_solve])
                        assert rmse < rmse0
                        assert mean_diff < mean_diff0
                        assert max_diff < max_diff0

                        rmse0, mean_diff0, max_diff0 = rmse, mean_diff, max_diff

    s = "\n" + title_line(text='SUMMARY', char='-', width=120, above=True, below=True)
    s += tabulate(tbl
        , headers=["u(r)", "stencil", "reduced", "M", "RMSE", "Mean diff", "Max diff", "Max rel", "cput_asmlb", 'cput_solve']
        , tablefmt="simple"
    )
    print(s)
    print(s, file=output)
