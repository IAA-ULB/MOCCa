import numpy as np
import scipy.sparse as sparse
import sympy
import pytest
from tabulate import tabulate

from MOCCaPy.mocca.mean_field.solvers.generalized_poisson import GeneralizedPoissonSolver, dia_entry_set, dia_entry_add
from mocca.mesh import LagrangeMesh, MeshQuantity


def gauss(dim, sigma):
    """A Gaussian function $g$ centered around the origin in `dim` dimensions with standard deviation `sigma`,
    and its Laplacian $f = \Delta g$

    Args:
         dim: number of spatial dimensions.
         sigma: standard deviation.
    Returns:
        ue_lambda, fe_lambda:
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

    return ue_lambda, fe_lambda


def test_1D_not_reduced():
    tbl = []

    reduced = False
    dim = 1
    sigma = 2
    ue_lambda, fe_lambda = gauss(dim=dim, sigma=sigma)

    M = 6
    h = 1.6
    for iter in range(10):
        M *= 2
        h *= .5

        mesh = LagrangeMesh(dim=dim, M=M, d=h, reduced=reduced, collect_boundary_points=True)
        f = MeshQuantity(mesh, name='f', n_components=1)
        ue = MeshQuantity(mesh, name='u', n_components=1)
        ue.data[:,0] = ue_lambda(mesh.g1D[0])

        for bc_scale in [None, 1e20]:
            f.data[:,0] = fe_lambda(mesh.g1D[0])
            pSolver = GeneralizedPoissonSolver(f,bc_scale=bc_scale)
            pSolver.assemble_rhs(ue_lambda)

            if iter == 0:
                # create a table to compare with poisson-bis
                t = []
                A = pSolver.L.toarray()
                b = pSolver.f.data[:,0]
                for i in range(M):
                    ti = [i]
                    for j in range(M):
                        ti.append(A[i,j])
                    ti.append(b[i])
                    t.append(ti)

            u = pSolver.solve()
            diff = np.abs(ue.data[:,0] - u)

            if iter == 0:
                for i in range(M):
                    t[i].extend([u[i], ue.data[i,0], diff[i]])
                headers = ['i']
                headers.extend([str(i) for i in range(M)])
                headers.extend(['f', 'u', 'ue', 'diff'])
                print(tabulate(t, tablefmt="simple", headers=headers))

            rmse = np.sqrt(np.sum(np.square(diff)) / M)
            mean_diff = float(np.mean(diff))
            max_diff = float(np.max(diff))
            print(f"{bc_scale=} {M=} : {rmse=} {mean_diff=} {max_diff=} ")
            tbl.append([bc_scale, M, rmse, mean_diff, max_diff])

    print(tabulate(tbl,  tablefmt="simple", headers=["bc_scale", "M", "RMSE", "Mean diff", "Max diff"]))


def test_2D_not_reduced():
    tbl = []

    reduced = False
    dim = 2
    sigma = 2
    ue_lambda, fe_lambda = gauss(dim=dim, sigma=sigma)

    M = 5
    h = 2
    for iter in range(5):
        M *=  2
        h *= .5
        assert M*h == 10.
        mesh = LagrangeMesh(dim=dim, M=M, d=h, reduced=reduced, collect_boundary_points=True)
        lx = M*h/2-h/2
        lx = [-lx, lx]
        for i in range(mesh.bp_l.size):
            assert (mesh.bp_xyz[i,0] in lx) or \
                   (mesh.bp_xyz[i,1] in lx)

        ue = MeshQuantity(mesh, name='u', n_components=1)
        ue.data[:,0] = ue_lambda(mesh.gridx.ravel(), mesh.gridy.ravel())

        f  = MeshQuantity(mesh, name='f', n_components=1)
        for bc_scale in [None, 1e20]:
            f .data[:,0] = fe_lambda(mesh.gridx.ravel(), mesh.gridy.ravel())
            pSolver = GeneralizedPoissonSolver(f,bc_scale=bc_scale)
            pSolver.assemble_rhs(ue_lambda)

            u = pSolver.solve()
            diff = np.abs(ue.data[:,0] - u)

            rmse = np.sqrt(np.sum(np.square(diff)) / M)
            mean_diff = float(np.mean(diff))
            max_diff = float(np.max(diff))
            print(f"{M=} : {rmse} {mean_diff=} {max_diff=} ")
            tbl.append([M, rmse, mean_diff, max_diff])

    print(tabulate(tbl,  tablefmt="simple", headers=["M", "RMSE", "Mean diff", "Max diff"]))


def test_3D_not_reduced():
    tbl = []

    reduced = False
    dim = 3
    sigma = 2
    ue_lambda, fe_lambda = gauss(dim=dim, sigma=sigma)

    M = 5
    h = 2
    for iter in range(5):
        M *=  2
        h *= .5
        assert M*h == 10.
        mesh = LagrangeMesh(dim=dim, M=M, d=h, reduced=reduced, collect_boundary_points=True)
        lx = M*h/2-h/2
        lx = [-lx, lx]
        for i in range(mesh.bp_l.size):
            assert (mesh.bp_xyz[i,0] in lx) or \
                   (mesh.bp_xyz[i,1] in lx) or \
                   (mesh.bp_xyz[i,2] in lx)

        ue = MeshQuantity(mesh, name='u', n_components=1)
        ue.data[:,0] = ue_lambda(mesh.gridx.ravel(), mesh.gridy.ravel(), mesh.gridz.ravel())

        f  = MeshQuantity(mesh, name='f', n_components=1)
        for bc_scale in [None, 1e20]:
            f.data[:,0] = fe_lambda(mesh.gridx.ravel(), mesh.gridy.ravel(), mesh.gridz.ravel())
            pSolver = GeneralizedPoissonSolver(f,bc_scale=bc_scale)
            pSolver.assemble_rhs(ue_lambda)

            u = pSolver.solve()
            diff = np.abs(ue.data[:,0] - u)

            rmse = np.sqrt(np.sum(np.square(diff)) / M)
            mean_diff = float(np.mean(diff))
            max_diff = float(np.max(diff))
            print(f"{M=} : {rmse} {mean_diff=} {max_diff=} ")
            tbl.append([M, rmse, mean_diff, max_diff])

    print(tabulate(tbl,  tablefmt="simple", headers=["M", "RMSE", "Mean diff", "Max diff"]))


def test_set_dia_entry():
    D = sparse.eye(4, k=1, format='dia')
    D44 = D.toarray()
    assert D44[0,1] == 1
    dia_entry_set(D, [(0, 1, 2)])
    D44 = D.toarray()
    assert D44[0,1] == 2
    with pytest.raises(ValueError):
        dia_entry_set(D, [(0, 2, 2)])

    D = sparse.eye(4, k=2, format='dia')
    D44 = D.toarray()
    assert D44[0, 2] == 1
    dia_entry_set(D, [(0, 2, 2)])
    D44 = D.toarray()
    assert D44[0, 2] == 2

    D = sparse.eye(4, k=-1, format='dia')
    D44 = D.toarray()
    assert D44[1, 0] == 1
    dia_entry_set(D, [(1, 0, 2),
                      (2, 1, 3)
                     ])
    D44 = D.toarray()
    assert D44[1,0] == 2
    assert D44[2,1] == 3

    dia_entry_add(D, [
        (1, 0, 2),
        (2, 1, 3),
    ])
    D44 = D.toarray()
    assert D44[1,0] == 4
    assert D44[2,1] == 6

    D = sparse.eye(4, k=-2, format='dia')
    D44 = D.toarray()
    assert D44[2, 0] == 1
    dia_entry_set(D, [
        (2, 0, 2),
    ])
    D44 = D.toarray()
    assert D44[2, 0] == 2


def test_collect_boundary_points_1D():
    dim = 1
    M = 10
    d = 1
    for reduced in [True, False]:
        mesh = LagrangeMesh(dim=dim, M=M, d=d, reduced=reduced)
        mesh.collect_boundary_points()
        if reduced:
            assert len(mesh.bp_l) == 1
            assert mesh.bp_l[0] == 4
            assert mesh.bp_xyz[0] == M/2 - d/2
        else:
            assert len(mesh.bp_l) == 2
            assert mesh.bp_l[0] == 0
            assert mesh.bp_xyz[0] == -(M/2 - d/2)
            assert mesh.bp_l[1] == M-1
            assert mesh.bp_xyz[1] ==  (M/2 - d / 2)
        pass


def test_collect_boundary_points_2D():
    dim = 2
    M = (4,6)
    d = 1
    xb = (M[0]-d)/2
    yb = (M[1]-d)/2
    x = [-xb, xb]
    y = [-yb, yb]

    for reduced in [
        (True, True),
        (False, False),
        (True, False),
        (False, True),
    ]:
        print(f"{reduced=}")
        mesh = LagrangeMesh(dim=dim, M=M, d=d, reduced=reduced)
        mesh.collect_boundary_points()
        nbp = mesh.bp_l.size
        if reduced == (True, True):
            assert nbp == 4
        elif reduced == (False, False):
            assert nbp == 16
        else:
            assert nbp == 8

        for i in range(nbp):
            assert mesh.bp_xyz[i,0] in x or mesh.bp_xyz[i,1] in y


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
        mesh.collect_boundary_points()
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
