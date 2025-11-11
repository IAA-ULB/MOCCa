import pytest
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

from mocca.mesh import LagrangeMesh
from mocca.mesh.observable import Observable

from pathlib import Path

from pandas.core.config_init import pc_width_doc

from MOCCaPy.mocca.mesh.lagrange import create_mesh

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

def bf(mesh, i, r, i_sign=1):
    two_pi_K = 2. * np.pi * ( mesh.gx[i] if not mesh.reduced[0] else
                              mesh.gx[i] * i_sign
                            ) / (mesh.box_width[0] * mesh.d[0])

    factor = np.sqrt(1. / mesh.box_width[0])

    arg = two_pi_K * r
    result = np.empty((r.shape[0], 2), dtype=float, order='F')
    result[:, 0] = np.cos(arg)
    result[:, 1] = np.sin(arg)
    result *= factor

    return result


def dbf_dx(mesh, i, r, i_sign=1):
    two_pi_K = 2. * np.pi * ( mesh.gx[i] if not mesh.reduced[0] else
                              mesh.gx[i] * i_sign
                            ) / (mesh.box_width[0] * mesh.d[0])

    factor = 1. / np.sqrt(mesh.box_width[0])

    arg = two_pi_K * r
    result = np.empty((r.shape[0], 2), dtype=float, order='F')
    result[:, 0] = -np.sin(arg)
    result[:, 1] =  np.cos(arg)
    result *= factor
    # By leaving out the factor two_pi_K, the derivative has the same
    # y_full-span as the basis function
    result *= -1 if two_pi_K < 0 else 1

    return result


def test_Observable_ctor():
    mesh = LagrangeMesh(dim=1, M=4, d=1., reduced=False)
    Q = Observable(mesh, n_components=1)

    data = np.zeros((mesh.linear_size, 2), order='F')


def test_differentiate_1D_x_non_reduced(debug=False):
    N = 3
    d = 1.

    r = np.linspace(-N * d, N * d, num=241)

    mesh = LagrangeMesh(dim=1, M=2*N, d=d, reduced=False, highest_derivative_order=1)
    D1 = mesh._get_D(axis=0, order=1)
    pw = mesh.basis_function(0, r)
    two_pi_K = abs(2. * np.pi * mesh.gx[0]) / (mesh.box_width[0] * mesh.d[0])

    # Validate D1 against a literal coding of the Ryssens et al 2015 eq 18
    # (correcting the sign error in that formula)
    for i in range(2*N):
        for j in range(2*N):
            if i == j:
                assert D1[j,i] == 0
            else:
                # d1ji = (-1)**(i-j) * np.pi / (2 * N * d * np.sin(np.pi*(i-j)/(2*N))) # (i-j) is wrong
                d1ji = (-1)**(i-j) * np.pi / (2 * N * d * np.sin(np.pi*(j-i)/(2*N)))   # (j-i) is correct
                print(f"D1[{j},{i}]: {D1[j,i]=}, {d1ji=}, diff={D1[j,i]-d1ji}")
                assert D1[j,i] == pytest.approx(d1ji), f"D1[{j},{i}]: {D1[j,i]=}, {d1ji=}, diff={D1[j,i]-d1ji}"

    fig, ax = plt.subplots()
    plt.title(f'{0}, x_i={mesh.gx[0]}')
    ax.plot([-N,N],[0,0],color='k')
    ax.plot([0,0],[-1,1],color='k')
    ax.plot(mesh.gx,  2*N*[0], 'rs')
    ax.plot(mesh.gx[0], [0], 'gs') # mark the mesh point to which the basis function corresponds
    ax.plot(r, pw[:,0], label='real')
    ax.plot(r, pw[:,1], label='imag')

    mbf = mesh.basis_function(ijk=0, r=mesh.gx)
    Q = Observable(mesh, data=mbf, symmetry=[1,-1])
    # check Q against bf() above
    Q_expected = bf(mesh, i=0, r=mesh.gx)
    ax.plot(mesh.gx, Q_expected[:,0], 'bo', label='bf real')
    assert Q.data == pytest.approx(Q_expected)
    # for i in range(2*N):
    #     print(f"{Q.data[i,0]} {Q_expected[i,0]} {Q.data[i,1]} {Q_expected[i,1]} ")
    dQdx = Q.differentiate(axes='x_full')

    # manually perform dQdx = D1 * Q.data - real part
    DQR = np.zeros((2*N,), dtype=float)
    for i in range(2*N): # loop over the rows of D1
        for j in range(2*N): # scalar product of D[i,:] and Q.data[:,0]
            DQR[i] += D1[i,j] * Q.data[j,0]

    for i in range(2 * N):
        print(f"{i=}: {DQR[i]} {dQdx[i,0]}")
        assert DQR[i] == pytest.approx(dQdx[i,0]), f"{i=}: {DQR[i]=}, {dQdx[i,0]=}"

    # manually perform dQdx = D1 * Q.data - imag part
    DQI = np.zeros((2*N,),dtype=float)
    for i in range(2*N): # loop over the rows of D1
        for j in range(2*N): # scalar product of D[i,:] and Q.data[:,1]
            DQI[i] += D1[i,j] * Q.data[j,1]

    for i in range(2 * N):
        print(f"{i=}: {DQI[i]} {dQdx[i,1]}")
        assert DQI[i] == pytest.approx(dQdx[i,1]), f"{i=}: {DQI[i]=}, {dQdx[i,1]=}"

    # check dQdx against dbf_dx
    dQdx_expected = dbf_dx(mesh, i=0, r=mesh.gx)
    #   the derivative apart from a factor that makes sure that the real part of dbf_dx conincides
    #   with the imaginary part of Q (because the latter is proportional to the derivative of the
    #   real part of Q).
    ax.plot(mesh.gx, dQdx_expected[:, 0], 'yo--', label='dbf/dx real')
    ax.plot(mesh.gx, dQdx_expected[:, 0]*two_pi_K, 'yo', label='dbf/dx real*')

    # ax.plot(mesh.gx, dQdx_expected[:,1], '--', label='dbf/dx imag')
    for i in range(2*N):
        print(f"{i=}: real {dQdx[i,0]} {dQdx_expected[i,0]}")
        assert dQdx[i, 0] == pytest.approx(two_pi_K*dQdx_expected[i, 0])

        print(f"{i=}: imag {dQdx[i,1]} {dQdx_expected[i,1]}")
        assert dQdx[i, 1] == pytest.approx(two_pi_K*dQdx_expected[i, 1])

    ratio = dQdx[:, 0] / Q.data[:, 1]
    print(f"dQdxR/QI={dQdx[:, 0] / Q.data[:, 1]}")
    assert ratio[0] > 0
    for rt in ratio:
        assert rt == pytest.approx(ratio[0])

    ratio = dQdx[:, 1] / Q.data[:, 0]
    print(f"dQdxI/QR={ratio}")
    assert ratio[0] < 0
    for rt in ratio:
        assert rt == pytest.approx(ratio[0])

    # ax.plot(mesh.gx, Q_expected[:,0], 'b--')
    # ax.plot(mesh.gx, Q_expected[:,1], 'r--')

    # ax.plot(mesh.gx, Q.data[:,0], 'co', label='Q real')
    # ax.plot(mesh.gx, Q.data[:,1], 'c*', label='Q imag')
    ax.plot(mesh.gx, dQdx[:,0], 'yx', label='dQ/dx real')
    # ax.plot(mesh.gx, dQdx[:,1], '*', label='dQ/dx imag')

    plt.legend()
    fig.savefig(png_folder/f"1D_basis_function_{0}_differentiation.png")
    plt.close(fig)

    # print(f"Q={Q.data}")
    # print(f"dQ/dx={dQdx}")

    print("test_differentiate_1D_x_non_reduced finished")


def test_differentiate_1D_x_reduced(debug=False):
    N = 3
    d = 1.

    mesh    = LagrangeMesh(dim=1, M=2*N, d=d, reduced=True , highest_derivative_order=1)
    mesh_nr = LagrangeMesh(dim=1, M=2*N, d=d, reduced=False, highest_derivative_order=1)

    bfr  = mesh   .basis_function(ijk=0, r=mesh   .gx)
    bfnr = mesh_nr.basis_function(ijk=3, r=mesh_nr.gx)
    for i in range(2*N):
        if i<N:
            print(f"bf{i} real: {bfnr[i,0]}")
        else:
            print(f"bf{i} real: {bfnr[i,0]} =? {bfr[i-N,0]}")
            assert bfnr[i, 0] == pytest.approx(bfr[i-N, 0])

    for i in range(2*N):
        if i<N:
            print(f"bf{i} imag: {bfnr[i,1]}")
        else:
            print(f"bf{i} imag: {bfnr[i,1]} =? {bfr[i-N,1]}")
            assert bfnr[i, 1] == pytest.approx(bfr[i-N, 1])

    Q = Observable(mesh   , data=bfr , symmetry=[1,-1] )
    R = Observable(mesh_nr, data=bfnr, symmetry=[1,-1]) # R for "R"eference.

    dQdx = Q.differentiate(axes='x_full')
    dRdx = R.differentiate(axes='x_full')

    for i in range(2*N):
        if i<N:
            print(f"dQ{i} real: {dRdx[i, 0]}")
        else:
            print(f"dQ{i} real: {dRdx[i, 0]} =? {dQdx[i-N, 0]}")
            assert dRdx[i, 0] == pytest.approx(dQdx[i-N, 0])

    for i in range(2 * N):
        if i < N:
            print(f"{i} imag: {dRdx[i, 1]}")
        else:
            print(f"{i} imag: {dRdx[i, 1]} =? {dQdx[i - N, 1]}")
            assert dRdx[i, 1] == pytest.approx(dQdx[i-N, 1])

    print("test_differentiate_1D_x_non_reduced finished")


def test_differentiate_1D_xx():
    N = 3
    d = 1.
    r = np.linspace(-N * d, N * d, num=241)

    # TODO: True case below
    for reduced in [
        False,
        # True,
    ]:
        mesh = LagrangeMesh(dim=1, M=2*N, d=d, reduced=reduced)

        if reduced == False:
            bfq = mesh.basis_function(ijk=0, r=mesh .gx)
            Q = Observable(mesh, data=bfq, symmetry=[1,-1])
            d2Qdx2 = Q.differentiate(axes='xx')

            # Re{d2Qdx2} is proportional to -Re{Q}
            # Im{d2Qdx2} is proportional to -Im{Q}
            for i in range(2 * N):
                print(f"real {i=}: {d2Qdx2[i, 0] / Q.data[i, 0]}")
                assert d2Qdx2[i,0] / Q.data[i,0] == pytest.approx(d2Qdx2[0, 0] / Q.data[0, 0]), f"real {i=}: {d2Qdx2[i, 0] / Q.data[i, 0]}"

            for i in range(2 * N):
                print(f"imag {i=}: {d2Qdx2[i, 1] / Q.data[i, 1]}")
                assert d2Qdx2[i,1] / Q.data[i,1] == pytest.approx(d2Qdx2[0, 1] / Q.data[0, 1]), f"imag {i=}: {d2Qdx2[i, 1] / Q.data[i, 1]}"

        else:
            # reduced == True case
            # Assert that a reduce mesh and a non-reduced mesh with the same parameters
            # yield identical D1 matrices.
            D1_reduced = mesh._get_D(axis=0, order=1)
            mesh_non_reduced = LagrangeMesh(dim=1, M=2*N, d=d, reduced=False, highest_derivative_order=1)
            D1_non_reduced = mesh_non_reduced._get_D(axis=0, order=1)
            for i in range(2*N):
                for j in range(2*N):
                    assert D1_reduced[i,j] == pytest.approx(D1_non_reduced[i,j])

    print("test_differentiate_1D finished")


def test_differentiate_1D_xxx():
    N = 3
    d = 1.
    # TODO: True case below
    for reduced in [
        False,
        # True,
    ]:
        mesh = LagrangeMesh(dim=1, M=2*N, d=d, reduced=reduced, highest_derivative_order=3)

        if reduced == False:
            bfq = mesh.basis_function(ijk=0, r=mesh.gx)
            Q = Observable(mesh, data=bfq, symmetry=[1,-1])
            d3Qdx3 = Q.differentiate(axes='xxx')

            # Re{d3Qdx3} is proportional to Im{Q}
            # Im{d3Qdx3} is proportional to Re{Q}
            for i in range(1,2 * N):
                print(f"real {i=}: {d3Qdx3[i, 0] / Q.data[i, 1]}")
                assert d3Qdx3[i, 0] / Q.data[i, 1] == pytest.approx(d3Qdx3[0, 0] / Q.data[0, 1]), f"real {i=}: {d3Qdx3[i, 0] / Q.data[i, 1]} {d3Qdx3[0, 0] / Q.data[0, 1]}"

            for i in range(1,2 * N):
                print(f"imag {i=}: {d3Qdx3[i, 1] / Q.data[i, 0]}")
                assert d3Qdx3[i, 1] / Q.data[i, 0] == pytest.approx(d3Qdx3[0, 1] / Q.data[0, 0]), f"imag {i=}: {d3Qdx3[i, 1] / Q.data[i, 1]} {d3Qdx3[0, 1] / Q.data[0, 0]}"

        else:
            # reduced == True case
            # Assert that a reduce mesh and a non-reduced mesh with the same parameters
            # yield identical D1 matrices.
            D1_reduced = mesh._get_D(axis=0, order=1)
            mesh_non_reduced = LagrangeMesh(dim=1, M=2*N, d=d, reduced=False, highest_derivative_order=1)
            D1_non_reduced = mesh_non_reduced._get_D(axis=0, order=1)
            for i in range(2*N):
                for j in range(2*N):
                    assert D1_reduced[i,j] == pytest.approx(D1_non_reduced[i,j])

    print("test_differentiate_1D finished")


def test_differentiate_2D_Hessian(debug=False):
    N = 3
    d = 1.
    dim = 2
    for reduced in [
        False,
        # True,
    ]:
        mesh = LagrangeMesh(dim=dim, M=2*N, d=1., reduced=reduced)
        bfq = mesh.basis_function(ijk=(0,0), r=mesh.grid)
        Q = Observable(mesh, data=bfq)
        H = Q.differentiate(axes='Hessian',debug=debug)
        print("test_differentiate_2D_Hessian finished")

def test_differentiate_3D_Hessian_non_reduced(debug=False):
    N = 3
    d = 1.
    dim = 3

    mesh = LagrangeMesh(dim=dim, M=2*N, d=1., reduced=False)
    bfq = mesh.basis_function(ijk=(0,0,0), r=mesh.grid)
    Q = Observable(mesh, data=bfq)
    Q.differentiate(axes='Hessian', debug=debug)

    with pytest.raises(ValueError):
        Q.differentiate(axes='Gradient', recompute=False, debug=debug)

    Q.differentiate(axes='Grad', recompute=False, debug=debug)

    print("test_differentiate_3D_Hessian_non_reduced finished")

def compare(desc, a, b, doassert=True):
    ok = (a == pytest.approx(b))
    s = f"{desc} : {a} =? {b} -> {ok}"

    # print(s)
    if doassert:
        assert a == pytest.approx(b), s

def test_differentiate_2D_Grad_reduced(debug=False):
    N = 3
    M = 2*N
    d = 1.
    dim = 2

    ############################################################################
    # This monster has bitten me once again!
    # Basis functions are NOT symmetric/skew-symmetric in 2D and 3D
    ############################################################################

    # We take a plane wave in the x_full-direction
    k = 1.5
    L = 2*N*d
    axes = 'x'

    mesh_2N = LagrangeMesh(dim=dim, M=M, d=1., reduced=False)
    rx = np.linspace(-L/2, L/2, 100)
    rxy = create_mesh(rx,mesh_2N.gy)

    # r = mesh_2N.gridx[:,0,0].reshape((M*M,), order='F')
    pw_full = mesh_2N.plane_wave(L=L, k=k, r=rxy[:,0])

    x_full = rxy[:,0].reshape((100,M), order='F')
    y_full = rxy[:,1].reshape((100,M), order='F')
    pw_full_real = pw_full[:, 0].reshape((100,M), order='F')
    pw_full_imag = pw_full[:, 1].reshape((100,M), order='F')

    if debug:
        fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
        plt.title(f'plane wave x_full real')
        ax.plot_surface(x_full, y_full, pw_full_real, cmap=cm.coolwarm, linewidth=0, antialiased=False)
        # fig.savefig(png_folder / f"2D_basis_function_{i}_real.png")
        plt.show()
        plt.close(fig)
    
        fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
        plt.title(f'plane wave x_full imag')
        ax.plot_surface(x_full, y_full, pw_full_imag, cmap=cm.coolwarm, linewidth=0, antialiased=False)
        # fig.savefig(png_folder / f"2D_basis_function_{i}_real.png")
        plt.show()
        plt.close(fig)

    mesh_N = LagrangeMesh(dim=dim, M=2*N, d=1., reduced=True, highest_derivative_order=1)
    pw_redu_real = pw_full_real[50:,3:]
    pw_redu_imag = pw_full_imag[50:,3:]
    x_redu = x_full[50:, 3:]
    y_redu = y_full[50:, 3:]

    pw_points = mesh_2N.plane_wave(L=L, k=k, r=mesh_2N.gridx.reshape((M*M,), order='F'))

    if debug:
        pw_points_real = pw_points[:, 0].reshape((M,M,), order='F')
        pw_points_imag = pw_points[:, 1].reshape((M,M,), order='F')
        fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
        plt.title(f'plane wave x_full real')
        ax.plot_surface(x_redu, y_redu, pw_redu_real, cmap=cm.coolwarm, linewidth=0, antialiased=False)

        for j in range(N,M):
            for i in range(N,M):
                ax.scatter(mesh_2N.gridx[i,j],mesh_2N.gridy[i,j],pw_points_real[i,j],marker='o',c='k')
        # # fig.savefig(png_folder / f"2D_basis_function_{i}_real.png")
        plt.show()
        plt.close(fig)

        fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
        plt.title(f'plane wave x_full imag')
        ax.plot_surface(x_redu, y_redu, pw_redu_imag, cmap=cm.coolwarm, linewidth=0, antialiased=False)
        for j in range(N,M):
            for i in range(N,M):
                ax.scatter(mesh_2N.gridx[i,j],mesh_2N.gridy[i,j],pw_points_imag[i,j],marker='o',c='k')
        # fig.savefig(png_folder / f"2D_basis_function_{i}_real.png")
        plt.show()
        plt.close(fig)

    QN = Observable(mesh_N , data=mesh_N.cast2linear(mesh_2N.cast2mesh(pw_points)[N:M, N:M, :]), symmetry=[1,-1], name='QN')
    Q2N = Observable(mesh_2N, data=mesh_2N.cast2linear(pw_points)                                  , symmetry=[1,-1], name='Q2N')

    QN.differentiate(axes=axes, debug=debug)
    Q2N.differentiate(axes=axes, debug=debug)

    # verify that the full and reduced mesh both yield the same derivative
    # verifying the x-direction
    print(f"Verifying {axes=}")
    dQN = mesh_N.cast2mesh(QN.derivatives['x'])
    dQ2N = mesh_2N.cast2mesh(Q2N.derivatives['x'])
    for i in range(N):
        for j in range(N):
                compare(f"Q real ({i},{j})",
                        dQN[  i,  j,0],
                        dQ2N[N+i,N+j,0],
                        doassert=True
                )
    for i in range(N):
        for j in range(N):
                compare(f"Q imag ({i},{j})",
                        dQN[  i,  j,1],
                        dQ2N[N+i,N+j,1],
                        doassert=True
                )

    # Now transpose the plane wave so that it goes in the y-direction and test the y-derivative:
    axes = 'y'
    pw_points = mesh_2N.cast2mesh(pw_points) # (M,M,2)
    pw_y = np.empty_like(pw_points)
    for i in range(M):
        for j in range(M):
            pw_y[i,j,:] = pw_points[j,i,:]

    print(f"Verifying {axes=}")
    QN = Observable(mesh_N, data=mesh_N.cast2linear(pw_y[N:M, N:M, :]), symmetry=[1,-1], name='QN')
    Q2N = Observable(mesh_2N, data=mesh_2N.cast2linear(pw_y)             , symmetry=[1,-1], name='Q2N')

    QN.differentiate(axes=axes, debug=debug)
    Q2N.differentiate(axes=axes, debug=debug)

    # verify that the full and reduced mesh both yield the same derivative
    dQN = mesh_N.cast2mesh(QN.derivatives['y'])
    dQ2N = mesh_2N.cast2mesh(Q2N.derivatives['y'])
    for i in range(N):
        for j in range(N):
                compare(f"Q real ({i},{j})",
                        dQN[  i,  j,0],
                        dQ2N[N+i,N+j,0],
                        doassert=True
                )
    for i in range(N):
        for j in range(N):
                compare(f"Q imag ({i},{j})",
                        dQN[  i,  j,1],
                        dQ2N[N+i,N+j,1],
                        doassert=True
                )

    print("test_differentiate_2D_Grad_non_reduced finished")


def test_differentiate_3D_Grad_reduced(debug=False):
    N = 3
    M = 2 * N
    d = 1.
    dim = 3

    ############################################################################
    # This monster has bitten me once again!
    # Basis functions are NOT symmetric/skew-symmetric in 2D and 3D
    ############################################################################

    # We take a plane wave in the x_full-direction
    k = 1.5
    L = 2 * N * d
    axes = 'x'

    mesh_2N = LagrangeMesh(dim=dim, M=M, d=1., reduced=False, highest_derivative_order=1)
    mesh_N  = LagrangeMesh(dim=dim, M=M, d=1., reduced=True , highest_derivative_order=1)

    pw_x = mesh_2N.plane_wave(L=L, k=k, r=mesh_2N.gridx.reshape((M * M * M,), order='F'))
    pw_mesh = mesh_2N.cast2mesh(pw_x)
    print(pw_mesh[N:M,0,0,0])
    print(pw_mesh[N:M,0,0,1])

    QN  = Observable(mesh_N , data=mesh_N .cast2linear(pw_mesh[N:M, N:M, N:M, :]), symmetry=[1, -1], name='QN')
    Q2N = Observable(mesh_2N, data=mesh_2N.cast2linear(pw_mesh)                  , symmetry=[1, -1], name='Q2N')

    QN .differentiate(axes=axes, debug=debug)
    Q2N.differentiate(axes=axes, debug=debug)

    # verify that the full and reduced mesh both yield the same derivative
    print(f"Verifying {axes=}")
    dQN  = mesh_N .cast2mesh(QN .derivatives[axes])
    dQ2N = mesh_2N.cast2mesh(Q2N.derivatives[axes])
    for iq in range(2):
        part = 'real' if iq == 0 else 'imag'
        for i in range(N):
            for j in range(N):
                for k in range(N):
                    compare(f"Q {part} ({i},{j},{k})",
                        dQN [  i,   j,   k, iq],
                        dQ2N[N+i, N+j, N+k, iq],
                        doassert=True
                        )

    # Now transpose the plane wave so that it goes in the y-direction and test the y-derivative:
    axes = 'y'

    pw_y = np.empty_like(pw_mesh)
    for i in range(M):
        for j in range(M):
            for j in range(M):
                for k in range(M):
                    pw_y[i, j, k, :] = pw_mesh[j, i, k, :]
    print(pw_y[0,N:M,0,0])
    print(pw_y[0,N:M,0,1])
    assert pw_mesh[N:M, 0, 0, 0] == pytest.approx(pw_y[0,N:M,0,0])
    assert pw_mesh[N:M, 0, 0, 1] == pytest.approx(pw_y[0,N:M,0,1])

    print(f"Verifying {axes=}")
    QN  = Observable(mesh_N , data=mesh_N .cast2linear(pw_y[N:M, N:M, N:M, :]), symmetry=[1, -1], name='QN')
    Q2N = Observable(mesh_2N, data=mesh_2N.cast2linear(pw_y)                  , symmetry=[1, -1], name='Q2N')

    QN .differentiate(axes=axes, debug=debug)
    Q2N.differentiate(axes=axes, debug=debug)

    # verify that the full and reduced mesh both yield the same derivative
    dQN  = mesh_N .cast2mesh(QN .derivatives['y'])
    dQ2N = mesh_2N.cast2mesh(Q2N.derivatives['y'])
    for iq in range(2):
        part = 'real' if iq == 0 else 'imag'
        for i in range(N):
            for j in range(N):
                for k in range(N):
                    compare(f"Q {part} ({i},{j},{k})",
                        dQN [  i,   j,   k, iq],
                        dQ2N[N+i, N+j, N+k, iq],
                        doassert=True
                        )

    # Now transpose the plane wave so that it goes in the z-direction and test the yz-derivative:
    axes = 'z'

    pw_z = np.empty_like(pw_mesh)
    for i in range(M):
        for j in range(M):
            for k in range(M):
                pw_z[i, j, k, :] = pw_mesh[k, j, i, :]
    print(pw_z[0,0,N:M,0])
    print(pw_z[0,0,N:M, 1])
    assert pw_mesh[N:M, 0, 0, 0] == pytest.approx(pw_z[0,0,N:M,0])
    assert pw_mesh[N:M, 0, 0, 1] == pytest.approx(pw_z[0,0,N:M,1])

    print(f"Verifying {axes=}")
    QN  = Observable(mesh_N , data=mesh_N .cast2linear(pw_z[N:M, N:M, N:M, :]), symmetry=[1, -1], name='QN')
    Q2N = Observable(mesh_2N, data=mesh_2N.cast2linear(pw_z)                  , symmetry=[1, -1], name='Q2N')

    QN .differentiate(axes=axes, debug=debug)
    Q2N.differentiate(axes=axes, debug=debug)

    # verify that the full and reduced mesh both yield the same derivative
    dQN  = mesh_N .cast2mesh(QN .derivatives['z'])
    dQ2N = mesh_2N.cast2mesh(Q2N.derivatives['z'])
    for iq in range(2):
        part = 'real' if iq == 0 else 'imag'
        for i in range(N):
            for j in range(N):
                for k in range(N):
                    compare(f"Q {part} ({i},{j},{k})",
                        dQN [  i,   j,   k, iq],
                        dQ2N[N+i, N+j, N+k, iq],
                        doassert=True
                        )

    print("test_differentiate_3D_Hessian_non_reduced finished")

