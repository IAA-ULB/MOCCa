import pytest
import numpy as np
import matplotlib.pyplot as plt

from mocca.mesh import LagrangeMesh
from mocca.mesh.observable import Observable

from pathlib import Path

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
    # y-span as the basis function
    result *= -1 if two_pi_K < 0 else 1

    return result


def test_Observable_ctor():
    mesh = LagrangeMesh(dim=1, M=4, d=1., reduced=False)
    Q = Observable(mesh, n_components=1)

    data = np.zeros((mesh.n_gridpoints(), 2), order='F')

def test_differentiate_1D_x(debug=False):
    N = 3
    d = 1.
    r = np.linspace(-N * d, N * d, num=241)

    for reduced in [
        False,
        # True,
    ]:
        mesh = LagrangeMesh(dim=1, M=2*N, d=d, reduced=reduced, highest_derivative_order=1)
        D1 = mesh.D[0,0]
        pw = mesh.basis_function(0, r)
        two_pi_K = abs(2. * np.pi * mesh.gx[0]) / (mesh.box_width[0] * mesh.d[0])

        # validate D1 against a literal coding of the Ryssens et al 2015 eq 18
        for i in range(2*N):
            for j in range(2*N):
                if i == j:
                    assert D1[j,i] == 0
                else:
                    # d1ji = (-1)**(i-j) * np.pi / (2 * N * d * np.sin(np.pi*(i-j)/(2*N)))
                    d1ji = (-1)**(i-j) * np.pi / (2 * N * d * np.sin(np.pi*(j-i)/(2*N)))
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

        Q = Observable(mesh, data = mesh.basis_function(ijk=0, r=mesh.gx))
        # check Q against bf() above
        Q_expected = bf(mesh, i=0, r=mesh.gx)
        ax.plot(mesh.gx, Q_expected[:,0], 'bo', label='bf real')
        assert Q.data == pytest.approx(Q_expected)
        # for i in range(2*N):
        #     print(f"{Q.data[i,0]} {Q_expected[i,0]} {Q.data[i,1]} {Q_expected[i,1]} ")
        dQdx = Q.differentiate(axes='x')

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

        print("test_differentiate_1D finished")

# def test_differentiate_2D_x(debug=False):
#     for reduced in [
#         False,
#         # True,
#     ]:
#         mesh = LagrangeMesh(dim=2, M=4, d=1., reduced=reduced)
#         Q = Observable(mesh, data=mesh.basis_function(ijk=(0,0), r=mesh.grid))
#         Q.differentiate(axes=['x', 'y'])
#         print("test_differentiate_2D finished")

def test_order_of_D():
    pass