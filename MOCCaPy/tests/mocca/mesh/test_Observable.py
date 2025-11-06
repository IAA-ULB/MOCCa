import pytest
import numpy as np
import matplotlib.pyplot as plt

from mocca.mesh import LagrangeMesh
from mocca.mesh.observable import Observable

from pathlib import Path
project_folder = Path(__file__)
while project_folder.name != 'tantalus_full':
    project_folder = project_folder.parent
test_folder = project_folder/f"MOCCaPy/tests/mocca/mesh"
png_folder = test_folder/"png"


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
        for i in range(2*N):
            for j in range(2*N):
                if i == j:
                    assert D1[j,i] == 0
                else:
                    d1ji = (-1)**(i-j) * np.pi / (2 * N * d * np.sin(np.pi*(i-j)/(2*N)))
                    print(f"D1[{j},{i}]: {D1[j,i]=}, {d1ji=}, diff={D1[j,i]-d1ji}")
                    assert D1[j,i] == pytest.approx(d1ji), f"D1[{j},{i}]: {D1[j,i]=}, {d1ji=}, diff={D1[j,i]-d1ji}"

        fig, ax = plt.subplots()
        plt.title(f'{0}, x_i={mesh.gx[0]}')
        ax.plot([-N,N],[0,0],color='k')
        ax.plot([0,0],[-1,1],color='k')
        ax.plot(mesh.gx,  2*N*[0], 'rs')
        ax.plot(mesh.gx[0], [0], 'gs')
        ax.plot(r, pw[:,0], label='real')
        ax.plot(r, pw[:,1], label='imag')

        Q = Observable(mesh, data = mesh.basis_function(ijk=0, r=mesh.gx))
        Qx = Q.differentiate(axes='x')

        DQ0 = np.zeros((6,),dtype=float)
        for i in range(2*N):
            for j in range(2*N):
                DQ0[i] += D1[i,j] * Q.data[j,0]
        # for i in range(6):
        #     DQ0[i] = np.sum(D1[i,:] * Q.data[:,0])

        for i in range(2 * N):
            print(f"{i=}: {DQ0[i]} {Qx[i,0]}")
            assert DQ0[i] == pytest.approx(Qx[i,0]), f"{i=}: {DQ0[i]=}, {Qx[i,0]=}"

        print(f"QxR/QI={Qx[:, 0] / Q.data[0, 1]}")
        print(f"Q={Qx[:, 1] / Q.data[0, 0]}")

        ax.plot(mesh.gx, Q.data[:,0], 'co', label='Q real')
        ax.plot(mesh.gx, Q.data[:,1], 'c*', label='Q imag')
        ax.plot(mesh.gx, Qx[:,0], 'yo', label='dQ/dx real')
        ax.plot(mesh.gx, Qx[:,1], 'y*', label='dQ/dx imag')

        plt.legend()
        fig.savefig(png_folder/f"1D_basis_function_{0}_differentiation.png")
        plt.close(fig)

        print(f"Q={Q.data}")
        print(f"dQ/dx={Qx}")

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