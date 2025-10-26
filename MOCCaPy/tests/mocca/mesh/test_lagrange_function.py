import pytest
import numpy as np
from matplotlib import pyplot as plt

from mocca.mesh import LagrangeMesh, lagrange_function

from pathlib import Path
project_folder = Path(__file__)
while project_folder.name != 'tantalus_full':
    project_folder = project_folder.parent


def test_lagrange_function():
    d = 0.5
    mesh  = LagrangeMesh(dim=1, n=6, d=d, reduced=False)
    print(mesh.gridx)
    N = len(mesh.gridx)//2
    for xi in mesh.gridx:
        assert np.isnan(lagrange_function(xi,xi,d,N))
        fgridx = lagrange_function(mesh.gridx,xi,d,N)
        for (f,x) in zip(fgridx,mesh.gridx):
            if x == xi:
                assert np.isnan(f) # corner case due to 0/0, which by l'Hopitals rule should be 1
            else:
                print(f"{f=}")
                assert f == pytest.approx(.0, abs=1e-15)
    xi = mesh.gridx[3]
    xi_close = mesh.gridx[3] + 1e-9
    f = lagrange_function(xi_close,xi,d,N)
    print(f"{xi_close=} -> {f}")
    assert f == pytest.approx(1.0, rel=1e-6)

def test_lagrange_function_plot():
    """Draw a plot of the lagrange functions on a 1D mesh"""
    d = 0.5
    mesh  = LagrangeMesh(dim=1, n=6, d=d, reduced=False)
    N = len(mesh.gridx)//2
    x = np.linspace(-N*d, N*d, num=61)
    print(mesh.gridx)
    print(x)
    for i in range(61):
        if i < 30:
            x[i] += 1e-9
        elif i > 30:
            x[i] -= 1e-9

    fig, ax = plt.subplots()
    sumf = np.zeros_like(x)
    for i in range(6):
        x_i = mesh.gridx[i]
        f = lagrange_function(x,x_i,d,N)
        sumf += f
        ax.plot(x, f)
    fig.savefig(project_folder/"MOCCaPy/tests/mocca/mesh/lagrange_function.png")
    fig.clear()

    fig, ax = plt.subplots()
    ax.plot(x, sumf)
    fig.savefig(project_folder/"MOCCaPy/tests/mocca/mesh/sum_lagrange_function.png")
    fig.clear()
