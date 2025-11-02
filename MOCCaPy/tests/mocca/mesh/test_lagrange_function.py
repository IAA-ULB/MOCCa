import pytest
import numpy as np
from matplotlib import pyplot as plt

from mocca.mesh.lagrange import LagrangeMesh
from mocca.mesh.lagrange_function import lagrange_function

from pathlib import Path
project_folder = Path(__file__)
while project_folder.name != 'tantalus_full':
    project_folder = project_folder.parent


def test_lagrange_function():
    d = 0.5
    mesh  = LagrangeMesh(dim=1, M=6, d=d, reduced=False)
    print(mesh.gridx)
    N = len(mesh.gridx)//2
    for xi in mesh.gridx:
        assert np.isnan(lagrange_function(xi,  xi, d, 2*N))
        fgridx = lagrange_function(mesh.gridx, xi, d, 2*N)
        for (f,x) in zip(fgridx,mesh.gridx):
            if x == xi:
                assert np.isnan(f) # corner case due to 0/0, which by l'Hopitals rule should be 1
            else:
                print(f"{f=}")
                assert f == pytest.approx(.0, abs=1e-14)
    xi = mesh.gridx[3]
    xi_close = mesh.gridx[3] + 1e-9
    f = lagrange_function(xi_close,xi,d,N)
    print(f"{xi_close=} -> {f}")
    assert f == pytest.approx(1.0, rel=1e-6)

def test_lagrange_function_plot():
    """Draw a plot of the lagrange functions on a 1D mesh"""
    d = 1.
    mesh  = LagrangeMesh(dim=1, M=6, d=d, reduced=False)
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
        figi, axi = plt.subplots()
        plt.title(f"lagrange function {i}")

        x_i = mesh.gridx[i]
        f = lagrange_function(x,x_i,d,2*N)
        sumf += f
        ax.plot(x, f)
        axi.plot(x, f)

        figi.savefig(project_folder/f"MOCCaPy/tests/mocca/mesh/lagrange_function_{i}.png")
        plt.close(figi)

    fig.savefig(project_folder/"MOCCaPy/tests/mocca/mesh/lagrange_functions.png")
    plt.close(fig)

    # A constant function cannot be interpolated by the Lagrange Functions because
    # they are antiperiodic, i.e. periodic, but with a sign change.  
    # See https://github.com/IAA-nuclear/tantalus_full/issues/52
    # I expected - wrongly, as it appears - that the sum of all Lagrange functions, 
    # which puts a one on each grid point would interpolate the constant function.
    fig, ax = plt.subplots()
    ax.plot(x, sumf)
    fig.savefig(project_folder/"MOCCaPy/tests/mocca/mesh/sum_lagrange_functions.png")
    plt.close(fig)

def test_lagrange_function_plot2():
    """Draw a plot of the lagrange functions on a 1D mesh, this time using mesh.lagrange_function"""
    d = 1.
    mesh  = LagrangeMesh(dim=1, M=6, d=d, reduced=False)
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
        figi, axi = plt.subplots()
        plt.title(f"lagrange function {i}")

        f = mesh.lagrange_function(x,i)
        sumf += f
        ax.plot(x, f)
        axi.plot(x, f)
        axi.plot(mesh.gridx, (2*N)*[0], 'sr')
        axi.plot(2*[mesh.gridx[i]],  [0,1], 'sg-')

        figi.savefig(project_folder/f"MOCCaPy/tests/mocca/mesh/lagrange_function_{i}.png")
        plt.close(figi)

    fig.savefig(project_folder/"MOCCaPy/tests/mocca/mesh/lagrange_functions.png")
    plt.close(fig)

    # A constant function cannot be interpolated by the Lagrange Functions because
    # they are antiperiodic, i.e. periodic, but with a sign change.
    # See https://github.com/IAA-nuclear/tantalus_full/issues/52
    # I expected - wrongly, as it appears - that the sum of all Lagrange functions,
    # which puts a one on each grid point would interpolate the constant function.
    fig, ax = plt.subplots()
    ax.plot(x, sumf)
    fig.savefig(project_folder/"MOCCaPy/tests/mocca/mesh/sum_lagrange_functions.png")
    plt.close(fig)
