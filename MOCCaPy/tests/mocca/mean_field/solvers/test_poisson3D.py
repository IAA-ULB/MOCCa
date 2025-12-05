import pytest
import numpy as np
import sympy

from mocca.mean_field.solvers.poisson3D import poisson3D
from mocca.mesh import LagrangeMesh, MeshQuantity
from mocca.util.timer import Timer

def test_gauss():
    """To be able to use shenfun to solve the Poisson equation arising in MOCCa, we must find out if
    we get a good solution by interpolating the MeshQuantity with the rhs of the Poisson equation on
    the quadrature points of shenfun. This is tested here by
    1. take an analytical potential `ue`
    2. take its laplacian `fe` which is the rhs of the Poisson equation
    3. evaluate `fe` on the LagrangeMesh -> fe_lm
    4. Call poisson3D on fe_lm
    """
    for reduced in [
        False,
        True,
    ]:
        print(f"{reduced=}")
        #----------------------------------------------------------------------------------------------
        # 1. take an analytical potential `ue`
        x, y, z = sympy.symbols("x,y,z", real=True)
        sigma = 3.0
        ue = sympy.exp(-0.5*(x**2 + y**2 + z**2)/sigma**2) # a gaussian curve as the potential

        # ----------------------------------------------------------------------------------------------
        # 2. take its laplacian `fe` which is the rhs of the Poisson equation
        fe = ue.diff(x, 2) + ue.diff(y, 2) + ue.diff(z, 2) # its laplacian serves as the rhs
        fe = fe.simplify()
        # ----------------------------------------------------------------------------------------------
        # 3. evaluate `fe` on the LagrangeMesh -> fe_lm
        d = 0.8
        M = 30
        mesh = LagrangeMesh(M=30, d=d, reduced=reduced)
        fe_lm = MeshQuantity(mesh, n_components=1, symmetry=1)
        # Convert to a NumPy-based function
        f = sympy.lambdify([x, y, z], fe, "numpy")
        # Transfer the rhs to the MeshQuantity
        gridx = mesh.cast2linear(mesh.gridx)
        gridy = mesh.cast2linear(mesh.gridy)
        gridz = mesh.cast2linear(mesh.gridz)
        fe_lm.data[:] = f(gridx, gridy, gridz)
        #         ^^^ mind this

        # ----------------------------------------------------------------------------------------------
        # 4. Call poisson3D on fe_lm
        poisson3D(fe_lm, out='rhs', verbose=False)

        # transfer also ue to the mesh
        ue_lm = MeshQuantity(mesh, n_components=1, symmetry=1)
        # Convert to a NumPy-based function
        u = sympy.lambdify([x, y, z], ue, "numpy")
        # Transfer the rhs to the MeshQuantity
        ue_lm.data[:] = f(gridx, gridy, gridz)

        delta = ue_lm.data - fe_lm.data
        print(f"RMSE = {np.sqrt((delta*delta).sum()/delta.size)}")

        Timer.report(sort=False)
        Timer.reset()