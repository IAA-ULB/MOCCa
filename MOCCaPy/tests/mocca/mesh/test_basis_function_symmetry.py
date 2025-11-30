import pytest
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib import cm

from mocca.mesh import LagrangeMesh
# from mocca.mesh.lagrange import create_mesh
# from mocca.mesh.mesh_quantity import MeshQuantity

project_folder = Path(__file__)
while project_folder.name != 'tantalus_full':
    project_folder = project_folder.parent
test_folder = project_folder/f"MOCCaPy/tests/mocca/mesh"

def test_symmetry():
    mesh = LagrangeMesh(dim=2, M=4, d=1, reduced=False)
    n_test_points = 10
    sign = np.empty((n_test_points,2),dtype=int,order='F')
    sign[:,0] = np.arange(n_test_points)
    for i in range(4):
        for j in range(4):
            print(f"bf({i},{j})")
            r = np.random.rand(mesh.dim*n_test_points).reshape((n_test_points, mesh.dim), order='F')
            bf_r = mesh.basis_function((i,j),r)
            rmx = r.copy()
            rmx[:,0] = -rmx[:,0]
            bf_rmx = mesh.basis_function((i,j), rmx)
            rmy = r.copy()
            rmy[:,1] = -rmy[:,1]
            bf_rmy = mesh.basis_function((i,j), rmy)
            for itp in range(n_test_points):
                for ic in range(2):
                    if bf_rmx[itp,ic] == bf_r[itp,ic]:
                        sign[itp,ic] = 1
                    elif bf_rmx[itp, ic] == -bf_r[itp, ic]:
                        sign[itp, ic] = -1
                    else:
                        sign[itp, ic] = 0
                print(f"{r[itp,:]} {bf_r[itp, :]} {sign[itp, :]}")
                print(f"{rmx[itp,:]} {bf_rmx[itp, :]}")
                print(f"{rmy[itp,:]} {bf_rmy[itp, :]}\n")
