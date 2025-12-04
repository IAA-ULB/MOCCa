"""
Unfinished code for solving the Poisson equation. This was a first start to the problem.
Instead we chose to rely on shenfun.
"""
import numpy as np

from .lagrange import create_mesh, Mesh


class PoissonMesh(Mesh):
    def __init__(self, M:int|tuple, d:int|float|tuple, reduced:tuple|bool=True, dim:int=0):
        """

        """
        super().init(M=M, d=d, reduced=reduced, dim=dim)
        self.initialize_gridpoints()

    def initialize_gridpoints(self, stencil=(-1,1)):
        ghost_zone_width = (stencil[1] - stencil[0] +1) // 5

        # initialize grid points:
        start = self.dim*[.0]
        n_reduced = self.dim*[0]
        g1D = self.dim*[0]

        for i, (n_i, reduced_i, d_i) in enumerate(zip(self.M, self.reduced, self.d)):
            if reduced_i:
                n_reduced[i] = n_i//2 + 2*ghost_zone_width
                start[i] = 0.5 - ghost_zone_width
            else:
                n_reduced[i] = n_i + 2*ghost_zone_width
                start[i] = -(n_i - 1)*0.5 - ghost_zone_width

            g1D[i] = np.linspace(start[i], start[i] + n_reduced[i] - 1, n_reduced[i])
            g1D[i] *= d_i

        self.gx = g1D[0]
        if self.dim > 1:
            self.gy = g1D[1]
        if self.dim > 2:
            self.gz = g1D[2]

        if self.dim == 3:
            self.grid = create_mesh(self.gx, self.gy, self.gz)

            self.gridx = np.empty(n_reduced, order='F')
            self.gridy = np.empty(n_reduced, order='F')
            self.gridz = np.empty(n_reduced, order='F')
            for k in range(n_reduced[2]):
                for j in range(n_reduced[1]):
                    self.gridx[:,j,k] = g1D[0]
            for k in range(n_reduced[2]):
                for i in range(n_reduced[0]):
                    self.gridy[i,:,k] = g1D[1]
            for j in range(n_reduced[1]):
                for i in range(n_reduced[0]):
                    self.gridz[i,j,:] = g1D[2]

        elif self.dim == 2:
            self.grid = create_mesh(self.gx, self.gy)

            self.gridx = np.empty(n_reduced, order='F')
            self.gridy = np.empty(n_reduced, order='F')
            for j in range(n_reduced[1]):
                self.gridx[:, j] = g1D[0]
            for i in range(n_reduced[0]):
                self.gridy[i, :] = g1D[1]

        else: # self.dim == 1
            self.grid = create_mesh(self.gx)
            self.gridx = g1D[0]
            # print(f"{self.gridx=}")

        self.mesh_shape   = self.gridx.shape
        self.ghost_zone_width = ghost_zone_width

        self.V = np.empty_like(self.gridx, order='F')
        self.f = np.empty_like(self.gridx, order='F')

    def apply_function(self, function):
        self.potential[:,:,:] = function(self.gridx, self.gridy, self.gridz)

    def __repr__(self):
        reduced = ''.join(['T' if r else 'F' for r in self.reduced])
        shape = 'x'.join([str(m) for m in self.mesh_shape])
        return (f"<LagrangeMesh(Mesh)[{shape}={self.linear_size}, {reduced=}]>")
