import numpy as np


class Observable:
    """Base class for observables."""
    def __init__(self, mesh, data=None, n_components=None, symmetry=1):
        """
        Args:
            data (np.array): Observable values on the mesh points. this array is reshaped as (n_gridpoints,n_components),
        """
        self.mesh = mesh

        if data is None:
            assert(n_components is not None)
            self.data = np.empty((mesh.n_gridpoints(), n_components), dtype=float, order='F')
        else:
            assert data.shape[0] == mesh.n_gridpoints()
            if len(data.shape) == 1:
                data = np.reshape(data, (mesh.n_gridpoints(), 1))
            if n_components is not None:
                assert(n_components == data.shape[1])
            self.data = data

        self.n_components = self.data.shape[1]

        self.symmetry = np.empty((mesh.dim, self.n_components), dtype=int, order='F')
        if isinstance(symmetry, int):
            self.symmetry[:,:] = symmetry
        elif isinstance(symmetry, (tuple,list)):
            assert(len(symmetry) == self.n_components)
            for iq in range(self.n_components):
                self.symmetry[:, iq] = symmetry[iq]
        else:
            assert(symmetry.shape == (mesh.dim, self.n_components))
            self.symmetry = symmetry
        for s in self.symmetry.ravel():
            assert s in [1,-1]

    @property
    def n_gridpoints(self):
        return self.mesh.n_gridpoints()

    def __getitem__(self, i) -> np.ndarray:
        """Get i-th component of the observable."""
        return self.data[:,i]

    def integrate(self):
        return self.mesh.integrate(self)

    def sign(self,iq):
        return self.symmetry[:,iq].ravel()