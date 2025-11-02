import numpy as np


class Observable:
    """Base class for observables."""
    def __init__(self, mesh, data=None, n_components=None, symmetry=1):
        """
        Args:
            mesh: a LagrangeMesh object
            data (np.array): Observable values on the mesh points. Array shape is (n_gridpoints,n_components),
                if None, an empty array is created with shape (mesh.n_gridpoints(), n_components).
            n_components (int): number of components to use. used if data is None. If both data and n_components
                are provided, then n_components must equal data.shape[1]
            symmetry: symmetry of the observable. If an integer is provided, then all components have this symmetry
                on all coordinate axes. If a list|tuple of integers is provided, then len(symmetry) == n_components
                must hold (if n_components is None, it is inferred from data), and symmtry[iq] prescribes the symmetry
                of component iq on all coordinate axes. If a np.ndarray is provided, then symmetry.shape must eaual
                (mesh.dim, n_components) and symmetry[i,iq] prescribes the symmetry of component iq on the i-th
                coordinate axis. All values must be either 1 (symmetric) or -1 (skew-symmetric).
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
        elif isinstance(symmetry, (tuple,list)) and not isinstance(symmetry[0], (tuple,list)):
            assert(len(symmetry) == self.n_components)
            for iq in range(self.n_components):
                self.symmetry[:, iq] = symmetry[iq]
        else:
            if not isinstance(symmetry, np.ndarray):
                symmetry = np.array(symmetry)
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