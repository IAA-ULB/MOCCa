import numpy as np


class Observable:
    """Base class for observables."""
    def __init__(self, data, symmetry=1):
        """
        Args:
            data (np.array): Observable values on the mesh points. this array is reshaped as (n_gridpoints,n_components),
        """
        self.data = data
        self.n_gridpoints = data.shape[0]

        if len(data.shape) == 1:
            self.n_components = 1
            self.data = data.reshape((self.n_gridpoints,1))
            self.shape = (1,) # the shape or the observable itself, e.g. a rank-2 tensor in 3D has shape (3,3)
            assert isinstance(symmetry, int)
            self.symmetry = np.array([symmetry], order='F')

        else:
            self.shape = data.shape[1:]
            self.n_components = int(np.prod(self.shape))
            self.data = data.reshape((self.n_gridpoints, self.n_components), order='F')
            if isinstance(symmetry,(int,float)):
                self.symmetry = symmetry * np.ones((self.n_components,), dtype=float, order='F')
            else:
                assert symmetry.shape == self.shape
                self.symmetry = np.array(symmetry, dtype=float, order='F')
                self.symmetry = self.symmetry.reshape((self.n_components,), order='F')

        for symm in self.symmetry:
            assert symm == 1. or symm == -1.

    def __getitem__(self, i) -> np.ndarray:
        """Get i-th component of the observable."""
        return self.data[:,i]

    def get_symmetry(self, i) -> int:
        return self.symmetry[i]