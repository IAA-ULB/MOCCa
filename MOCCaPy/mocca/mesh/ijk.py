import numpy as np
# TODO Maybe class IJK isn't needed after all. Remove it?

class IJK:
    """IJK is an index object to iterate over the grid points of a rectangular mesh, like the LagrangeMesh. It
    provides a `ig` attribute running from 1 to n_gridpoints, and attributes `i`, `j` and `k`, referring to the
    spatial indices of the grid points. The `ig` attribute corresponds to flattened arrays, and the `i`, `j` and `k`
    attributes correspond to the unflattened arrays.
    """
    def __init__(self, shape):
        """IJK constructor.
        Args:
            shape (tuple): shape of the mesh, can be 1D, 2D, or 3D.
        """
        self.dim = len(shape)
        # self.shape is a 3D tuple, even if self.dim<3, self.shape[i] == 1 for i>self.dim-1
        self.shape = list(shape)
        self.shape.extend([1,1])
        self.shape = self.shape[0:3]
        self.shape = tuple(self.shape)
        self.ig_end = int(np.prod(np.array(shape)))
        self.reset()


    def inc(self):
        """Increment the index to point to the next grid point in memory.

        Returns:
            The spatial index of the grid point pointed to after being incremented. Three values (`i`, `j` and `k`) are
            returned even if the mesh is 1D or 2D. In the latter case, only `i`, resp. `i` and `j` are relevant.
        """
        self.ig += 1
        self.index[0] += 1
        if self.index[0] == self.shape[0]:
            self.index[0] = 0
            if self.dim > 1:
                self.index[1] += 1
                if self.index[1] == self.shape[1]:
                    self.index[1] = 0
                    if self.dim > 2:
                        self.index[2] += 1
        return self.index

    def done(self) -> bool:
        """Returns `True` if the last `inc` call incremented past the end of the grid points, `False` otherwise."""
        return self.ig >= self.ig_end

    def ijk2ig(self, ijk):
        """ijk->ig"""
        return ijk[0] + ijk[1] * self.shape[0] + ijk[2] * (self.shape[0] * self.shape[1])

    def ig2ijk(self, ig):
        """ig->ijk
        returns:
            [i,j,k] corresponding to ig. Always returns 3 indices, even if `self.dim < 3`. Indices for non-existing
            dimensions are zero. This behavior is consistent with `self.index`.
        """
        ijk =[0,0,0]
        if self.dim > 2:
            ijk[2] = ig // (self.shape[0]*self.shape[1])
            ig %= (self.shape[0]*self.shape[1])
        if self.dim > 1:
            ijk[1] = ig // self.shape[0]
            ig %= self.shape[0]
        ijk[0] = ig
        return ijk


    def reset(self, ig: int = None, ijk=None):
        """Resets the index to point to the first grid point in the mesh."""
        # TODO implement ig and ijk

        if ig is None and ijk is None:
            self.index = [0,0,0]
            self.ig = 0
        elif ig is not None:
            if self.dim == 3:
                self.index[2] = ig//self.shape[0]*self.shape[1]
        elif ijk is not None:
            if self.dim == 3:
                self.ig = ijk[0] + ijk[1] * self.shape[0] + ijk[2] * (self.shape[0] * self.shape[1])
