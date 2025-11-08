import numpy as np


class Observable:
    """Base class for observables."""
    def __init__(self, mesh, data=None, n_components=None, symmetry=1, name=''):
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
            name: optional name of the observable.
        """
        self.name = name

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

        self.derivatives = {} # A dictionary where derivatives will be stored. Keys are `str` combining the characters
            # 'x', 'y', 'z', e.g. 'xyz' corresponds to d^3/dxdydz, Accummulated derivatives, as e.g. the 'Laplacian'
            # (d^2/dx^2 + d^2/dy^2 + d^2/dz^2) can be keys too.
        self._derivative_is_uptodate = {}

    @property
    def dim(self):
        return self.mesh.dim

    @property
    def shape(self):
        return self.data.shape

    @property
    def n_components(self):
        return self.data.shape[1]

    @property
    def n_gridpoints(self):
        return self.mesh.n_gridpoints()

    def __getitem__(self, i) -> np.ndarray:
        """Get i-th component of the observable."""
        return self.data[:,i]

    def sign(self,iq):
        """Return the symmetry signs of component iq for the different coordinate axes. A 1 implies symmetric and -1 is
        skew-symmetric."""
        return self.symmetry[:,iq].ravel()

    # Differentiation
    #---------------------------------------------------------------------------
    def differentiate(self, axes:str|list[str], recompute:bool=True):
        """Compute some spatial derivative(s) of this observable. all components are differentiated

         Args:
             axes: See doc-string of  LagrangeMesh.differentiate
             recompute: Recompute all the derivatives specified in axes. If False, previously computed derivatives are
                recycled, and the method behaves as a getter. Derivatives that hadn't been computed yet, are computed.
                After modifying `self.data`, `self.differentiate` should, obviously, be called with 'recompute=True'.
                Ideally, you request all needed derivatives at once in a single call. You may want to split the
                >>> Q = Observable(...)
                >>> Q.data = ... # modify the observable's data, derivatives are now outdated
                Recompute the derivative for the new `Q.data`
                >>> Q.differentiate(axes=['x','y','z']) # recompute=True by default
                >>> Q.differentiate(axes=['Laplacian'], recompute=False) # The Laplacian can reuse 1st orde derivatives
        Returns:
            the derivative if a single derivative or derivative tensor (e.g. 'Grad') was requested. If a list of
            derivatives was requested, `None` is returned and the user must access the derivatives via the dict
            `Q.derivatives`.
        """
        # TODO: find a way of automatically calling invalidate_derivatives() after updating the Observable's data member?
        #       that would avoid specifying recompute.

        # design constraints:
        # -V store results internally (dict self.derivatives)
        # -V honour symmetry of derivatives 'xyx' == 'xxy' must be computed and stored only once
        # -V enable reusing previous computations: if we need 'xxy' and 'xx' is known, compute as D1y * 'xx', if 'y' is
        #    known, compute as D2x * 'y', otherwise compute (from scratch) as D2x * D1y * q

        if recompute:
            self.invalidate_derivatives()

        result = self.mesh.differentiate(self, axes=axes)
        return result

    def invalidate_derivatives(self):
        for axes in self._derivative_is_uptodate.keys():
            self._derivative_is_uptodate[axes] = False
        for axes in self.derivatives.keys():
            self._derivative_is_uptodate[axes] = False

    def derivative_set_uptodate(self, axes, value=True):
        """Indicate that the derivative wrt axes was computed after modifying the observable's data,
        and, thus, that the derivative uptodate relative to the observable's data."""
        self._derivative_is_uptodate[axes] = value

    def derivative_is_uptodate(self, axes):
        """Is the derivative wrt axes upto date?"""
        self._derivative_is_uptodate.get(axes, False)

    # Forwarding methods: Since the observable stores (a reference to) the mesh on which it is defined, we can call
    # LagrangeMesh methods directly on the Observable.
    def integrate(self):
        return self.mesh.integrate(self)
