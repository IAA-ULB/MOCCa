import numpy as np


def is_composite(axes:str) -> bool:
    """By convention composite derivative are represented by str object starting with a capital."""
    return axes[0].isupper()


def sort_axes(axes:str) -> str:
    return f"{len(axes)}{axes}" if not is_composite(axes) else \
           axes


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
    def _get_result_array(self, axes, shape):
        if not axes in self.derivatives:
            self.derivatives[axes] = np.empty(shape, dtype=np.ndarray)
        return self.derivatives[axes]

    def differentiate(self, axes:str|list[str], recompute:bool=True, debug=False):
        """Compute some spatial derivative(s) of this observable. all components are differentiated

         Args:
             axes:
                str: a str of 'x'|'y'|'z' characters, indicating the differentiation order and axes.
                    E.g. 'xxz' requests differentiation d^3/dxdxdz.
                    Alternatetively, the name of a composite derivative:
                    - 'Grad': all first order derivatives
                    - 'Hessian': all second order derivatives
                    - 'Laplacian' d^2/dx^2 + d^2/dy^2 + d^2/dz^2)
                    - 'Tensor3': all third order derivatives
                    - 'Tensor4': all fourth order derivatives
                    These return an numpy array with the corresponding tensor of partial derivatives.
                    'Laplacian', is an exception because it is a scalar differentiation operator, and
                    therefor the result of d^2/dxdx + d^2/dy^2 + d^2/dz^2)Q is returned.
                list: A list of the above strings is also accepted, requesting several derivatives
                    at once. The list can internally be manipulated to allow storing intermediate
                    derivatives that can be reused to speed up computation. E.g. when requesting 'xz'
                    and 'yz', it is advantageous to compute d/dz first and apply d/dx and d/dy to it,
                    thereby saving one matrix application.
             recompute: If true (=default) all the derivatives needed by the axes request are recomputed.
                If False, derivatives computed in previous calls to `Observable.differentiate()` can be
                reused as a starting point for the requested derivatives.
                After modifying the Observable, `differentiate` should, obviously, be called with
                'recompute=True' (=default). It may be practical to request all needed derivatives
                in a single differentate() call. Somtimes it may be more practical to split the
                request over several calls where the first call uses `recompute=True` and succeeding
                calls use `recompute=False`. The succeeding call could e.g. correspond to increasingly
                higher order derivatives
                >>> Q = Observable(...)
                >>> Q.data = ... # modify the observable's data, derivatives are now outdated
                >>> Q.differentiate(axes=['Grad'])
                >>> Q.differentiate(axes=['Hessian'], recompute=False)
                The first call uses `recompute=True` and marks all previously computed derivatives as
                not uptodate, then computes all 1st order derivatives. The second call proceeds to comppute
                all 2nd order derivatives and reuse 'y' and 'z' computed in the first call in the computation
                of the cross derivatives 'xy', 'xx' and 'yz'.
        Returns:
            If axes is a str, returns the result in the form of a numpy array of floats. For composite str
            representing tensor differentiaton operators ('Grad, 'Hessian', ...) the returned result is a
            numpy array (the tensor) of numpy arrays of floats.
            In the case of a list, None is returned and the user must access the individual derivatives as
            `self.derivatives[axes:str]`
        """
        # TODO: find a way of automatically calling invalidate_derivatives() after updating the Observable's data member?
        #       that would avoid specifying recompute.

        # design constraints:
        # -V store results internally (dict self.derivatives)
        # -V honour symmetry of derivatives 'xyx' == 'xxy' must be computed and stored only once
        # -V enable reusing previous computations: if we need 'xxy' and 'xx' is known, compute as D1y * 'xx', if 'y' is
        #    known, compute as D2x * 'y', otherwise compute (from scratch) as D2x * D1y * q
        if debug:
            print(f"Debug log>  {axes=}")

        if recompute:
            self.invalidate_derivatives()
            # Cast self.data in mesh shape which is required for the einsum calls.
            # The data structures for the derivatives are then automatically in the right shape too.
            self.mesh.cast_observable_in_mesh_shape(self)
            self._composite_done = set()

        if isinstance(axes, str):
            if self.derivative_is_uptodate(axes):
                # Uptodate derivative already available. This method can be used as a getter.
                if debug:
                    print(f"Debug log>  reusing {axes=}")
                return self.derivatives[axes]

            if is_composite(axes):
                # axes is a multi-component derivative. Hence, self.dim >= 2 must hold.
                assert self.dim >= 2
                if not axes in ['Grad', 'Hessian', 'Laplacian', 'Tensor3', 'Tensor4']:
                    if debug:
                        print(f"Debug log>  {axes=} unknown composite derivative.")
                    raise ValueError(f"Unknown composite derivatve {axes}, allowed={['Grad', 'Hessian', 'Laplacian', 'Tensor3', 'Tensor4']}")

                # Wrap axes in a list to allow manipulations for reusing intermediate results
                if not axes in self._composite_done:
                    self._composite_done.add(axes)
                    self.differentiate(axes=[axes], recompute=False, debug=debug)

                # Use the above computed partial derivatives to compute the result
                if axes == 'Grad':
                    result = self._get_result_array(axes, (self.dim,))

                    result[0] = self.derivatives['x']
                    result[1] = self.derivatives['y']
                    if self.dim == 3:
                        result[2] = self.derivatives['z']

                elif axes == 'Hessian':
                    result = self._get_result_array(axes, (self.dim, self.dim))

                    result[0, 0] = self.derivatives['xx']
                    result[0, 1] = self.derivatives['xy']
                    result[1, 0] = self.derivatives['xy']
                    result[1, 1] = self.derivatives['yy']
                    if self.dim == 3:
                        result[0, 2] = self.derivatives['xz']
                        result[2, 0] = self.derivatives['xz']
                        result[1, 2] = self.derivatives['yz']
                        result[2, 1] = self.derivatives['yz']
                        result[2, 2] = self.derivatives['zz']

                elif axes == 'Laplacian':
                    if axes in self.derivatives:
                        result = self.derivatives[axes]
                    else:
                        result = np.zeros_like(self.data)
                        self.derivatives[axes] = result

                    result += self.derivatives['xx']
                    result += self.derivatives['yy']
                    if self.dim >= 2:
                        result += self.derivatives['zz']

                elif axes == 'Tensor3':
                    raise NotImplementedError(f"{axes=} is not implemented.")
                    result = self._get_result_array(axes, 3*(self.dim, ))
                    # TODO: implement

                elif axes == 'Tensor4':
                    raise NotImplementedError(f"{axes=} is not implemented.")
                    result = self._get_result_array(axes, 4*(self.dim, ))
                    # TODO: implement

                else:
                    raise NotImplementedError(f"{axes=} is not implemented.")

                return result

            else:  # not composite
                # All simple derivatives. `axes` is composed as a sequence of 'x'|'y'|'z' characters.
                nx, ny, nz = axes.count('x'), axes.count('y'), axes.count('z')
                assert nx + ny + nz == len(axes), \
                       f"Extraneous characters in '{axes}', only 'x', 'y', 'and 'z' are allowed"
                # sort the `axes` str, as the order of differentiation is immaterial
                axes = nx*'x' + ny*'y' + nz*'z' # E.g. 'xyzx' -> 'xxyz', which is  evaluated as Dx2*Dy*Dz*Q

                if not axes in self.derivatives:
                    # allocate memory
                    self.derivatives[axes] = np.empty_like(self.data)
                out = self.derivatives[axes]

                # This is where the responsibility of Observable ends and the responsibility of
                # the mesh object (typically, LagrangeMesh) begins.
                if debug:
                    print(f"Debug log>  mesh.differentiate(Q=self, axes='{axes}', out=out)")
                self.mesh.differentiate(Q=self, axes=axes, out=out)

                self.derivative_set_uptodate(axes)

                return out

        elif isinstance(axes, list):
            # Handle lists of derivatives
            # Add components to allow reuse of derivatives:
            if 'Grad' in axes:
                assert self.dim > 1, f"1D LagrangeMesh objects do not support composite derivatives: '{axes}'."
                axes = ['x', 'y', 'z'] + axes

            if 'Laplacian' in axes:
                assert self.dim > 1, f"1D LagrangeMesh objects do not support composite derivatives: '{axes}'."
                axes = ['xx', 'yy', 'zz'] + axes

            if 'Hessian' in axes:
                assert self.dim > 1, f"1D LagrangeMesh objects do not support composite derivatives: '{axes}'."
                axes = ['y', 'z',
                        'xx', 'xy', 'xz',
                        'yy', 'yz',
                        'zz'
                        ] + axes
                # As D2x*Q is computed equally efficient as D1x*dQdx and the the same holds for y and z,
                # there is no point in reusing 'x' for 'xx'. However, reusing 'y' for 'xy' and 'z' for 'xz'
                # and 'yz' is beneficial.
                # The genaral rule is that everything starting with 'x' can be dropped, except the one needed
                # by the composite, c.q 'Hessian'.

            if 'Tensor3' in axes:
                assert self.dim > 1, f"1D LagrangeMesh objects do not support composite derivatives: '{axes}'."
                axes = ['y', 'z',
                        'yy', 'yz', 'zz',
                        'xxx', 'xxy', 'xxz',
                        'xyy', 'xyz',
                        'xzz',
                        'yyy', 'yyz',
                        'yzz',
                        'zzz',
                        ] + axes
                # 'x', 'xx' and and 'xy' are dropped for the same reason as above.

            if 'Tensor4' in axes:
                assert self.dim > 1, f"1D LagrangeMesh objects do not support composite derivatives: '{axes}'."
                axes = ['y', 'z',
                        'yy', 'xz', 'yz', 'zz',
                        'yyy', 'yyz', 'yzz', 'zzz',
                        'xxxx', 'xxxy', 'xxxz',
                        'xxyy', 'xxyz',
                        'xxzz',
                        'xyyy', 'xyyz',
                        'xyzz',
                        'xzzz',
                        'yyyy', 'yyyz',
                        'yyzz',
                        'yzzz',
                        'zzzz',
                        ] + axes
                # all starting with 'x' and order < 4 are dropped.

            if self.dim == 2:
                # Remove entries containing 'z':
                axes = [ax for ax in axes if 'z' not in ax]
            # Remove duplicate entries:
            axes = list(set(axes))
            # Sort the list in-place (the sorting key ensures that low order derivatives are
            # computed first, enabling optimal reuse):
            axes.sort(key=sort_axes)
            if debug:
                print(f"Debug log>  {axes=}")
            # Process the list:
            for ax in axes:
                self.differentiate(axes=ax, recompute=False, debug=debug)

            return None  # returning a list would make no sense, the user must access the requested derivatives via
                         # `self.derivatives`

        else: # Axes should be eiter str or list
            raise ValueError(f"Axes of type {type(axes)} not supported ({axes=}).")


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
        return self._derivative_is_uptodate.get(axes, False)

    # Forwarding methods: Since the observable stores (a reference to) the mesh on which it is defined, we can call
    # LagrangeMesh methods directly on the Observable.
    def integrate(self):
        return self.mesh.integrate(self)
