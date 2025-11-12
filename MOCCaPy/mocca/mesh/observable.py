import numpy as np


def is_composite(axes:str) -> bool:
    """By convention composite derivative are represented by str object starting with a capital."""
    return axes[0].isupper()


def sort_axes(axes:str) -> str:
    """A sort key for sorting axes lists according to these criteria.

      - single partial derivatives before composite: e.g. `xyz` < 'Hessian'
      - single partial derivatives in alphabetical order and increasing differentiation order: 'x' < 'y' < 'yy'
      - composite derivatives from low to high differentiation order: e.g. `Grad` < 'Hessian' < 'Laplacian' < 'Tensor3'
        < 'Tensor4'. (Happens to be alphabetical too).
    """
    return f"{len(axes)}{axes}" if not is_composite(axes) else \
           axes


class Observable:
    """Base class for observables."""
    def __init__(self, mesh, n_components=None, data=None, symmetry=None, name=''):
        """
        Args:
            mesh: a mesh object - Currently, only LagrangeMesh objects are supported.
            data (np.array): Observable values on the mesh points. Array shape is (n_gridpoints, n_components),
                if None, an empty array is created with shape (mesh.linear_size, n_components).
            n_components (int): number of components to use. used if data is None. If both data and n_components
                are provided, then n_components must equal data.shape[1]
            symmetry: symmetry behavior of the observable's components wrt the coordinate axes of the mesh.
                The argument is converted to a numpy array of shape (mesh.dim, n_components). Possible values
                are 0, 1, or -1 implying, resp. no symmetry, symmetric or skew-symmetric behavior of the selected
                component on the selected coordinate axis. The value is critical for computing derivatives on
                reduced axes.
            name: optional name of the observable.
        Raises:
            RuntimeWarning if symmetry is not set for reduced coordinate axes.
        """
        self.name = name

        self.mesh = mesh

        if data is None:
            assert(n_components is not None)
            self.data = np.empty((mesh.linear_size, n_components), dtype=float, order='F')
            self.n_components = n_components
        else:
            assert data.shape[0] == mesh.linear_size, f"{data.shape[0]=} <> {mesh.linear_size=}"
            if len(data.shape) == 1:
                data = np.reshape(data, (mesh.linear_size, 1), order='F')
            if n_components is not None:
                assert(n_components == data.shape[1])
            self.data = data
        self.n_components = self.data.size // mesh.linear_size

        if symmetry is None:
            self.symmetry = None
            if any(self.mesh.reduced):
                raise UserWarning(f"Observable {self.name}: {symmetry=} was specified.\n"
                                  f"\tThis will yield ValueErrors when taking derivatives or interpolating."
                                  )

        else:
            self.symmetry = np.zeros((mesh.dim, self.n_components), dtype=int, order='F')
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
               assert s in [1,0,-1]

            # Verify that symmetry is specified for reduced axes.
            for idim in range(self.mesh.dim):
                if mesh.reduced[idim]:
                    for iq in range(self.n_components):
                        if self.symmetry[idim,iq] == 0:
                            raise UserWarning(f"Observable {self.name}: No symmetry specified for component {iq}.\n"
                                              f"\tThis will yield ValueErrors when taking derivatives or interpolating."
                                             )
        # Grid-based access
        self.dataG = self.mesh.cast2grid(self.data)
        
        self.derivatives = {} # A dictionary where derivatives will be stored. Keys are `str` combining the characters
            # 'x', 'y', 'z', e.g. 'xyz' corresponds to d^3/dxdydz, Accummulated derivatives, as e.g. the 'Laplacian'
            # (d^2/dx^2 + d^2/dy^2 + d^2/dz^2) can be keys too.
        self.derivativesG = {} # A dictionary where derivatives will be stored. Keys are `str` combining the characters
        self._derivative_is_uptodate = {}
        

    def __repr__(self):
        return f"Observable {self.name} {self.data.shape}"

    @property
    def shape(self):
        return self.data.shape

    @property
    def n_gridpoints(self):
        return self.mesh.linear_size

    def __getitem__(self, index:tuple) -> np.ndarray:
        """Access data and derivatives in linear or grid based way.

        >> Q['L', l, iq] # linear access by linear index l and component index iq (slices work too!)
        >> Q['G', i, j, k, iq] # Grid based access by grid index (i,j,k) and component index iq (slices work too!)
        >> Q['xy', 'L', l, iq] # linear access of d2Q/dxdy by linear index l and component index iq (slices work too!)
        >> Q['Laplacian', 'G', i, j, k, iq] # Grid based access of Laplacian derivative by grid index (i,j,k) and component index iq (slices work too!)

        Args:
            index (tuple): ([axes:str,] 'L'|'G', l | i[,i[,k]], iq), respectively
                - an optional axes str indicating a simple derivative, e.g. 'xy', 'z', 'Laplacian', ... Composite derivatives
                  like 'Grad', 'Hessian', 'Tensor3' and 'Tensor4' are not supported. The 'Laplacian' is an exception
                  because it is a scalar.
                - a str 'L' or 'G' requesting linear, resp. grid-based access.
                - 'L' is followed by a single linear index l (int), 'G' is followed by a grid index i[,j[,k]] (ints).
                  The number of ints must equal `self.mesh.dim`. Each int can be a slice too
                - iq (int) component index or slice.

        Remark:
            Looping through data/derivatives using this mechanism will be slow.
        """
        if isinstance(index[1], str):
            # Accessing derivatives
            axes = index[0]
            if is_composite(axes) and axes != 'Laplacian':
                raise ValueError(f"Accessing composite derivatives, like '{axes}', of Observable ({name}) is forbidden"
                                 f" (except for 'Laplacian'.")
            assert index[1] in 'LG', f"Access identifier must be 'L' or 'G', got {index[1]}."
            return self.derivatives [axes][index[2:]] if index[0] == 'L' else \
                   self.derivativesG[axes][index[2:]]

        else:
            # Accessing data
            assert index[0] in 'LG', f"Access identifier must be 'L' or 'G', got {index[1]}."
            return self.data [axes][index[1:]] if index[0] == 'L' else \
                   self.dataG[axes][index[1:]]
        
    def sign(self, iq:int):
        """Return the symmetry signs of component iq for the different coordinate axes.
        On reduced axes 1 implies symmetric and -1 skew-symmetric behavior.
        On non-reduced axes a 1 is returned by default.
        """
        result = np.ones((self.mesh.dim,), dtype=float, order='F')
        for idim in range(self.mesh.dim):
            if self.mesh.reduced[idim]:
                if self.symmetry is None:
                    raise ValueError(f"Observable {self.name}: No symmetry behavior specified for all components on reduced axis {'xyz'[idim]}.\n"
                                     f"\tInterpolation not possible.")
                if self.symmetry[idim, iq] == 0:
                    raise ValueError(f"Observable {self.name}: No symmetry behavior specified for component {iq} on reduced axis {'xyz'[idim]}.\n"
                                     f"\tInterpolation not possible.")
                result[idim] = self.symmetry[idim,iq]
            else:
                if self.symmetry is None:
                    result[idim] = 0
                else:
                    result[idim] = self.symmetry[idim,iq]
        return result

    # Differentiation
    #---------------------------------------------------------------------------
    def _get_tensor_of_derivatives(self, axes, shape):
        """Return a numpy array for a tensor derivative (Grad, Hessian, Tensor3, Tensor3),
        the elemenst of which are numpy arrays themselves containing the corresponding partial
        derivatives on the mesh points.
        """
        if not axes in self.derivativesG:
            self.derivativesG[axes] = np.empty(shape, dtype=np.ndarray)
        return self.derivativesG[axes]

    def grad(self, access='L'):
        """Compose Grad = [d/dx d/dy d/dz] from current derivatives. As the individual
        components the underlying data are automatically updated.

        Args:
            access: 'L' for linear access, 'G' grid based access.
        Returns:
            a 1-tensor represented by a numpy array of numpy arrays.
        """
        derivatives = self.derivatives if access == 'L' else \
                      self.derivativesG

        result = self._get_tensor_of_derivatives(axes, (self.mesh.dim,))
        result[0] = derivatives['x']
        result[1] = derivatives['y']
        if self.mesh.dim == 3:
            result[2] = self.derivatives['z']

        return result

    def hessian(self, access='L'):
        """Compose the Hessian from current derivatives (=all 2nd order derivatives).
        As the individual components the underlying data are automatically updated.

        Args:
            access: 'L' for linear access, 'G' grid based access.
        Returns:
            a 2-tensor represented by a numpy array of numpy arrays.
        """
        derivatives = self.derivatives if access == 'L' else \
                      self.derivativesG

        result = self._get_tensor_of_derivatives(axes, 2*(self.mesh.dim,))
        for i in range(self.mesh.dim):
            for j in range(self.mesh.dim):
                    axes = ('xyz'[i] +
                            'xyz'[j])
                    axes = ''.join(sorted(axes))
                    result[i,j] = derivatives[axes]

        return result

    def tensor3(self, access='L'):
        """Compose the Tensor3 from current derivatives (=all 3rd order derivatives).
        As the individual components the underlying data are automatically updated.

        Args:
            access: 'L' for linear access, 'G' grid based access.
        Returns:
            a 3-tensor represented by a numpy array of numpy arrays.
        """
        derivatives = self.derivatives if access == 'L' else \
                      self.derivativesG

        result = self._get_tensor_of_derivatives(axes, 3*(self.mesh.dim, ))
        for i in range(self.mesh.dim):
            for j in range(self.mesh.dim):
                for k in range(self.mesh.dim):
                    axes = ('xyz'[i] +
                            'xyz'[j] +
                            'xyz'[k])
                    axes = ''.join(sorted(axes))
                    result[i,j,k] = derivatives[axes]

        return result

    def tensor3(self, access='L'):
        """Compose the Tensor4 from current derivatives (=all 4th order derivatives).
        As the individual components the underlying data are automatically updated.

        Args:
            access: 'L' for linear access, 'G' grid based access.
        Returns:
            a 4-tensor represented by a numpy array of numpy arrays.
        """
        derivatives = self.derivatives if access == 'L' else \
                      self.derivativesG

        result = self._get_tensor_of_derivatives(axes, 4*(self.mesh.dim, ))
        for i in range(self.mesh.dim):
            for j in range(self.mesh.dim):
                for k in range(self.mesh.dim):
                    for l in range(self.mesh.dim):
                        axes = ('xyz'[i] +
                                'xyz'[j] +
                                'xyz'[k] +
                                'xyz'[l])
                        axes = ''.join(sorted(axes))
                        result[i,j,k,l] = self.derivativesG[axes]
        return result

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
            `self.derivatives[axes:str]` or `self.derivativesG[axes:str]`
        """
        # TODO: find a way of automatically calling invalidate_derivatives() after updating the Observable's data member?
        #       that would avoid specifying recompute.

        # design constraints:
        # -V store results internally (dict self.derivatives and self.derivativesG)
        # -V honour symmetry of derivatives 'xyx' == 'xxy' must be computed and stored only once
        # -V enable reusing previous computations: if we need 'xxy' and 'xx' is known, compute as D1y * 'xx', if 'y' is
        #    known, compute as D2x * 'y', otherwise compute (from scratch) as D2x * D1y * q
        if debug:
            print(f"Debug log>  {axes=}")

        if recompute:
            self.invalidate_derivatives()
            self._composite_done = set()

        if isinstance(axes, str):
            if self.derivative_is_uptodate(axes):
                # Uptodate derivative already available. This method can be used as a getter.
                # The grid based accessor is returned.
                if debug:
                    print(f"Debug log>  reusing {axes=}")
                return self.derivativesG[axes]

            if is_composite(axes):
                # axes is a multi-component derivative. Hence, self.mesh.dim >= 2 must hold.
                assert self.mesh.dim >= 2
                if not axes in ['Grad', 'Hessian', 'Laplacian', 'Tensor3', 'Tensor4']:
                    if debug:
                        print(f"Debug log>  {axes=} unknown composite derivative.")
                    raise ValueError(f"Unknown composite derivatve {axes}, allowed={['Grad', 'Hessian', 'Laplacian', 'Tensor3', 'Tensor4']}")

                # Wrap axes in a list to allow manipulations for reusing intermediate results
                if not axes in self._composite_done:
                    self._composite_done.add(axes)
                    self.differentiate(axes=[axes], recompute=False, debug=debug)

                if axes == 'Laplacian':
                    # Reuse or allocate memory
                    if axes in self.derivativesG:
                        result = self.derivativesG[axes]
                    else:
                        result = np.empty_like(self.dataG)
                        self.derivativesG[axes] = result
                        self.derivatives [axes] = self.mesh.cast2linear(result)

                    # Compute
                    if self.mesh.dim == 2:
                        result[:,:] = self.derivativesG['xx'] + \
                                      self.derivativesG['yy']
                    else:
                        result[:,:,:] = self.derivativesG['xx'] + \
                                        self.derivativesG['yy'] + \
                                        self.derivativesG['zz']
                    return result

            else:  # not composite
                # All simple derivatives. `axes` is composed as a sequence of 'x'|'y'|'z' characters.
                nx, ny, nz = axes.count('x'), axes.count('y'), axes.count('z')
                assert nx + ny + nz == len(axes), \
                       f"Extraneous characters in '{axes}', only 'x', 'y', 'and 'z' are allowed"
                # sort the `axes` str, as the order of differentiation is immaterial
                axes = nx*'x' + ny*'y' + nz*'z' # E.g. 'xyzx' -> 'xxyz', which is  evaluated as Dx2*Dy*Dz*Q

                if not axes in self.derivativesG:
                    # allocate memory
                    a = np.empty_like(self.dataG)
                    self.derivativesG[axes] = a
                    self.derivatives [axes] = self.mesh.cast2linear(a)

                out = self.derivativesG[axes]

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
                assert self.mesh.dim > 1, f"1D LagrangeMesh objects do not support composite derivatives: '{axes}'."
                assert self.mesh.highest_derivative_order >= 1, f"Required by axes='Grad'."
                axes = ['x', 'y', 'z'] + axes

            if 'Laplacian' in axes:
                assert self.mesh.dim > 1, f"1D LagrangeMesh objects do not support composite derivatives: '{axes}'."
                assert self.mesh.highest_derivative_order >= 2, f"Required by axes='Laplacian'."
                axes = ['xx', 'yy', 'zz'] + axes

            if 'Hessian' in axes:
                assert self.mesh.dim > 1, f"1D LagrangeMesh objects do not support composite derivatives: '{axes}'."
                assert self.mesh.highest_derivative_order >= 2, f"Required by axes='Hessian'."
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
                assert self.mesh.dim > 1, f"1D LagrangeMesh objects do not support composite derivatives: '{axes}'."
                assert self.mesh.highest_derivative_order >= 3, f"Required by axes='Tensor3'."
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
                assert self.mesh.dim > 1, f"1D LagrangeMesh objects do not support composite derivatives: '{axes}'."
                assert self.mesh.highest_derivative_order >= 4, f"Required by axes='Tensor4'."
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

            if self.mesh.dim == 2:
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

            return None  # returning a list is impractical

        else: # Axes should be eiter str or list
            raise ValueError(f"Axes of type {type(axes)} not supported ({axes=}).")

    def derivative_symmetry(self, axes):
        """Return the symmetry of the derivative of this observable wrt axes .
        """
        symmetry = self.symmetry.copy()

        # an even number of differentiations keeps the symmetry sign
        # an odd number of differentiations flips the symmetry sign
        # x-axis
        if self.mesh.reduced[0]:
            if axes.count('x') % 2:
                symmetry[0, :] *= -1
        # y-axis
        if self.mesh.dim > 1 and self.mesh.reduced[1]:
            if axes.count('y') % 2:
                symmetry[1, :] *= -1
        # z-axis
        if self.mesh.dim > 2 and self.mesh.reduced[2]:
            if axes.count('z') % 2:
                symmetry[2, :] *= -1

        return symmetry

    def invalidate_derivatives(self):
        """Mark all derivatives as outdated."""
        for axes in self._derivative_is_uptodate.keys():
            self._derivative_is_uptodate[axes] = False
        for axes in self.derivativesG.keys():
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
        """Integrate the observable over the simulation volume."""
        return self.mesh.integrate(self)
