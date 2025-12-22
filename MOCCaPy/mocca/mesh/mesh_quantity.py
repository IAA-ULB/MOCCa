import numpy as np
from numpy import ndarray


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


class MeshQuantity:
    """Base class for observables."""
    def __init__(self, mesh, n_components=None, data=None, symmetry=None, name=''):
        """
        Args:
            mesh: a mesh object - Currently, only LagrangeMesh objects are supported.
            data (np.array): MeshQuantity values on the mesh points. Array shape is (n_gridpoints, n_components),
                if None, an empty array is created with shape (mesh.linear_size, n_components).
            n_components (int): number of components to use. used if data is None. If both data and n_components
                are provided, then n_components must equal data.shape[1]
            symmetry: symmetry behavior of the MeshQuantity's components wrt the coordinate axes of the mesh.
                The argument is converted to a numpy array of shape (mesh.dim, n_components). Possible values
                are 0, 1, or -1 implying, resp. no symmetry, symmetric or skew-symmetric behavior of the selected
                component on the selected coordinate axis. The value is critical for computing derivatives on
                reduced axes.
            name: optional name of the MeshQuantity.
        Raises:
            ValueError: if n_components is None and data is None.
            RuntimeWarning if symmetry is not set for reduced coordinate axes.
        """
        self.name = name

        self.mesh = mesh

        if data is None:
            if n_components is None:
                raise ValueError(
                    f"Parameter n_components must be specified if data is None (got {n_components=}).")
            self.data = np.empty((mesh.linear_size, n_components), dtype=np.float64, order='F')
            self.n_components = n_components
        else:
            if not (data.shape[0] == mesh.linear_size):
                raise ValueError(f"{data.shape[0]=} <> {mesh.linear_size=}")

            if len(data.shape) == 1:
                data = np.reshape(data, (mesh.linear_size, 1), order='F')
            elif len(data.shape) > 2:
                # The component index is not 1D, e.g. for HFPsi, it is (4,n_total_wf)
                # Flatten it out
                data = np.reshape(data, (mesh.linear_size, np.prod(data.shape[1:])), order='F')
            if n_components is not None:
                if not (n_components == data.shape[1]):
                    raise ValueError(
                        f"If specified together with data, parameter n_components must "
                        f"be equal to data.shape[1] (got {n_components=})."
                    )
            self.data = data
        self.n_components = self.data.size // mesh.linear_size

        if symmetry is None:
            self.symmetry = None
            if any(self.mesh.reduced):
                raise UserWarning(f"MeshQuantity {self.name}: {symmetry=} was specified.\n"
                                  f"\tThis will yield ValueErrors when taking derivatives, interpolating, solving "
                                  f"generalized Poisson equations."
                                  )

        else:
            self.set_symmetry(symmetry)


            # Verify that symmetry is specified for reduced axes.
            for idim in range(self.mesh.dim):
                if mesh.reduced[idim]:
                    for iq in range(self.n_components):
                        if self.symmetry[iq,idim] == 0:
                            raise UserWarning(f"MeshQuantity {self.name}: No symmetry specified for component {iq}.\n"
                                              f"\tThis will yield ValueErrors when taking derivatives or interpolating."
                                             )
        # Grid-based access
        self.dataG = self.mesh.cast2grid(self.data)
        
        # Dictionaries where derivatives will be stored. Keys are `str` combining the characters
        # 'x', 'y', 'z', e.g. 'xyz' corresponds to d^3/dxdydz, Accummulated derivatives, as e.g. the 'Laplacian'
        # (d^2/dx^2 + d^2/dy^2 + d^2/dz^2) can be keys too.
        # self.derivatives returns arrays with linear access
        # self.derivativesG returns arrays with grid-based access
        self.derivatives = {}
        self.derivativesG = {} # A dictionary where derivatives will be stored. Keys are `str` combining the characters
        self._derivative_is_uptodate = {}

    def set_symmetry(self, symmetry):
        """"""
        self.symmetry = np.empty((self.n_components, self.mesh.dim), dtype=np.int32, order='F')
        if isinstance(symmetry, int):
            self.symmetry[:, :] = symmetry

        elif isinstance(symmetry, (tuple, list)):
            symmetry = np.array(symmetry)
            # symmetry.shape must be
            # (1) either (self.mesh.dim,)               -> all components identical
            # (2) or (self.n_components,1)              -> all dimensions identical, per component
            # (3) or (self.n_components, self.mesh.dim) -> all components and dimensions specified separately
            if symmetry.shape == (self.mesh.dim,):  # (1)
                for iq in range(self.n_components):
                    self.symmetry[iq, :] = symmetry[:]

            elif symmetry.shape == (self.n_components, 1):  # (2)
                for idim in range(self.mesh.dim):
                    self.symmetry[:, idim] = symmetry[:, 0]

            elif symmetry.shape == (self.n_components, self.mesh.dim):  # (3)
                # may raise "ValueError: could not broadcast input array ..."
                self.symmetry[:, :] = symmetry

            else:
                raise ValueError(f"Symmetry {symmetry} not understood. Either specify:\n"
                                 f"  - s                use int s for all components and all dimensions \n"
                                 f"  - [sx<,sy<,sz>>]   use this for all components\n"
                                 f"  - [[s],            use this for component 0 (all dimensions)\n"
                                 f"     [s],            use this for component 1 (all dimensions)\n"
                                 f"     ...]"
                                 f"  - [[sx<,sy<,sz>>], use this for component 0\n"
                                 f"     [sx<,sy<,sz>>], use this for component 1\n"
                                 f"     ...           ]\n"
                                 f"(Tuples may be used instead of lists).")

        elif isinstance(symmetry, np.ndarray):
            # may raise BroadcastError
            self.symmetry[:, :] = symmetry

        else:
            raise TypeError(
                f"Symmetry specifications must be of type int|list|tuple|list[list]|tuple[tuple]|np.ndarray, not '{type(symmetry)}'.'")

        assert len(self.symmetry.shape) == 2

        for s in self.symmetry.ravel():
            if not s in [1, 0, -1]:
                raise ValueError(
                    f"Symmetry values must be +1, -1 or 0  (got {s})."
                )

    def __repr__(self):
        return f"{self.name}:{self.data.shape}"

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
                raise ValueError(f"Accessing composite derivatives, like '{axes}', of MeshQuantity ({name}) is forbidden"
                                 f" (except for 'Laplacian'.")
            if not (index[1] in 'LG'):
                raise ValueError("Access identifier must be 'L' or 'G', got {index[1]}.")

            return self.derivatives [axes][index[2:]] if index[0] == 'L' else \
                   self.derivativesG[axes][index[2:]]

        else:
            # Accessing data
            if not (index[0] in 'LG'):
                raise ValueError(f"Access identifier must be 'L' or 'G', got {index[1]}.")

            return self.data [axes][index[1:]] if index[0] == 'L' else \
                   self.dataG[axes][index[1:]]
        
    def sign(self, iq:int):
        """Return the symmetry signs of component iq for the different coordinate axes.
        On reduced axes 1 implies symmetric and -1 skew-symmetric behavior.
        On non-reduced axes a 1 is returned by default.
        """
        result = np.ones((self.mesh.dim,), dtype=np.float64, order='F')
        for idim in range(self.mesh.dim):
            if self.mesh.reduced[idim]:
                if self.symmetry is None:
                    raise ValueError(f"MeshQuantity {self.name}: No symmetry behavior specified for all components on reduced axis {'xyz'[idim]}.\n"
                                     f"\tInterpolation not possible.")
                if self.symmetry[iq,idim] == 0:
                    raise ValueError(f"MeshQuantity {self.name}: No symmetry behavior specified for component {iq} on reduced axis {'xyz'[idim]}.\n"
                                     f"\tInterpolation not possible.")
                result[idim] = self.symmetry[iq,idim]
            else:
                if self.symmetry is None:
                    result[idim] = 0
                else:
                    result[idim] = self.symmetry[iq,idim]
        return result

    # Differentiation
    #---------------------------------------------------------------------------
    def grad(self, access='L'):
        """Compose Grad = [d/dx d/dy d/dz] from current derivatives. 
        As the individual components are references, the underlying data are automatically 
        updated.

        Args:
            access: 'L' for linear access, 'G' grid based access.
        Returns:
            a 1-tensor represented by a numpy array of numpy arrays.
        """
        derivatives = self.derivatives if access == 'L' else \
                      self.derivativesG

        result = np.empty((self.mesh.dim,), dtype=np.ndarray)
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

        result = np.empty(2*(self.mesh.dim,), dtype=np.ndarray)
        for i in range(self.mesh.dim):
            for j in range(self.mesh.dim):
                    axes = ('xyz'[i] +
                            'xyz'[j])
                    axes = ''.join(sorted(axes))
                    result[i,j] = derivatives[axes]

        return result

    def tensor3(self, access='L'):
        """Compose the Tensor3 from current derivatives (=all 3rd order derivatives).
        As the individual components are references, the underlying data are automatically 
        updated.

        Args:
            access: 'L' for linear access, 'G' grid based access.
        Returns:
            a 3-tensor represented by a numpy array of numpy arrays.
        """
        derivatives = self.derivatives if access == 'L' else \
                      self.derivativesG

        result = np.empty(3*(self.mesh.dim,), dtype=np.ndarray)
        for i in range(self.mesh.dim):
            for j in range(self.mesh.dim):
                for k in range(self.mesh.dim):
                    axes = ('xyz'[i] +
                            'xyz'[j] +
                            'xyz'[k])
                    axes = ''.join(sorted(axes))
                    result[i,j,k] = derivatives[axes]

        return result

    def tensor4(self, access='L'):
        """Compose the Tensor4 from current derivatives (=all 4th order derivatives).
        As the individual components are references, the underlying data are automatically 
        updated.

        Args:
            access: 'L' for linear access, 'G' grid based access.
        Returns:
            a 4-tensor represented by a numpy array of numpy arrays.
        """
        derivatives = self.derivatives if access == 'L' else \
                      self.derivativesG

        result = np.empty(4*(self.mesh.dim,), dtype=np.ndarray)
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

    def _differentiate1(self, axes:str, debug:bool):
        """Private method for computing a single derivative."""
        if is_composite(axes):
            if axes == 'Laplacian':
                # Reuse or allocate memory
                if axes in self.derivativesG:
                    result = self.derivativesG[axes]
                else:
                    result = np.empty_like(self.dataG)
                    self.derivativesG[axes] = result
                    self.derivatives[axes] = self.mesh.cast2linear(result)

                # Compute
                if self.mesh.dim == 2:
                    result[:, :] = self.derivativesG['xx'] + \
                                   self.derivativesG['yy']
                else:
                    result[:, :, :] = self.derivativesG['xx'] + \
                                      self.derivativesG['yy'] + \
                                      self.derivativesG['zz']

            # For composite derivatives other than 'Laplacian' there is nothing to do,
            # as their components have been computed already, and they must not be combined.

        else:
            # All simple derivatives. `axes` is composed as a sequence of 'x'|'y'|'z' characters.
            nx, ny, nz = axes.count('x'), axes.count('y'), axes.count('z')
            if not (nx + ny + nz == len(axes)): raise ValueError(
                f"Extraneous characters in '{axes}', only 'x', 'y', 'and 'z' are allowed"
            )
            # sort the `axes` str, as the order of differentiation is immaterial
            axes = nx * 'x' + ny * 'y' + nz * 'z'  # E.g. 'xyzx' -> 'xxyz', which is  evaluated as Dx2*Dy*Dz*Q

            if not axes in self.derivativesG:
                # allocate memory
                a = np.empty_like(self.dataG)
                self.derivativesG[axes] = a
                self.derivatives[axes] = self.mesh.cast2linear(a)

            out = self.derivativesG[axes]

            # This is where the responsibility of MeshQuantity ends and the responsibility of
            # the mesh object (typically, LagrangeMesh) begins.
            if debug:
                print(f"Debug log>  mesh.differentiate(Q=self, axes='{axes}', out=out)")

            self.mesh.differentiate(Q=self, axes=axes, out=out)

        self.derivative_set_uptodate(axes)


    def differentiate(self, axes:str|list[str], access='L', recompute:bool=True, debug=False):
        """Compute some spatial derivative(s) of this MeshQuantity's components.

        If composite derivatives are requested ('Grad', 'Hessian', 'Laplacian', ...) the `axes` list
        is first completed with all the needed simple derivatives.

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
                list: A list of the above strings is also accepted, requesting several derivatives
                    at once. The list can internally be manipulated to allow storing intermediate
                    derivatives that can be reused to speed up computation. E.g. when requesting 'xz'
                    and 'yz', it is advantageous to compute d/dz first and apply d/dx and d/dy to it,
                    thereby saving one matrix application.
            access: (Optional) access method for the returned array, if `isinstance(axes,str)`,
                ignored otherwise. 'L'=linear, 'G'=grid-based. This does not influence the
                computation. Default is linear access in agreement with standard MOCCa data
                structures.
             recompute: If `True` (=default) all the derivatives needed by the `axes` request are recomputed.
                If `False`, derivatives computed in previous calls to `MeshQuantity.differentiate()` can be
                reused as a starting point for the requested derivatives.
                After modifying the MeshQuantity, `differentiate` should, obviously, be called with
                'recompute=True' (=default). Sometimes it is more practical to split the derivatives
                request over several calls where the first call uses `recompute=True` and succeeding
                calls use `recompute=False`. The succeeding calls can e.g. request increasingly
                higher order derivatives.
                >>> Q = MeshQuantity(...)
                >>> Q.data = ... # modify the MeshQuantity's data, derivatives are now outdated
                >>> Q.differentiate(axes=['Grad'])
                >>> Q.differentiate(axes=['Hessian'], recompute=False)
                The first call uses `recompute=True` and requires all previously computed derivatives to be
                recomputed, then proceeds computing all 1st order derivatives. The second call proceeds to comppute
                all 2nd order derivatives and reuses 'y' and 'z' computed in the first call in the computation
                of the cross derivatives 'xy', 'xz' and 'yz'.
        Returns:
            A numpy array with the requested derivative if `axes` is a `str` referring to a scalar derivative
            (such as 'x', `yz`, or 'Laplacian`). Otherwise, if `axes` is a `list` or refers to a tensor-like
            derivative (such as `Grad`, 'Hessian', ...), `None` is returned.
            The shape of the returned array depends on the `access` parameter.

        Raises:
            ValueError
              - if one requests higher order derivatives than allowed by
               `LagrangeMesh.highest_derivative_order`,
              - in case of an invalid `axes` parameter,
              - in case of a composite derivative on a 1D mesh.

        All computed derivatives are stored internally and can be accessed by the Observables as
        `MeshQuantity.derivatives[axes:str]` (linear access) or `MeshQuantity.derivativesG[axes:str]` (grid-based
        access).

        For tensor-like derivatives, tensor-shaped structures can be created using  the `grad`, `hessian`,
        `tensor3` and `tensor4` member functions.
        """
        if isinstance(axes, str):
            self.differentiate([axes], access=access, recompute=recompute, debug=debug)
            if is_composite(axes) and axes != "Laplacian":
                return
            else:
                return self.derivatives [axes] if access == 'L' else \
                       self.derivativesG[axes]

        if recompute:
            self.invalidate_derivatives()

        if self.mesh.dim == 1:
            # No completion needed
            for ax in axes:
                if is_composite(ax):
                    raise ValueError(
                        f"1D LagrangeMesh objects do not support composite derivatives: '{axes}'."
                    )
        else:
            for ax in axes:
                if is_composite(ax) and ax not in ['Grad', 'Hessian', 'Laplacian', 'Tensor3', 'Tensor4']:
                    raise ValueError(f"Unrecognized composite derivative: {ax}")

            # Add components to allow reuse of derivatives:
            if 'Grad' in axes:
                if not self.mesh.highest_derivative_order >=1:
                    raise ValueError(f"Derivative '{axes}' incompatible with highest derivative order: {self.mesh.highest_derivative_order}.")
                axes = ['x', 'y', 'z'] + axes

            if 'Laplacian' in axes:
                if not self.mesh.highest_derivative_order >=2:
                    raise ValueError(f"Derivative '{axes}' incompatible with highest derivative order: {self.mesh.highest_derivative_order}.")
                axes = ['xx', 'yy', 'zz'] + axes

            if 'Hessian' in axes:
                if not self.mesh.highest_derivative_order >=2:
                    raise ValueError(f"Derivative '{axes}' incompatible with highest derivative order: {self.mesh.highest_derivative_order}.")
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
                if not self.mesh.highest_derivative_order >=3:
                    raise ValueError(f"Derivative '{axes}' incompatible with highest derivative order: {self.mesh.highest_derivative_order}.")
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
                if not self.mesh.highest_derivative_order >=4:
                    raise ValueError(f"Derivative '{axes}' incompatible with highest derivative order: {self.mesh.highest_derivative_order}.")
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
            self._differentiate1(axes=ax, debug=debug)

    def derivative_symmetry(self, axes):
        """Return the symmetry of the derivative of this MeshQuantity wrt axes .
        """
        symmetry = self.symmetry.copy()

        # an even number of differentiations keeps the symmetry sign
        # an odd number of differentiations flips the symmetry sign
        # x-axis
        if self.mesh.reduced[0]:
            if axes.count('x') % 2:
                symmetry[:,0] *= -1
        # y-axis
        if self.mesh.dim > 1 and self.mesh.reduced[1]:
            if axes.count('y') % 2:
                symmetry[:,1] *= -1
        # z-axis
        if self.mesh.dim > 2 and self.mesh.reduced[2]:
            if axes.count('z') % 2:
                symmetry[:,2] *= -1

        return symmetry

    def invalidate_derivatives(self):
        """Mark all derivatives as outdated."""
        for axes in self._derivative_is_uptodate.keys():
            self._derivative_is_uptodate[axes] = False
        for axes in self.derivativesG.keys():
            self._derivative_is_uptodate[axes] = False

    def derivative_set_uptodate(self, axes, value=True):
        """Indicate that the derivative wrt axes was computed after modifying the MeshQuantity's data,
        and, thus, that the derivative uptodate relative to the MeshQuantity's data."""
        self._derivative_is_uptodate[axes] = value

    def derivative_is_uptodate(self, axes):
        """Is the derivative wrt axes upto date?"""
        return self._derivative_is_uptodate.get(axes, False)

    # Forwarding methods: Since the MeshQuantity stores (a reference to) the mesh on which it is defined, we can call
    # LagrangeMesh methods directly on the MeshQuantity.
    def integrate(self):
        """Integrate the MeshQuantity over the simulation volume."""
        return self.mesh.integrate(self)

    def dbg_assert(self):
        """Assert some conditions that may indicate bugs when `False`."""

        print(f"\nMeshQuantity.dbg_assert() called on instance `{self}`")
        # Test that access methods still correctly share memory. See issues/51.
        assert np.shares_memory(self.data, self.dataG), \
            (f"`data` and `dataG` are expected to share memory with different shapes. "
             f"See https://github.com/IAA-nuclear/tantalus_full/issues/51.")
        for ax,derivative in self.derivatives.items():
            assert np.shares_memory(derivative, self.derivatives[ax]), \
                (f"`data` and `dataG` are expected to share memory with different shapes. "
                 f"See https://github.com/IAA-nuclear/tantalus_full/issues/51.")
