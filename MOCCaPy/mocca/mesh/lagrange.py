import numpy as np
import numpy.typing as npt

from mocca.mesh.observable import Observable
from mocca.mesh.lagrange_function import lagrange_function

def create_mesh(gx, gy=None, gz=None):
    """Create a rectangular mesh from 1D arrays with coordinates.

    Args:
        gx (ndarray): 1D array
        gy (ndarray): 1D array
        gz (ndarray): 1D array
            gx, gy, gz may refer to the same array.
    Returns:
        an array (order='F)
        if only gx is provided, gx is returned
        if only gx and gy are provided, gxy of shape (nx*ny,2) is returned
        if gx, gy and gz are provided, gxzy of shape (nx*ny*nz,3) is returned
    """
    if gz is not None:
        assert gy is not None, "Argument gy is required"
        nx = len(gx)
        ny = len(gy)
        nz = len(gz)
        gxyz = np.empty((nx*ny*nz,3), dtype=float, order='F')
        gxy = create_mesh(gx, gy)
        nxy = nx*ny
        for iz in range(nz):
            gxyz[iz * nxy:(iz + 1) * nxy, 0:2] = gxy
            gxyz[iz * nxy:(iz + 1) * nxy, 2] = gz[iz]
        return gxyz

    elif gy is not None:
        nx = len(gx)
        ny = len(gy)
        gxy = np.empty((nx*ny,2), dtype=float, order='F')
        for iy in range(ny):
            gxy[iy * nx:(iy + 1) * nx, 0] = gx
            gxy[iy * nx:(iy + 1) * nx, 1] = gy[iy]
        return gxy
    else:
        return gx


class LagrangeMesh:
    def __init__(self, M: int|tuple,
                       d: int|float|tuple,
                       dim: int=0,
                       bc:str='antiperiodic',
                       reduced:tuple|bool=True,
                       shift:tuple|float=.0,
                 ) -> None:
        """Construct a Lagrange mesh in 1, 2 or 3 dimensions.

        Args:
            dim (int): dimension of mesh. if not specified, dim is guessed as len(M), where M must be a tuple.
            M: tuple with number of grid points in the respective dimensions. Each component of M=2N must be even,
                yielding N grid points at both sides of the origin. If a single M is specified, M is the same in every
                direction.
                The distinction between M and N is consistent with  Ryssens et al, PHYSICAL REVIEW C 92, 064318 (2015)
                section III.B Lagrange-mesh representation.
            d: spacing of points in the respective dimensions. If d is a float or an int d is the same in each direction.
            bc: boundary condition type. 'antiperiodic' or 'periodic'.
            reduced: restrict the mesh to the positive half-axis.
            shift: subtract shift from the grid points. If non-zero, corresponding `reduced` entry must be `False`.

        Raises:
            AssertionError: in case of invalid choices

        Remark:
            The `reduced` parameter is derived from symmetry considerations and may at some point - when the complexity
             of Hephaestos is taken into account - be replaced with a `Symmetry` object. For the time being, however,
             we content with explicitly indicating which coordinate axes must be 'reduced'
        """
        # initialize M and d as dim-tuples
        if dim == 0 and isinstance(M, tuple):
            self.dim = len(M)
        else:
            self.dim = dim

        assert 1 <= self.dim <= 3

        if isinstance(M, int):
            assert M > 0, "M must be strictly positive."
            assert M % 2 == 0, "M must be an even number."
            self.M = tuple(self.dim*[M])
        else:
            assert isinstance(M, tuple)
            assert len(M) == self.dim
            for mi in M:
                assert mi >= 0, "M must be strictly positive."
                assert mi % 2 == 0, "M must be an even number."
            for mi in M:
                assert isinstance(mi, int)
            self.M = M

        # box spacing (ints are converted to floats)
        if isinstance(d,tuple):
            assert len(d) == self.dim
            for di in d:
                assert di > 0, "d must be strictly positive."
            self.d = np.array([float(di) for di in d])
        else:
            assert d > 0, "d must be strictly positive."
            self.d = np.array(self.dim*[float(d)])

        # validate boundary condition
        assert bc in ['antiperiodic', 'periodic']
        self.bc = bc
        # convenience attributes
        self.antiperiodic = bc == 'antiperiodic'
        self.periodic = not self.antiperiodic # since there are only 2 options.

        # validate reduced
        if isinstance(reduced, bool):
            self.reduced = tuple(self.dim*[reduced])
        else:
            assert isinstance(reduced, tuple)
            assert len(reduced) == self.dim
            self.reduced = reduced

        # Compute dv (for integration)
        self.dv = np.prod(self.d) * 2 ** self.reduced.count(True)

        # Compute the unreduced (!) box widths
        # The full (unreduced box width is needed by the plane wave base functions
        self.box_width = np.array([M*d for (M,d) in zip(self.M, self.d)])

        # validate shift
        if isinstance(shift, float):
            self.shift = tuple(self.dim*[shift])
        else:
            assert isinstance(shift, tuple)
            assert len(shift) == self.dim
            self.shift = shift
        for (r,s) in zip(self.reduced, self.shift):
            if r:
                assert s == 0., "A nonzero shift cannot be applied when reduced is True."

        # initialize grid points:
        start = self.dim*[.0]
        n_reduced = self.dim*[0]
        g1D = self.dim*[0]
        for i, (n_i, shift_i, reduced_i, d_i) in enumerate(zip(self.M, self.shift, self.reduced, self.d)):
            if reduced_i:
                n_reduced[i] = n_i//2
                start[i] = 0.5
            else:
                n_reduced[i] = n_i
                start[i] = -(n_i - 1)*0.5 - shift_i/d_i

            g1D[i] = np.linspace(start[i], start[i] + n_reduced[i] - 1, n_reduced[i])
            g1D[i] *= d_i
            # print(f"{i=} {g1D[i]=}")

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
            # for k in range(n_reduced[2]):
            #     for j in range(n_reduced[1]):
            #         for i in range(n_reduced[0]):
            #             self.gridx[i, j, k] = g1D[0][i]
            #             self.gridy[i, j, k] = g1D[1][j]
            #             self.gridz[i, j, k] = g1D[2][k]
            for k in range(n_reduced[2]):
                for j in range(n_reduced[1]):
                    self.gridx[:,j,k] = g1D[0]
            for k in range(n_reduced[2]):
                for i in range(n_reduced[0]):
                    self.gridy[i,:,k] = g1D[1]
            for j in range(n_reduced[1]):
                for i in range(n_reduced[0]):
                    self.gridz[i,j,:] = g1D[2]
            # print(f"{self.gridx=}")
            # print(f"{self.gridy=}")
            # print(f"{self.gridz=}")
            # print(f"{self.gridx.ravel(order='F')=}")
            # print(f"{self.gridy.ravel(order='F')=}")
            # print(f"{self.gridz.ravel(order='F')=}")

        elif self.dim == 2:
            self.grid = create_mesh(self.gx, self.gy)

            self.gridx = np.empty(n_reduced, order='F')
            self.gridy = np.empty(n_reduced, order='F')
            for j in range(n_reduced[1]):
                self.gridx[:, j] = g1D[0]
            for i in range(n_reduced[0]):
                self.gridy[i, :] = g1D[1]
            # print(f"{self.gridx=}")
            # print(f"{self.gridy=}")
            # print(f"{self.gridx.ravel(order='F')=}")
            # print(f"{self.gridy.ravel(order='F')=}")

        else: # self.dim == 1
            self.grid = create_mesh(self.gx)
            self.gridx = g1D[0]
            # print(f"{self.gridx=}")

        self.shape = self.gridx.shape
        self.flat_shape = (int(np.prod(self.shape)),)

    def flatten(self, array=None):
        """Flatten the mesh array.

        Args:
            array: if None, self.gridx|y|z are flattened. Otherwise, array is flattened.

        Returns:
            None or the flattened array.

        Raises:
            AssertionError: if array.shape != self.shape.
        """
        if array is None:
            # self._reshape(flat=True)
            if len(self.gridx.shape) == 1:
                # gridx/gridy/gridz already flat
                # (We assume that either all or none of them are flat or unflattened at the same time)
                return None
            else:
                try:
                    self.gridx = self.flatten(self.gridx)
                    # AttributeError raised if one of gridy/z does not exist
                    self.gridy = self.flatten(self.gridy)
                    self.gridz = self.flatten(self.gridz)
                except AttributeError:
                    pass

        else:
            if array.shape[0] == self.flat_shape[0]:
                # already flat
                return array
            else:
                # collapse the first dim dimensions into one and copy the remaining dimensions:
                flat_shape = list(self.flat_shape) # a single element
                for d in array.shape[self.dim:]:
                    flat_shape.append(d)
                flat_shape = tuple(flat_shape)
                return array.reshape(flat_shape, order='F')

    def unflatten(self, array=None) -> None:
        """Unflatten the mesh array. For most computations

        Args:
            array: if None, self.gridx|y|z are unflattened. Otherwise, array is unflattened.

        Returns:
            None or the unflattened array.

        Raises:
            AssertionError: if array.shape != self.flat_shape.
        """
        if array is None:
            if len(self.gridx.shape) == self.dim:
                # gridx/gridy/gridz already unflattened
                # (We assume that either all or none of them are flat or unflattened at the same time)
                return None
            else:
                try:
                    self.gridx = self.unflatten(self.gridx)
                    # AttributeError raised if one of gridy/z does not exist
                    self.gridy = self.unflatten(self.gridy)
                    self.gridz = self.unflatten(self.gridz)
                except AttributeError:
                    pass
        else:
            if array.shape[0:self.dim] == self.shape[0]:
                # already unflattened
                return array
            else:
                # Take the unflattened shape and append the multicomponent dimensions of array
                shape = list(self.shape)
                for d in array.shape[1:]:
                    shape.append(d)
                shape = tuple(shape)
                return array.reshape(shape, order='F')

    def apply(self, function):
        """Apply a function on the mesh, i.e. compute the function value on every grid point.

        Args:
            function: function f(x,y,z) to be applied on the mesh array. If f accepts additional parameters, a closure
                must be defined.

        Returns:
            a flat mesh array (self.flat_shape).

        Raises:
            AssertionError: if `not q.shape in [self.shape, self.flat_shape]`.
        """
        # TODO use arguments for output variables?
        self.flatten()
        return function(self.gridx, self.gridy, self.gridz)


    def integrate(self, Q):
        """Compute the integral of a scalar quantity `q` on the mesh.

        Args:
            Q (Observable) : discretised on the grid.

        Returns:
            a scalar:

        Raises:
            AssertionError: if `not q.shape in [self.shape, self.flat_shape]`.
        """
        n_components = Q.n_components
        result = np.zeros((n_components,), dtype=float)
        for ic in range(n_components):
            q = Q[ic]
            result[ic] = q.sum() * self.dv
        return result if n_components > 1 else result[0]


    def plane_wave_1D(self, L:float, k:float, r:np.array):
        """Evaluate the plane wave basis function at with wave vector `k` at position `r`.

        According to Ryssens et al, PHYSICAL REVIEW C 92, 064318 (2015), eq(14)
        According to [eq 5.1] in https://github.com/IAA-nuclear/tantalus_full/blob/MOCCaPy/MOCCaPy/mocca/mesh/lagrange.pdf

        Args:
            k: wave vector, shape is `(self.dim,)`. values must be odd half integers ±1/2, ±3/2, ...
            r: position vector, shape is `(self.dim, nr), nr being the number of evaluation points.
            as_complex: if True, return a complex array. otherwise return two arrays with the real and imaginary parts,
                resp.

        Returns:
            array of shape `(nr,2)` with complex and imaginary parts
        """
        arg = (2*np.pi * k / L) * r
        result = np.empty((len(r),2), dtype=float, order='F')
        result[:,0] = np.cos(arg)
        result[:,1] = np.sin(arg)
        result *= np.sqrt(1/L)
        return result


    def basis_function(self, ijk, r, negative_axis=False):
        """Evaluate the plane wave basis function corresponding to grid point `x_ijk = [x_i,y_j,z_k]` 
        at position `r`. The wave vector `k` is related to `x_i` as `k = x_i/dx`.

        According to [eq 5.1] in https://github.com/IAA-nuclear/tantalus_full/blob/MOCCaPy/MOCCaPy/mocca/mesh/lagrange.pdf

        Args:
            ijk: grid point index. A D-tuple of grid indices with `D == self.dim'.
        """
        if self.dim == 3:
            i = ijk[0]
            j = ijk[1]
            k = ijk[2]
            two_pi_K = np.array([ 2. * np.pi * ( self.gx[i] if not self.reduced[0] else
                                                 self.gx[i] if i>0 else
                                                -self.gx[i] ) / (self.box_width[0] * self.d[0])
                                , 2. * np.pi * ( self.gy[j] if not self.reduced[0] else
                                                 self.gy[j] if i>0 else
                                                -self.gy[j] ) / (self.box_width[1] * self.d[1])
                                , 2. * np.pi * ( self.gz[k] if not self.reduced[0] else
                                                 self.gz[k] if i>0 else
                                                -self.gz[k] ) / (self.box_width[2] * self.d[2])
                                ])
            arg = r @ two_pi_K
            factor = np.sqrt(1/(self.box_width[0] * self.box_width[1] * self.box_width[2]))

        elif self.dim == 2:
            i = ijk[0]
            j = ijk[1]
            two_pi_K = np.array([ 2. * np.pi * ( self.gx[i] if not self.reduced[0] else
                                                 self.gx[i] if i>0 else
                                                -self.gx[i] ) / (self.box_width[0] * self.d[0])
                                , 2. * np.pi * ( self.gy[j] if not self.reduced[0] else
                                                 self.gy[j] if i>0 else
                                                -self.gy[j] ) / (self.box_width[1] * self.d[1])
                                ])
            arg = r @ two_pi_K
            factor = np.sqrt(1/(self.box_width[0] * self.box_width[1]))

        elif self.dim == 1:
            # ijk == i
            i = ijk
            two_pi_K = 2. * np.pi * self.gx[i]  / (self.box_width[0] * self.d[0])
            if self.reduced[0] and negative_axis:
                arg = r * (-two_pi_K)
            else:
                arg = r * two_pi_K

            factor = np.sqrt(1 / self.box_width[0])

        nr = r.shape[0]
        result = np.empty((nr,2), dtype=float, order='F')
        result[:,0] = np.cos(arg)
        result[:,1] = np.sin(arg)
        result *= factor
        return result

    def derive1(self, Q:npt.NDArray):
        """Compute the 1st order derivative of a scalar quantity `q` on the mesh.

        Args:
            q: scalar quantity discretised on the grid. Thus `q.shape in [self.shape, self.flat_shape]` evaluates
                to True

        Returns:


        Raises:
            AssertionError: if `not q.shape in [self.shape, self.flat_shape]`.
        """
        # TODO : implement

    def derive2(self, Q:npt.NDArray):
        """Compute the 2nd order derivative of a scalar quantity `q` on the mesh.

        Args:
            q: scalar quantity discretised on the grid. Thus `q.shape in [self.shape, self.flat_shape]` evaluates
                to True

        Returns:


        Raises:
            AssertionError: if `not q.shape in [self.shape, self.flat_shape]`.
        """
        # TODO : implement
        
    
    def is_flat(self, Q:Observable) -> bool:
        """Determine if `q` is flat.
        Args:
            Q: a quantity discretised on the grid.
        """
        return Q.data.shape[0] == self.flat_shape[0]
    
    
    def is_unflattened(self, Q:npt.NDArray) -> bool:
        """Determine if `q` is unflattened.
        Args:
            Q: a quantity discretised on the grid.
        """
        return Q.shape[0:self.dim] == self.shape        


    def is_grid_quantity(self, Q:npt.NDArray) -> bool:
        """Determine if `Q` is a grid quantity.
        Args:
            Q: a quantity discretised on the grid.
        """
        return self.is_flat(Q) or self.is_unflattened(Q)


    def n_gridpoints(self):
        return np.prod(self.shape)
    
    def get_number_of_components(self, Q:npt.NDArray):
        """compute the number of components in agrid quantity `Q`.
        Args:
            Q: a quantity discretised on the grid.
        """
        return len(Q)//self.flat_shape[0]
    
    
    def get_component(self, Q:npt.NDArray, i:int):
        """Get the i-th component of `Q`."""
        return Q[i*self.flat_shape[0]:(i+1)*self.flat_shape[0]]


    def interpolate(self, Q:Observable, r:npt.NDArray) -> npt.NDArray:
        """Interpolate a quantity `Q` on the mesh.

        Args:
            Q: representation of a scalar quantity on the grid.
            r: array of `p` points at which to interpolate the scalar quantity `q`. 'r.shape == (p, self.dim)'
        """
        self.flatten()
        
        if self.dim == 3:
            return self._interpolate3D(Q, r)
        elif self.dim == 2:
            return self._interpolate2D(Q, r)
        else:
            return self._interpolate1D(Q, r)


    def _interpolate1D(self, Q:Observable, r:npt.NDArray) -> npt.NDArray:
        """Interpolate `Q` on a 1D  mesh."""
        nr = r.shape[0]
        nq = Q.n_components
        Qr = np.zeros(shape=(nr,nq), dtype=float, order='F')
        N = len(self.gridx) if self.reduced[0] else len(self.gridx)//2
        for ig in range(self.n_gridpoints()):
            for iq in range(nq):
                q = Q[iq]
                sign = Q.symmetry[iq]
                fxa = self.lagrange_function(r, ig, sign)
                if len(fxa) == 1:
                    Qr[:,iq] += q[ig] * fxa[0]
                else:
                    Qr[:, iq] += q[ig] * (fxa[0] + fxa[1])
        return Qr


    def _interpolate2D(self, Q:Observable, r:npt.NDArray) -> npt.NDArray:
        """Interpolate `Q` on a 2D mesh.

        Args:
            Q: An observable with values specified on all grid points
            r: array of points at which to interpolate Q. 'r.shape == (nr, self.dim)'
        """
        nr = r.shape[0]
        nq = Q.n_components
        Qr = np.zeros(shape=(nr,nq), dtype=float, order='F')
        N = len(self.gridx) if self.reduced[0] else len(self.gridx)//2
        ig = 0
        for iy in range(self.M[1]):
            for ix in range(self.M[0]):
                for iq in range(nq):
                    q = Q[iq]
                    sign = Q.symmetry[iq]
                    fxya = self.lagrange_function (r, (ix,iy), sign)
                    if len(fxya) == 1:
                        Qr[:,iq] += q[ig] * fxya[0]
                    else:
                        Qr[:, iq] += q[ig] * (fxya[0] + fxya[1])
                ig += 1
        return Qr
        # nr = r.shape[0]
        # nq = Q.n_components
        # Qr = np.zeros(shape=(nr,nq), dtype=float, order='F') # Q interpolated at r
        # if self.reduced[0]:
        #     if self.reduced[1]:
        #         # both axes reduced
        #         nx = self.M[0] // 2
        #         ny = self.M[1] // 2
        #         ig = 0
        #         for ix in range(nx):
        #             for iy in range(ny):
        #                 for iq in range(nq):
        #                     q = Q[iq]
        #                     sign = Q.symmetry[iq]
        #                     lf__x  = lagrange_function(r,  self.gridx[ix], self.d[0], self.M[0])
        #                     lf_mx = lagrange_function(r, -self.gridx[ix], self.d[0], self.M[0])
        #                     lf__y  = lagrange_function(r,  self.gridy[iy], self.d[1], self.M[1])
        #                     lf_my = lagrange_function(r, -self.gridy[iy], self.d[1], self.M[1])
        #
        #                     Qr[:,iq] += q[ig] * ( (        lf__x        * lf__y ) +
        #                                           (        lf__x * sign * lf_my ) +
        #                                           ( sign * lf_mx        * lf__y ) +
        #                                           (        lf_mx        * lf_my ) ) # sign * sign == 1
        #                 ig +=1
        #
        #     else:
        #         # only x-axis reduced
        #         nx = self.M[0] // 2
        #         ny = self.M[1]
        #         ig = 0
        #         for ix in range(nx):
        #             for iy in range(ny):
        #                 for iq in range(nq):
        #                     q = Q[iq]
        #                     sign = Q.symmetry[iq]
        #                     lf__x = lagrange_function(r,  self.gridx[ix], self.d[0], self.M[0])
        #                     lf_mx = lagrange_function(r, -self.gridx[ix], self.d[0], self.M[0])
        #                     lf__y = lagrange_function(r,  self.gridy[iy], self.d[1], self.M[1])
        #
        #                     Qr[:,iq] += q[ig] * ( lf__x + sign * lf_mx ) * lf__y
        #                 ig += 1
        #
        # else:
        #     if self.reduced[1]:
        #         # only y-axis reduced
        #         nx = self.M[0]
        #         ny = self.M[1] // 2
        #         ig = 0
        #         for ix in range(nx):
        #             for iy in range(ny):
        #                 for iq in range(nq):
        #                     q = Q[iq]
        #                     sign = Q.symmetry[iq]
        #                     lf__x  = lagrange_function(r,  self.gridx[ix], self.d[0], self.M[0])
        #                     lf__y  = lagrange_function(r,  self.gridy[iy], self.d[1], self.M[1])
        #                     lf_my = lagrange_function(r, -self.gridy[iy], self.d[1], self.M[1])
        #
        #                     Qr[:,iq] += q[ig] * lf__x * ( lf__y + sign * lf_my)
        #                 ig += 1
        #     else:
        #         # none of the axes reduced
        #         nx = self.M[0]
        #         ny = self.M[1]
        #         ig = 0
        #         for ix in range(nx):
        #             for iy in range(ny):
        #                 for iq in range(nq):
        #                     q = Q[iq]
        #                     lf__x  = lagrange_function(r,  self.gridx[ix], self.d[0], self.M[0])
        #                     lf__y  = lagrange_function(r,  self.gridy[iy], self.d[1], self.M[1])
        #                     Qr[:,iq] += q[ig] * lf__x * lf__y
        #                 ig +=1
        #
        # return Qr


    def _interpolate3D(self, Q:Observable, r:npt.NDArray) -> npt.NDArray:
        """Interpolate `Q` on a 2D  mesh."""
        nr = r.shape[0]
        nq = Q.n_components
        Qr = np.zeros(shape=(nr,nq), dtype=float, order='F') # Q interpolated at r

        if self.reduced[0]:
            if self.reduced[1]:
                if self.reduced[2]: # True True True
                    nx = self.M[0] // 2
                    ny = self.M[1] // 2
                    nz = self.M[2] // 2
                    for ix in range(nx):
                        for iy in range(ny):
                            for iz in range(nz):
                                for iq in range(nq):
                                    q = Q[iq]
                                    sign = Q.symmetry[iq]
                                    lf__x = lagrange_function(r,  self.gridx[ix], self.d[0], self.M[0])
                                    lf_mx = lagrange_function(r, -self.gridx[ix], self.d[0], self.M[0])
                                    lf__y = lagrange_function(r,  self.gridy[iy], self.d[1], self.M[1])
                                    lf_my = lagrange_function(r, -self.gridy[iy], self.d[1], self.M[1])
                                    lf__z = lagrange_function(r,  self.gridz[iz], self.d[2], self.M[2])
                                    lf_mz = lagrange_function(r, -self.gridz[iz], self.d[2], self.M[2])

                                    Qr[:,iq] += q[ig] * ( (        lf__x        * lf__y        * lf__z ) + # 0
                                                          (        lf__x        * lf__y * sign * lf_mz ) + # 1
                                                          (        lf__x * sign * lf_my        * lf__z ) + # 1
                                                          (        lf__x        * lf_my        * lf_mz ) + # 2  sign * sign == 1
                                                          ( sign * lf_mx        * lf__y        * lf__z ) + # 1
                                                          (        lf_mx        * lf__y        * lf_mz ) + # 2  sign * sign == 1
                                                          (        lf_mx        * lf_my        * lf__z ) + # 2  sign * sign == 1
                                                          ( sign * lf_mx        * lf_my        * lf_mz ) ) # 3  sign * sign * sign == sign

                else:               # True True False
                    nx = self.M[0] // 2
                    ny = self.M[1] // 2
                    nz = self.M[2]
                    for ix in range(nx):
                        for iy in range(ny):
                            for iz in range(nz):
                                for iq in range(nq):
                                    q = Q[iq]
                                    sign = Q.symmetry[iq]
                                    lf__x = lagrange_function(r,  self.gridx[ix], self.d[0], self.M[0])
                                    lf_mx = lagrange_function(r, -self.gridx[ix], self.d[0], self.M[0])
                                    lf__y = lagrange_function(r,  self.gridy[iy], self.d[1], self.M[1])
                                    lf_my = lagrange_function(r, -self.gridy[iy], self.d[1], self.M[1])
                                    lf__z = lagrange_function(r,  self.gridz[iz], self.d[2], self.M[2])

                                    Qr[:,iq] += q[ig] * ( (        lf__x        * lf__y ) + # 0
                                                          (        lf__x * sign * lf_my ) + # 1
                                                          ( sign * lf_mx        * lf__y ) + # 1
                                                          (        lf_mx        * lf_my ) ) * lf__z # 2  sign * sign == 1

            else:
                if self.reduced[2]: # True False True
                    nx = self.M[0] // 2
                    ny = self.M[1]
                    nz = self.M[2] // 2
                    for ix in range(nx):
                        for iy in range(ny):
                            for iz in range(nz):
                                for iq in range(nq):
                                    q = Q[iq]
                                    sign = Q.symmetry[iq]
                                    lf__x = lagrange_function(r,  self.gridx[ix], self.d[0], self.M[0])
                                    lf_mx = lagrange_function(r, -self.gridx[ix], self.d[0], self.M[0])
                                    lf__y = lagrange_function(r,  self.gridy[iy], self.d[1], self.M[1])
                                    lf__z = lagrange_function(r,  self.gridz[iz], self.d[2], self.M[2])
                                    lf_mz = lagrange_function(r, -self.gridz[iz], self.d[2], self.M[2])

                                    Qr[:,iq] += q[ig] * ( (        lf__x         * lf__z ) + # 0
                                                          (        lf__x  * sign * lf_mz ) + # 1
                                                          ( sign * lf_mx         * lf__z ) + # 1
                                                          (        lf_mx         * lf_mz ) ) * lf__y # 2  sign * sign == 1

                else:               # True False False
                    nx = self.M[0] // 2
                    ny = self.M[1]
                    nz = self.M[2]
                    for ix in range(nx):
                        for iy in range(ny):
                            for iz in range(nz):
                                for iq in range(nq):
                                    q = Q[iq]
                                    sign = Q.symmetry[iq]
                                    lf__x = lagrange_function(r,  self.gridx[ix], self.d[0], self.M[0])
                                    lf_mx = lagrange_function(r, -self.gridx[ix], self.d[0], self.M[0])
                                    lf__y = lagrange_function(r,  self.gridy[iy], self.d[1], self.M[1])
                                    lf__z = lagrange_function(r,  self.gridz[iz], self.d[2], self.M[2])

                                    Qr[:,iq] += q[ig] * ( ( lf__x + sign * lf_mx ) ) * lf__y * lf__z  # 1

        else:
            if self.reduced[1]:
                if self.reduced[2]: # False True True
                    nx = self.M[0]
                    ny = self.M[1] // 2
                    nz = self.M[2] // 2
                    for ix in range(nx):
                        for iy in range(ny):
                            for iz in range(nz):
                                for iq in range(nq):
                                    q = Q[iq]
                                    sign = Q.symmetry[iq]
                                    lf__x = lagrange_function(r,  self.gridx[ix], self.d[0], self.M[0])
                                    lf__y = lagrange_function(r,  self.gridy[iy], self.d[1], self.M[1])
                                    lf_my = lagrange_function(r, -self.gridy[iy], self.d[1], self.M[1])
                                    lf__z = lagrange_function(r,  self.gridz[iz], self.d[2], self.M[2])
                                    lf_mz = lagrange_function(r, -self.gridz[iz], self.d[2], self.M[2])

                                    Qr[:,iq] += q[ig] * ( lf__x * (        lf__y        * lf__z ) + # 0
                                                                  (        lf__y * sign * lf_mz ) + # 1
                                                                  ( sign * lf_my        * lf__z ) + # 1
                                                                  (        lf_my        * lf_mz ) ) # 2  sign * sign == 1

                else:               # False True False
                    nx = self.M[0]
                    ny = self.M[1] // 2
                    nz = self.M[2]
                    for ix in range(nx):
                        for iy in range(ny):
                            for iz in range(nz):
                                for iq in range(nq):
                                    q = Q[iq]
                                    sign = Q.symmetry[iq]
                                    lf__x = lagrange_function(r,  self.gridx[ix], self.d[0], self.M[0])
                                    lf__y = lagrange_function(r,  self.gridy[iy], self.d[1], self.M[1])
                                    lf_my = lagrange_function(r, -self.gridy[iy], self.d[1], self.M[1])
                                    lf__z = lagrange_function(r,  self.gridz[iz], self.d[2], self.M[2])

                                    Qr[:,iq] += q[ig] * lf__x * ( lf__y + sign * lf_my ) * lf__z  # 1

            else:                   # False False True
                if self.reduced[2]:
                    nx = self.M[0]
                    ny = self.M[1]
                    nz = self.M[2] // 2
                    for ix in range(nx):
                        for iy in range(ny):
                            for iz in range(nz):
                                for iq in range(nq):
                                    q = Q[iq]
                                    sign = Q.symmetry[iq]
                                    lf__x = lagrange_function(r,  self.gridx[ix], self.d[0], self.M[0])
                                    lf__y = lagrange_function(r,  self.gridy[iy], self.d[1], self.M[1])
                                    lf__z = lagrange_function(r,  self.gridz[iz], self.d[2], self.M[2])
                                    lf_mz = lagrange_function(r, -self.gridz[iz], self.d[2], self.M[2])

                                    Qr[:,iq] += q[ig] * lf__x * lf__y * ( lf__z + sign * lf_mz ) # 1

                else:               # False False False
                    nx = self.M[0]
                    ny = self.M[1]
                    nz = self.M[2]
                    for ix in range(nx):
                        for iy in range(ny):
                            for iz in range(nz):
                                for iq in range(nq):
                                    q = Q[iq]
                                    sign = Q.symmetry[iq]
                                    lf__x = lagrange_function(r,  self.gridx[ix], self.d[0], self.M[0])
                                    lf__y = lagrange_function(r,  self.gridy[iy], self.d[1], self.M[1])
                                    lf__z = lagrange_function(r,  self.gridz[iz], self.d[2], self.M[2])

                                    Qr[:,iq] += q[ig] * ( lf__x * lf__y * lf__z )

        return Qr

    def lagrange_function(self, x:npt.NDArray|float, ijk:tuple|int, sign=1):
        """Evaluate the Lagrange function corresponding to the `ijk` grid point at `x`.
        Args:
            x : array of points at which to evaluate the Lagrange function.
            ijk : index of the grid point selecting the gridpoint of the Lagrange function.
            sign : the symmetry ±1. For reduced axes the lagrange function of negative grid points
                is multplied by `sign`.
        Returns:
            1D Array with values of the Lagrange function corresponding to the `ijk` grid point at `x`.
            If an axis is reduced, several arrays are returned, one for each combination of positeve
            and negative grid points.
            reduced == [True,True,True] -> + + +   [True,True] -> + +
                                           + + -                  + -
                                           + - +                  - +
                                           + - -                  - -
                                           - + +        [True] -> +
                                           - + -                  -
                                           - - +
                                           - - -
            If one of the axes is not reduced the rows with a - on that axes are removed
            re
        """
        # The original idea to select grid points on the negative axes with negative indices doesn't work because
        # the first index is 0, and its sign is lost, as Python does not discriminate between -0 and +0.

        if self.dim == 3:
            i,j,k = ijk
            x_i, y_j, z_k = self.gx[i], self.gy[j], self.gz[k]

            if self.reduced[0]:
                fxa = [ lagrange_function(x[:, 0],  x_i, self.d[0], self.M[0]),
                        lagrange_function(x[:, 0], -x_i, self.d[0], self.M[0]) ]
                if sign == -1:
                    fxa[1] *= sign

                if self.reduced[1]:
                    fya  = [ lagrange_function(x[:, 1],  y_j, self.d[1], self.M[1]),
                             lagrange_function(x[:, 1], -y_j, self.d[1], self.M[1]) ]
                    if sign == -1:
                        fya[1] *= sign
                    if self.reduced[2]:
                        fza = [ lagrange_function(x[:, 2],  z_k, self.d[2], self.M[2]),
                                lagrange_function(x[:, 2], -z_k, self.d[2], self.M[2]) ]
                        if sign == -1:
                            fza[1] *= sign
                    else:
                        fza = [ lagrange_function(x[:, 2], z_k, self.d[2], self.M[2]) ]
                else:
                    fya = [ lagrange_function(x[:, 1],  y_j, self.d[1], self.M[1]) ]
                    if self.reduced[2]:
                        fza = [ lagrange_function(x[:, 2],  z_k, self.d[2], self.M[2]),
                                lagrange_function(x[:, 2], -z_k, self.d[2], self.M[2]) ]
                        if sign == -1:
                            fza[1] *= sign
                    else:
                        fza = [ lagrange_function(x[:, 2], z_k, self.d[2], self.M[2]) ]
            else:
                fxa = [ lagrange_function(x[:, 0],  x_i, self.d[0], self.M[0]) ]
                if self.reduced[1]:
                    fya = [ lagrange_function(x[:, 1],  y_j, self.d[1], self.M[1]),
                            lagrange_function(x[:, 1], -y_j, self.d[1], self.M[1]) ]
                    if sign == -1:
                        fya[1] *= sign
                    if self.reduced[2]:
                        fza = [ lagrange_function(x[:, 2],  z_k, self.d[2], self.M[2]),
                                lagrange_function(x[:, 2], -z_k, self.d[2], self.M[2]) ]
                        if sign == -1:
                            fza[1] *= sign
                    else:
                        fza = [ lagrange_function(x[:, 2], z_k, self.d[2], self.M[2]) ]
                else:
                    fya = [ lagrange_function(x[:, 1],  y_j, self.d[1], self.M[1]) ]
                    if self.reduced[2]:
                        fza = [ lagrange_function(x[:, 2],  z_k, self.d[2], self.M[2]),
                                lagrange_function(x[:, 2], -z_k, self.d[2], self.M[2]) ]
                        if sign == -1:
                            fza[1] *= sign
                    else:
                        fza = [ lagrange_function(x[:, 2], z_k, self.d[2], self.M[2]) ]

            result = []
            for fx in fxa:
                for fy in fya:
                    for fz in fza:
                        result.append( fx * fy * fz )
            return result

        elif self.dim == 2:
            i, j = ijk
            x_i, y_j = self.gx[i], self.gy[j]

            if self.reduced[0]:
                fxa = [ lagrange_function(x[:, 0],  x_i, self.d[0], self.M[0]),
                        lagrange_function(x[:, 0], -x_i, self.d[0], self.M[0]) ]
                if sign == -1:
                    fxa[1] *= sign

                if self.reduced[1]:
                    fya = [ lagrange_function(x[:, 1],  y_j, self.d[1], self.M[1]),
                            lagrange_function(x[:, 1], -y_j, self.d[1], self.M[1]) ]
                    if sign == -1:
                        fya[1] *= sign
                else:
                    fya = [ lagrange_function(x[:, 1],  y_j, self.d[1], self.M[1]) ]
            else:
                fxa = [lagrange_function(x[:, 0], x_i, self.d[0], self.M[0])]
                if self.reduced[1]:
                    fya = [ lagrange_function(x[:,1],  y_j, self.d[1], self.M[1]),
                            lagrange_function(x[:,1], -y_j, self.d[1], self.M[1]) ]
                    if sign == -1:
                        fya[1] *= sign
                else:
                    fya = [ lagrange_function(x[:,1],  y_j, self.d[1], self.M[1]) ]

            result = []
            for fx in fxa:
                for fy in fya:
                    result.append(fx * fy )
            return result

        else:
            i = ijk
            x_i = self.gridx[i]

            if self.reduced[0]:
                fxa = [ lagrange_function(x,  x_i, self.d[0], self.M[0]),
                        lagrange_function(x, -x_i, self.d[0], self.M[0]) ]
                if sign == -1:
                    fxa[1] *= sign
            else:
                fxa = [ lagrange_function(x, x_i, self.d[0], self.M[0]) ]

            return fxa
