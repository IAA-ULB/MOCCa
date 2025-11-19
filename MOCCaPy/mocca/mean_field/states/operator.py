import numpy as np

from .hfpsi import HFPsi


class Operator:
    """Base class for operators.

    Attributes:
        ket: ket object based on hfpsi
        bra: bra object based on hfpsi
        O_ket: the result of applying the operator to ket. this must be created in one of the methods
            compute_derivatives, add_local_terms, add_non_local_terms.
    Methods:
        compute_derivatives, add_local_terms, add_non_local_terms: these methods can be overridden to
            define the action of the operator to ket resulting in an observable O_ket. In principle,
            there is only need for one method to be overridden, one that implements the action of the
            operator on the ket. However, the names of the three methods serve as a hint to the things
            that must be computed.
        compute_matrix_representation: compute <bra|O_ket >. This must not be overridden
    """
    def __init__(self, hfpsi):
        """Initialize the operator with a wave function `hfpsi` to operate on."""
        self.ket = hfpsi
        self.bra = hfpsi
        # Allocate space for a matrix represention of this operator relative to hfpsi.
        # It has the same structure as hfblocks
        self.matrix = hfpsi.allocate_matrix_representation()

    def compute_derivatives(self):
        """Override this to compute derivatives."""
        pass

    def add_local_terms(self):
        """Override this to compute derivatives.
        Usually, this is where O_ket is created as
            `self.O_ket = HFPsi.like(self.ket)`
        This creates an empty HFPsi object with the same shape as self.ket.
        """
        pass

    def add_non_local_terms(self):
        """Override this to compute derivatives."""
        pass

    def compute_matrix_representation(self):
        """Compute the matrix representation of the operator wrt `self.ket`.

        Raises:
            AttributeError: if self.O_ket is not created in one of the override methods.

        """
        self.compute_derivatives()
        self.add_local_terms()
        self.add_non_local_terms()

        try:
            O_ket = self.O_ket
        except AttributeError:
            raise AttributeError(f"Attribute 'self.O_ket' is missing in class {self.__class__.__name__}. \n"
                                 f"One of the derived methods compute_derivatives, add_local_terms, "
                                 f"add_non_local_terms must define `self.O_ket`, where the action of "
                                 f"the operator on `self.ket` is stored.")
        bra = self.bra
        for ib in range(8):
            block = self.matrix[ib]
            for j in range(*self.bra.hfblockrange[ib]):
                for i in range(*self.ket.hfblockrange[ib]):
                    # TODO: use out= ? for performance
                    block[i-self.ket.hfblockrange[ib][0],j-self.bra.hfblockrange[ib][0]] = (
                        bra.data[:,   4*j] * O_ket.data[:,   4*i] +
                        bra.data[:, 1+4*j] * O_ket.data[:, 1+4*i] +
                        bra.data[:, 2+4*j] * O_ket.data[:, 2+4*i] +
                        bra.data[:, 3+4*j] * O_ket.data[:, 3+4*i]
                    ).sum(0)
        return self.matrix


class Overlap(Operator):
    """Compute the matrix representation of the overlap operator. <hfpsi|hfpsi>."""
    def __init__(self, hfpsi):
        super().__init__(hfpsi)

    def add_local_terms(self):
        self.O_ket = self.ket


class HamiltonianWoodsSaxon(Operator):
    """Compute the matrix representation of the hamiltonian h with a Woods-Saxon potential. <hfpsi|h|hfpsi>."""
    def __init__(self, hfpsi, hbm=20.73553000, V0=-50.0, r0=1.25, a=0.5):
        """Initialize the operator with a wave function `hfpsi` to operate on and set the parameters for the
        operator.

        Args:
            hbm: prefactor of the kinetic energy operator.
            V0, r0, a:  parameters of the Woods-Saxon potential
                (cfr https://en.wikipedia.org/wiki/Woods–Saxon_potential).
        """
        super().__init__(hfpsi)
        self.hbm = hbm
        self.V0 = V0
        self.ainv = 1/a
        self.R = r0 * np.pow(hfpsi.n_neutrons + hfpsi.n_protons, 1/3)

    def compute_derivatives(self):
        """Compute the derivatives needed for the kinetic energy operator."""
        axes = 'Laplacian' if self.ket.mesh.dim >= 2 else \
               'xx'
        self.ket.differentiate(axes)

    def add_local_terms(self):
        """Create self.O_ket, fill it with the Woods-Saxon potential, and add the kinetic energy."""
        self.O_ket = HFPsi.like(self.ket)

        # TODO: speed up with numba decorators?

        def VWoodsSaxon1D(r):
            return self.V0 / (1. + np.exp(self.ainv * (r - self.R)))

        def VWoodsSaxon2D(x,y):
            return self.V0 / (1. + np.exp(self.ainv * (np.sqrt(x*x + y*y) - self.R)))

        def VWoodsSaxon3D(x,y,z):
            return self.V0 / (1. + np.exp(self.ainv * (np.sqrt(x*x + y*y + z*z) - self.R)))

        dim = self.ket.mesh.dim
        if dim == 3:
            VWoodsSaxon = VWoodsSaxon3D
        elif dim == 2:
            VWoodsSaxon = VWoodsSaxon2D
        else: # dim == 1:
            VWoodsSaxon = VWoodsSaxon1D

        Vr = self.ket.mesh.apply(VWoodsSaxon)
        self.O_ket.data = Vr * self.O_ket.data
        axes = 'Laplacian' if self.ket.mesh.dim >= 2 else \
               'xx'
        self.O_ket.data += self.hbm * self.ket.derivatives[axes]


