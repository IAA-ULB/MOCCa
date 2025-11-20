import numpy as np


# ==============================================================================
# Base classes
# ==============================================================================
class Operator:
    """Base class for operators.

    Attributes:
        ket: ket object based on hfpsi
        bra: bra object based on hfpsi
        O_ket: the result of applying the operator to ket. this must be created in one of the methods
            add_local_terms, add_non_local_terms.
        mesh: mesh on which hfpsi is discretised.
        derivatives: list of derivatives needed by the operator.
    Methods:
        add_local_terms, add_non_local_terms: these methods can be overridden to
            define the action of the operator to ket resulting in an observable O_ket. In principle,
            there is only need for one method to be overridden, one that implements the action of the
            operator on the ket. However, the names of the three methods serve as a hint to the things
            that must be computed.
        compute_matrix_representation: compute <bra|O_ket >. This must not be overridden
    """
    def __init__(self, hfpsi, derivatives=None):
        """Initialize the operator with a wave function `hfpsi` to operate on.
        Args:
            derivatives: str|list of derivatives needed by the operator.
        """
        self.ket = hfpsi
        self.bra = hfpsi
        self.mesh = hfpsi.mesh # convenient
        # Allocate space for a matrix represention of this operator relative to hfpsi.
        # It has the same structure as hfblocks
        self.hfblockrange = hfpsi.hfblockrange
        self.derivatives = [] if derivatives is None else \
                           [derivatives] if isinstance(derivatives, str) else \
                            derivatives
        if self.mesh.dim == 1:
            # replace 'laplacian' with 'xx'
            self.derivatives = [ 'xx' if axes == 'Laplacian' else axes for axes in self.derivatives]

    # --------------------------------------------------------------------------
    # Methods to be overridden by derived classes (at least one of them)
    # --------------------------------------------------------------------------
    def apply_base_operator(self):
        """In case the derived class is to be used as a base class itself.
        See KineticEnergyOperator for an example.
        """
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

    # --------------------------------------------------------------------------
    # Methods NOT to be overridden by derived classes
    # --------------------------------------------------------------------------
    def compute_derivatives(self):
        """Compute all derivatives needed by the operator."""
        self.ket.differentiate(self.derivatives)

    def compute_action(self):
        """Compute the action of the operator on the wave function.

        Returns:
            self.O_ket: O|ket>
        """
        self.compute_derivatives()
        self.apply_base_operator()
        self.add_local_terms()
        self.add_non_local_terms()

        if not hasattr(self, 'O_ket'):
            raise AttributeError(
                f"Calling {self.__class__.__name__}.compute_action() failed to create "
                f"Attribute 'self.O_ket'. This is probably an implementation error."
            )

        return self.O_ket

    def compute_matrix_representation(self):
        """Compute the matrix representation of the operator wrt `self.ket`.

        Raises:
            AttributeError: if self.O_ket is not created in one of the override methods or if they
                were not called.
        """
        self.compute_action()
        if not hasattr(self, 'O_ket'):
            raise AttributeError(f"Attribute 'self.O_ket' is missing in class {self.__class__.__name__}. \n"
                                 f"One of the derived methods compute_derivatives, add_local_terms, "
                                 f"add_non_local_terms must define `self.O_ket`, where the action of "
                                 f"the operator on `self.ket` is stored.")

        if not hasattr(self, 'matrix'):
            self.matrix = self.ket.allocate_matrix_representation()

        bra_data   = self.  bra.data.reshape((self.mesh.linear_size, 4, self.  bra.n_total_wf), order='F')
        O_ket_data = self.O_ket.data.reshape((self.mesh.linear_size, 4, self.O_ket.n_total_wf), order='F')
        for ib in range(8):
            block = self.matrix[ib]
            # Alternative formulation
            blockstart, blockstop =  self.hfblockrange[ib][0], self.hfblockrange[ib][1]
            bra_data_ib   =   bra_data[:,:,blockstart:blockstop]
            O_ket_data_ib = O_ket_data[:,:,blockstart:blockstop]
            np.einsum("ijk,ijl->kl", bra_data_ib, O_ket_data_ib, out=block) * self.mesh.dv
        return self.matrix

    def compute_diagonal_elements(self):
        """Compute only the diagonal elements of the operator wrt `self.ket`.

        Raises:
            AttributeError: if self.O_ket is not created in one of the override methods or if they
                were not called.

        """
        self.compute_action()

        if not hasattr(self, 'matrix'):
            self.diagonal = np.empty(self.ket.n_total_wf, dtype=np.float64)

        bra_data   = self.  bra.data.reshape((self.mesh.linear_size, 4, self.bra  .n_total_wf), order='F')
        O_ket_data = self.O_ket.data.reshape((self.mesh.linear_size, 4, self.O_ket.n_total_wf), order='F')
        np.einsum("ijk,ijk->k", bra_data, O_ket_data, out=self.diagonal) * self.mesh.dv
        return self.diagonal
