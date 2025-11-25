import numpy as np


# ==============================================================================
# Base classes
# ==============================================================================
class Operator:
    """Base class for operators.

    Attributes:
        ket: ket object based on hfpsi
        O_ket: the result of applying the operator to ket. this must be created in one of the methods
            add_local_terms, add_non_local_terms.
        mesh: mesh on which hfpsi is discretized.
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
        if not hasattr(hfpsi, 'd3'):
            hfpsi.d3 = hfpsi.data.reshape(hfpsi.spwf_shape, order='F')

        self.ket = hfpsi

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

        return self.O_ket

    def compute_matrix_representation(self, recompute_action=True):
        """Compute the matrix representation of the operator wrt `self.ket`.

        Args:
            recompute_action: bool, If True, the action will be recomputed, i.e. `self.compute_action()`
                will be called. This parameter `recompute_action` must be true (default), if the `ket`
                member was modified after the last call to `compute_action()`. Otherwise, some computation
                can be saved by passing `recompute_action=False`.

        Raises:
            AttributeError: if self.O_ket is not created in one of the override methods or if they
                were not called.
        """
        if recompute_action:
            self.compute_action()

        if not hasattr(self, 'matrix'):
            self.matrix = self.ket.allocate_matrix_representation()

        for ib in range(8):
            block_ib = self.matrix[ib]
            blockstart, blockstop =  self.hfblockrange[ib][0], self.hfblockrange[ib][1]
            if blockstart < blockstop:
                ket_d3_ib   = self.  ket.d3[:,:,blockstart:blockstop]
                O_ket_d3_ib = self.O_ket.d3[:,:,blockstart:blockstop]
                np.einsum("ijk,ijl->kl", ket_d3_ib, O_ket_d3_ib, out=block_ib, order='F', optimize=True)
                block_ib *= self.mesh.dv

        return self.matrix

    def compute_diagonal_elements(self, recompute_action=True):
        """Compute only the diagonal elements of the operator wrt `self.ket`.

        Args:
            recompute_action: bool, If True, the action will be recomputed, i.e. `self.compute_action()`
                will be called. This parameter `recompute_action` must be true (default), if the `ket`
                member was modified after the last call to `compute_action()`. Otherwise, some computation
                can be saved by passing `recompute_action=False`.

        Raises:
            AttributeError: if self.O_ket is not created in one of the override methods or if they
                were not called.

        """
        if recompute_action:
            self.compute_action()

        if not hasattr(self, 'diagonal'):
            self.diagonal = np.empty(self.ket.n_total_wf, dtype=np.float64)

        np.einsum("ijk,ijk->k", self.ket.d3, self.O_ket.d3, out=self.diagonal, order='F', optimize=True)
        self.diagonal *= self.mesh.dv

        return self.diagonal

    def compute_dispersion(self, recompute_action=True):
        """Compute the dispersion and the diagonal of this operator: <bra_i|O^dagger O|ket_i> - <bra_i|O|ket_i>**2

        Args:
            recompute_action: bool, If True, the action will be recomputed, i.e. `self.compute_action()`
                will be called. This parameter `recompute_action` must be true (default), if the `ket`
                member was modified after the last call to `compute_action()`. Otherwise, some computation
                can be saved by passing `recompute_action=False`.
        Returns:
            diagonal, dispersion: both
        """
        if recompute_action:
            self.compute_action()

        if not hasattr(self, 'diagonal'):
            self.diagonal = np.empty(self.ket.n_total_wf, dtype=np.float64)

        if not hasattr(self, 'dispersion'):
            self.dispersion = np.empty(self.ket.n_total_wf, dtype=np.float64)

        # compute the diagonal:
        np.einsum("ijk,ijk->k", self.ket.d3, self.O_ket.d3, out=self.diagonal, order='F', optimize=True)
        self.diagonal *= self.mesh.dv

        # compute the dispersion:
        np.einsum("ijk,ijk->k", self.O_ket.d3, self.O_ket.d3, out=self.dispersion, order='F', optimize=True)
        self.dispersion *= self.mesh.dv
        self.dispersion -= (self.diagonal)**2

        return self.diagonal, self.dispersion
