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
            define the action of the operator to ket resulting in a MeshQuantity O_ket. In principle,
            there is only need for one method to be overridden, one that implements the action of the
            operator on the ket. However, the names of the three methods serve as a hint to the things
            that must be computed.
        compute_matrix_representation: compute <bra|O_ket >. This must not be overridden
    """
    def __init__(self, mfs, derivatives=None):
        """Initialize the operator with a mean-field state object. If you need to apply the operator to
        several mean-field state objects, you must instantiate one Operator object for each mean-field
        state. The Operator objects themselves consumes little memory.

        Args:
            mfs: mean-field state object to which the operator will be applied.
            derivatives: str|list of derivatives needed by the operator. See the `axes` parameter
                of MeshQuantity.differentiate() for details
        """
        # if not hasattr(hfpsi, 'd3'):
        #     hfpsi.d3 = hfpsi.data.reshape(hfpsi.spwf_shape, order='F')

        self.operand_data = OperandData(mfs)

        self.derivatives = [] if derivatives is None else \
                           [derivatives] if isinstance(derivatives, str) else \
                            derivatives

        if self.operand_data.mesh.dim == 1:
            # replace 'laplacian' with 'xx'
            self.derivatives = [ 'xx' if axes == 'Laplacian' else axes for axes in self.derivatives]

    def __repr__(self):
        return f"<Operator[{self.operand_data.mfs}]>"

    # --------------------------------------------------------------------------
    # Methods to be overridden by derived classes (at least one of them)
    # --------------------------------------------------------------------------
    def apply_base_operator(self):
        """In case the derived class is to be used as a base class itself.
        See KineticEnergyOperator for an example.
        """
        raise NotImplementedError(f"Method apply_base_operator() must be overridden.")

    def add_local_terms(self):
        """Override this to compute local terms of the operator."""
        pass

    def add_non_local_terms(self):
        """Override this to compute nonlocal terms of the operator."""
        pass

    # --------------------------------------------------------------------------
    # Methods NOT to be overridden by derived classes
    # --------------------------------------------------------------------------
    def invalidate(self):
        """
        Args:
            mfs: mean-field state object
        """
        self.operand_data.invalidate()

    def compute_derivatives(self):
        """Compute all derivatives needed by the operator."""
        self.operand_data.ket.differentiate(self.derivatives)

    def compute_action(self):
        """Compute the action of the operator on the wave function.

        Returns:
            operand_data.O_ket: O|ket>
        """
        if not self.operand_data.actn_uptodate:
            self.compute_derivatives()
            self.apply_base_operator()
            self.add_local_terms    ()
            self.add_non_local_terms()

            self.operand_data.actn_uptodate = True
            self.operand_data.mtrx_uptodate = False
            self.operand_data.diag_uptodate = False
            self.operand_data.disp_uptodate = False

        return self.operand_data.O_ket

    def compute_matrix_representation(self):
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
        if not self.operand_data.actn_uptodate:
            self.compute_action()

        if not self.operand_data.mtrx_uptodate:
            if self.operand_data.matrix is None:
                self.operand_data.matrix = self.operand_data.ket.allocate_matrix_representation()

            for ib in range(8):
                block_ib = self.operand_data.matrix[ib]
                blockstart, blockstop =  self.operand_data.hfblockrange[ib][0], self.operand_data.hfblockrange[ib][1]
                if blockstart < blockstop:
                    ket_d3_ib   = self.operand_data.  ket.d3[:,:,blockstart:blockstop]
                    O_ket_d3_ib = self.operand_data.O_ket.d3[:,:,blockstart:blockstop]
                    np.einsum("ijk,ijl->kl", ket_d3_ib, O_ket_d3_ib, out=block_ib, order='F', optimize=True)
                    block_ib *= self.operand_data.mesh.dv

            self.operand_data.mtrx_uptodate = True

        return self.operand_data.matrix

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
        if not self.operand_data.actn_uptodate:
            self.compute_action()

        if not self.operand_data.diag_uptodate:
            if self.operand_data.diagonal is None:
                self.operand_data.diagonal = np.empty(self.operand_data.ket.n_total_wf, dtype=np.float64)

            np.einsum("ijk,ijk->k", self.operand_data.ket.d3, self.operand_data.O_ket.d3, out=self.operand_data.diagonal, order='F', optimize=True)
            self.operand_data.diagonal *= self.operand_data.mesh.dv
            self.operand_data.diag_uptodate = True

        return self.operand_data.diagonal

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
        if not self.operand_data.actn_uptodate:
            self.compute_action()

        if not self.operand_data.diag_uptodate:
            self.compute_diagonal_elements()

        if not self.operand_data.disp_uptodate:
            if self.operand_data.dispersion is None:
                self.operand_data.dispersion = np.empty(self.operand_data.ket.n_total_wf, dtype=np.float64)

            # compute the dispersion:
            np.einsum("ijk,ijk->k", self.operand_data.O_ket.d3, self.operand_data.O_ket.d3, out=self.operand_data.dispersion, order='F', optimize=True)
            self.operand_data.dispersion *= self.operand_data.mesh.dv
            self.operand_data.dispersion -= (self.operand_data.diagonal)**2
            self.operand_data.disp_uptodate = True

        return self.operand_data.diagonal, self.operand_data.dispersion


class OperandData:
    def __init__(self, mfs):
        """
        Args:
            mfs: mean field state
        """
        self.mfs = mfs
        self.mesh = mfs.hfpsi.mesh
        self.hfblocks = mfs.hfpsi.hfblocks
        self.hfblockrange = mfs.hfpsi.hfblockrange
        self.ket = mfs.hfpsi
        self.O_ket = None  # where the action of the operator is stored
        self.diagonal = None            # where the diagonal elements of the operator are stored
        self.dispersion = None          # where the dispersion elements of the operator are stored
        self.matrix = None              # where the matrix representation of the operator are stored
        self.actn_uptodate = False
        self.diag_uptodate = False
        self.disp_uptodate = False
        self.mtrx_uptodate = False

    def invalidate(self):
        """Mark all intermediate and end results as not uptodate."""
        self.actn_uptodate = False
