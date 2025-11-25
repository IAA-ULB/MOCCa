# The Diagonalization subproblem
# see Ryssens et al (2019) Eur. Phys. J. A (2019) 55:93 section 3.

import numpy as np
import pytest

from mocca.mean_field.operators.overlap import (Overlap)
from mocca.util.timer import Timer


def _gradient_descent_step(ket_d3, h_ket_d3, epsilon, alpha):
    """Take a single gradient descent step (eq 63 ine see Ryssens et al (2019)
    Eur. Phys. J. A (2019) 55:93 section 3.

    Args:
        ket_data: HFPsi.data member of ket.
        h_ket_data: HFPsi.data member of h_ket containing the hamiltonian applied
            to ket.
        epsilon: contains the diagonal elements of the hamiltonian wrt ket:
            diag(<ket|h|ket>).
    """
    # DO NOT overwrite epsilon! it is needed by gramm_schmidt.
    f = epsilon * alpha + 1.0
    ket_d3_copy = ket_d3.copy()
    for i in range(epsilon.size):
        ket_d3[:, :, i] *= f[i] # *= (1 + alph*epsilon+ii)

    h_ket_d3 *= alpha
    # if we do not update ket, nothing should change. That turned out to be correct, so the hamiltionian's
    # ket member is not accidently changed, and its O_ket member is correctly recomputed.
    ket_d3 -= h_ket_d3
    pass

class DSP:
    """A class for evolving the single particle wave functions as part of the
    Diagonalization SubProblem, either by gradient descent step or heavy ball
    dynamics stepping. The aim is to diagonalize the Hamiltonian wrt HFPsi.
    """
    def __init__(self, hamiltonian, alpha, mu=0.0):
        """Constructor for the DSP class

        Args:
            hamiltonian: operator representing the Hamiltonian to be diagonalized.
                It contains a HFPsi object.
            alpha: scalar parameter for step size of gradient descent. (see Ryssens et
                al. (2019) Eur. Phys. J. A (2019) 55:93 section 3.)
            mu: scalar parameter for momentum step size of heavy ball dynamics. If mu
                is zero, a gradient descent step is performed. (see Ryssens et al (2019)
                Eur. Phys. J. A (2019) 55:93 section 3.)
        """
        if not (0 < alpha < 1):
            raise ValueError(f"`alpha` parameter must be in ]0,1[ (got {alpha=}).")
        self.alpha = alpha

        if not (0 <= mu < 1):
            raise ValueError(f"`mu` parameter must be in [0,1[ (got {mu=}).")
        self.mu = mu

        self.hamiltonian = hamiltonian

    # @Timer("DSP.step()")
    def step(self, check=False, debug=False):
        """Apply a single step (gradient descent or heavy ball dynamics, iff self.mu>0).

        1. Compute the diagonal elements of the hamiltonian
        2. Update the single particle wave functions (HFPsi)
        3. Orthogonalize and normalize the single particle wave functions (HFPsi)
        """

        # 1. Compute the diagonal elements of the hamiltonian
        if debug or self._counter == 0:
            epsilon = self.hamiltonian.compute_diagonal_elements()
            print(f"\niter = {self._counter}: h_ii = {self.hamiltonian.diagonal}")
            # print(  f"iter = {self._counter}: d_ii = {self.hamiltonian.dispersion}")
            # epsilon, dispersion = self.hamiltonian.compute_dispersion()
            # print(f"\niter = {self._counter}: h_ii = {self.hamiltonian.diagonal}")
            # print(  f"iter = {self._counter}: d_ii = {self.hamiltonian.dispersion}")
        else:
            epsilon = self.hamiltonian.compute_diagonal_elements()

        # 2. Update the single particle wave functions (HFPsi)
        #    This requires ket, epsilon and h_ket (both of which have been computed
        #    in the call to self.hamiltonian.compute_diagonal_elements())
        ket_d3   = self.hamiltonian.  ket.d3
        h_ket_d3 = self.hamiltonian.O_ket.d3
        if self.mu == 0:
            _gradient_descent_step(ket_d3, h_ket_d3, epsilon, self.alpha)
        else:
            # heavy ball dynamics step
            if hasattr(self, '_mu_ket_nextprev_data'):
                self._mu_ket_nextprev_d3[:,:] = ket_d3
                self._mu_ket_nextprev_d3[:,:] *= self.mu

                epsilon *= -self.alpha
                epsilon += (1 + self.mu)

                for i in range(epsilon.size):
                    ket_d3[:, :, i] *= epsilon[i]

                h_ket_d3 *= self.alpha
                h_ket_d3 += self._mu_ket_prev_d3
                ket_d3 -= h_ket_d3

                # now we can overwrite _mu_ket_prev_d3
                self._mu_ket_prev_d3[:,:] = self._mu_ket_nextprev_d3
            else:
                self._mu_ket_nextprev_d3 = np.empty_like(ket_d3)
                self._mu_ket_prev_d3     = np.empty_like(ket_d3)
                self._mu_ket_prev_d3[:,:] = ket_d3
                self._mu_ket_prev_d3[:,:] *= self.mu
                _gradient_descent_step(ket_d3, h_ket_d3, epsilon, self.alpha)

        # 3. Orthogonalize and normalize the single particle wave functions (HFPsi)
        gramm_schmidt(self.hamiltonian.ket, self.hamiltonian.diagonal, normalize=True, check=check)

    def evolve(self, nsteps=1, check=True, debug=False):
        """"""
        if not hasattr(self, '_step'):
            self._counter = 0

        for i in range(nsteps):
            self.step(check=check, debug=debug)
            self._counter +=1


        self.hamiltonian.compute_dispersion()
        print(f"iter = {self._counter}: h_ii = {self.hamiltonian.diagonal}")
        print(f"iter = {self._counter}: d_ii = {self.hamiltonian.dispersion}")

# @Timer("_projector()")
def _projector(hfpsi_d3_ib, i, j):
    """Used by gramm_schmidt."""
    Overlap_ij = np.einsum("hk,hk", hfpsi_d3_ib[:, :, i], hfpsi_d3_ib[:, :, j], order='F', optimize=True)
    Overlap_jj = np.einsum("hk,hk", hfpsi_d3_ib[:, :, j], hfpsi_d3_ib[:, :, j], order='F', optimize=True)
    # both are missing a factor mesh.dv but that doesn's matter because of the quotient
    return Overlap_ij / Overlap_jj

# @Timer('gramm_schmidt()')
def gramm_schmidt(hfpsi, order=None, normalize=True, check=False):
    """Gramm-Schmidt orthogonalization of a wave function.

    Args:
        hfpsi: wave function to orthogonalize
        order: ndarray of length hfpsi.n_total_wf. This array determines the selection order of the
            single particle wave functions uin the orthogonalization process. The lowest spwf is
            selected first. If None, the order is the order of the single particle wave functions
            in hfpsi.data.
        normalize: whether to normalize the single particle wave functions after orthogonalization.
        check: if True, verifies that the overlap matrix is unity. (For debugging purposes only).
    """
    hfpsi_d3 = hfpsi.d3
    for ib in range(8):
        range_ib = hfpsi.hfblockrange[ib]
        if range_ib[1] > range_ib[0]:
            # block is not empty
            # Restrict all data structures to symmetry block ib
            hfpsi_d3_ib = hfpsi_d3[:,:,range_ib[0]:range_ib[1]]
            if not order is None:
                order_ib = np.argsort(order[range_ib[0]:range_ib[1]])
            n = range_ib[1]-range_ib[0]
            for I in range(1,n):
                i = order_ib[I] if (order is not None) else I
                for J in range(I):
                    j = order_ib[J] if (order is not None) else J
                    hfpsi_d3_ib[:,:,i] -= _projector(hfpsi_d3_ib,i,j) * hfpsi_d3_ib[:,:,j]
                # for J in range(0,I):
                #     j = order_ib[J]
                #     Oij = np.einsum("hk,hk", hfpsi_d33_ib[:,:,i], hfpsi_d33_ib[:,:,j])
                #     print(f"{ib=} ({i},{j}) {Oij=}")

    if normalize:
        hfpsi.normalize()

    if check:
        overlap = Overlap(hfpsi)
        overlap.compute_matrix_representation()
        for ib in range(8):
            n_ib = hfpsi.hfblocks[ib]
            block_ib = overlap.matrix[ib]
            for i in range(0, n_ib):
                for j in range(i, n_ib):
                    if i == j:
                        if normalize:
                            assert block_ib[i,i] == pytest.approx(1.0), f"block[{ib}][{i},{j}] {block_ib[i, j]}"
                    else:
                        assert block_ib[i,j] == pytest.approx(0.0), f"block[{ib}][{i},{j}] {block_ib[i, j]}"