# The Diagonalisation subproblem
# see Ryssens et al (2019) Eur. Phys. J. A (2019) 55:93 section 3.

import numpy as np

from mocca.mean_field.operators import Overlap

class DSP:
    """A class for evolving the spwfs as part of the Diagonalisation SubProblem,
    either gradient descent step or heavy ball dynamics stepping.
    """
    def __init__(self, hamiltonian, alpha, mu=0.0):
        """A class for evolving the spwfs as part of the Diagonalisation SubProblem,
         either gradient descent step or heavy ball dynamics stepping.

        Args:
            hamiltonian: operator representing the Hamiltonian to be diagonalised.
            alpha: scalar parameter for step size of gradient descent.
            mu: scalar parameter for momentum step size of heavy ball dynamics. If mu
                is zero, a gradient descent step is performed.
        """
        self.alpha = alpha
        assert 0 < alpha < 1
        self.mu = mu
        assert 0 <= mu < 1
        self.hamiltonian = hamiltonian

    def step(self):
        """Apply one step (gradient descent or heavy ball dynamics, iff self.mu>0)."""

        def gradient_descent_step(ket_data, h_ket_data, epsilon):
            # assumes that epsilon is already compputed in self.epsilon
            epsilon *= -self.alpha
            epsilon += 1.
            kd3 = ket_data.reshape((ket_data.shape[0], 4, ket_data.shape[1]//4), order='F')
            for i in range(epsilon.size):
                kd3[:,:,i] *= epsilon[i]

            h_ket_data *= self.alpha
            ket_data -= h_ket_data

        epsilon = self.hamiltonian.compute_diagonal_elements()
        ket_data   = self.hamiltonian.  ket.data
        h_ket_data = self.hamiltonian.O_ket.data

        if self.mu == 0:
            gradient_descent_step(ket_data, h_ket_data, epsilon)
        else:
            # heavy ball dynamics step
            if hasattr(self, 'mu_ket_nextprev_data'):
                self.mu_ket_nextprev_data[:,:] = ket_data
                self.mu_ket_nextprev_data[:,:] *= self.mu

                epsilon *= -self.alpha
                epsilon += (1 + self.mu)
                kd3 = ket_data.reshape((ket_data.shape[0], 4, ket_data.shape[1] // 4), order='F')
                for i in range(epsilon.size):
                    kd3[:, :, i] *= epsilon[i]

                h_ket_data *= self.alpha
                h_ket_data += self.mu_ket_prev_data
                ket_data -= h_ket_data

                # now we can overwrite mu_ket_prev_data
                self.mu_ket_prev_data[:,:] = self.mu_ket_nextprev_data
            else:
                self.mu_ket_nextprev_data = np.empty_like(ket_data)
                self.mu_ket_prev_data     = np.empty_like(ket_data)
                self.mu_ket_prev_data[:,:] = ket_data
                self.mu_ket_prev_data[:,:] *= self.mu
                gradient_descent_step(ket_data, h_ket_data, epsilon)

    def evolve(self, nsteps=None):
        """"""
        gs = GrammSchmidt(self.hamiltonian.ket)
        for i in range(nsteps):
            self.step(n)
            self.ket.normalize()


def gramm_schmidt(hfpsi, epsilon=None, normalize=True):
    """Gramm-Schmidt orthogonalisaton of a wave function.

    Args:
        hfpsi: wave function to orthogonalize
        epsilon: diagonal elements of the hamiltonian. determines the selection of the
            next spwf to orthogonalize (lowest first).
        normalize: whether to normalize the hfpsi.
    """
    def projector(hfpsi_data3_ib, i, j):
        Oij = np.einsum("hk,hk", hfpsi_data3_ib[:,:,i], hfpsi_data3_ib[:,:,j])
        Ojj = np.einsum("hk,hk", hfpsi_data3_ib[:,:,j], hfpsi_data3_ib[:,:,j])
        # both are missing a factor mesh.dv but that doesn's matter because of the quotient
        return Oij / Ojj

    # overlap = Overlap(self.hfpsi)
    hfpsi_data3 = hfpsi.data.reshape((hfpsi.data.shape[0], 4, hfpsi.data.shape[1]//4), order='F')
    for ib in range(8):
        range_ib = hfpsi.hfblockrange[ib]
        if range_ib[1] > range_ib[0]:
            # block is not empty
            # Restrict all data structures to symmetry block ib
            hfpsi_data3_ib = hfpsi_data3[:,:,range_ib[0]:range_ib[1]]
            epsilon_ib     = epsilon    [    range_ib[0]:range_ib[1]]
            if not epsilon is None:
                hfpsi_data3_ib[:,:,:] *= epsilon_ib
            order_ib = np.argsort(epsilon_ib)
            n = range_ib[1]-range_ib[0]
            for I in range(1,n):
                i = order_ib[I] if (epsilon is not None) else I
                for J in range(I):
                    j = order_ib[J] if (epsilon is not None) else J
                    hfpsi_data3_ib[:,:,i] -= projector(hfpsi_data3_ib,i,j) * hfpsi_data3_ib[:,:,j]
                # for J in range(0,I):
                #     j = order_ib[J]
                #     Oij = np.einsum("hk,hk", hfpsi_data3_ib[:,:,i], hfpsi_data3_ib[:,:,j])
                #     print(f"{ib=} ({i},{j}) {Oij=}")

    if normalize:
        hfpsi.normalize()

