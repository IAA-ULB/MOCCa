# The Diagonalisation subproblem
# see Ryssens et al (2019) Eur. Phys. J. A (2019) 55:93 section 3.

import numpy as np

from mocca.mean_field.states import Overlap

class DSP:
    """"""
    def __init__(self, hamiltonian, alpha, mu=0.0):
        """A class for evolving the spwfs as part of the Diagonalisation SubProblem,
         either gradient descent step or heavy ball dynamics stepping.

        Args:
            hamiltonian: operator representing the Hamiltonian to be diagonalised.
            alpha: scalar parameter for step size of gradient descent.
            mu : scalar parameter for momentum step size of heavy ball dynamics. If mu
                is zero, a gradient descent step is performed.
        """
        self.alpha = alpha
        assert 0 < alpha < 1
        self.mu = mu
        assert 0 < mu < 1
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


class GrammSchmidt:
    """Class for Gramm-Schmidt orthogonalisaton of HFPsi objects"""
    def __init__(self, hfpsi, epsilon):
        """
        Args:
            hfpsi: wave function
            epsilon: di
        """
        self.hfpsi = hfpsi
        self.epsilon = epsilon

    def ortogonalize(self):
        """

        """
        overlap = Overlap(self.hfpsi)
        self.overlap_matrix = overlap.compute_matrix_representation()
        hfpsi_data3 = hfpsi.data.reshape((hfpsi.data.shape[0], 4, hfpsi.data.shape[1]//4), order='F')
        for ib in range(8):
            epsilon_ib = epsilon[*self.hfpsi.hfblockrange[ib]]
            if epsilon_ib.size > 0:
                order = np.argsort(epsilon_ib) + self.hfpsi.hfblockrange[ib][0]
                for i in order[1:]:
                    for j in order[0:i]
                        hfpsi_data3[:,:,i] -= (self.overlap_matrix[i,j] / self.overlap_matrix[j,j]) *hfpsi_data3[:,:,j]

        self.normalize()

    def normalize(self):
        """

        """
        overlap = Overlap(self.hfpsi)
        self.overlap_matrix = overlap.compute_matrix_representation()
        hfpsi_data3 = hfpsi.data.reshape((hfpsi.data.shape[0], 4, hfpsi.data.shape[1]//4), order='F')
        for i_spwf in range(self.hfpsi.n_total_wf):
            hfpsi_data3[:,:,i_spwf] /= self.overlap_matrix[i_spwf,i_spwf]

