# The Diagonalisation subproblem
# see Ryssens et al (2019) Eur. Phys. J. A (2019) 55:93 section 3.

import numpy as np



class DSP:
    """"""
    def __init__(self, hamiltonian, alpha, mu=0.0):
        """Apply a gradient descent step or a heavy ball dynamics step to
        hamiltonian.

        Args:
            hamiltonian: operator representing the Hamiltonian.
            alpha: scalar parameter for step size of gradient descent.
            mu : scalar parameter for step size of heavy ball dynamics. If mu
                is zero, a gradient descent step is performed.
        """
        self.alpha = alpha
        self.mu = mu
        self.hamiltonian = hamiltonian

    def step(self):
        """Apply one gradient descent step or a heavy ball dynamics step."""

        def gradient_descent_step(self, ket_data, h_ket_data, epsilon):
            # assumes that epsilon is already compputed in self.epsilon
            epsilon *= -self.alpha
            epsilon += 1.
            ket_data *= epsilon

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
                ket_data *= epsilon

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
