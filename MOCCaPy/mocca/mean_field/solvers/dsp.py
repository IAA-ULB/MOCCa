# The Diagonalization subproblem
# see Ryssens et al (2019) Eur. Phys. J. A (2019) 55:93 section 3.

import numpy as np
import pytest

from mocca.mean_field.states import gramm_schmidt
from mocca.util.timer import Timer


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

    def __repr__(self):
        return f"DSP(H={self.hamiltonian}, alpha={self.alpha}, mu={self.mu})" if (self.mu > 0) else \
               f"DSP(H={self.hamiltonian}, alpha={self.alpha})"

    # @Timer("DSP.step()")
    def step(self, check=False, verbose=False):
        """Apply a single step (gradient descent or heavy ball dynamics, iff self.mu>0).

        1. Compute the diagonal elements of the hamiltonian
        2. Update the single particle wave functions (HFPsi)
        3. Orthogonalize and normalize the single particle wave functions (HFPsi)
        4. Invalidate the operand data of the hamiltonian operator.
        """

        # 1. Compute the diagonal elements of the hamiltonian
        epsilon = self.hamiltonian.compute_diagonal_elements()
        if verbose:
            print(f"iter = {self._counter}: h_ii = {epsilon}")

        # 2. Update the single particle wave functions (HFPsi)
        #    This requires ket, epsilon and h_ket (both of which have been computed
        #    in the call to self.hamiltonian.compute_diagonal_elements())
        ket_d3   = self.hamiltonian.operand_data.  ket.d3
        h_ket_d3 = self.hamiltonian.operand_data.O_ket.d3
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
        gramm_schmidt(self.hamiltonian.operand_data.ket, self.hamiltonian.operand_data.diagonal, normalize=True, check=check)

        # 4. Invalidate the operand data of the hamiltonian operator.
        self.hamiltonian.invalidate()

    def evolve(self, nsteps=100, residual=1e-8, check=False, verbose=False):
        """
        Args:
            nsteps (int): maximum number of steps to evolve
            residual: convergence is assumed when the average dispersions of occupied
                levels becomes <= residual.
            check (bool): whether to check ortohonality in gramm_schmidt
            verbose:
        """
        if not hasattr(self, '_step'):
            self._counter = 0

        for i in range(nsteps):
            self.step(check=check, verbose=verbose)
            self._counter +=1
            if self._counter % 10 == 0:
                h_ii, d2_ii = self.hamiltonian.compute_dispersion()
                rho_ii, _, _ = self.hamiltonian.operand_data.mfs.occupancies(h_ii=h_ii)
                res = np.einsum('i,i', d2_ii, rho_ii) / (0.5*(self.hamiltonian.operand_data.mfs.n_neutrons + self.hamiltonian.operand_data.mfs.n_protons))
                #   this is basically the average of the dispersions of occupied energy levels (rho_ii == 1.0)
                converged = res <= residual
                print(f"iter={self._counter}: res = {res:.2e} {'<=' if converged else '>'} {residual}{' : CONVERGED' if converged else ''}")
                if converged:
                    print(self.hamiltonian)
                    break


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
        ket_d3[:, :, i] *= f[i] # *= (1 + alph*epsilon_ii)

    h_ket_d3 *= alpha
    # if we do not update ket, nothing should change. That turned out to be correct, so the hamiltionian's
    # ket member is not accidently changed, and its O_ket member is correctly recomputed.
    ket_d3 -= h_ket_d3

