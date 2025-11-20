import numpy as np

from .operator import Operator

# ==============================================================================
# Base classes
# ==============================================================================
class KineticEnergyOperator(Operator):
    """Implements the kinetic energy operator. Intended as a base class for
    Hamiltonian operators, derived classes must add the potential operator.
    """

    def __init__(self, hfpsi, hbm=20.73553000, extra_derivatives=None):
        """Initialize the operator with a wave function `hfpsi` to operate on and set the
        parameters for the kinetic energy operator.
        Derived classes must implement add_local_terms and/or add_non_local_terms.

        Args:
            hbm: prefactor of the kinetic energy operator. A single value considers the masses of neutrons
                and protons to be identical, a tuple of 2 values (hbar^2/2m_n, hbar^2/2m_p) considers
                the masses to differ.
            extra_derivatives: extra derivatives needed by the derived Operator. The Laplacian operator for the
                kinetic energy operator is automatically added.
        """
        super().__init__(hfpsi, derivatives='Laplacian')

        if extra_derivatives is not None:
            if isinstance(extra_derivatives, str):
                self.derivatives.append(extra_derivatives)
            elif isinstance(extra_derivatives, list):
                self.derivatives.extend(extra_derivatives)
            else:
                raise ValueError(f"{extra_derivatives=}, must be a string or a list, not {type(extra_derivatives)}.")

        self.nabla = 'xx' if self.mesh.dim == 1 else 'Laplacian'
        self.hbm = hbm

    def apply_base_operator(self):
        """apply the kinetic energy operator."""
        self.compute_derivatives()

        self.O_ket = self.ket.clone()

        if isinstance(self.hbm, float):
            # Neutrons and protons are treated equally (mass)
            self.O_ket.data[:, :] = self.hbm * self.ket.derivatives[self.nabla]

        else:
            # Neutrons and protons are treated differently (mass)
            hbm_n = self.hbm[0]
            hbm_p = self.hbm[1]
            # Blocks[0:4] are for neutrons
            # Blocks[4:8] are for protons
            n = self.O_ket.hfblockrange[3][1]  # end of neutron range in the spwfs and begin of proton range
            self.O_ket.data[:, :, :n] = hbm_n * self.ket.derivatives[nabla][:, :, :n]
            self.O_ket.data[:, :, n:] = hbm_p * self.ket.derivatives[nabla][:, :, n:]


# ==============================================================================
# Derived classes
# ==============================================================================
class HamiltonianWoodsSaxon(Operator):
    """Compute the matrix representation of the hamiltonian h with a Woods-Saxon potential. <hfpsi|h|hfpsi>."""

    def __init__(self, hfpsi, hbm=20.73553000, V0=-50.0, r0=1.25, a=0.5):
        """Initialize the operator with a wave function `hfpsi` to operate on and set the parameters for the
        operator.

        Args:
            hbm: prefactor of the kinetic energy operator. A single value considers the masses of neutrons
                and protons to be identical, a tuple of 2 values (hbar^2/2m_n, hbar^2/2m_p) considers
                the masses to differ.
            V0, r0, a:  parameters of the Woods-Saxon potential
                (cfr https://en.wikipedia.org/wiki/Woods–Saxon_potential).
        """
        super().__init__(hfpsi)
        self.hbm = hbm
        self.V0 = V0
        self.ainv = 1 / a
        self.R = r0 * np.pow(hfpsi.n_neutrons + hfpsi.n_protons, 1 / 3)
        self.derivatives = ['Laplacian' if self.mesh.dim >= 2 else 'xx']

    def add_local_terms(self):
        """Create self.O_ket, fill it with the Woods-Saxon potential, and add the kinetic energy."""
        self.O_ket = self.ket.clone()

        # TODO: speed up with numba decorators?

        def V_WoodsSaxon1D(r):
            return self.V0 / (1. + np.exp(self.ainv * (r - self.R)))

        def V_WoodsSaxon2D(x, y):
            return self.V0 / (1. + np.exp(self.ainv * (np.sqrt(x * x + y * y) - self.R)))

        def V_WoodsSaxon3D(x, y, z):
            return self.V0 / (1. + np.exp(self.ainv * (np.sqrt(x * x + y * y + z * z) - self.R)))

        dim = self.mesh.dim
        if dim == 3:
            V_WoodsSaxon = V_WoodsSaxon3D
        elif dim == 2:
            V_WoodsSaxon = V_WoodsSaxon2D
        else:  # dim == 1:
            V_WoodsSaxon = V_WoodsSaxon1D

        # aply the Woods-Saxon potential
        Vr = self.mesh.apply(V_WoodsSaxon)
        self.O_ket.data = Vr * self.ket.data

        nabla = 'Laplacian' if self.mesh.dim >= 2 else \
            'xx'
        if isinstance(self.hbm, float):
            # Neutrons and protons are treated equally (mass)
            self.O_ket.data += self.hbm * self.ket.derivatives[nabla]
        else:
            # Neutrons and protons are treated differently (mass)
            hbm_n = self.hbm[0]
            hbm_p = self.hbm[1]
            # Blocks[0:4] are for neutrons
            # Blocks[4:8] are for protons
            n = self.O_ket.hfblockrange[3][1]  # end of neutron range in the spwfs and begin of proton range
            self.O_ket.data[:, :, :n] += hbm_n * self.ket.derivatives[nabla][:, :, :n]
            self.O_ket.data[:, :, n:] += hbm_p * self.ket.derivatives[nabla][:, :, n:]


class HamiltonianWoodsSaxon2(KineticEnergyOperator):
    """Compute the matrix representation of the hamiltonian h with a Woods-Saxon potential. <hfpsi|h|hfpsi>.

    Does the same as the class above, but derives from KineticEnergyOperator.
    """

    def __init__(self, hfpsi, hbm=20.73553000, V0=-50.0, r0=1.25, a=0.5):
        """Initialize the operator with a wave function `hfpsi` to operate on and set the parameters for the
        operator.

        Args:
            hbm: prefactor of the kinetic energy operator. A single value considers the masses of neutrons
                and protons to be identical, a tuple of 2 values (hbar^2/2m_n, hbar^2/2m_p) considers
                the masses to differ.
            V0, r0, a:  parameters of the Woods-Saxon potential
                (cfr https://en.wikipedia.org/wiki/Woods–Saxon_potential).
        """
        super().__init__(hfpsi, hbm)
        self.V0 = V0
        self.ainv = 1 / a
        self.R = r0 * np.pow(hfpsi.n_neutrons + hfpsi.n_protons, 1 / 3)

    def add_local_terms(self):
        """Create self.O_ket, fill it with the Woods-Saxon potential, and add the kinetic energy."""

        # TODO: speed up with numba decorators?
        # Define Woods-Saxon potential functions
        def V_WoodsSaxon1D(r):
            return self.V0 / (1. + np.exp(self.ainv * (r - self.R)))

        def V_WoodsSaxon2D(x, y):
            return self.V0 / (1. + np.exp(self.ainv * (np.sqrt(x * x + y * y) - self.R)))

        def V_WoodsSaxon3D(x, y, z):
            return self.V0 / (1. + np.exp(self.ainv * (np.sqrt(x * x + y * y + z * z) - self.R)))

        # Select Woods-Saxon potential functions
        V_WoodsSaxon = V_WoodsSaxon3D if self.mesh.dim == 3 else \
            V_WoodsSaxon2D if self.mesh.dim == 2 else \
                V_WoodsSaxon1D

        # Apply the Woods-Saxon potential to the mesh points
        Vr = self.mesh.apply(V_WoodsSaxon)

        # add to O_ket
        self.O_ket.data += Vr * self.O_ket.data

