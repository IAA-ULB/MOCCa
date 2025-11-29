import numpy as np

from .operator import Operator

# ==============================================================================
# Base classes
# ==============================================================================
class KineticEnergyOperator(Operator):
    """Implements the kinetic energy operator. Intended as a base class for
    Hamiltonian operators, derived classes must add the potential operator.
    """

    def __init__(self, mfs, hbm=20.73553000, extra_derivatives=None):
        """Initialize the operator with a wave function `hfpsi` to operate on and set the
        parameters for the kinetic energy operator.
        Derived classes must implement add_local_terms and/or add_non_local_terms.

        Args:
            mfs: mean-field state object
            hbm: prefactor of the kinetic energy operator. A single value considers the masses of neutrons
                and protons to be identical, a tuple of 2 values (hbar^2/2m_n, hbar^2/2m_p) considers
                the masses to differ.
            extra_derivatives: extra derivatives needed by the derived Operator. The Laplacian operator for the
                kinetic energy operator is automatically added.
        """
        super().__init__(mfs, derivatives='Laplacian')

        if extra_derivatives is not None:
            if isinstance(extra_derivatives, str):
                self.derivatives.append(extra_derivatives)
            elif isinstance(extra_derivatives, list):
                self.derivatives.extend(extra_derivatives)
            else:
                raise ValueError(f"{extra_derivatives=}, must be a string or a list, not {type(extra_derivatives)}.")

        self.Delta = 'xx' if self.operand_data.mesh.dim == 1 else 'Laplacian'
        # The kinetic energy term has a minus sign which we incorporate in hbm.
        if isinstance(hbm, float):
            if hbm > 0:
                hbm = -hbm
        else:
            hbm = [ -h if h > 0 else h for h in hbm]
        self.hbm = hbm


    def apply_base_operator(self):
        """apply the kinetic energy operator."""

        if self.operand_data.O_ket is None:
            self.operand_data.O_ket = self.operand_data.ket.clone()

        if isinstance(self.hbm, float):
            # Neutrons and protons are treated equally (mass)
            assert self.hbm < 0, "hbm must incorporate the minus sign in the Kinetic Energy Operator"
            self.operand_data.O_ket.data[:, :] = self.hbm * self.operand_data.ket.derivatives[self.Delta]

        else:
            # Neutrons and protons are treated differently (mass)
            hbm_n = self.hbm[0]
            hbm_p = self.hbm[1]
            assert (hbm_n < 0) and (hbm_p < 0), "hbm must incorporate the minus sign in the Kinetic Energy Operator"

            # Blocks[0:4] are for neutrons
            # Blocks[4:8] are for protons
            n = self.operand_data.O_ket.hfblockrange[3][1]  # end of neutron range in the spwfs and begin of proton range
            Delta3 = self.operand_data.ket.derivatives[self.Delta].reshape(self.operand_data.ket.spwf_shape, order='F')
            self.operand_data.O_ket.data[:, :, :n] = hbm_n * Delta3[:, :, :]
            self.operand_data.O_ket.data[:, :, n:] = hbm_p * Delta3[:, :, :]

    def __str__(self):
        """
        cfr https://github.com/IAA-nuclear/tantalus_full/issues/66
        "index" -> the wavefunction index in memory
        "shell" -> the number of particles you can fit in the levels up to and including the current one; i.e. 2 * the position of the level in the energ-ordering. It is the column "n" in the example.
        "parity" -> parity of the spwf; +1 in blocks 1,2,5,6; -1 in blocks 3,4,7,8 (fortran indices)
        "signature" -> +1 for now
        "occupation" -> the occupation of the spwfs; I'm not sure if you already construct this?
        "energy" -> the single-particle energy, or rather h_ii
        "MPIrank" -> the rank that stores this particular spwf; 0 for now.
        "Dispersion" -> the dispersion of the single-particle energy,

        Blok (0): neutrons with positive parity en signature +i
        Blok (1): neutrons with positive parity en signature -i
        Blok (2): neutrons with negative parity en signature +i
        Blok (3): neutrons with negative parity en signature -i
        Blok (4): protons  with positive parity en signature +i
        Blok (5): protons  with positive parity en signature -i
        Blok (6): protons  with negative parity en signature +i
        Blok (7): protons  with negative parity en signature -i
        """
        col_parity = np.empty(self.O_ket.data.shape[1], dtype=int)

def V_WoodsSaxon1D(r, V0, ainv, R):
    return V0 / (1. + np.exp(ainv * (r - R)))

def V_WoodsSaxon2D(x, y, V0, ainv, R):
    return V0 / (1. + np.exp(ainv * (np.sqrt(x * x + y * y) - R)))

def V_WoodsSaxon3D(x, y, z, V0, ainv, R):
    return V0 / (1. + np.exp(ainv * (np.sqrt(x * x + y * y + z * z) - R)))


#=======================================================================================================================
# Derived classes
#=======================================================================================================================
class HamiltonianWoodsSaxon(KineticEnergyOperator):
    """A hamiltonian with a Woods-Saxon potential. <hfpsi|h|hfpsi>."""

    def __init__(self, mfs, hbm=20.73553000, V0=-50.0, r0=1.25, a=0.5):
        """Initialize the operator with a wave function `hfpsi` to operate on and set the parameters for the
        operator.

        Args:
            hbm: prefactor of the kinetic energy operator. A single value considers the masses of neutrons
                and protons to be identical, a tuple of 2 values (hbar^2/2m_n, hbar^2/2m_p) considers
                the masses to differ.
            V0, r0, a:  parameters of the Woods-Saxon potential
                (cfr https://en.wikipedia.org/wiki/Woods–Saxon_potential).
        """
        super().__init__(mfs, hbm)
        self.V0 = V0 if (V0 < 0) else -V0 # incorporate the minus sign of the potential
        self.ainv = 1 / a
        # Note that this line makes the hamiltionian implicitly dependent on the system
        self.R = r0 * np.pow(mfs.hfpsi.n_neutrons + mfs.hfpsi.n_protons, 1 / 3)
        self.ws_parms = {
            'V0'  : self.V0,
            'ainv': self.ainv,
            'R'   : self.R,
        }

    def add_local_terms(self):
        """Create self.O_ket, fill it with the Woods-Saxon potential, and add the kinetic energy."""

        # Select Woods-Saxon potential functions
        V_WoodsSaxon = V_WoodsSaxon3D if self.operand_data.mesh.dim == 3 else \
            V_WoodsSaxon2D if self.operand_data.mesh.dim == 2 else \
                V_WoodsSaxon1D

        # Apply the Woods-Saxon potential to the mesh points
        Vr = self.operand_data.mesh.apply(V_WoodsSaxon, self.ws_parms)

        # add to O_ket
        self.operand_data.O_ket.data += Vr * self.operand_data.ket.data

