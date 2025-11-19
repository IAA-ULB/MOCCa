import numpy as np

from .hfpsi import HFPsi

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
        self.matrix = hfpsi.allocate_matrix_representation()
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

    def compute_matrix_representation(self):
        """Compute the matrix representation of the operator wrt `self.ket`.

        Raises:
            AttributeError: if self.O_ket is not created in one of the override methods.

        """
        self.compute_derivatives()
        self.apply_base_operator()
        self.add_local_terms()
        self.add_non_local_terms()

        try:
            O_ket = self.O_ket
        except AttributeError:
            raise AttributeError(f"Attribute 'self.O_ket' is missing in class {self.__class__.__name__}. \n"
                                 f"One of the derived methods compute_derivatives, add_local_terms, "
                                 f"add_non_local_terms must define `self.O_ket`, where the action of "
                                 f"the operator on `self.ket` is stored.")
        bra = self.bra
        for ib in range(8):
            block = self.matrix[ib]
            for j in range(*self.bra.hfblockrange[ib]):
                for i in range(*self.ket.hfblockrange[ib]):
                    # TODO: use out= ? for performance
                    block[i-self.ket.hfblockrange[ib][0],j-self.bra.hfblockrange[ib][0]] = (
                        bra.data[:,   4*j] * O_ket.data[:,   4*i] +
                        bra.data[:, 1+4*j] * O_ket.data[:, 1+4*i] +
                        bra.data[:, 2+4*j] * O_ket.data[:, 2+4*i] +
                        bra.data[:, 3+4*j] * O_ket.data[:, 3+4*i]
                    ).sum(0)
        return self.matrix


class KineticEnergyOperator(Operator):
    """Implements the kinetic energy operator. Intended as a base class for Hamiltonian operators, 
    derived classes must add the potential operator."""
    
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

        self.O_ket = HFPsi.like(self.ket)

        if isinstance(self.hbm, float):
            # Neutrons and protons are treated equally (mass)
            self.O_ket.data[:,:] = self.hbm * self.ket.derivatives[self.nabla]

        else:
            # Neutrons and protons are treated differently (mass)
            hbm_n = self.hbm[0]
            hbm_p = self.hbm[1]
            # Blocks[0:4] are for neutrons
            # Blocks[4:8] are for protons
            n = self.O_ket.hfblockrange[3][1] # end of neutron range in the spwfs and begin of proton range
            self.O_ket.data[:,:,:n] = hbm_n * self.ket.derivatives[nabla][:,:,:n]
            self.O_ket.data[:,:,n:] = hbm_p * self.ket.derivatives[nabla][:,:,n:]


# ==============================================================================
# Derived classes
# ==============================================================================
class Overlap(Operator):
    """Compute the matrix representation of the overlap operator. <hfpsi|hfpsi>."""
    def __init__(self, hfpsi):
        super().__init__(hfpsi)

    def add_local_terms(self):
        self.O_ket = self.ket


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
        self.ainv = 1/a
        self.R = r0 * np.pow(hfpsi.n_neutrons + hfpsi.n_protons, 1/3)
        self.derivatives = ['Laplacian' if self.mesh.dim >= 2 else 'xx']

    def add_local_terms(self):
        """Create self.O_ket, fill it with the Woods-Saxon potential, and add the kinetic energy."""
        self.O_ket = HFPsi.like(self.ket)

        # TODO: speed up with numba decorators?

        def V_WoodsSaxon1D(r):
            return self.V0 / (1. + np.exp(self.ainv * (r - self.R)))

        def V_WoodsSaxon2D(x,y):
            return self.V0 / (1. + np.exp(self.ainv * (np.sqrt(x*x + y*y) - self.R)))

        def V_WoodsSaxon3D(x,y,z):
            return self.V0 / (1. + np.exp(self.ainv * (np.sqrt(x*x + y*y + z*z) - self.R)))

        dim = self.mesh.dim
        if dim == 3:
            V_WoodsSaxon = V_WoodsSaxon3D
        elif dim == 2:
            V_WoodsSaxon = V_WoodsSaxon2D
        else: # dim == 1:
            V_WoodsSaxon = V_WoodsSaxon1D

        # aply the Woods-Saxon potential
        Vr = self.mesh.apply(V_WoodsSaxon)
        self.O_ket.data = Vr * self.O_ket.data

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
            n = self.O_ket.hfblockrange[3][1] # end of neutron range in the spwfs and begin of proton range
            self.O_ket.data[:,:,:n] += hbm_n * self.ket.derivatives[nabla][:,:,:n]
            self.O_ket.data[:,:,n:] += hbm_p * self.ket.derivatives[nabla][:,:,n:]

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
        super().__init__(hfpsi,hbm)
        self.V0 = V0
        self.ainv = 1/a
        self.R = r0 * np.pow(hfpsi.n_neutrons + hfpsi.n_protons, 1/3)
        
    def add_local_terms(self):
        """Create self.O_ket, fill it with the Woods-Saxon potential, and add the kinetic energy."""
    
        # TODO: speed up with numba decorators?
        # Define Woods-Saxon potential functions
        def V_WoodsSaxon1D(r):
            return self.V0 / (1. + np.exp(self.ainv * (r - self.R)))

        def V_WoodsSaxon2D(x,y):
            return self.V0 / (1. + np.exp(self.ainv * (np.sqrt(x*x + y*y) - self.R)))

        def V_WoodsSaxon3D(x,y,z):
            return self.V0 / (1. + np.exp(self.ainv * (np.sqrt(x*x + y*y + z*z) - self.R)))

        # Select Woods-Saxon potential functions
        V_WoodsSaxon = V_WoodsSaxon3D if self.mesh.dim == 3 else \
                       V_WoodsSaxon2D if self.mesh.dim == 2 else \
                       V_WoodsSaxon1D

        # Apply the Woods-Saxon potential to the mesh points
        Vr = self.mesh.apply(V_WoodsSaxon)

        # add to O_ket
        self.O_ket.data += Vr * self.O_ket.data

