import numpy as np
from tabulate import tabulate

from .operator import Operator
from mocca.util import title_line
from mocca.edf import EDF

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
        self.set_property_list()

    def __repr__(self):
        return f"<KineticEnergyOperator(Operator)[{self.operand_data.mfs}]>"

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

    def set_property_list(self, property_list=None):
        self.property_list = property_list

    def __str__(self):
        """Create a string representation of the mean-field state in the form of a table.

        Feature requested in issues/66.
        """
        # TODO: add customization of the table.
        mfs_name = self.operand_data.mfs.__class__.__name__
        h_ii, d2h = self.compute_dispersion()
        # self.operand_data.mfs.rho[:], shell, order_n, order_p = self.operand_data.mfs.occupancies(h_ii=h_ii)
        self.operand_data.mfs.rho[:], order_n, order_p = self.operand_data.mfs.occupancies(h_ii=h_ii)
        P = np.ones(self.operand_data.mfs.n_total_wf, dtype=np.int32)
        S = np.ones(self.operand_data.mfs.n_total_wf, dtype=np.int32)
        for ib in range(8):
            if ib % 4 >= 2:
                r = self.operand_data.hfblockrange[ib]
                for i in range(*r):
                    P[i] = -1
            if ib % 2 == 1:
                r = self.operand_data.hfblockrange[ib]
                for i in range(*r):
                    S[i] = -1
        tbl = SpwfTable(order_n, order_p, mfs_name)
        tbl.add_column('i', np.arange(self.operand_data.ket.n_total_wf)
                                 , property_list=self.property_list)
        tbl.add_column('p', P    , property_list=self.property_list)
        tbl.add_column('s', S    , property_list=self.property_list)
        tbl.add_column('occ', self.operand_data.mfs.rho
                                 , property_list=self.property_list)
        tbl.add_column('E', h_ii , property_list=self.property_list)
        tbl.add_column('d2h', d2h, property_list=self.property_list)
        return str(tbl)

# helper class
class SpwfTable:
    """Create spwf tables like
    ```
    ------------- SlaterDeterminant -----------------
    -------------------------------------------------
    --- Neutron functions ---------------------------
      n    i    p    s    occ          E          d2h
    ---  ---  ---  ---  -----  ---------  -----------
      2    0    1    1      1  -39.0633   1.70414e-08
      4    7   -1    1      1  -29.0719   2.67243e-08
      6    8   -1    1      1  -29.0719   2.67228e-08
      8    9   -1    1      1  -29.0719   2.67255e-08
     10    3    1    1      1  -17.6487   9.24864e-07
     ...
    -------------------------------------------------
    --- Proton functions ----------------------------
      n    i    p    s    occ          E          d2h
    ---  ---  ---  ---  -----  ---------  -----------
      2   15    1    1      1  -39.0633   1.97081e-08
      4   22   -1    1      1  -29.0719   2.91585e-08
      6   24   -1    1      1  -29.0719   2.91673e-08
      8   23   -1    1      1  -29.0719   2.91703e-08
     10   18    1    1      1  -17.6487   1.01297e-06
     ...
    -------------------------------------------------
    ```
    Useful for monitoring progress.
    As requested in https://github.com/IAA-nuclear/tantalus_full/issues/66
    """
    def __init__(self, order_n, order_p, name):
        """
        Args:
            order_n: indices of the spwf's when sorted from low to high energy
                (diagonal elements of the hamiltonian)
            order_p: indices of the protons when sorted from low to high energy
                (diagonal elements of the hamiltonian)
            name: class name of the mean-field state object, Appears in the title line.
        """
        self.order_n = order_n
        self.n = len(order_n)
        self.order_p = order_p
        self.p = len(order_p)
        self.order = np.concatenate((order_n, order_p))
        self.n_columns = []
        self.p_columns = []
        self.headers = []
        self.name = name
        n_spwfs = len(self.order)
        shell = 2*np.ones(n_spwfs, dtype=np.int32)
        self.add_column('n', shell, property_list=None)

    def add_column(self, name, data, property_list):
        """Add a column to the table, The order of adding is also the print order.

        Args:
            name: column header
            data: np.ndarray. Its length is the total number of spwfs in the mean-field state.
                The order of the elements is the order of the spwfs in the mean-field state.
                They are automatically sorted from low to high energy.
            property_list: list of property names to include in the table. None includes all
                known properties.
        """
        if (property_list is not None) and \
           (name not in property_list):
            return

        self.headers.append(name)
        if name == 'n':
            self.n_columns.append(np.cumsum(data[self.order[:self.n]]))
            self.p_columns.append(np.cumsum(data[self.order[self.n:]]))
        else:
            self.n_columns.append(data[self.order[:self.n]])
            self.p_columns.append(data[self.order[self.n:]])

    def transpose(self, columns):
        """Transpose the columns into rows so that tabulate can handle it."""
        ncols = len(columns)
        nrows = len(columns[0])
        rows = [ ]
        for irow in range(nrows):
            rows.append([])
            for icol in range(ncols):
                rows[irow].append(columns[icol][irow])
        return rows

    def __str__(self):
        """Create a string representation of the mean-field state."""
        rows = self.transpose(self.n_columns)
        sn = tabulate(rows, headers=self.headers, tablefmt="simple")
        sn +='\n'

        w = sn.index('\n')
        s = title_line(text=self.name, char='-', width=w, start=-1, below=True)
        s += title_line(text='Neutron functions', char='-', width=w)
        s += sn
        s += title_line(text='Proton functions', char='-', width=w, above=True)
        rows = self.transpose(self.p_columns)
        s += tabulate(rows, headers=self.headers, tablefmt="simple")
        s += '\n'
        s += title_line(char='-', width=w)
        return s

#=======================================================================================================================
# Derived classes
#=======================================================================================================================
def V_WoodsSaxon1D(r, V0, ainv, R):
    return V0 / (1. + np.exp(ainv * (r - R)))

def V_WoodsSaxon2D(x, y, V0, ainv, R):
    return V0 / (1. + np.exp(ainv * (np.sqrt(x * x + y * y) - R)))

def V_WoodsSaxon3D(x, y, z, V0, ainv, R):
    return V0 / (1. + np.exp(ainv * (np.sqrt(x * x + y * y + z * z) - R)))

class HamiltonianWoodsSaxon(KineticEnergyOperator):
    """A hamiltonian with a Woods-Saxon potential. <hfpsi|h|hfpsi>."""

    def __init__(self, mfs, hbm=20.73553000, V0=-50.0, r0=1.25, a=0.5):
        """Initialize the operator with a mean-field state instance `mfs` to operate on
        and set the parameters for the operator.

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

    def __repr__(self):
        return f"<HamiltonianWoodsSaxon(KineticEnergyOperator)[{self.operand_data.mfs}]>"

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

class EdfHamiltonian(KineticEnergyOperator):
    """Hamiltonian based an energy density functional."""

    def __init__(self, mfs, hbm=20.73553000, edf=None):
        """Initialize the operator with a wave function `hfpsi` to operate on and set the parameters for the
        operator.

        Args:
            hbm: prefactor of the kinetic energy operator. A single value considers the masses of neutrons
                and protons to be identical, a tuple of 2 values (hbar^2/2m_n, hbar^2/2m_p) considers
                the masses to differ.
            V0, r0, a:  parameters of the Woods-Saxon potential
                (cfr https://en.wikipedia.org/wiki/Woods–Saxon_potential).
        """
        if not issubclass(edf.__class__, EDF):
            raise ValueError("Parameter edf=`{edf=}`, must be a subclass of mocca.edf.EDF")

        super().__init__(mfs, hbm)
        self.edf = edf

    # TODO: override base methods.