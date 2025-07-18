import numpy as np
import pyscalapack

scalapack = pyscalapack(
    "/apps/antwerpen/zen2/rocky8/impi/2021.13.0-intel-compilers-2024.2.0/mpi/2021.13/lib/libmpi.so",
    "/apps/antwerpen/zen2/rocky8/imkl/2024.2.0/mkl/2024.2/lib/intel64/libmkl_core.so",
    "/apps/antwerpen/zen2/rocky8/imkl/2024.2.0/mkl/2024.2/lib/intel64/libmkl_sequential.so",
    "/apps/antwerpen/zen2/rocky8/imkl/2024.2.0/mkl/2024.2/lib/intel64/libmkl_intel_lp64.so",
    "/apps/antwerpen/zen2/rocky8/imkl/2024.2.0/mkl/2024.2/lib/intel64/libmkl_blacs_intelmpi_lp64.so",
    "/apps/antwerpen/zen2/rocky8/imkl/2024.2.0/mkl/2024.2/lib/intel64/libmkl_scalapack_lp64.so",
    )  


class Elpa:
    def __init__(self):
        """An Elpa-like object for pyev.
        Not all pyelpa functionality is provided (yet)
        """

    @classmethod
    def from_distributed_matrix(cls, a):
        """Initialize an Elpa objectc with values from a distributed matrix `a`

        Args:
            a (DistributedMatrix): matrix to get values from
        """
        self = cls()
        # Set parameters the matrix and it's MPI distribution
        self.processor_layout = a.processor_layout
        self.na = a.na
        self.nev = a.nev
        self.local_nrows = a.na_rows
        self.local_ncols = a.na_cols
        self.nblk = a.nblk
        
        # Setup
        # self.setup()
        # if desired, set tunable run-time options
        # self.set_integer("solver", ELPA_SOLVER_2STAGE)
        return self
    
    def eigenvectors(self, a_data, ev, q_data):
        """
        Args:
            a_data : (input) symmetric matrix of which the eigenvalues/eigenvectors must be computed.
                distributed matrix data pointer, typically the data member of a numpy array
            ev : (output) computed eigenvalues. 1D numpy array with length na
            q_data : (output) computed eigenvectors. 
                distributed matrix data pointer, typically the data member of a numpy array
        """
        with self._get_context() as context:
            a_sc = context.array(self.na, self.na, self.nblk, self.nblk, data=a_data)
            q_sc = context.array(self.na, self.na, self.nblk, self.nblk, data=q_data)
            return a_sc.pdsyev(ev,q_sc)
        # caveat: here it is assumed that nev == na

    # --------------------------------------------
    # member functions for internal use (I think)

    def _get_context(self):
        """Get a pyscalapack Context"""
        return scalapack(b'C', self.processor_layout.np_rows, self.processor_layout.np_cols)
    
    def _get_context0(self):
        """Get a pyscalapack Context on rank 0 """
        return scalapack(b'C', 1, 1)
    




