import numpy as np
import sys
import pyscalapack
from mpi4py import MPI
from elpa import Elpa


scalapack = pyscalapack(
    "/apps/antwerpen/skylake/rocky9/ELPA/2024.05.001-intel-2024a/lib/libelpa.so",
    "/data/antwerpen/201/vsc20170/tantalus_full/experiments/elpa/pyelpa/libelpa_object_interface.so"
)  
with Elpa(scalapack) as e:
    na = 8
    nev= 8
    na_rows = 2
    na_cols = 2
    nblk = 2
    mpi_comm_world = MPI.COMM_WORLD.py2f()
    my_prow = 1
    my_pcol = 1

    e.set(b'na' , na)
    e.set(b'nev', nev)
    e.set(b"local_nrows", na_rows)
    e.set(b"local_ncols", na_cols)
    e.set(b"nblk", nblk)
    e.set(b"mpi_comm_parent", mpi_comm_world)
    e.set(b"process_row", my_prow)
    e.set(b"process_col", my_pcol)

    