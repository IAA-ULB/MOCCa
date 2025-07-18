from pyev.helpers import *
import inspect

def test_import():
    print()
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    header = f"{inspect.currentframe().f_code.co_name}()-{comm2str(comm)}"
    print(header, f': backend={use_pyelpa}')

def test_barrier():
    print()
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    comm.Barrier()
    header = f"{inspect.currentframe().f_code.co_name}()-{comm2str(comm)}"
    print(header, f': backend={use_pyelpa}')

def test_size():
    print()
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    assert comm.size > 0
    header = f"{inspect.currentframe().f_code.co_name}()-{comm2str(comm)}"
    print(header, f': backend={use_pyelpa}')
 