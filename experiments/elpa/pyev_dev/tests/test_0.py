import pytest
from pyev.helpers import *
import inspect
import sys

def test_0():

    # print(f"{use_elpa()=}")
    if use_pyelpa:
        import pyelpa as pyev
    from pyev import ProcessorLayout, DistributedMatrix, Elpa

    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    header = f"{inspect.currentframe().f_code.co_name}()-{comm2str(comm)}"
    print(header, f': backend={use_pyelpa}')


def test_size():
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    assert comm.size > 0
    header = f"{inspect.currentframe().f_code.co_name}()-{comm2str(comm)}"
    print(header, f': backend={use_pyelpa}')



# ==============================================================================
# The code below is for debugging a particular test
# (normally all tests are run with pytest)
# ==============================================================================
if __name__ == "__main__":
    the_test_you_want_to_debug = test_0

    print("__main__ running", the_test_you_want_to_debug)

    # modify use_elpa if a command line argument is given:
    use_pyelpa.from_argv()

    the_test_you_want_to_debug()

    print("-*# finished #*-")
    exit(0)
# ==============================================================================
