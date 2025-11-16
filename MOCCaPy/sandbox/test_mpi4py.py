"""
Test mpi4py framework
https://mpi4py.readthedocs.io/en/stable/
run this test suite as
mpirun -np 2 python -m pytest test_mpi4py.py
"""
# TODO: find a solution for this problem:
# The `is_using_mpi` method seems to detect `MPI` only if it is defined in the same scope as where
# `from mpi4py import MPI` is included.
# - I watched https://www.youtube.com/watch?v=p_UQ7tzUFLo from ArjanCodes, on "The Real Reason the
#   Singleton Pattern Exists". Interestingly emough, there is a Python construct that is actually
#   a singleton: the module.
#  - https://stackoverflow.com/questions/1109870/python-singleton-into-multiprocessing
#
# We want to be able to test whether mpi4py.MPI is available, and has more than one process, in a
# simple way, without importing mpi4py if we do not need it,
"""
# these tests work fine, but the mechanism does not allow what we want.

def test_without_mpi4py():

    def is_using_mpi():
        try:
            return MPI.COMM_WORLD.Get_size() > 1 # raises NameError if MPI not defined.
        except Exception as e:
            print(type(e), e) # sandbox/test_mpi4py.py <class 'NameError'> name 'MPI' is not defined
            return False

    assert not is_using_mpi()

def test_with_mpi4py():
    from mpi4py import MPI

    def is_using_mpi():
        try:
            return MPI.COMM_WORLD.Get_size() > 0
        except Exception as e:
            return False

    assert is_using_mpi() == True
    rank = MPI.COMM_WORLD.Get_rank()
    size = MPI.COMM_WORLD.Get_size()
    print(f"mpi rank: {rank}/{size}")
    assert 0 <= rank < size
"""

