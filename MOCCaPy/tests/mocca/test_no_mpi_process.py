import pytest

def test_process():
    import mocca.no_mpi_process as process
    assert process.using_mpi == False
    assert process.size == 1
    assert process.rank == 0
