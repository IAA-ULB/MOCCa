import pytest

def test_process():
    import mocca.util.no_mpi_process as process
    assert process.using_mpi == False
    assert process.size == 1
    assert process.rank == 0
    # test code selection
    if process.using_mpi:
        assert False
    else:
        assert True
