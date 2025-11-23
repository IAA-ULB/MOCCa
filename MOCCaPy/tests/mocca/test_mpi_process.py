import pytest

def test_process():
    try:
        import mocca.util.mpi_process as process
    except UserWarning as w:
        print(f"\n{w.__class__.__name__}:\n{w}")
    else:
        assert process.size >1
        # test code selection
        if process.using_mpi:
            assert True
        else:
            assert False
