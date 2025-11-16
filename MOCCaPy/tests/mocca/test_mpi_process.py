import pytest

def test_process():
    try:
        import mocca.mpi_process as process
    except UserWarning as w:
        print(f"\n{w.__class__.__name__}:\n{w}")
    else:
        assert process.size >1

