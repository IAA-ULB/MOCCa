import pytest
from pathlib import Path

from mocca.edf import EDF
from mocca.param import Param

project_folder = Path(__file__).parent.parent.parent.parent.parent
assert project_folder.name == 'tantalus_full'

def test_EDF_ctor():
    """Construct a Param object and assert that the EDF object is raising an AssertionError. """
    param = Param(project_folder / 'parameterizations/BSkG1.param')
    with pytest.raises(AssertionError):
        # This must always raise AssertionError because the base class is always incompatible with a parameterization object
        edf = EDF(param)