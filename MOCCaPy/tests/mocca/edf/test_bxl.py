import pytest
from pathlib import Path

from mocca.edf.bxl import BXL
from mocca.param import Param

project_folder = Path(__file__).parent.parent.parent.parent.parent
assert project_folder.name == 'tantalus_full'

def test_BXL_ctor():
    # test the the BXL ctor which basically tests the read_param method
    param = Param(project_folder/"parameterizations/BSkG1.param")
    bxl = BXL(param)
    with pytest.raises(AssertionError):
        param = Param(project_folder / "parameterizations/MSk7.param")
        bxl = BXL(param)

