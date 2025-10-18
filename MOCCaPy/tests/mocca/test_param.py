import pytest
from pathlib import Path

from mocca.param import Param
from mocca.edf.bxl import BXL

project_folder = Path(__file__).parent.parent.parent.parent
assert project_folder.name == 'tantalus_full'


def test_Param():
    # Try all .param files (without assert_config) to make sure that all corner cases are handeled.
    for p in (Path(project_folder)/'parameterizations').glob('*.param'):
        if p.name in ['forces.param']:
            # Following Wouter's advice, we do not support forces.param
            # (because the latter contains many .param files in one.)
            pass
        else:
            print(p, end='...')
            param = Param(p)
            print(' OK')

def test_BXL_from_param():
    param = Param(project_folder / "parameterizations/BSkG1.param")
    bxl = param.create_EDF()
    assert isinstance(bxl, BXL)

def test_FunctionalNotImplemented_from_param():
    param = Param(Path(__file__).parent/"FunctionalNotImplemented.param")
    with pytest.raises(ModuleNotFoundError):
        nlo = param.create_EDF()