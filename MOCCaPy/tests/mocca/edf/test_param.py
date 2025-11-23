import pytest
from pathlib import Path

from mocca.edf import Param
from mocca.edf.bxl import BXL

project_folder = Path(__file__).parent.parent.parent.parent.parent
assert project_folder.name == 'tantalus_full'


def test_Param0():
    """Test setting Param object's attributes."""
    param = Param(Path(__file__).parent/"FunctionalNotImplemented.param")
    assert param.name == "FunctionalNotImplemented"
    assert param.func_file == "FunctionalNotImplemented.func"
    assert param.p1  == 1
    assert param.p2  == 2.0
    assert param.p3  == 3e-9
    assert param.p4  == 'str'
    assert param.p5  == "str"
    assert param.p6  == [6, 6]
    assert param.p7  == 7
    assert param.p8  == 8
    assert param.p11 == 11
    assert param.p12 == 12.0
    assert param.p13 == 13e-9
    assert param.p14 == 'str'
    assert param.p15 == "str"
    assert param.p16 == [16, 16]
    assert param.p17 == 17
    assert param.p18 == 18
    assert param.p19 == True
    assert param.p20 == False


def test_Param():
    """Try all .param files to make sure that all corner cases are handled well."""
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
    """Create a BXL object from a Param object that has `func_file=='BXL'`."""
    param = Param(project_folder / "parameterizations/BSkG1.param")
    bxl = param.create_EDF()
    assert isinstance(bxl, BXL)

    param = Param("BSkG1")

def test_FunctionalNotImplemented_from_param():
    """Create a BXL object from a Param object that specifies a functional that is
    not implemented."""
    param = Param(Path(__file__).parent/"FunctionalNotImplemented.param")
    with pytest.raises(ModuleNotFoundError):
        nlo = param.create_EDF()