import pytest
from pathlib import Path

from mocca.edf.BXL import BXL
project_folder = Path(__file__).parent.parent.parent.parent.parent
print(project_folder)

def test_ctor():
    # test the the BXL ctor which basically tests the read_param method
    print(project_folder)
    bxl = BXL(param_file=project_folder/"parameterizations/BSkG1.param")
    with pytest.raises(AssertionError):
        bxl = BXL(param_file=project_folder/"parameterizations/MSk7.param")

def test_read_param():
    # Try all .param files (without assert_config) to make sure that all corner cases are handeled.
    for p in (Path(project_folder)/'parameterizations').glob('*.param'):
        print(p)
        if p.name in ['BSkG1.param', 'forces.param']:
            # BSkG1.param is already tested in test_ctor, and we do not support forces.param
            # (because the latter contains many .param files in one.)
            pass
        else:
            bxl = BXL(param_file=p, assert_config=False)
