import pytest
from pathlib import Path

from mocca.edf.BXL import BXL
project_folder = Path(__file__).parent.parent.parent.parent.parent
print(project_folder)

def test_ctor():
    print(project_folder)
    bxl = BXL(param_file=project_folder/"parameterizations/BSkG1.param")
    with pytest.raises(AssertionError):
        bxl = BXL(param_file=project_folder/"parameterizations/MSk7.param")

def test_read_param():
    for p in (Path(project_folder)/'parameterizations').glob('*.param'):
        print(p)
        if p.name in ['BSkG1.param', 'forces.param']:
            pass
        else:
            bxl = BXL(param_file=p, assert_config=False)
