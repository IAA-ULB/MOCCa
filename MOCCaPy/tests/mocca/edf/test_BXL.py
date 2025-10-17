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