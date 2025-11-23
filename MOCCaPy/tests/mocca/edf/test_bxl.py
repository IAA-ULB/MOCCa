import pytest
from pathlib import Path

from mocca.edf import Param
from mocca.edf.bxl import BXL

project_folder = Path(__file__).parent.parent.parent.parent.parent
assert project_folder.name == 'tantalus_full'

def test_BXL_ctor():
    """Construct a BXL object from a Param object."""
    param = Param(project_folder/"parameterizations/BSkG1.param")
    bxl = BXL(param)

