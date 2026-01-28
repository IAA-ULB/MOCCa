"""
    PyTest-compatible tests for the hephaestos submodule
"""
import glob
import pytest
from mocca.hephaestos.hephaestos import read_functional_from_file

func_names = glob.glob('mocca/hephaestos/func_files/*.func')

@pytest.mark.parametrize("func_file", func_names)
def test_read_functional_from_file(func_file):
    """
        Execute the read_functional_from_file function for a given 

        Args:
            func_file :  string 
                file containing the functional specification
    """
    dic = read_functional_from_file(func_file)
    print(dic['coupling_constants'])

if __name__ == "__main__":
    test_read_functional_from_file("NLO.func")
