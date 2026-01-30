"""
    PyTest-compatible tests for the hephaestos submodule
"""
from pathlib import Path
import pytest
from pytest import approx
import sys
from importlib import import_module

from mocca.edf.param import Param
from mocca.hephaestos import read_functional_from_file, generate_EDF_class

path2MOCCaPy = Path(__file__).parent
while not path2MOCCaPy.name == 'MOCCaPy':
    path2MOCCaPy = path2MOCCaPy.parent

func_files = (path2MOCCaPy/'mocca/hephaestos/func_files/').glob('*.func')


tmp_folder = path2MOCCaPy/"tests/mocca/hephaestos/tmp"
sys.path.insert(0, str(tmp_folder))

@pytest.mark.parametrize("func_file", func_files)
def test_read_functional_from_file(func_file):
    """
        (1) Execute the read_functional_from_file function for a given .func file
        (2) generate a corresponding EDF file
        (3) Import it to check for syntax errors

        Args:
            func_file :  string
                file containing the functional specification
    """
    edf_specification = read_functional_from_file(func_file)

    generate_EDF_class(edf_specification, location=tmp_folder)

    # For now, we skip all functional files with '-' in the name -> this generates invalid Python code
    # ET: that is fixed by replacing '-' with '_"
    print('Loading the new EDF class')
    new_edf = import_module(edf_specification['module_name'])
    print(new_edf.__dir__())

def test_calculate_coupling_constants_NLO():
    """
        Compare the calculation of the time-even coupling constants with known values
        for the SLy4 parameterization of the NLO EDF.

        TODO: extend this to include the time-odd coupling constants

        Args:
            None
    """

    # Go to the top folder
    project_folder = Path(__file__).parent.parent.parent.parent.parent

    edf_specification = read_functional_from_file(path2MOCCaPy / 'mocca/hephaestos/func_files/NLO.func')
    generate_EDF_class(edf_specification, location=tmp_folder)
    nlo = import_module(edf_specification['module_name'])

    SLy4_param = Param(path2MOCCaPy / 'mocca/data/parameterizations/SLy4.param')
    SLy4       = nlo.NLO(SLy4_param)

    # Compare to known coupling constants for SLy4
    """
    Raw MOCCa output
        ______________________________________________________________________________
                                                Particle-hole terms
        Term                                   Isospin     # #G      Value
        ______________________________________________________________________________
        E_D_I_I_D_I_I                         | 0 0     |  1  1|     -933.342375000000
        E_D_I_I_D_I_I                         | 1 1     |  2  2|      830.052485500000
        E_D_I_I_D_I_I_D_I_I                   | 0 0 0   |  3  1|      861.062525661662
        E_D_I_I_D_I_I_D_I_I                   | 0 1 1   |  4  2|    -1064.273281717815
        E_D_I_I_D_I_I_D_I_I                   | 0 0 0   |  5  3|        0.000000000000
        E_D_I_I_D_I_I_D_I_I                   | 0 1 1   |  6  4|       -0.000000000000
        E_D_I_Sm_D_I_Sm                       | 0 0     |  7  1|     -207.824235500000
        E_D_I_Sm_D_I_Sm                       | 1 1     |  8  2|      311.114125000000
        E_D_I_I_D_I_Sm_D_I_Sm                 | 0 0 0   |  9  1|      490.231597943373
        E_D_I_I_D_I_Sm_D_I_Sm                 | 0 1 1   | 10  2|     -287.020841887221
        E_D_I_I_D_I_Sm_D_I_Sm                 | 0 0 0   | 11  3|        0.000000000000
        E_D_I_I_D_I_Sm_D_I_Sm                 | 0 1 1   | 12  4|       -0.000000000000
        E_D_I_I_D_Nm_Nm                       | 0 0     | 13  1|       57.128687500000
        E_D_I_I_D_Nm_Nm                       | 1 1     | 14  2|       24.656736500000
        E_D_I_I_Lap_D_I_I                     | 0 0     | 15  1|      -76.996203125000
        E_D_I_I_Lap_D_I_I                     | 1 1     | 16  2|       15.657135125000
        E_C_I_NmSk_C_I_NmSk                   | 0 0     | 17  1|       -0.000000000000
        E_C_I_NmSk_C_I_NmSk                   | 1 1     | 18  2|       -0.000000000000
        E_C_I_Nm_C_I_Nm                       | 0 0     | 19  1|      -57.128687500000
        E_C_I_Nm_C_I_Nm                       | 1 1     | 20  2|      -24.656736500000
        E_D_I_Sm_Lap_D_I_Sm                   | 0 0     | 21  1|        0.000000000000
        E_D_I_Sm_Lap_D_I_Sm                   | 1 1     | 22  2|        0.000000000000
        E_D_I_Sm_D_Nk_NkSm                    | 0 0     | 23  1|        0.000000000000
        E_D_I_Sm_D_Nk_NkSm                    | 1 1     | 24  2|        0.000000000000
        E_D_I_I_Derm_C_I_NxmSxm               | 0 0     | 25  1|      -92.250000000000
        E_D_I_I_Derm_C_I_NxmSxm               | 1 1     | 26  2|      -30.750000000000
        E_C_I_Nm_Derxm_D_I_Sxm                | 0 0     | 27  1|      -92.250000000000
        E_C_I_Nm_Derxm_D_I_Sxm                | 1 1     | 28  2|      -30.750000000000
        ______________________________________________________________________________
                                                Pairing terms
        ______________________________________________________________________________
        E_DP_I_I_DP_I_I                       | n n     | 29  1|     -312.500000000000
        E_DP_I_I_DP_I_I                       | p p     | 30  2|     -312.500000000000
        E_D_I_I_DP_I_I_DP_I_I                 | 0 n n   | 31  1|     1953.125043655747
        E_D_I_I_DP_I_I_DP_I_I                 | 0 p p   | 32  2|     1953.125043655747
        E_CP_I_I_CP_I_I                       | n n     | 33  1|     -312.500000000000
        E_CP_I_I_CP_I_I                       | p p     | 34  2|     -312.500000000000
        E_D_I_I_CP_I_I_CP_I_I                 | 0 n n   | 35  1|     1953.125043655747
        E_D_I_I_CP_I_I_CP_I_I                 | 0 p p   | 36  2|     1953.125043655747
        ______________________________________________________________________________
    """

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Particle-hole terms
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # LO central time-even
    assert SLy4.coupling_constants[('E_D_I_I_D_I_I'          , ('0', '0'), '1')]            == approx (  -933.342375000000)
    assert SLy4.coupling_constants[('E_D_I_I_D_I_I'          , ('1', '1'), '1')]            == approx (   830.052485500000)
    assert SLy4.coupling_constants[('E_D_I_I_D_I_I_D_I_I'    , ('0', '0', '0'), 'sigma')]   == approx (   861.062525661662)
    assert SLy4.coupling_constants[('E_D_I_I_D_I_I_D_I_I'    , ('0', '1', '1'), 'sigma')]   == approx ( -1064.273281717815)
    assert SLy4.coupling_constants[('E_D_I_I_D_I_I_D_I_I'    , ('0', '0', '0'), 'sigma_b')] == approx (     0.0)
    assert SLy4.coupling_constants[('E_D_I_I_D_I_I_D_I_I'    , ('0', '1', '1'), 'sigma_b')] == approx (     0.0)

    # NLO central time-even
    assert SLy4.coupling_constants[('E_D_I_I_D_Nm_Nm'        , ('0', '0'), '1')]            == approx (    57.128687500000)
    assert SLy4.coupling_constants[('E_D_I_I_D_Nm_Nm'        , ('1', '1'), '1')]            == approx (    24.656736500000)
    assert SLy4.coupling_constants[('E_D_I_I_Lap_D_I_I'      , ('0', '0'), '1')]            == approx (   -76.996203125000)
    assert SLy4.coupling_constants[('E_D_I_I_Lap_D_I_I'      , ('1', '1'), '1')]            == approx (    15.657135125000)
    assert SLy4.coupling_constants[('E_C_I_NmSk_C_I_NmSk'    , ('0', '0'), '1')]            == approx (     0.0)
    assert SLy4.coupling_constants[('E_C_I_NmSk_C_I_NmSk'    , ('1', '1'), '1')]            == approx (     0.0)

    # NLO spin-orbit time-even
    assert SLy4.coupling_constants[('E_D_I_I_Derm_C_I_NxmSxm', ('0', '0'), '1')]            == approx (   -92.250000000000)
    assert SLy4.coupling_constants[('E_D_I_I_Derm_C_I_NxmSxm', ('1', '1'), '1')]            == approx (   -30.750000000000)

    # Pairing terms
    assert SLy4.coupling_constants[('E_DP_I_I_DP_I_I', ('n', 'n'), '1')]                    == approx (  -312.500000000000)
    assert SLy4.coupling_constants[('E_DP_I_I_DP_I_I', ('p', 'p'), '1')]                    == approx (  -312.500000000000)
    assert SLy4.coupling_constants[('E_D_I_I_DP_I_I_DP_I_I', ('0', 'n', 'n'), 'sigmap')]    == approx (  1953.125043655747)
    assert SLy4.coupling_constants[('E_D_I_I_DP_I_I_DP_I_I', ('0', 'p', 'p'), 'sigmap')]    == approx (  1953.125043655747)


if __name__ == "__main__":
    test_read_functional_from_file('mocca/hephaestos/func_files/NLO.func')
    test_calculate_coupling_constants()
