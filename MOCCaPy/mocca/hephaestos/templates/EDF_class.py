"""
$DESCRIPTION
"""
from mocca.edf.param import Param
from mocca.edf import EDF
from mocca.edf.ccfunctions import skyrme_cc, skyrme_ct

class $FUNC_NAME(EDF):
    """ A class implementing a functional of the $FUNC_NAME type"""

    def __init__(self, param: Param):
        """ Constructor

        Args:
            param : a Param object (parameterization)

        """
        super().__init__(param)
        self.calculate_coupling_constants()

    def calculate_coupling_constants(self):
        """
            Calculate the coupling constants from the Skyrme parameters
        """

        self.coupling_constants={}
        #---------------------------------------------------------------------
        # HEPHAESTOS inserts code that calculates the coupling constants here
        # <<<<< HEPHAESTOS
$CC_CALCULATION
        # >>>>> HEPHAESTOS
        #---------------------------------------------------------------------



