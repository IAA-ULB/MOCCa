from mocca.param import Param
from mocca.edf import EDF


class BXL(EDF):
    """Concrete Energy Density Functional class."""

    def __init__(self, param: Param):
        """BXL constructor.

        Args:
            param : a Param object (parameterization)

        """
        super().__init__(param)


