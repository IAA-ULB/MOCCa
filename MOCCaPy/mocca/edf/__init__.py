from .param import Param

class EDF:
    """Abstract Energy Density Functional class."""
    def __init__(self, param:Param):
        """
        Args:
            param :a parameterization object
            edf_name : name of the derived EDF class. if provided asserts that the param argument
                is compatible with the derived EDF class. E.g. if `self.param.func_file == `BXL.func``
                the derived EDF class must be `BXL`.
        """
        self.param = param

        CONFIG = param.func_file.replace('.func', '')
        if not  (CONFIG == self.__class__.__name__):
            raise ValueError(
                f"\n  EDF class `{self.__class__.__name__}` is incompatible with .param file `{self.param.filepath}`"\
                f"\n  .param file requires EDF class `{CONFIG}`."
            )




