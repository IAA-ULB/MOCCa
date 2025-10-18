from mocca.param import Param

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

        class_name = self.__class__.__name__
        CONFIG = param.func_file.replace('.func', '')
        assert CONFIG == class_name, f"\n  EDF class `{class_name}` is incompatible with .param file `{self.param.filepath}`"\
                                     f"\n  .param file requires EDF class `{CONFIG}`."




