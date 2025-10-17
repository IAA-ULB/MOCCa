from pathlib import Path
from mocca.edf.read_param import read_param


class BXL:
    """Concrete Energy Density Functional class."""

    def __init__(self, param_file: Path):
        """BXL constructor.

        Args:
            param_file (Path): Path to the .param file.

        """
        read_param(param_file, self)


if __name__ == "__main__":
    bxl = BXL(Path("../../../parameterizations/MSk7.param"))
    print(bxl)

