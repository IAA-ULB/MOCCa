from pathlib import Path
import re
from importlib import import_module

from src_heph.fortran_templates.GenTermExpression_templates import edent_fam

_array_element_pattern = re.compile(r"^(\w+)\(d+\)$")

class Param:
    def __init__(self, filepath: Path):
        """Read EDF parameters from file and store them as attributes of self.

        Args:
            filepath (Path): Path to the .param file.
        Raises:
            RuntimeError: when reading the .param file fails.
        """
        self.filepath = filepath
        with open(filepath) as f:
            lines = f.readlines()
            for line in lines:
                if line.startswith("!"):                 # comment line
                    pass
                elif line.startswith("&skf"):            # start of fortran name list
                    pass
                elif line.startswith("/"):               # end of fortran name list
                    pass # next block we're done
                elif line == "\n":                         # empty line
                    pass
                else:                                    # key-value pair
                    if '!' in line:
                        line = line.split('!')[0]      # remove trailing comment and strip
                    if '#' in line:
                        line = line.split('#')[0]      # remove trailing comment and strip
                    n = line.count('=')
                    if n == 1:
                        self._handle_single_entry(line)
                    else:
                        # more than one entry on this line
                        sub_lines = line.split(',')
                        for line in sub_lines:
                            self._handle_single_entry(line)

    def _handle_single_entry(self, line):
        key, value = line.split("=")
        key = key.strip()
        value = value.strip()
        try:
            if key == 'name':
                key = 'param_name'
            else:
                m = _array_element_pattern.match(key)          # key has form 'name(1)'
                if m:
                    key = m.group(1)

            if value.endswith(","):             # remove trailing comma
                value = value[:-1]
            if value.startswith("'") or \
               value.startswith('"'):           # str
                value = value[1:-1]
            elif value == ".true.":             # bool
                value = True
            elif value == ".false.":            # bool
                value = False
            elif '.' in value or \
                 'e' in value or \
                 'E' in value or \
                 'd' in value:                  # e.g. 1.5, 1e-8, 1E-8
                if 'd' in value:
                    value = value.replace('d', 'e')     # 1d-8 not supported by python
                if ',' in value:                # array of floats
                    values = value.split(',')
                    value = [float(value.strip()) for value in values]
                else:
                    value = float(value)
            else:                       # int
                value = int(value)

            setattr(self, key, value)

        except Exception as e:
            raise RuntimeError(f"Error reading file {self.filepath} at entry {line}.\n"
                               f"Error message: {e}")

    def create_EDF(self):
        """Create an EDF object defined by `self.func_file`.

        Raises:
            ModuleNotFoundError: If the functional `self.func_file` is not found in mocca.edf.
        """

        edf_name = self.func_file.replace('.func', '')
        edf_module_name = f"mocca.edf.{edf_name.lower()}"
        edf_module = import_module(edf_module_name)

        create_cmd = f"edf_module.{edf_name}(self)"
        edf = eval(create_cmd)
        return edf


