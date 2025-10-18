from pathlib import Path
import re

pattern = re.compile(r"^(\w+)\(d+\)$")

def read_param(filepath: Path, edf_object, assert_config=True):
    """Read EDF parameters from file and store them in the edf_object (as attributes).

    Args:
        filepath (Path): Path to the .param file.
        edf_object: a concrete EDF object, e.g. a BXL object, which will store the parameters.
        assert_config (bool): whether or not to raise an assertion error when func_file does not
            correspond to the class of edf_object. assert_config=False is only used in testing.
    Raises:
        AssertionError: when func_file does not correspond to the class of edf_object. e.g class BXL
            requires func_file = 'BXL.func'.
    """
    edf_object.filepath = filepath
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
                    _handle_line(line, edf_object,assert_config=assert_config)
                else:
                    # more than one entry on this line
                    sub_lines = line.split(',')
                    for line in sub_lines:
                        _handle_line(line, edf_object, assert_config=assert_config)

def _handle_line(line, edf_object, assert_config):
                key, value = line.split("=")
                key = key.strip()
                value = value.strip()
                try:
                    if key == 'name':
                        key = 'param_name'
                    else:
                        m = pattern.match(key)          # key has form 'name(1)'
                        if m:
                            key = m.group(1)

                    if key == 'func_file':              # Assert that the .param file is compatible with edf_object
                        CONFIG = value[1:-1].replace('.func', '')
                        if assert_config:
                            assert CONFIG == edf_object.__class__.__name__, f"EDF class `{edf_object.__class__.__name__}` is incompatible with .param file `{edf_object.filepath}`"

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

                    setattr(edf_object, key, value)

                except AssertionError:
                    raise
                except Exception as e:
                    print(f"{e} Error in {read_param}: {key=}={value=} ")
                    print(f"{line=}")
                    raise




