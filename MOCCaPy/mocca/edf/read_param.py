from pathlib import Path
import re

def read_param(filepath: Path, edf_object):
    """Read EDF parameters from file and store them in the edf_object.

    Args:
        filepath (Path): Path to the .param file.
        edf_object: a concrete EDF object, e.g. a BXL object, which will store the parameters.
    """
    pattern = re.compile(r"^(\w+)\(d+\)$")
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
            elif line == "":                         # empty line
                pass
            else:                                    # key-value pair
                if '!' in line:
                    line = line.split('!')[0]      # remove trailing comment and strip
                key, value = line.split("=")
                key = key.strip()
                value = value.strip()
                try:
                    if key == 'name':
                        key = 'param_name'
                    else:
                        m = pattern.match(key)
                        if m:
                            key = m.group(1)

                    if key == 'func_file':              # Assert that the .param file is compatible with edf_object
                        CONFIG = value[1:-1].replace('.func', '')
                        assert CONFIG == edf_object.__class__.__name__, f"EDF class {edf_object.__class__.__name__} is incompatible with .param file {filepath}"

                    if value.endswith(","):             # remove trailing comma
                        value = value[:-1]
                    if value.startswith("'"):           # str
                        value = value[1:-1]
                    elif value == ".true.":             # bool
                        value = True
                    elif value == ".false.":            # bool
                        value = False
                    elif '.' in value:                  # float
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




