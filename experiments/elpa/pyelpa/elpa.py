# from pdb import set_trace as BREAKPOINT

import ctypes

_DBG = True
class Elpa:
    """
    Python context manager object holding a (Fortran) Elpa object.
    """
    elpa_api_version = None
    def __init__(self, scalapack, elpa_api_version=20240501):
        """
        """
        self.scalapack = scalapack
        if _DBG:
            print("Elpa.__init__")
        scalapack.elpa_initialize()
        if _DBG:
            print("Elpa.__init__ done")
    
    def __enter__(self):
        return self
    
    def __exit__(self, type, value, traceback):
        """
        """
        if _DBG:
            print("Elpa.__exit__")
        self.scalapack.elpa_finalize()
        if _DBG:
            print("Elpa.__exit__ done")

    def set(self, name:str, val):
        if _DBG:
            print(f"Elpa.set({name=},{val=})")
        
        c_error = ctypes.c_int(len(name))

        if isinstance(val, int):
            c_val = ctypes.c_int(val)
        elif isinstance(val, float):
            c_val = ctypes.c_double(val)
        elif isinstance(val, np.float32):
            c_val = ctypes.c_float(val)
        
        self.scalapack.set_integer(name, c_val, c_error)
        
        if c_error.value != 0:
            raise RuntimeError(f"ERROR : e%set({name=}, {val=}, error={c_error.value})")
