import ctypes
import os
import sys
import numpy as np


_DBG = True
class Wrapper():

    @staticmethod
    def _default_dll_loader(lib):
        return ctypes.CDLL(lib, mode=os.RTLD_LAZY | os.RTLD_GLOBAL)

    def __init__(self, *libs, loader=None):
        """
        Create scalapack library handle from the given pathes as scalapack dynamic linked lbiraries.

        Parameters
        ----------
        *libs : list[str]
            The dynamics linked libraries full pathes or just the file names if in the default directory.
        loader : Callable
            The loader for the dynamic shared libraries.
        """
        if loader is None:
            loader = self._default_dll_loader
        self.libs = [loader(lib) for lib in libs]
        self.function_database = {}

    def __getattr__(self, name, suffix='_'):
        """
        Get a function from the scalapack library, and wrap it by fortran function wrapper.

        Parameters
        ----------
        name : str
            The function name, without suffix "_" in fortran function.

        Returns
        -------
            The function wrapped by fortran function wrapper.
        """
        # All function will be cached in function_database
        if name not in self.function_database:
            # Get the real function name since fortran function has a suffix "_" in its name.
            real_name = name + suffix
            # Look it up in all libraries
            for lib in self.libs:
                # If it is in a library
                if hasattr(lib, real_name):
                    # Get it, wrapper it, save it, and return.
                    self.function_database[name] = self._fortran_function(getattr(lib, real_name))
                    break
            else:
                raise AttributeError(f"No function named '{real_name}' in the libraries")
            print(f"function {real_name} found.")
        return self.function_database[name]

    class Val():
        """
        An auxiliary type to indicate an argument should be passed by value to a fortran function.
        """

        def __init__(self, value):
            """
            Create Val wrapper

            Parameters
            ----------
            value
                The value should be passed by value.
            """
            self.value = value

    @classmethod
    def _resolve_arg(cls, arg):
        """
        Resolve argument as the fortran style: pass by reference by default, except indicated by `Val`.
        """
        if _DBG:
            print(f"_resolve_arg : {type(arg)=} {arg=}") # for debugging
        if isinstance(arg, cls.Val):
            if _DBG:
                print(f"_resolve_arg : {type(arg)=} {arg=} cls.Val") # for debugging
            # This argument is specified to pass by value
            return arg.value
        elif isinstance(arg, int):
            if _DBG:
                print(f"_resolve_arg : {type(arg)=} {arg=} int") # for debugging
            # This is a python int, wrap it in c_int
            arg = ctypes.c_int(arg)
            return ctypes.byref(arg)
        elif isinstance(arg, bytes):
            if _DBG:
                print(f"_resolve_arg : {type(arg)=} {arg=} bytes") # for debugging
            # This is a python bytes, wrap it in c_char_p
            arg = ctypes.c_char_p(arg)
            return arg
        elif isinstance(arg, np.ndarray):
            if _DBG:
                print(f"_resolve_arg : {type(arg)=} {arg=} np.ndarray") # for debugging
            return arg.ctypes.data_as(ctypes.c_void_p)
        elif isinstance(arg, ctypes.c_void_p):
            if _DBG:
                print(f"_resolve_arg : {type(arg)=} {arg=} ctypes.c_void_p") # for debugging
            return arg.value
        else:
            if _DBG:
                print(f"_resolve_arg : {type(arg)=} {arg=} else") # for debugging
            # This must be already a ctypes object, get the reference.
            return ctypes.byref(arg)

    @classmethod
    def _fortran_function(cls, function):
        """
        Wrapper for fortran function in dynamic linked library.

        See Also
        --------
        Scalapack._resolve_arg : Resolve argument to fit the fortran calling style.
        """

        def result(*args):
            if _DBG:
                print(f"{len(args)} arguments")
            return function(*(cls._resolve_arg(arg) for arg in args))

        result.__doc__ = function.__doc__
        return result
