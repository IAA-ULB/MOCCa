import ctypes
import numpy as np

n_bytes_f64 = 8
nrows = 2
ncols = 5

clib = ctypes.cdll.LoadLibrary("libc.so.6")

clib.memcpy.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.c_size_t]
clib.memcpy.restype = ctypes.c_void_p

arr_from = np.arange(nrows * ncols).astype(np.float64)
arr_to = np.empty(shape=(nrows, ncols), dtype=np.float64)

print('arr_from:', arr_from)
print('arr_to:', arr_to)

print('\ncalling clib.memcpy ...\n')
clib.memcpy(arr_to, arr_from, nrows * ncols * n_bytes_f64)

print('arr_from:', arr_from)
print('arr_to:', arr_to)