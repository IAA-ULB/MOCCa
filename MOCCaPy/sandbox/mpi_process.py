"""
Import this module to include code that is only executed when there are >1 MPI processes.

```
import mpi_process as process

if process.using_mpi:
    # Code below only executed when there are >1 MPI processes.
    ...
```

This module imports and exposes`mpi4py.MPI`.
A UserWarning is raised when `mpi4py.MPI` was imported, but only 1 process was requested.
"""
from mpi4py import MPI

using_mpi = True
rank = MPI.COMM_WORLD.Get_rank()
size = MPI.COMM_WORLD.Get_size()

if  size == 1:
    raise UserWarning(
        "Running using mpi4py with only 1 process.\n"
        "Maybe you forgot `mpirun -n <n_processes>`?"
    )