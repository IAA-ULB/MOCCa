"""
Import this module to exclude code that is only executed when there are >1 MPI processes.

```
import mpi_process as process

if process.using_mpi:
    # Code below only executed when there are >1 MPI processes.
    ...
```
"""
using_mpi = False
rank = 0
size = 1
