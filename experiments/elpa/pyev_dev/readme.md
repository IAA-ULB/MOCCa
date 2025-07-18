# pyev

Wrapper for PyscaLAPACK and pyelpa. 

> ***What?*** 
> Choose a backend (PyscaLAPACK|pyelpa) and run exactly the same code.

> ***How?*** 
> pyev is either an alias for pyelpa and mimicks pyelpa if PyScaLAPACK is choosen as a backend.

Let's make it even simpler. If you want to use pyelpa:

``` python
import pyelpa

# code using pyelpa
...
```
if you want to use PyScaLAPACK:

``` python
import pyev as pyelpa

# same code as above, now using PyScaLAPACK as a backend
...
```
So, pyev provides the same functionality as pyelpa, but but uses PyScaLAPACK.

## testing

Each test - obviously - has two cases, one with pyscalapack as backend, and one with pyelpa as backend. 

Using pytest with mpi is generally supported with plugins, e.g.:

- [pytest-mpi](https://pytest-mpi.readthedocs.io/)
- [pytest-isolate-mpi](https://pytest-isolate-mpi.readthedocs.io/en/stable/
  )

This works fine for the pyscalapack backend. 

For the pyelpa backend, however, python seems to import mpi4py before the individual tests are entered. That yields

```shell
NotImplementedError: Please load the pyelpa module before mpi4py, otherwise there will be MPI problems.
```

None of the pytest MPI plugins seems to be ok with this.

However, without the plugins (!) everything works just fine:

``` shell
> USE_PYELPA=0 mpirun -n 4 python -m pytest -s tests/test_1.py 
```
Above, the backend is selected via the environment variable USE_PYELPA, which if set to `1` uses pyelpa as backend and pyev/pyscalapack otherwise.

In a python script or pytest script the backend can be selected based on the `USE_PYELPA` environment variable:
s
``` python
from pyev.helpers import use_pyelpa
if use_pyelpa:
    from pyelpa import ProcessorLayout, DistributedMatrix, Elpa
else:
    from pyev   import ProcessorLayout, DistributedMatrix, Elpa
from mpi4py import MPI

# your code using either pyev or pyelpa should work without change
```

> [!Note]
> As the classes `ProcessorLayout` and `DistributedMatrix` are defined in file `distributedmatrix.py` which for pyev a copy of that in pyelpa, it does not matter from which package they are imported. The `Elpa` class, however is different in pyev (see `wrapper.py`) and pyelpa (see `wrapper.pyx`). 

Alternatively, you can set `use_pyelpa` from `sys.argv`:

``` python
# file my_script.py
from pyev.helpers import use_pyelpa
import sys
use_pyelpa.from_argv()
if use_pyelpa:
    from pyelpa import ProcessorLayout, DistributedMatrix, Elpa
else:
    from pyev   import ProcessorLayout, DistributedMatrix, Elpa
from mpi4py import MPI
...
```
Here, the following command line options are valid:

``` shell
> python my_script.py --use-pyelpa=0|1 # -> False|True (pyev/pyscalapack|pyelpa)
> python my_script.py --use-pyelpa 0|1 # -> False|True (pyev/pyscalapack|pyelpa)
> python my_script.py -s               # -> False      (pyev/pyscalapack       ) 
> python my_script.py -e               # -> True       (                 pyelpa) 
```

> [!Note]
> Not all Elpa functionality has (yet) been implemented in pyev with the `pyscalapack` backend.

Al
``` shell
vsc20170@login9@breniac [1027] ~/workspace/tantalus_full/experiments/elpa/pyev_dev
> USE_PYELPA=1 mpirun -n 4 python -m pytest -x tests/ 
============================= test session starts ==============================
platform linux -- Python 3.12.3, pytest-8.2.2, pluggy-1.5.0
rootdir: /data/antwerpen/201/vsc20170/tantalus_full/experiments/elpa/pyev_dev
plugins: xdist-3.6.1
collected 181 items

tests/test_0.py ============================= test session starts ==============================
platform linux -- Python 3.12.3, pytest-8.2.2, pluggy-1.5.0
rootdir: /data/antwerpen/201/vsc20170/tantalus_full/experiments/elpa/pyev_dev
plugins: xdist-3.6.1
collected 181 items

tests/test_0.py ============================= test session starts ==============================
platform linux -- Python 3.12.3, pytest-8.2.2, pluggy-1.5.0
rootdir: /data/antwerpen/201/vsc20170/tantalus_full/experiments/elpa/pyev_dev
plugins: xdist-3.6.1
collected 181 items

tests/test_0.py ============================= test session starts ==============================
platform linux -- Python 3.12.3, pytest-8.2.2, pluggy-1.5.0
rootdir: /data/antwerpen/201/vsc20170/tantalus_full/experiments/elpa/pyev_dev
plugins: xdist-3.6.1
collected 181 items

tests/test_0.py .....                                                    [  1%].
tests/test_1.py ..                                                       [  1%]
tests/test_1.py                                                          [  1%]
tests/test_1.py                                                          [  1%]
tests/test_1.py ...............                                          [  3%]
tests/test_mpi4py.py .                                                   [  3%]
tests/test_mpi4py.py 
tests/test_mpi4py.py ...                                                 [  3%]
tests/test_mpi4py.py .........                                           [  4%]
tests/test_numroc.py 
tests/test_numroc.py                                                     [  4%]
tests/test_numroc.py                                                     [  4%]
tests/test_numroc.py .                                                   [  5%]
tests/test_pyev_import.py 
tests/test_pyev_import.py                                                [  5%]
tests/test_pyev_import.py ..                                             [  5%]
tests/test_pyev_import.py ....                                           [  6%]
tests/test_with_mpi.py                                                   [  6%]
tests/test_with_mpi.py 
tests/test_with_mpi.py .                                                 [  6%].
tests/test_with_mpi.py ................................................
.......................................................................
.......................................................................
...                                                                      [ 33%]
.......................................................................
.......................................................................
.......................................................................
.......................................................................
...                                                                      [ 73%]
.......................................................................
.......................................................................
..................................                                       [100%]

============================= 181 passed in 38.35s =============================
.......                                                                  [100%]

============================= 181 passed in 39.35s =============================
....                                                                     [100%]

============================= 181 passed in 39.73s =============================
.....                                                                   [100%]

============================= 181 passed in 42.27s =============================
```

With `USE_PYELPA=0` several tests in `tests/test_with_mpi,py` fail, due to missing functionality in `pyev.Elpa`.