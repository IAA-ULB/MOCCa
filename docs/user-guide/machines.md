# Machine-specific remarks
## IAA machines at ULB

Compilation 

TO BE DESCRIBED

 - mention ifort launching
 - HDF5 installation?

## Lucia - the TIER-0 of CENAERO

The recommended `make.inc` files for computations on LUCIA are

- `make.inc.cray-serial` : for non-MPI calculations
- `make.inc.cray-parallel` : for MPI calculations

These rely on the CRAY compiler suite, 

Please refer to the original documentation of this machine [here](https://doc.lucia.cenaero.be/). 



 - recommended compiler = cray

 - module load Prg-env Cray
 - module load Cray
 - ...

## LUMI

Compilers provided by CRAY are the recommended option for LUMI. If you want to
use these, then you will need to load the following modules:

- cray-python
- PrgEnv-cray
- craype-hugepages2M
- cray-hdf5 OR cray-hdf5-parallel (optional)

Only the first two are strictly necessary for compilation. The third avoids
segfaults that might occur in large calculations. The loading of a HDF5 module
is only necessary if you enable it at compiletime.

The recommended make.inc files are make.inc.cray-serial and make.inc.cray-parallel;
note that these enable HDF5 support by default.

## MareNostrum
 We don't have acces yet.
