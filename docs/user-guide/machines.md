# Machine-specific remarks

## IAA machines at ULB

The astropcs at the IAA are a bit of a heterogeneous set; nevertheless they all have access to intel compilers 
so the recommended make.inc is `make.inc.intel-serial`. Remember that modern ifort installations on linux machines 
require the explicit setting of environment variables through something like 

    source /opt/intel/oneapi/setvars.sh --force &> /dev/null

which I recommend you add to your `.bashrc`.


Given that HDF5 libraries can be installed in multiple ways, there is not much help I can systematically get you 
beyond the linking options that work on astropc19, included in `make.inc.intel-serial-hdf5`.


## CECI machines

On all CECI machines, you should load the relevant Python libraries to run Hephaestos. The relevant commands are:
```
module load Python 
module load Scipy-bundle
```
which I recommend you add to your `.bashrc` file.

Most CECI machines offer intel compilers; on those I would advocate to also execute
```
module load intel
```
and to compile with `make.inc.intel-serial`.

On Lyra, there is currently no support for intel compilers. There, I recommend compiling with `make.inc.serial`.


## Lucia - the TIER-0 of CENAERO

The recommendation is to use CRAY compilers on LUCIA, i.e. use the following `make.inc` files 

- `make.inc.cray-serial`   : for non-MPI calculations
- `make.inc.cray-parallel` : for MPI calculations

To correctly compile the code, you will need to load:
```
module load Cray
module load CPE
module load PrgEnv-cray
module load cray-python
```
The ordering of these commands is not arbitrary.

You can optionally enable HDF5: 
```
module load cray-hdf5 OR cray-hdf5-parallel 
```

Please do not forget to refer to the original documentation of this machine [here](https://doc.lucia.cenaero.be/). 

## MareNostrum
 We don't have acces yet to this machine.
 
## LUMI

The recommendation is to use CRAY compilers on LUMI, i.e. use the following `make.inc` files 

- `make.inc.cray-serial`   : for non-MPI calculations
- `make.inc.cray-parallel` : for MPI calculations

To correctly compile the code, you will need to load 
```
module load cray-python
module load PrgEnv-cray
module load craype-hugepages2M
```

and optionally
```
module load cray-hdf5 OR cray-hdf5-parallel 
```

Only the first two are strictly necessary for compilation. The third (`craype-hugeâges2M`) avoids
segfaults that might occur in large calculations. The loading of a HDF5 module
is only necessary if you enable it at compiletime.

The recommended make.inc files are make.inc.cray-serial and make.inc.cray-parallel;
note that these enable HDF5 support by default.

