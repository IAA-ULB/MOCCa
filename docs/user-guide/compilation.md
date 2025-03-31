# Compiling MOCCa

## Prerequisites

MOCCa has several dependencies; to succesfully compiler

1. **Python3**: a working installation that includes the numpy and scipy libraries.
2. **Linear algebra packages**: BLAS in all cases and either LAPACK (no MPI) or SCALAPACK (MPI).
3. **HDF5**: this is only required if you request HDF5 I/O.
4. [**make.inc**](make.inc.md): a compilation file that is configured correctly for your system.
5. [a configuration file](config.md): this describes thes physics you need MOCCa to handle.

Our recommendation is to use prepackaged make.inc and configuration files as much as possible.
These can be found in the `make_includes/` and `configs`, respectively. If you need more information,
please refer to the pages on these subjects in the navigation menu.

## Basic compilation

Supposing you have taken care of the prerequisites, you can now proceed to compile the
code by invoking make in the main directory as follows:

    make CONFIG=$configname CALCTYPE=$calctype

 where

 - configname : the name of the configuration file in the `configs/` folder that you wish to use.

    For example, type `CONFIG=BXL` to use `BXL.py` located in the configs folder.

 - calctype   : the type of calculation that you want to do
     - NUCLEI: (default) impose antiperiodic boundary conditions to model finite nuclei.
     - PASTA: impose periodic boundary conditions to model crystalline lattices.

If nothing fails, the compilation process will create an executable named `Tantalus.$configname.exe`
inside the `exec/` directory.

## Advanced compilation

 TO BE DESCRIBED
 - Makefile command line options
 - EXENAME specification
 - make.inc command line option
