# Compiling MOCCa

## Prerequisites

MOCCa has several dependencies; to succesfully compiler

1. **Python3**: a working installation that includes the numpy and scipy libraries.
2. **Linear algebra packages**: BLAS in all cases and either LAPACK (no MPI) or SCALAPACK (MPI).
3. **HDF5**: this is only required if you request HDF5 I/O -- which is highly recommended!
4. [**make.inc**](make.inc.md): a compilation file that is configured correctly for your system.
5. [a configuration file](config.md): this describes thes physics you need MOCCa to handle.

Our recommendation is to use prepackaged make.inc and configuration files as much as possible.
These can be found in the `make_includes/` and `configs`, respectively. If you need more information,
please refer to the pages on these subjects in the navigation menu.

## Basic compilation

Supposing you have taken care of the prerequisites, you can now proceed to compile the
code by invoking make in the main directory as follows:

    make $MC  CONFIG=$configname CALCTYPE=$calctype

 where

 - MC : the type of many-body calculation you would like to perform; either `mf` for static 
        mean-field or `fam` for linear response calculations. 

 - configname : the name of the configuration file in the `configs/` folder that you wish to use.

    For example, type `CONFIG=BXL` to use `BXL.py` located in the configs folder.

 - calctype   : the type of calculation that you want to do
     - NUCLEI: (default) impose antiperiodic boundary conditions to model finite nuclei.
     - PASTA: impose periodic boundary conditions to model crystalline lattices.

If nothing fails, the compilation process will either create `Tantalus.$configname.exe` (`MC = mf`)
or `fam.$configname.exe` (`MC=fam`) inside the `exec/` directory.

## Advanced compilation

There are many tricks you can play with the MOCCa compilation process since you can override any variable
in make.inc or the Makefile from the command line. I here list only a few typical commands that have been useful in the past, 
you can surely come up with more to suit your workflow.

- `make EXENAME=$exe`:
   Manually specify the name of the final executable as `$exe`, overriding the default convention.
