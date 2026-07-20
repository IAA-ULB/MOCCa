## The make.inc file

The make.inc file specifies various compilation options that are not related to
physics: from directory names to compiler flags and library linking.

It is **not recommended** to write a make.inc file from scratch. Instead, the code
comes with a predefined set of files that are stored in the make_include/ directory.
Their names are variations on the following theme

      make.inc.gnu-serial-hdf5

in which
      gnu    => indicates compilers provided by gnu (gfortran,mpifort,....)
      serial => indicates that MPI parallelisation is disabled
      hdf5   => indicates that HDF5 in- and output is enabled

At the time of writing, make_includes/ has predefined files for nearly all
possible combinations of

      compilers      : gnu, cray, intel
      parallelisation: serial, MPI or OpenMP
      HDF5           : yes or no

as well as separate sets of compilation options more suited to debugging.

!!! warning 
    Only rather specific calculations can benefit from parallelisation. 
    Specifically: large-scale pasta calculations should leverage MPI, 
    linear response calculations should leverage OpenMP. All other 
    calculations will not benefit from either.

To use a predefined file from the make_includes/ folder, it suffices to copy
the file into the main MOCCa directory and rename it. For example, execute
the following comand from inside the main directory:

      cp make_include/make.inc.gnu-serial-hdf5 make.inc

If you insist on writing your own make.inc file or want to modify an existing
one for use, the following table summarizes all the options you are required to set.

| Variable           | Description
|:-------------------|:---------------------------------------------------------
|EXECDIR             | storage of the final executables
|MF_SRC_DIR          | storage of the source code for the mean-field code as processed by Hephaestos
|MF_OBJ_DIR          | storage for intermediate object files  for the mean-field code
|MF_MOD_DIR          | storage for final module files for the mean-field code
|                    | 
|FAM_SRC_DIR         | storage of the source code for the linear response code as processed by Hephaestos
|FAM_OBJ_DIR         | storage for intermediate object files for the linear response code
|FAM_MOD_DIR         | storage for final module files for the linear response code
|                    |
|PYTHON_CMD          | Python interpreter used to invoke Hephaestos
|                    | 
|CXX                 | compiler invokation to be used; e.g. 'gfortran' or 'mpiifort'
|OPTFLAGS            | compiler flags related to optimisation
|CXXFLAGS            | compiler flags not related to optimisation
|PREPFLAG            | the relevant syntax to invoke the preprocessor for your choice of CXX
|DEBUG_LEVEL         | debugging flags to pass to the compiler
|                    |
| USE_MPI            | whether (1) or not (0) to compile with MPI support 
|                    | 
| USE_HDF5           | whether to offer HDF5 support, yes(1) or no (0)
| HDF5_LIB           | linking statements for the HDF5 library
| 
| LINEAR_ALGEBRA_LIB | linking statements for (Sca)LAPACK and BLAS
