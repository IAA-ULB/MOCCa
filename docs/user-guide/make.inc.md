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
      parallelisation: serial, parallel
      HDF5           : yes or no

as well as separate sets of compilation options more suited to debugging.

To use a predefined file from the make_includes/ folder, it suffices to copy
the file into the main Tantalus directory and rename it. For example, execute
the following comand from inside the main directory:

      cp make_include/make.inc.gnu-serial-hdf5 make.inc

If you insist on writing your own make.inc file or want to modify an existing
one for use, the following table summarizes all the options you are required to set.

| Variable           | Description
|:-------------------|:---------------------------------------------------------
|EXECDIR             | storage of the final executables
|SRCDIR              | storage of the source code as processed by Hephaestos
|OBJDIR              | storage for intermediate object files
|MODDIR              | storage for final module files
|                    |
|CXX                 | compiler invokation to be used; e.g. 'gfortran' or 'mpiifort'
|OPTFLAGS            | compiler flags related to optimisation
|CXXFLAGS            | compiler flags not related to optimisation
|PREPFLAG            | the relevant syntax to invoke the preprocessor for your choice of CXX
|PRE                 | steps to do before compilation.
|USE_MPI             | set to 0/1 to disable/enable MPI parallelisation
|                    |
| USE_HDF5           | whether to offer HDF5 support, yes(1) or no (0)
| HDF5_LIB           | linking statemetns for the HDF5 library
| LINEAR_ALGEBRA_LIB | linking statements for (Sca)LAPACK and BLAS
|--------------------|---------------------------------------------------------
