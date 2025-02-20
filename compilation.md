# Compiling MOCCa

## Prerequisites

To compile MOCCa, you will need the following

1. a working installation of some version of Python3, including basic Python libraries such as numpy and scipy.
2. an appropriate linear algebra package. This means BLAS in all cases and either LAPACK or SCALAPACK depending on whether you choose to enable MPI.
3. If you (optionally) enable HDF5 in- and output, you will need a working HDF5 installation that matches your choice regarding serial or MPI compilation.
4. a make.inc file that is configured correctly for your system.
5. a configuration file that reflects the physics you want to describe.

## The make.inc file

The make.inc file specifies various compilation options that are not related to
physics: from directory names to compiler flags and library linking.

It is not recommended to write a make.inc file from scratch. Instead, the code
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
one for use, these are the options that it should set

* Directory structure. 
    - EXECDIR  : directory for storage of the final executables
    - SRCDIR   : source code as processed by Hephaestos
    - OBJDIR   : storage for intermediate object files
    - MODDIR   : storage for final module files

* Compiler options
    - CXX      : compiler invokation to be used; e.g. 'gfortran' or 'mpiifort'
    - OPTFLAGS : compiler flags related to optimisation
    - CXXFLAGS : compiler flags not related to optimisation
    - PREPFLAG : the relevant syntax to invoke the preprocessor for your choice of CXX
    - PRE      : steps to do before compilation.  
    - USE_MPI  : set to 0/1 to dis/enable MPI parallelisation 

* Linking options
    - USE_HDF5           : whether to offer HDF5 support, yes(1) or no (0)
    - HDF5_LIB           : linking statemetns for the HDF5 library
    - LINEAR_ALGEBRA_LIB : linking statements for (Sca)LAPACK and BLAS

## The configuration file 

 TO BE DESCRIBED

## Basic compilation

With all of the prequisites taken care of, you can now proceed to compile the 
code by invoking make in the main directory as follows:
 
    make CONFIG=$configname CALCTYPE=$calctype

 where
 - configname : the name of the configuration file (without the .py extension) from the configs/ folder that you wish to use. 
     For example, type "CONFIG=BXL" to use BXL.py. 
 - calctype   : the type of calculation that you want to do, either "NUCLEI" or "PASTA". 
 Note that "CALCTYPE" defaults to "NUCLEI" such that you can omit that specification for executables targetting isolated nuclei.

If nothing fails, the compilation process will create an executable named

   Tantalus.$configname.exe 
   
inside the exec/ directory.
 
## Advanced compilation

 TO BE DESCRIBED
 - Makefile command line options
 - EXENAME specification 
 - make.inc command line option

## Machine-specific remarks
### astropc's
 TO BE DESCRIBED
 
 - mention ifort launching 
 - HDF5 installation? 
### Lucia
 TO BE DESCRIBED

  relevant make.inc = make.inc.cray-serial
 
 - recommended compiler = cray

 - module load Prg-env Cray
 - module load Cray 
 - ...

### LUMI

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

### MareNostrum
 We don't have acces yet.
