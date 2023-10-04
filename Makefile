#-------------------------------------------------------------------------------
#
# Makefile for the succesfull compilation of different Tantalus executables.
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# For succesful compilation, one needs
#
# * a complete copy of the Tantalus repository, including
#   - a working Hephaestos version
#   - a correct set of Tantalus source code template files
#
# * a working installation of some version of Python3 to run Hephaestos
#   - which should have access to basic Python libraries and numpy in particular
#
# * a configuration file to run Hephaestos with; a ton are provided in configs/.
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
# Basic useage:
#
#   > make
#
# which will compile a standard executable to be placed in the exec/ folder
# using gfortran (provided it is installed).
#
# For more control, specify additional options either in this Makefile itself
# or on the command line. For example:
#
#   > make CONFIG=BXL COMPILER=ifort
#
# will compile a Tantalus executable based on the BXL.py configuration file
# (look in the configs/ folder) using the Intel ifort compiler.
#
# The option to specify compiler and CONFIG file should be sufficient for most
# users; changing any other options is at your own risk.
#
# Executables at the end will be named Tantalus.$(CONFIG).exe and be placed
# in the $(EXECDIR) configured below.
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# OPTIONS
# - - - - -
#  CXX      : compiler to use
#  CONFIG   : name of configutation file in the config/ folder
#            (without trailing .py)
#  OPTFLAGS : optimisation compiler flags
#  CXXFLAGS : other compiler flags
#  PRE      : steps to do before compilation
#  DEBUG    : 0 => no debugging options
#             1 => full debugging options
#  USE_MPI  : 0 => no MPI
#             1 => MPI (note: your compiler should be MPI-capable to use this)
#  EXECDIR  : directory for storage of the final executables
#  SRCDIR   : source code as processed by Hephaestos
#  OBJDIR   : storage for intermediate object files
#  MODDIR   : storage for final module files
#
# Notes
# - - - -
# 1. this is a Makefile, so you can essentially override ANY AND ALL
#    variables from the command line by simply passing them as argument to
#    make. The options I document above are only the ones that I think are
#    relevant to normal useage.
# 2. All relevant directories will be created if they don't exist already.
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Compilers that are currently "pre-configured" with appropriate optimisation
# flags etc.
#
# 1. gfortran
#    optimisation: -O3, -O3 -ffast-math, -Ofast, -O3 -funroll-loops,
#                  -Ofast -funroll-loops
#    versions tested: 8.5.0, 9.4.0
#
# 2. ifort
#    optimisation: -Ofast, -O3, -O3 -xhost
#    versions tested: 2021.1
#
# 3. cray compilers
#    optimisation: -O2, -O3, -O3 -hfp3
#    versions tested: 14.0.3
#-------------------------------------------------------------------------------
# Acknowledgment:
#   the organisation of this Makefile as well as a bunch of options are
#   inspired by the Makefile of the HFBTHO v4 code, see the repository of
#   P. Marević et al., Computer Physics Communications 276, 108367 (2022).
#-------------------------------------------------------------------------------

################################################################################
# Physics details (modify as you want)
################################################################################
CONFIG  := default
EXENAME := Tantalus.$(CONFIG).exe

################################################################################
# Compilation details (this section should be modified as you see fit)
################################################################################

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Compiler type 
COMPILER      :=  gfortran
#
# Note: this is NOT the compiler wrapper that will get invoked for execution
#       it is rather the "type" of compiler (or the organisation behind it)
#       this variable is used to set different compilation options for which
#       syntax and/or linking might not be identical
#       A set of default options are present for
#         - gfortran  by GNU
#         - ifort     by Intel
#         - cray      by Cray
#
#       The primary reason that COMPILER and CXX are different is because 
#       vendors have different compiler wrappers for different modes 
#       (i.e. with or without MPI) but which nevertheless have similar options.
#       A bonus reason is that this can easily account for versions, i.e. 
#       compilation with gfortran-9.3 will get the right options set.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# USE_MPI
#  => 0 if inactive
#  => 1 if active
USE_MPI := 0
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Actual compiler wrapper that gets invoked
#  I provide default options based on the USE_MPI and COMPILER options
#  but it is up to the user to make sure that CXX and COMPILER match. The 
#  compilation will obviously fail if, say, an INTEL compiler gets invoked with 
#  GNU options. If you want to specify a specific wrapper, fill the following 
#  line
CXX :=
ifeq ($(CXX), )
ifeq ($(COMPILER),gfortran)
  ifeq ($(USE_MPI),1)
    CXX := mpifort 
    # on the systems available to me, this is the wrapper for
    # MPI-enabled GFORTRAN
  else
    CXX := gfortran
endif
else ifeq ($(COMPILER),ifort)
  ifeq ($(USE_MPI),1)
    CXX := mpiifort # on the systems available to me, this is the wrapper for
                   # MPI-enabled IFORT
  else
    CXX := ifort
  endif
else ifeq ($(COMPILER), cray)
  CXX := ftn
  # I have yet to figure out MPI with CRAY compilers
endif
endif
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Directory names  (will be created if they don't exist)
#
# - EXECDIR = directory for storage of the final executables
# - SRCDIR  = source code as processed by Hephaestos
# - OBJDIR  = storage for intermediate object files
# - MODDIR  = storage for final module files
EXECDIR :=  exec
SRCDIR  :=   src
OBJDIR  :=   obj
MODDIR  :=   mod


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# DEBUG
# => 0 : compile without debugging options
# => 1 : compile with debugging options for each compiler
DEBUG   := 0

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Libraries for linear algebra
# This can be specified on the command line, but is in practice compiler based
ifeq ($(COMPILER),gfortran)
	# versions of gfortran should link to OPENBLAS
	LIBS := -lopenblas
else ifeq ($(COMPILER),ifort)
  # ifort compiler should link to the new Intel math library
	LIBS := -mkl
else ifeq ($(COMPILER), cray)
  # Cray compilers don't need specific linking to my knowledge
	LIBS :=
endif

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Python interpreter with which to invoke Hephaestos
PYTHON_CMD := python3

################################################################################
# Setting compiler flags based on information provided above.
# (This section should NOT be modified in principle)
################################################################################

# 1. set some compiler-specific options concerning storage etc.
ifeq ($(COMPILER),gfortran)
	CXXFLAGS := -J$(MODDIR)
else ifeq ($(COMPILER),ifort)
	CXXFLAGS := -module $(MODDIR) -assume realloc-lhs -assume byterecl -no-wrap-margin
else ifeq ($(CXX),ftn)
	CXXFLAGS := -J$(MODDIR)
endif

# 2. set compiler-specific optimisation level
# .... when in production mode
ifeq ($(DEBUG),0)
  ifeq ($(COMPILER),gfortran)
	  OPTFLAGS := -O3
  else ifeq ($(COMPILER),ifort)
	  OPTFLAGS := -Ofast
  else ifeq ($(COMPILER),ftn)
	  OPTFLAGS := -O3 # to be tested if optimal
  endif
else
  ifeq ($(COMPILER),gfortran)
	  OPTFLAGS := -O0 -g -Wall -Wno-uninitialized -fbacktrace
  else ifeq ($(COMPILER),ifort)
	  OPTFLAGS := -g -traceback
  else ifeq ($(COMPILER),ftn)
	  OPTFLAGS := -e c -e D # to be tested if optimal
  endif
endif
################################################################################
# Sanity checks (should not be modified)
################################################################################

# 1. Check if compiler exists
ifeq (, $(shell which $(CXX)))
 $(error "The compiler wrapper CXX=$(CXX) cannot be invoked. Is it installed?")
endif

# 2. Check if config file exists
CONFIGPATH := configs/$(CONFIG).py
ifeq ("$(wildcard $(CONFIGPATH))", "")
 $(error "The configuration file '$(CONFIG).py' does not seem to exist in the configs/ directory.")
endif

# 2. Check if python3 interpreter exists
ifeq (, $(shell which $(PYTHON_CMD)))
 $(error "The python interpreter $(PYTHON_CMD) does not seem to be installed.")
endif

################################################################################
# Files to be compiled (This section should NOT be modified in principle)
################################################################################

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Tantalus source files
TARGET :=   Tantalus.exe
SRC    :=   compilation.f90 geninfo.f90 timing.f90 constants.f90 
SRC    +=   sphericalharmonics.f90 folding.f90   
SRC    +=   nil8.f90 derivatives.f90 precondition.f90 wavefunctions.f90
SRC    +=   pairingcutoffs.f90 parameterization.f90
SRC    +=   pairing_strengths.f90 basis_transform.f90 hartree-fock.f90 BCS.f90 
SRC    +=   HFB_gradient.f90 HFB_direct.f90 HFB.f90
SRC    +=   pairing.f90 densities.f90 moments.f90 
SRC    +=   coulomb.f90 cranking.f90 momentsofinertia.f90 transform.f90  
SRC    +=   functional.f90 fission_MOI.f90  evolution.f90 scfiteration.f90 
SRC    +=   IO.f90 temperature_projection.f90 convergence.f90 printing.f90 
SRC    +=   tantalus.version.f90
SINGLE_SRC = $(SRC) run_single.f90
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Nilsson source files
NIL_SRC := compilation.f90 timing.f90 geninfo.f90 derivatives.f90 nil8.f90   
NIL_SRC += wavefunctions.f90 gennilsson.f90
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Translate source files into object files
SINGLE_OBJ  :=  $(patsubst %.f90,$(OBJDIR)/%.o,$(SINGLE_SRC))
NIL_OBJ     :=  $(patsubst %.f90,$(OBJDIR)/%.o,$(NIL_SRC))

################################################################################
# Explicit precompilation steps 
#
# 1) Run Hephaestos to preprocess the entire code
# 2) Get version information from git
# 3) Get compiler information 
# 4) set the version and compiler info in the source code
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
PRE         :=  run_heph getgitinfo getcompilerinfo setversioninfo 
PRE_NIL     :=  cp_nil 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Internal (to the compiler) preprocessing directives
#    -cpp      => explicitly enable preprocessing
#    -DUSE_MPI => enable (1) or disable (0) MPI (see above) 
PREPROCESSOR :=  -cpp -DUSE_MPI=$(USE_MPI)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

################################################################################
# Recipes (This section should NOT be modified in principle)
################################################################################

all: single

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Creation of required directories 
$(EXECDIR)/:
	mkdir -p $(EXECDIR)/

$(OBJDIR)/:
	mkdir -p  $(OBJDIR)/

$(MODDIR)/:
	mkdir -p  $(MODDIR)/

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
single: $(PRE) $(SINGLE_OBJ)
	$(CXX) $(OPTFLAGS) $(CXXFLAGS) $(PREPROCESSOR) -o $@ $(SINGLE_OBJ) $(LIBS)
	mv single exec/$(EXENAME)

run_heph:
  # Run Hephaestos with the correct configuration file
	python3 Hephaestos.py $(CONFIG) 

gen_nilsson: $(PRE_NIL) $(NIL_OBJ)
	$(CXX) $(OPTFLAGS) $(CXXFLAGS) $(PREPROCESSOR) -o $@ $(NIL_OBJ) $(LIBS)
	mv gen_nilsson exec/$(EXENAME)

clean:
	rm  -f $(OBJDIR)/*.o
	rm  -f $(MODDIR)/*.mod

$(OBJDIR)/%.o : $(SRCDIR)/%.f90 | $(OBJDIR)/ $(MODDIR)/ exec/ 
	$(CXX)  $(OPTFLAGS) $(CXXFLAGS) $(PREPROCESSOR) -c  $< -o $@ 

setversioninfo:
# Copy the git information into the main code, so it can be printed
	@cp $(SRCDIR)/tantalus.f90 $(SRCDIR)/tantalus.version.f90
	@sed -i.bak 's/VERSION1/"${GIT_INFO1}"/' $(SRCDIR)/tantalus.version.f90
	@sed -i.bak 's/VERSION2/"${GIT_INFO2}"/' $(SRCDIR)/tantalus.version.f90
	@sed -i.bak 's/VERSION3/"${GIT_INFO3}"/' $(SRCDIR)/tantalus.version.f90
	@sed -i.bak 's/VERSION4/"${GIT_INFO4}"/' $(SRCDIR)/tantalus.version.f90
#Copy the compiler information
	@sed -i.bak 's/COMPCOMP/"${COMPVERSION}"/' $(SRCDIR)/tantalus.version.f90
	@sed -i.bak 's/CFLAGS/"${CXXFLAGS}"/'      $(SRCDIR)/tantalus.version.f90
	@sed -i.bak 's/OPTFLAGS/"${OPTFLAGS}"/'    $(SRCDIR)/tantalus.version.f90
	@rm $(SRCDIR)/tantalus.version.f90.bak

getgitinfo:
# Get information from 'git show', to see what kind of build this is.
	$(eval GIT_INFO1=$(shell git show   | grep 'commit ' | head -1))
	$(eval GIT_INFO2=$(shell git show   | grep 'Author:' | head -1))
	$(eval GIT_INFO3=$(shell git show   | grep 'Date:'   | head -1))
	$(eval GIT_INFO4=$(shell git branch | grep '*'       | head -1 | cut -c2- ))

getcompilerinfo:
  # Get information from 'CXX --version'
	$(eval COMPVERSION=$(shell $(CXX) --version | head -1))

cp_nil:
	cp src_orig/gennilsson.f90 $(SRCDIR)/gennilsson.f90

################################################################################

