#------------------------------------------------------------------------------------
#
# Makefile for the succesfull compilation of different Tantalus executables.
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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
# * a make.inc file configured for your system; a few are provided in make_include/.
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
# After having verified that make.inc exists (and is correct!), you can
# compile the code with:
#
#    > make
#
# which will compile a standard executable to be placed in the exec/ folder.
# For more control, specify additional options either in make.inc or on the
# command line. For example:
#
#   > make CONFIG=BXL
#
# will compile a Tantalus executable based on the BXL.py configuration file
# (look in the configs/ folder). By default, executables at the end will be
# named Tantalus.$(CONFIG).exe and be placed in the $(EXECDIR) configured
# in make.inc.
#
# An alternative is to not use a make.inc file but directly pass the Makefile the
# extension of the make.inc file in the make_include folder as in the following
# example:
#
#   > make INCLUDE=gnu-serial
#
# which uses the options specified in make_include/make.inc.gnu-serial.
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# The Makefile requires the following to be set. This can either be done inside
# make.inc or from the command line.
#
# Directory structure (will get created if they don't exist)
# EXECDIR  : directory for storage of the final executables
# SRCDIR   : source code for the mean-field executable as processed by Hephaestos
# FAMDIR   : source code for the FAM executable as processed by Hephaestos
# OBJDIR   : storage for intermediate object files
# MODDIR   : storage for final module files
#
# PHYSICS options
# CONFIG   : name of configutation file in the config/ folder
# CALCTYPE : NUCLEI or PASTA
# FAM      : 0/1, depending on whether you want to compile the FAM executable
#
# Compiler options
# CXX      : compiler invokation to be used.
# OPTFLAGS : optimisation compiler flags
# CXXFLAGS : other compiler flags
# MODFLAG  :
# PREPFLAG : syntax to invoke the preprocessor
# PRE      : steps to do before compilation, leave empty to NOT run Hephaestos
# USE_MPI  : 0 => no MPI
#            1 => MPI (note: your compiler should be MPI-capable to use this)
#
# Linking options
# USE_HDF5           : whether to offer HDF5 support, yes(1) or no (0)
# HDF5_LIB           : linking statemetns for the HDF5 library
# LINEAR_ALGEBRA_LIB : linking statements for (Sca)LAPACK and BLAS
#
#
# There is are two specific options that are not set in the example make.inc 
# files, EXENAME and FAMNAME. These allow you to override the default naming 
# scheme of the executables.
#-------------------------------------------------------------------------------
# Acknowledgment:
#   the organisation of this Makefile as well as a bunch of options are
#   inspired by the Makefile of the HFBTHO v4 code, see the repository of
#   P. Marević et al., Computer Physics Communications 276, 108367 (2022).
#-------------------------------------------------------------------------------

INCLUDE=make.inc
INCLUDEFILE=$(INCLUDE)
################################################################################
# Go get the compilation settings
ifeq ("$(wildcard $(INCLUDEFILE))","")
  INCLUDE_ALT := make_include/make.inc.$(INCLUDEFILE)
ifeq ("$(wildcard $(INCLUDE_ALT))","")
  $(error $(INCLUDE_ALT) was not found; Tantalus cannot be compiled.)
else
  INCLUDEFILE=$(INCLUDE_ALT)
endif
endif
include $(INCLUDEFILE)
################################################################################

################################################################################
# Physics details (modify as you want)
################################################################################
CONFIG  := default
EXENAME := Tantalus.$(CONFIG).exe
FAMNAME := fam.$(CONFIG).exe
################################################################################
# Compilation details (this section should be modified as you see fit)
################################################################################

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Type of calculation aimed at: 
#    'NUCLEI':  finite nuclei
#    'PASTA' : nuclear pasta
#
# These categories are the main ones and their choice determines
#    PASTA        = 0/1 
#    USE_Periodic = 0/1 
#    DENSUM       = 0/1 
# but it is possible that the user might want to set these flags differently...
CALCTYPE :=NUCLEI
ifeq ($(CALCTYPE),PASTA)
  PASTA        := 1
  USE_Periodic := 1
else
ifeq ($(CALCTYPE),NUCLEI)
  PASTA        := 0
  USE_Periodic := 0
else
  $(error "Invalid value of CALCTYPE. $(CALCTYPE)")
endif
endif
ifeq ($(USE_Periodic),1)
  DENSUM       := 1
else
  DENSUM       := 0 
endif
# .... but not all combination are meaningful!
ifeq ($(CALCTYPE),PASTA)
ifeq ($(USE_Periodic),0)
    $(error "Periodic boundary conditions should be enforced when attempting pasta calculations.")
endif
endif
ifeq ($(USE_Periodic),1)
ifeq ($(DENSUM),0)
    $(error "Periodic boundary conditions require setting DENSUM = 1")
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
MF_SRC  :=   compilation.f90 geninfo.f90 timing.f90 constants.f90 diaglib.f90
ifeq ($(USE_HDF5),1)
MF_SRC    +=   hdf5_auxiliary.f90
endif
MF_SRC  +=   sphericalharmonics.f90 folding.f90 nil8.f90
MF_SRC  +=   derivatives.f90 vectors.f90 precondition.f90 wavefunctions.f90
MF_SRC  +=   pairingcutoffs.f90 parameterization.f90
MF_SRC  +=   pairing_strengths.f90 basis_transform.f90 hartree-fock.f90 BCS.f90
MF_SRC  +=   HFB_gradient.f90 HFB_direct.f90 HFB.f90
MF_SRC  +=   pairing.f90 densities.f90 moments.f90
MF_SRC  +=   coulomb.f90 cranking.f90 momentsofinertia.f90 transform.f90
MF_SRC  +=   fission_MOI.f90 functional.f90 evolution.f90 scfiteration.f90
MF_SRC  +=   IO_aux.f90 IO_wf.f90
MF_SRC  +=   IO.f90 convergence.f90 printing.f90 version.f90 tantalus.f90
MF_SRC  +=   run_single.f90
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# FAM source files
FAM_SRC :=   compilation.f90 geninfo.f90 timing.f90 constants.f90 diaglib.f90
ifeq ($(USE_HDF5),1)
FAM_SRC    +=   hdf5_auxiliary.f90
endif
FAM_SRC +=   sphericalharmonics.f90 folding.f90 nil8.f90
FAM_SRC +=   derivatives.f90 vectors.f90 precondition.f90 wavefunctions.f90
FAM_SRC +=   pairingcutoffs.f90 parameterization.f90
FAM_SRC +=   pairing_strengths.f90 basis_transform.f90 hartree-fock.f90 BCS.f90
FAM_SRC +=   HFB_gradient.f90 HFB_direct.f90 HFB.f90
FAM_SRC +=   pairing.f90 densities.f90 moments.f90
FAM_SRC +=   coulomb.f90 cranking.f90 momentsofinertia.f90 transform.f90
FAM_SRC +=   fission_MOI.f90 functional.f90 evolution.f90 scfiteration.f90
FAM_SRC +=   fam_gmres.f90 fam.f90 fam_testing.f90
FAM_SRC  +=  IO_aux.f90 IO_wf.f90
FAM_SRC +=   IO.f90 convergence.f90 printing.f90 version.f90 tantalus.f90
FAM_SRC +=   fam_run.f90
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Nilsson source files
NIL_SRC := compilation.f90 geninfo.f90 timing.f90 derivatives.f90 nil8.f90
NIL_SRC += wavefunctions.f90 gennilsson.f90

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Translate source files into object files
MF_OBJ  :=  $(patsubst %.f90,$(MF_OBJ_DIR)/%.o ,$(MF_SRC))
FAM_OBJ :=  $(patsubst %.f90,$(FAM_OBJ_DIR)/%.o,$(FAM_SRC))
NIL_OBJ :=  $(patsubst %.f90,$(MF_OBJ_DIR)/%.o ,$(NIL_SRC))

################################################################################
# Explicit precompilation steps:
#
# 1) Check for the existence of all directories
# 2) Run Hephaestos to preprocess the entire code
# 3) Get version information from git
# 4) Get compiler information
# 5) set the version and compiler info in the source code
#
# Note that these steps depend somewhat on whether the mf or fam executable
# is being built.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
PRE_MF  :=  $(MF_SRC_DIR)/  $(MF_OBJ_DIR)/  $(MF_MOD_DIR)/ $(EXEC_DIR)/
PRE_MF  +=  run_heph_mf  getgitinfo getcompilerinfo setversioninfo_mf
PRE_FAM :=  $(FAM_SRC_DIR)/ $(FAM_OBJ_DIR)/ $(FAM_MOD_DIR)/ $(EXEC_DIR)/
PRE_FAM +=  run_heph_fam getgitinfo getcompilerinfo setversioninfo_fam
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Internal (to the compiler) preprocessing directives
#    -cpp      => explicitly enable preprocessing
#    -DUSE_MPI => enable (1) or disable (0) MPI (see above)

DIRECTIVES :=-DUSE_MPI=$(USE_MPI) -DUSE_Periodic=$(USE_Periodic) -DPASTA=$(PASTA)
DIRECTIVES +=-DDENSUM=$(DENSUM) -DDEBUG_LEVEL=$(DEBUG_LEVEL) -DUSE_HDF5=$(USE_HDF5)
PREPROCESSOR :=  $(PREPFLAG) $(DIRECTIVES)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

################################################################################
# Recipes (This section should NOT be modified in principle)
################################################################################
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
all: mf fam

mf: $(PRE_MF) $(MF_OBJ)
	$(CXX) $(OPTFLAGS) $(CXXFLAGS) $(PREPROCESSOR) -o $@ $(MF_OBJ) $(LINEAR_ALGEBRA_LIB) $(HDF5_LIB)
	mv mf exec/$(EXENAME)

fam: $(PRE_FAM) $(FAM_OBJ)
	$(CXX) $(OPTFLAGS) $(CXXFLAGS) $(PREPROCESSOR) -o $@ $(FAM_OBJ) $(LINEAR_ALGEBRA_LIB) $(HDF5_LIB)
	mv fam exec/$(FAMNAME)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Creation of required directories
$(EXEC_DIR)/:
	mkdir -p $(EXEC_DIR)/

# mean-field code
$(MF_OBJ_DIR)/:
	mkdir -p $(MF_OBJ_DIR)/

$(MF_MOD_DIR)/:
	mkdir -p $(MF_MOD_DIR)/

$(MF_SRC_DIR)/:
	mkdir -p $(MF_SRC_DIR)/

# FAM code
$(FAM_OBJ_DIR)/:
	mkdir -p $(FAM_OBJ_DIR)/

$(FAM_MOD_DIR)/:
	mkdir -p $(FAM_MOD_DIR)/

$(FAM_SRC_DIR)/:
	mkdir -p $(FAM_SRC_DIR)/

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Hephaestos runs
run_heph_mf:
# Run Hephaestos with the correct configuration file and information from
# the Makefile for the mean-field source code
	python3 Hephaestos.py $(CONFIG) mf $(DENSUM)

run_heph_fam:
# Run Hephaestos with the correct configuration file and information from
# the Makefile for the FAM source code
	python3 Hephaestos.py $(CONFIG) fam $(DENSUM)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
gen_nilsson: $(PRE) $(PRE_NIL) $(NIL_OBJ)
	$(CXX) $(OPTFLAGS) $(CXXFLAGS) $(PREPROCESSOR) -o $@ $(NIL_OBJ) $(LINEAR_ALGEBRA_LIB) $(HDF5_LIB)
	mv gen_nilsson exec/gen_nilsson.exe

$(MF_OBJ_DIR)/%.o : $(MF_SRC_DIR)/%.f90 | $(MF_MOD_DIR)/ $(MF_MOD_DIR)/ exec/
	$(CXX)  $(OPTFLAGS) $(CXXFLAGS) $(PREPROCESSOR) $(MODFLAG) $(MF_MOD_DIR) -c  $< -o $@ $(HDF5_LIB)

$(FAM_OBJ_DIR)/%.o : $(FAM_SRC_DIR)/%.f90 | $(FAM_OBJ_DIR)/ $(FAM_MOD_DIR)/ exec/
	$(CXX)  $(OPTFLAGS) $(CXXFLAGS) $(PREPROCESSOR) $(MODFLAG) $(FAM_MOD_DIR) -c  $< -o $@ $(HDF5_LIB)

setversioninfo_mf:
# Copy the git information into the main code, so it can be printed
	@sed -i.bak 's~VTAG~"${GIT_INFO5}"~'     $(MF_SRC_DIR)/version.f90
	@sed -i.bak 's~VERSION1~"${GIT_INFO1}"~' $(MF_SRC_DIR)/version.f90
	@sed -i.bak 's~VERSION2~"${GIT_INFO2}"~' $(MF_SRC_DIR)/version.f90
	@sed -i.bak 's~VERSION3~"${GIT_INFO3}"~' $(MF_SRC_DIR)/version.f90
# Copy the compiler information
	@sed -i.bak 's!COMPCOMP!"${COMPVERSION}"!' $(MF_SRC_DIR)/version.f90
# The above command uses '!' as sed delimiter, because Ubuntu sometimes uses ~ for kernel versions
	@sed -i.bak 's/CFLAGS/"${CXXFLAGS}"/'      $(MF_SRC_DIR)/version.f90
	@sed -i.bak 's/OPTFLAGS/"${OPTFLAGS}"/'    $(MF_SRC_DIR)/version.f90
	@rm $(MF_SRC_DIR)/version.f90.bak

setversioninfo_fam:
# Copy the git information into the main code, so it can be printed
	@sed -i.bak 's~VTAG~"${GIT_INFO5}"~'     $(FAM_SRC_DIR)/version.f90
	@sed -i.bak 's~VERSION1~"${GIT_INFO1}"~' $(FAM_SRC_DIR)/version.f90
	@sed -i.bak 's~VERSION2~"${GIT_INFO2}"~' $(FAM_SRC_DIR)/version.f90
	@sed -i.bak 's~VERSION3~"${GIT_INFO3}"~' $(FAM_SRC_DIR)/version.f90
# Copy the compiler information
	@sed -i.bak 's!COMPCOMP!"${COMPVERSION}"!' $(FAM_SRC_DIR)/version.f90
# The above command uses '!' as sed delimiter, because Ubuntu sometimes uses ~ for kernel versions
	@sed -i.bak 's/CFLAGS/"${CXXFLAGS}"/'      $(FAM_SRC_DIR)/version.f90
	@sed -i.bak 's/OPTFLAGS/"${OPTFLAGS}"/'    $(FAM_SRC_DIR)/version.f90
	@rm $(FAM_SRC_DIR)/version.f90.bak


getgitinfo:
# Get information from 'git show', to see what kind of build this is.
	$(eval GIT_INFO1=$(shell git show   | grep 'commit ' | head -1))
	$(eval GIT_INFO2=$(shell git show   | grep 'Author:' | head -1))
	$(eval GIT_INFO3=$(shell git show   | grep 'Date:'   | head -1))
# 	$(eval GIT_INFO4=$(shell git branch | grep '*'       | head -1 | cut -c2- ))
	$(eval GIT_INFO5=$(shell git describe --tags --always ))
# 	echo $(GIT_INFO5)
getcompilerinfo:
# Get information from 'CXX --version'
	$(eval COMPVERSION=$(shell $(CXX) --version | head -1))

cp_nil:
	cp src_orig/gennilsson.f90 $(SRCDIR)/gennilsson.f90

cp_fam:	
	cp src_orig/fam.f90 $(SRCDIR)/fam.f90

clean:
	rm  -f $(MF_OBJ_DIR)/*.o
	rm  -f $(MF_MOD_DIR)/*.mod
	rm  -f $(FAM_OBJ_DIR)/*.o
	rm  -f $(FAM_MOD_DIR)/*.mod
################################################################################

