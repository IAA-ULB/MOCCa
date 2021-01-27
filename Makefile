#-------------------------------------------------------------------------------
# Make options included at this point
#
#   all:    Build both single and mpi executables.
#           Note that Hephaestos is only used once to generate the code.
#   single: Ordinary mode, one calculation. Corresponds to run_single.f90.
#   mpi:    (Possibly) multiple runs. Corresponds to run_mpi.f90.
#
#-------------------------------------------------------------------------------
#
# OPTIONS:
#                                                    DEFAULT
#  CXX      : compiler to use                        gfortran
#  CXXFLAGS : compiler flags                         -O3 -J$(MODDIR) -Wall
#  PRE      : steps to do before compilation         run Hephaestos
#  CONFIG   : name of configutation file in the      default
#             config/ folder
#  
# Note that config file needs to supply all the necessary information to 
# Hephaestos: a functional file as well as the information on the symmetries
# to conserve and the symmetries that can be expected on the input file.
#
#-------------------------------------------------------------------------------
# Executables at the end will be named
#  Tantalus.$(CONFIG).exe        => single
#  Tantalus.$(CONFIG).mpi.exe    => mpi
#-------------------------------------------------------------------------------

OBJDIR :=   obj
SRCDIR :=   src
MODDIR :=   mod

# Public src directory
SRCPUBLIC := ../tantalus_public/src/

TARGET :=   Tantalus.exe
SRC    :=   compilation.f90 geninfo.f90 timing.f90 constants.f90 
SRC    +=   sphericalharmonics.f90 folding.f90   
SRC    +=   nil8.f90 derivatives.f90 precondition.f90 wavefunctions.f90
SRC    +=   pairingcutoffs.f90 parameterization.f90 hartree-fock.f90 BCS.f90 
SRC    +=   HFB.f90 pairing.f90 densities.f90 moments.f90 coulomb.f90 
SRC    +=   cranking.f90 momentsofinertia.f90 transform.f90 functional.f90 
SRC    +=   evolution.f90  scfiteration.f90 
SRC    +=   IO.f90 temperature_projection.f90 convergence.f90 printing.f90 
SRC    +=   tantalus.version.f90

SINGLE_SRC = $(SRC) run_single.f90
MPI_SRC    = $(SRC) multirun_example.f90

################################################################################
# Compiler details
CXX      :=  gfortran
OPENMP   :=  

ifneq (,$(findstring gfortran,$(CXX)))
#	CXXFLAGS := -J$(MODDIR) -Wall -fbacktrace -g3
  OPENMP   := 
	CXXFLAGS := -O3 -J$(MODDIR) -Wall -Wno-uninitialized $(OPENMP) 
    LIBS   := -llapack -lblas
else ifeq ($(CXX),ifort)
  OPENMP   := 
	CXXFLAGS := -O3  $(OPENMP) -assume realloc-lhs -assume byterecl -no-wrap-margin -module $(MODDIR)
    LIBS   := -mkl
endif

################################################################################
# Precompilation instructions
PRE         :=  run_heph getgitinfo setversioninfo 
SINGLE_OBJ  :=  $(patsubst %.f90,$(OBJDIR)/%.o,$(SINGLE_SRC))
MPI_OBJ     :=  $(patsubst %.f90,$(OBJDIR)/%.o,$(MPI_SRC))

# Default configuration
CONFIG   :=  default

mpi:    EXENAME:= Tantalus.$(CONFIG).mpi.exe
single: EXENAME:= Tantalus.$(CONFIG).exe


################################################################################
# Recipes
all: single mpi

single: $(PRE) $(SINGLE_OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $(SINGLE_OBJ) $(LIBS)
	mv single exec/$(EXENAME)

mpi: $(PRE) $(MPI_OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $(MPI_OBJ) $(LIBS)
	mv mpi exec/$(EXENAME)

run_heph:
  # Run Hephaestos with the correct configuration file
	python3 Hephaestos.py $(CONFIG) 

clean:
	rm  -f $(OBJDIR)/*.o
	rm  -f $(MODDIR)/*.mod

$(OBJDIR)/%.o : $(SRCDIR)/%.f90
	$(CXX) $(CXXFLAGS) -c  $< -o $@ 

setversioninfo:
  # Copy the git information into the main code, so it can be printed
	cp $(SRCDIR)/tantalus.f90 $(SRCDIR)/tantalus.version.f90
	sed -i.bak 's/VERSION1/"${GIT_INFO1}"/' $(SRCDIR)/tantalus.version.f90 
	sed -i.bak 's/VERSION2/"${GIT_INFO2}"/' $(SRCDIR)/tantalus.version.f90 
	sed -i.bak 's/VERSION3/"${GIT_INFO3}"/' $(SRCDIR)/tantalus.version.f90 
	rm $(SRCDIR)/tantalus.version.f90.bak

getgitinfo:
	# Get information from 'git show', to see what kind of build this is.
	$(eval GIT_INFO1=$(shell git show | grep 'commit ' | head -1))
	$(eval GIT_INFO2=$(shell git show | grep 'Author:' | head -1))
	$(eval GIT_INFO3=$(shell git show | grep 'Date:'   | head -1))

################################################################################

