OBJDIR :=   obj
SRCDIR :=   src
MODDIR :=   mod

# Public src directory
SRCPUBLIC := ../tantalus_public/src/

TARGET :=   Tantalus.exe
SRC    :=   compilation.f90 geninfo.f90 constants.f90 sphericalharmonics.f90
SRC    +=   folding.f90   
SRC    +=   diag.f90 nil8.f90 derivatives.f90 precondition.f90 wavefunctions.f90
SRC    +=   pairingcutoffs.f90 parameterization.f90 hartree-fock.f90 BCS.f90 
SRC    +=   HFB.f90 pairing.f90 densities.f90 moments.f90 coulomb.f90 
SRC    +=   momentsofinertia.f90  functional.f90 evolution.f90 scfiteration.f90 
SRC    +=   IO.f90 temperature_projection.f90 printing.f90 tantalus.version.f90

SINGLE_SRC = $(SRC) run_single.f90
MPI_SRC    = $(SRC) run_mpi.f90

################################################################################
# Compiler details
CXX      :=  gfortran

ifeq ($(CXX),gfortran)
#	CXXFLAGS := -J$(MODDIR) -Wall -fbacktrace -g3
	CXXFLAGS := -O3 -J$(MODDIR) -Wall 
else ifeq ($(CXX),ifort)
	CXXFLAGS := -O3 -assume realloc-lhs -assume byterecl -no-wrap-margin -module $(MODDIR)
endif

################################################################################
# Precompilation instructions
PRE         :=  run_heph getgitinfo setversioninfo 
SINGLE_OBJ  :=  $(patsubst %.f90,$(OBJDIR)/%.o,$(SINGLE_SRC))
MPI_OBJ     :=  $(patsubst %.f90,$(OBJDIR)/%.o,$(MPI_SRC))

# Default functional is an NLO one
FUNC   :=  NLO.func

mpi:    EXENAME:= Tantalus.$(FUNC).mpi.exe
single: EXENAME:= Tantalus.$(FUNC).exe

LIBS   := -llapack -lblas

################################################################################
# Recipes
single: $(PRE) $(SINGLE_OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $(SINGLE_OBJ) $(LIBS)
	mv single exec/$(EXENAME)

mpi: $(PRE) $(MPI_OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $(MPI_OBJ) $(LIBS)
	mv mpi exec/$(EXENAME)

run_heph:
  # Run Hephaestos with the correct .func file to generate the source code in 
  # src folder.
	python Hephaestos.py $(FUNC) 

clean:
	rm  -f $(OBJDIR)/*.o
	rm  -f $(MODDIR)/*.mod

$(OBJDIR)/%.o : $(SRCDIR)/%.f90
	$(CXX) $(CXXFLAGS) -c  $< -o $@ 

genpublic: run_heph getgitinfo setversioninfo
  # Generate a public version of Tantalus with Hephaestos in src_public
	cp $(SRCDIR)/*.f90 $(SRCPUBLIC)/

setversioninfo:
  # Copy the git information into the main code, so it can be printed
	cp $(SRCDIR)/tantalus.f90 $(SRCDIR)/tantalus.version.f90
	sed -i.bak 's/VERSION1/"${GIT_INFO1}"/' $(SRCDIR)/tantalus.version.f90 
	sed -i.bak 's/VERSION2/"${GIT_INFO2}"/' $(SRCDIR)/tantalus.version.f90 
	sed -i.bak 's/VERSION3/"${GIT_INFO3}"/' $(SRCDIR)/tantalus.version.f90 
	rm $(SRCDIR)/tantalus.version.f90.bak

getgitinfo:
	# Get information from 'git show', to see what kind of build this is.
	$(eval GIT_INFO1=$(shell git show | grep 'commit '))
	$(eval GIT_INFO2=$(shell git show | grep 'Author:'))
	$(eval GIT_INFO3=$(shell git show | grep 'Date:'))

