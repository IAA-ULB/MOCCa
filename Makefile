OBJDIR :=   obj
SRCDIR :=   src
MODDIR :=   mod

TARGET :=   Tantalus.exe
SRC    :=   compilation.f90 geninfo.f90 constants.f90 sphericalharmonics.f90    
SRC    +=   diag.f90 nil8.f90 derivatives.f90 precondition.f90 wavefunctions.f90
SRC    +=   pairingcutoffs.f90 hartree-fock.f90 BCS.f90 HFB.f90 pairing.f90
SRC    +=   densities.f90 moments.f90 parameterization.f90 coulomb.f90 
SRC    +=   functional.f90 evolution.f90 scfiteration.f90 IO.f90 
SRC    +=   printing.f90 tantalus.version.f90

################################################################################
# Compiler details
CXX      :=  gfortran
CXXFLAGS := -O3 -J$(MODDIR) -Wall

################################################################################
# Precompilation instructions
PRE    :=  run_heph getgitinfo setversioninfo 
OBJ    :=  $(patsubst %.f90,$(OBJDIR)/%.o,$(SRC))
# Default functional is an NLO one
FUNC   :=  NLO.func

EXENAME:= Tantalus.$(FUNC).exe

LIBS = -llapack  -lblas

################################################################################
# Recipes
$(TARGET): $(PRE) $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJ) $(LIBS)
	mv $(TARGET) exec/$(EXENAME)

run_heph:
  # Run Hephaestos with the correct .func file to generate the source code in 
  # src folder.
	python Hephaestos.py $(FUNC) 

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
	$(eval GIT_INFO1=$(shell git show | grep 'commit '))
	$(eval GIT_INFO2=$(shell git show | grep 'Author:'))
	$(eval GIT_INFO3=$(shell git show | grep 'Date:'))

