OBJDIR :=   obj
SRCDIR :=   src
MODDIR :=   mod

TARGET :=   Tantalus.exe
SRC    :=   compilation.f90 geninfo.f90 constants.f90 sphericalharmonics.f90 diag.f90 nil8.f90 derivatives.f90 precondition.f90 wavefunctions.f90 hartree-fock.f90 densities.f90 coulomb.f90 functional.f90 evolution.f90 IO.f90 tantalus.f90

################################################################################
# Compiler details
CXX      :=  gfortran
CXXFLAGS := -O3 -J$(MODDIR) -Wall

################################################################################
# Precompilation instructions
PRE    :=  run_heph
OBJ    :=  $(patsubst %.f90,$(OBJDIR)/%.o,$(SRC))
# Default functional is an NLO one
FUNC   :=  NLO.func

EXENAME:= Tantalus.$(FUNC).exe

LIBS = -llapack  -lblas


$(TARGET): $(PRE) $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJ) $(LIBS)
	mv $(TARGET) exec/$(EXENAME)

run_heph:
  # Run Hephaestos with the correct .func file
	python Hephaestos.py $(FUNC) #> hephaestos.out

clean:
	rm  -f $(OBJDIR)/*
	rm  -f $(MODDIR)/*

$(OBJDIR)/%.o : $(SRCDIR)/%.f90
	$(CXX) $(CXXFLAGS) -c  $< -o $@ 
