OBJDIR :=   obj
SRCDIR :=   src
MODDIR :=   mod

TARGET :=   Tantalus.exe
SRC    :=   compilation.f90 geninfo.f90 constants.f90 diag.f90 nil8.f90 derivatives.f90 wavefunctions.f90 hartree-fock.f90 densities.f90 functional.f90 IO.f90 tantalus.f90

CXX :=      gfortran
CXXFLAGS := -O3 -J$(MODDIR) -Wall

PRE    :=  run_heph
OBJ    :=  $(patsubst %.f90,$(OBJDIR)/%.o,$(SRC))


$(TARGET): $(PRE) $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJ) $(LIBS)

run_heph:
	python Hephaestos.py #> hephaestos.out

clean:
	rm  -f $(OBJDIR)/*
	rm  -f $(MODDIR)/*

$(OBJDIR)/%.o : $(SRCDIR)/%.f90
	$(CXX) $(CXXFLAGS) -c  $< -o $@ 
