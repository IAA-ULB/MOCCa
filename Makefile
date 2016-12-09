OBJDIR :=   obj
SRCDIR :=   src
MODDIR :=   mod

TARGET :=   Tantalus.exe
SRC    :=   compilation.f90 geninfo.f90 derivatives.f90 wavefunctions.f90 tantalus.f90

CXX :=      gfortran
CXXFLAGS := -O3 -J$(MODDIR)

PRE    :=  run_heph
OBJ    :=  $(patsubst %.f90,$(OBJDIR)/%.o,$(SRC))


$(TARGET): $(PRE) $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJ) $(LIBS)

run_heph:
	python Hephaestos.py

clean:
	rm  -f $(OBJDIR)/*
	rm  -f $(MODDIR)/*

$(OBJDIR)/%.o : $(SRCDIR)/%.f90
	$(CXX) $(CXXFLAGS) -c  $< -o $@ 
