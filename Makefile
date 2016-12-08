OBJDIR :=   obj
SRCDIR :=   src
MODDIR :=   mod

TARGET :=   Tantalus.exe
SRC    :=   compilation.f90 geninfo.f90 wavefunctions.f90 tantalus.f90

CXX :=      gfortran
CXXFLAGS := -O3 -J$(MODDIR)

OBJ    :=  $(patsubst %.f90,$(OBJDIR)/%.o,$(SRC))


$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)
	
clean:
	rm  -f $(OBJDIR)/*
	rm  -f $(MODDIR)/*

$(OBJDIR)/%.o : $(SRCDIR)/%.f90
	$(CXX) $(CXXFLAGS) -c  $< -o $@ 
