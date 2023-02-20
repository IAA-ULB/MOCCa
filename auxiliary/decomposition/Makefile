
OBJDIR :=   obj
SRCDIR :=   src
MODDIR :=   mod
EXENAME:= decomposition.exe
SRC    :=   compilation.f90 sphericalharmonics.f90 interpolation.f90 decomposition.f90
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

OBJ  :=  $(patsubst %.f90,$(OBJDIR)/%.o,$(SRC))

single: $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJ) $(LIBS)
	mv single exec/$(EXENAME)

exec/:
	mkdir -p exec/

$(OBJDIR)/:
	mkdir -p  obj/

$(MODDIR)/:
	mkdir -p  mod/

$(OBJDIR)/%.o : $(SRCDIR)/%.f90 | $(OBJDIR)/ $(MODDIR)/ 
	$(CXX) $(CXXFLAGS) -c  $< -o $@ 

clean:
	rm  -f $(OBJDIR)/*.o
	rm  -f $(MODDIR)/*.mod
