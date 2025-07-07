
exe: src.F90
	mpiifx -I/data/antwerpen/201/vsc20170/buildingElpa/elpa/_install-intel/include/elpa-2025.01.002/modules \
		-L/data/antwerpen/201/vsc20170/buildingElpa/elpa/_install-intel/lib \
		-lmkl_blacs_intelmpi_lp64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lelpa -o exe \
		src.F90 

