
exe: test_real_example_nonrandom.f90 
	mpiifx -I$(EBROOTELPA)/include/elpa-2024.05.001/modules/ \
	       -lmkl_blacs_intelmpi_lp64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lelpa -o exe \
		   test_real_example_nonrandom.f90

clean:
	rm ./exe