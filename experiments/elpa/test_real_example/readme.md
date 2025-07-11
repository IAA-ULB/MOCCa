This example was found [here](https://docs.rcc.fsu.edu/software/elpa/).

Fhe original source file from the web page above is saved as `test_real_example.orig.f90`. 

We made small changes in `test_real_example.f90`. Mainly, comments as to our understanding of the program, print statements, variable renamings, set `nev = na` (computing all eigenvalues), ...

Finally, in  `test_real_example_nonrandom.f90`, we adapt the code to use a non-random matrix where `a(i,i) = 1.5` and `a(i,j) = 1/sqrt(i+j)`, in order to verify the correctness of the solution.



