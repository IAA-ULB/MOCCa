This example was found [here](https://docs.rcc.fsu.edu/software/elpa/).

Fhe original source file from the web page above is saved as `test_real_example.orig.f90`. 

We made small changes in `test_real_example.f90`. Mainly, comments as to our understanding of the program, print statements, variable renamings, set `nev = na` (computing all eigenvalues), ...

Finally, in  `test_real_example_nonrandom.f90`, we adapt the code to use a non-random matrix where `a(i,i) = 1.5` and `a(i,j) = 1/sqrt(|i-j|)`, in order to verify the correctness of the solution. ()

That seems to work too:

ELPA solution:
``` 
> ./cl.sh 8 2
mpiifx -I/apps/antwerpen/skylake/rocky9/ELPA/2024.05.001-intel-2024a/include/elpa-2024.05.001/modules/ \
       -lmkl_blacs_intelmpi_lp64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lelpa -o exe \
           test_real_example_nonrandom.f90
Source file : experiments/elpa/test_real_example/test_real_example.f90
matrix size     :       8
blocking factor :       2
#eigenvalues.   :       8
[0/4] -> 2x2 (0,0)
| Past BLACS_Gridinfo.
| Past scalapack descriptor setup.
[1/4] -> 2x2 (1,0)
[2/4] -> 2x2 (0,1)
[3/4] -> 2x2 (1,1)
0  1.50  1.00  0.50  0.45
0  1.00  1.50  0.58  0.50
0  0.50  0.58  1.50  1.00
0  0.45  0.50  1.00  1.50
| rank 0 : Symmetric matrix block has been set up.
1  0.71  1.00  0.71  0.58
1  0.58  0.71  1.00  0.71
1  0.41  0.45  0.71  1.00
1  0.38  0.41  0.58  0.71
| rank 1 : Symmetric matrix block has been set up.
2  0.71  0.58  0.41  0.38
2  1.00  0.71  0.45  0.41
2  0.71  1.00  0.71  0.58
2  0.58  0.71  1.00  0.71
| rank 2 : Symmetric matrix block has been set up.
3  1.50  1.00  0.50  0.45
3  1.00  1.50  0.58  0.50
3  0.50  0.58  1.50  1.00
3  0.45  0.50  1.00  1.50
| rank 3 : Symmetric matrix block has been set up.
| Entering one-step ELPA solver 
| ... 
 Eigenvalues:
           1  0.306609541843836     
           2  0.358272781626069     
           3  0.450363891062035     
           4  0.601606913192053     
           5  0.832799343069226     
           6   1.24365754057420     
           7   2.00769896784700     
           8   6.19899102078558     
| One-step ELPA solver complete.
| cleaning up.
```
numpy.linalg.eigh solution
```
> ./test_real_example_nonrandom.py 8 2
[1.5        1.         0.70710678 0.57735027 0.5        0.4472136  0.40824829 0.37796447]
[1.         1.5        1.         0.70710678 0.57735027 0.5        0.4472136  0.40824829]
[0.70710678 1.         1.5        1.         0.70710678 0.57735027 0.5        0.4472136 ]
[0.57735027 0.70710678 1.         1.5        1.         0.70710678 0.57735027 0.5       ]
[0.5        0.57735027 0.70710678 1.         1.5        1.         0.70710678 0.57735027]
[0.4472136  0.5        0.57735027 0.70710678 1.         1.5        1.         0.70710678]
[0.40824829 0.4472136  0.5        0.57735027 0.70710678 1.         1.5        1.        ]
[0.37796447 0.40824829 0.4472136  0.5        0.57735027 0.70710678 1.         1.5       ]
Eigenvalues:
0 0.3066095515041723
1 0.3582726976185496
2 0.4503637822206513
3 0.6016068130970025
4 0.832799273606104
5 1.243657532882371
6 2.0076990549702614
7 6.198991294100891
```

pyscalapack.pdsyev solution (experiments/elpa/test_real_example/test_real_example_nonrandom_pyscalapack.py)
```
Eigenvalues:
0 0.30660955150417296
1 0.35827269761854946
2 0.45036378222065193
3 0.6016068130970024
4 0.8327992736061045
5 1.2436575328823696
6 2.0076990549702622
7 6.198991294100887
norm_diff=4.041654044610298e-15
```

All three approaches give (approximately) the same solution, so we're good.

After installing `pyelpa` that too gives the same solution:

```
*** elpa.eigenvectors solution ***
Eigenvalues:
0 0.30660955150417274
1 0.3582726976185492
2 0.4503637822206516
3 0.6016068130970021
4 0.8327992736061037
5 1.2436575328823705
6 2.007699054970261
7 6.19899129410089
norm_diff=1.3299526536404359e-15
```

Perhaps it would be nice to have a common frontend for both pyscalapack and pyelpa, so that we can run the same code with both backends and easily produce timings.