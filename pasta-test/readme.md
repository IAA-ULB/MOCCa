### a test-case for the ELPA development. 

The original files are in `pasta-test.tar.gz`

Idn bijlage een voorbeeldje: 3x3x3 spherische clusters op relatief hoge
dichtheid. Er zijn vier bestanden

  * run*sh: het runscriptje; geconfigureerd voor de setup van Nikolai,
    maar het zou je duidelijk moeten maken hoe de andere bestandjes aan
    MOCCa gevoederd worden.
  * data* : de input data voor MOCCa
  * out*  : dit is de STDOUT voor de code; die bevat eigenlijk alle
    fysieke informatie.
  * inp*pot: dit bevat onze initialisatie voor de code, bepaald door een
    semiklassieke benadering van het probleem.

Nog wat opmerkingen:

  * Zoals je kan zien in het runscript, kan dit voorbeeld runnen op 8
    nodes; Nikolai heeft het ook getest op 6 & 10 nodes.
  * De rekentijd is wel serieus als je evenveel iteraties wilt doen: ~9u
    voor 400 iteraties. Je kan echter zoveel iteraties doen als je
    budget hebt door met het keyword "maxiter" in het data* bestandje te
    spelen. Schatting: ~1.5 minuten per iteratie op 6 nodes.

The test needs to be run with Tantalus.BXL.exe, compile with HDF5 support.

We zetten `maxiter=10` (ipv 400 -> 9h rekentijd)

| nodes | output file                 | walltime[s] | Spwf evolution | Subspace rotation: Spwf transformation | Subspace rotation: Matrix diagonalisation |  Matrix elements of h |
| ----: | :-------------------------- | ----------: | -------------: | ------------------: | ---------------------: | --------------------: | 
|     4 | [link](./slurm-2144074.out) |   3933.2006 |.        10.27% |              60.82% |                 18.41% |                 5.17%
|     8 | [link](./slurm-2143960.out) |   2107.8961 |         12.98% |              51.46% |                 21.03% |                 7.66%
|    16 | [link](./slurm-2143961.out) |    986.5521 |         11.63% |              49.82% |                 21.19% |                 3.94%

1. subroutine `Cholesky_orthonormalisation` in [wavefunctions.f90](./) => orthonormalisatie voor MPI berekeningen; daar zitten grote PDGEMM, PDPQTRF en PDTRSM calls in.
2. subroutine apply_subspace_rotation in evolution.f90 => diagonalisatie van de single-particle hamiltoniaan in de subspace gegenereerd door de golffuncties in het geheugen; grote PDSYEV en PDGEMM calls.
3. function calc_sphamil in evolution.f90 => berekening van de matrix elementen van de single-particle hamiltoniaan in de subspace gegenereerd door de golffuncties in het geheugen; grote PDGEMM call.

ScaLAPACK calls to 
- [PDGEMM](https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-fortran/2023-0/p-gemm.html) Computes a scalar-matrix-matrix product and adds the result to a scalar-matrix product for distributed matrices. (MM stands for matrix-matrix)
- [PDPOTRF](https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-fortran/2023-0/p-potrf.html) Computes the Cholesky factorization of a symmetric (Hermitian) positive-definite distributed matrix.
- [PDTRSM](https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-fortran/2023-0/p-trsm.html) Solves a distributed matrix equation (one matrix operand is triangular).
- [PDSYEV](https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-fortran/2023-0/p-syev.html) Computes selected eigenvalues and eigenvectors of a symmetric matrix.


The `P` inf position 1 is for Parallel, I guess.
The `D` in position 2 stands for DOUBLE PRECISION.

#meshpoints : 30x30x30 = 27000
#spwfs neutron 