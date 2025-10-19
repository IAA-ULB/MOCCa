# The Lagrange-mesh method
| heading                                                           | page |
|:------------------------------------------------------------------|-----:|
| Abstract                                                          |    1 |
| 1. **Introduction**                                               |    1 |
| 2 **Presentation of the Lagrange-mesh method**                    |    4 |
| 2.1 Gauss quadrature                                              |    4 |
| 2.2 Lagrange functions                                            |    6 |
| 2.3 Lagrange-mesh equations for bound states                      |    8 |
| 2.4 Comments                                                      |   10 |
| 2.5 Exact and approximate matrix elements                         |   11 |
| 2.6 Translation, scaling, mapping, parity projection              |   12 |
| 2.7 Regularization of a singularity                               |   14 |
| 2.8 Multidimensional Lagrange meshes                              |   17 |
| 2.9 Mean values and transition matrix elements                    |   18 |
| 3 **Explicit expressions for various Lagrange meshes**            |   19 |
| 3.1 Lagrange meshes from orthogonal polynomials                   |   19 |
| 3.1.1 Summary of properties of orthogonal polynomials             |   19 |
| 3.1.2 Lagrange functions                                          |   20 |
| 3.1.3 Exact matrix elements                                       |   21 |
| 3.1.4 Classical orthogonal polynomials                            |   21 |
| 3.2 Lagrange-Hermite mesh                                         |   22 |
| 3.3.2 Lagrange-Laguerre functions                                 |   25 |
| 3.3.4 Laguerre mesh regularized by x for α≥ 0                     |   28 |
| 3.3.5 Laguerre mesh regularized by x3/2 for α≥ 0                  |   29 |
| 3.3.6 Modified Laguerre mesh in x2                                |   30 |
| 3.3.7 Modified Laguerre mesh in x2 regularized by x               |   31 |
| 3.4 Lagrange-Legendre meshes                                      |   32 |
| 3.4.1 Generalized Legendre polynomials                            |   32 |
| 3.4.2 Lagrange-Legendre functions                                 |   33 |
| 3.4.3 Values and matrix elements                                  |   33 |
| 3.4.4 Legendre mesh regularized by √1− x2                         |   34 |
| 3.4.5 Shifted Legendre mesh on (0,1) regularized by x             |   35 |
| 3.4.6 Shifted Legendre mesh on (0,1) regularized by x3/2          |   37 |
| 3.4.7 Shifted Legendre mesh on (0,1) regularized by x(1− x)       |   38 |
| 3.5 Lagrange-Jacobi meshes                                        |   39 |
| 3.5.1 Jacobi polynomials                                          |   39 |
| 3.5.2 Lagrange-Jacobi functions                                   |   40 |
| 3.5.3 Values and matrix elements                                  |   40 |
| 3.5.4 α= 0 Lagrange-Jacobi mesh over (0,1)                        |   41 |
| 3.5.5 Regularized α= 0 Lagrange-Jacobi mesh over (0,1)            |   43 |
| 3.6 Lagrange meshes based on non-classical orthogonal polynomials |   44 |
| 3.6.1 Non-classical orthogonal polynomials                        |   44 |
| 3.6.2 Symmetric Lagrange-Gaussian mesh on (−∞,+∞)                 |   45 |
| 3.6.3 Asymmetric Lagrange-Gaussian mesh on (−∞,+∞)                |   48 |
| 3.6.4 Lagrange-Gaussian mesh on (0,∞)                             |   50 |
| 3.7 Lagrange meshes based on periodic functions                   |   51 |                   
| 3.7.1 Lagrange-Fourier mesh                                       |   51 |
| 3.7.2 First sine meshes                                           |   53 |
| 3.7.4 Cosine meshes                                               |   58 |
| 3.7.5 Cardinal sine (sinc) mesh                                   |   60 |
| 3.8 Numerical considerations                                      |   62 |
| 4 **Comparison with related methods**                             |   62 |
| 4.1 Historical aspects                                            |   62 |
| 4.2 Discrete-variable representation                              |   64 |
| 4.3 Quadrature discretization method                              |   66 |
| 4.4 Sinc and Fourier methods                                      |   67 |
| 4.5 Miscellaneous                                                 |   68 |
| 5 **Potential and two-body bound states**                         |   69 |
| 5.1 On the choice of a Lagrange mesh                              |   69 |
| 5.2 Harmonic oscillators                                          |   69 |
| 5.3 Morse potential                                               |   73 |
| 5.4 Hydrogen atom                                                 |   76 |
| 5.4.1 Energies and wave functions                                 |   76 |
| 5.4.2 Exactness in the regularized Laguerre-mesh method           |   79 |
| 5.4.3 Static polarizabilities                                     |   80 |
| 5.5 Confined hydrogen atom                                        |   82 |
| 5.6 Hydrogen atom in a strong magnetic field                      |   85 |
| 5.7 Lagrange-mesh method in momentum space                        |   87 |
| 5.8 One-dimensional Dirac oscillator                              |   90 |
| 6 **Two-body continuum**                                          |  101 | 
| 6.1 R-matrix method on a Lagrange mesh                            |  101 |
| 6.2 Eﬀective-range expansion                                      |  103 |
| 6.3 R matrix at high orbital momenta                              |  104 |
| 6.4 Multichannel R matrix                                         |  105 |
| 6.5 Strength functions                                            |  107 |
| 6.6 Complex scaling                                               |  109 |
| 7 **Three-body bound states**                                     |  113 |
| 7.1 Choice of coordinate system                                   |  113 |
| 7.2 Lagrange-mesh method in perimetric coordinates                |  114 |
| 7.3 Unconfined and confined helium atom                           |  116 |
| 7.4 Hydrogen molecular ion                                        |  119 |
| 7.5 Antiprotonic helium atom                                      |  120 |
| 7.6 Lagrange-mesh method in hyperspherical coordinates            |  123 |
| 7.7 Particles in a hyperradial potential                          |  125 |
| 7.8 Helium trimer                                                 |  127 |
| 7.9 Two-neutron halo nuclei                                       |  129 |
| 8 **Miscellaneous applications**                                  |  132 |
| 8.1 Real and imaginary time propagation                           |  132 |
| 8.2 Forced harmonic oscillator                                    |  135 |
| 8.3 Bose-Einstein condensates                                     |  136 |
| 8.4 Hartree-Fock calculations with a finite-range nuclear force   |  139 |
| 8.5 Translations and rotations on a Lagrange mesh                 |  142 |
| 9 **Conclusion and outlook**                                      |  146 |
| Acknowledgments                                                   |  148 |
| Appendix A: Proof of relation (2.68)                              |  148 |
| Appendix B: Summation formulas                                    |  149 |
| Appendix C: Polarizabilities                                      |  150 |
| Appendix D: Dirac hydrogenic atom with two mesh points            |  151 |
| Appendix E: System with separable imaginary part                  |  153 |
| References                                                        |  153 |

