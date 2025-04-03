# The nucleus namelist

 This namelist contains the parameters determining 
 
 1. the proton and neutron numbers in our many-body wavefunction.
 2. the treatment of the chemical potentials of both species.
 3. the convergence criteria.
 4. a few miscellaneous options to organise the whole of the calculation.



| Parameter               | Type      | Default     | Explanation                                                 |
|-------------------------|-----------|-------------|-------------------------------------------------------------|
| **Nucleus**             |           |             |                                                             |
| `neutrons`              | real*8    | 12.0        |  The targetted number of neutrons                           |
| `protons`               | real*8    | 12.0        |  The targetted number of protons                            |
|                         |           |             |                                                             |
| **Chemical potentials** |           |             |                                                             |
| `fixfermi`              | logical   | .false.     |  If .true., freeze the chemical potential.                  |
| `mun`                   | real*8    | -10e8       |  The neutron chemical potential if `fixfermi` is set.       |
| `mup`                   | real*8    | -10e8       |  The proton chemical potential if `fixfermi` is set.        |
|                         |           |             |                                                             |
| **Convergence criteria**|           |             |  [More information](../convergence.md).             |
| `energy_prec`           | real*8    | 1e-9        |  Energy criterion                                           | 
| `moment_prec`           | real*8    | 1e-3        |  Multipole moment criterion                                 | 
| `disp_prec`             | real*8    | 1e-5        |  Spwf dispersion criterion                                  |
| `gradient_prec`         | real*8    | 1e+0        |  Energy gradient criterion                                  |
| `fermi_prec`            | real*8    | 1e-3        |  Fermi energy criterion                                     | 
| `angmom_prec`           | real*8    | 1e-3        |  Angular momentum criterion                                 |
|                         |           |             |                                                             |
| **Miscellaneous**       |           |             |                                                             |
| `balancing_strategy`    | integer   |  1          |  Not operational currently.                                 | 
| `store_derivatives`     | logical   |  .true.     |  If .true., store the derivatives of the spwfs.             |
|                         |           |             |  If .false., recalculate them on the fly.                   |



Some further remarks of general interest:

 1. The `neutrons` and `protons` are declared as reals since users might want to enforce non-vanishing pairing by setting them to non-integer values.
 2. The meaning of the parameters related to convergence criteria is described in detail [here](../convergence.md).
 3. `balancing_strategy` is a parameter that will determine the (future) comportment of the code when utilising MPI parallelism.
 4. The value of `store_derivatives` has big implications for the codes resource useage. If `.false.`, the code has minimal memory footprint but requires some additional CPU time.
    If set to `.true.`, the codes memory footpritn will be much larger - up to factor 5!. Particularly for production runs, the CPU time cost of setting `store_derivatives = .false.` is not worth the increased memory cost.
