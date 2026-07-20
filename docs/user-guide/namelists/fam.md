# The FAM namelist

 This namelist contains the parameters for solving linear response via FAM. 


| Parameter               | Type      | Default     | Explanation                                                  |
|-------------------------|-----------|------------:|--------------------------------------------------------------|
| **FAM frequency**       |           |             |                                                              |
| `omega`                 | real*8    | 0.0   [MeV] | Single omega frequency of the external perturbation          |
| `omega_min`             | real*8    | 0.0   [MeV] | Minimal value of frequency range                             |
| `omega_max`             | real*8    | 30.0  [MeV] | Maximal value of frequency range                             |
| `omega_step`            | real*8    | 1.0   [MeV] | Step size used within the frequency range                    |
| `smear`                 | real*8    | 1.0   [MeV] | Complex smearing of the frequecy                             |
| **External field**      |           |             |                                                              |
| `operator_type`         | string    | 'multipole' | Type of the perturbing external field                        |
| `l`                     | integer   | -1          | Multipolarity of the external field                          |
| `m`                     | integer   | -1          | Projection (orientation) of the external field               |
| `eff_charge_n`          | real*8    | +1.0   [e]  | Effective charge for neutrons                                |
| `eff_charge_p`          | real*8    | +1.0   [e]  | Effective charge for protons                                 |
| `remove_spurious`       | boolean   | true        | Subtract the spurious mode ?                                 |
| **Convergence**         |           |             |                                                              |
| `mixingscheme`          | integer   | 0           | Convergence strategy : 0 = GMRES, 1 = linear mixing          |
| `maxiter`               | integer   | 100         | maximal number of FAM iterations                             |
| `maxhist`               | integer   | 30          | maximal history size of GMRES                                |
| `fam_precision`         | real*8    | 1e-5        | convergence criterion of FAM iterations                      |
| `fam_lin_mix`           | real*8    | 0.3         | linear mixing coefficient  (lower=slower)                    |



Some further remarks of general interest:

 1. When passing `omega_min` `omega_max` and `omega_step`, the code will ignore `omega` and solve FAM for a range of frequencies omega = `omega_min`  + k * `omega_step` < `omega_max` for k=0,...
 2. `operator_type` sets the external field. Current options are : 
    - `multipole`        : multipole moment Q_{lm} [default], i.e. solid harmonics f_lm = r^l Y_lm, except for the monopole operator f_00 = r^2
    - `particle number`  : particle number operator N
    - `Zcom`             : center-of mass z-coordinate
    - `Zmomentum`        : center-of mass z momentum
 3. By setting `eff_charge_n` = `eff_charge_p`, the external field corresponds to an isoscalar operator, while setting `eff_charge_n` = - `eff_charge_p` defines an isovector perturbation. Note that one often takes neutron and proton charges as Z/A and N/A respectively. 
 4. The convergence measure `fam_precision` means different things depending on the mixing scheme. In case one uses GMRES (`mixingscheme` = 0 [default]) `fam_precision` sets the convergence criterion for the normalised GMRES residual ||b - A(x)||/||b||. When using linear mixing 
(`mixingscheme` = 1), then convergence is met once ||FAM(X) - X|| / ||X|| < `fam_precision` and ||FAM(Y) - Y|| / ||Y|| < `fam_precision`
 5. Note that the filename used to output the strength function **must** be passed as `famfile` in the `IO` namelist. 