# Changelog of the development branch 

## v1.0.3-dev 

  -  Added spurious mode subtraction for translational modes (Zcom, Zmomentum)
  -  Refactored external field handling: get_external_field → get_qpme_op
  -  Enabled FAM calculations with BSkG3, BSkG4, BSkG5 parameterisation
  -  Added validation: effective charges now restricted to multipole/particle-number operators only
  -  Improved EWSR (Energy-Weighted Sum Rule) support for dipole (L=1) perturbations
  -  Added Pbreak_HFB_vector and Pbreak_HFB_matrix functions for HFB matrix restructuring
  -  New documentation: docs/user-guide/breaking-symmetries.md with comprehensive guidance
  -  Re-enabled full functionality of write_timeodd_densities and write_densities
  -  Makefile: better error handling (Hephaestos failures now halt compilation)
  -  Increased filename length limit from 80 to 120 characters
  -  Enhanced fam_pairing_zero_mode test with robust fitting strategy
  -  Added new fam_spurious_mode_subtraction test case
  -  Extended CI workflows to test BXL and BXL-N2LO functionals
  -  New integration test for spurious mode subtraction
  -  New config.md with detailed configuration file reference
  -  Added MathJAX support for LaTeX rendering in documentation
  -  12 new symmetry visualization figures (Boxes_*.png/.eps)
  -  Updated contributor list and dependencies in README
  - Corrected deprecated scipy.sph_harm → scipy.sph_harm_y usage 
  - Resolved OpenMP parallelization conflicts in FAMQRPA
  -  Bug fix in functional files: replaced Density% with R% for proper external vmicro calls
