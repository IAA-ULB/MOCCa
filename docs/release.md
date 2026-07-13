# Release notes

## Rationale of the versioning 

 We aim to adhere to the [semantic versioning](https://www.semver.org) standard: code 
 versions get labels as 
 ```
     vA.B.C
 ```
where

 - A indicates the **major version number**, which gets incremented on large-scale changes happen 
   that (perhaps) affect backwards compatibility.
 - B is the **minor version number**, which gets incremented when new features release.
 - C is the **patch version number** that gets incremented on small implementations, documentation changes, small bugfixes, ... 

This page will list release notes for all releases that increment at least the minor version, but 
not the patch version.

## v1.0.0 MOCCa's debut
  With a minor cleaning, some documentation improvements and an opensource licence (AGPLv3), MOCCa is now ready to face the world!

  - v0.4.2: Improved FAMQRPA, new functionality and several bugfixed (also on the mean-field level!)
  - v0.4.1: Restored efficiency for constrained calculations that took a hit in v0.3.0 (W. Ryssens)


## v0.4.0 - FAMQRPA 
  Major new feature: full implementation of linear response functionality with pairing (P. Demol) 
  Notable changes from v0.3.1-0.3.9
  - v0.3.9: exact folding with Gaussian charge factors (L. Gonzalez)
  - v0.3.8: new building framework, significantly enhanced compilation speed when developing (W. Ryssens)
  - v0.3.7: a collection of small changes (multiple authors)
  - v0.3.6: sum rules infrastructure (P. Demol)
  - v0.3.5: a collection of small bugfixes (W. Ryssens)
  - v0.3.4: official jump off for the MOCCaPy project 
  - v0.3.3: FAMRPA operational (P. Demol)
  - v0.3.1: Implementation of vibrational correction from a collective Hamiltonian (L. Gonzalez)

## v0.3.0 - long-awaited merges 
  - merging of the `difficult-merge' pasta AND fam branches into master
  - support for reading and writing HDF5 files
  - calculation of the moments of inertia is now performed by default
  - automated testing expanded
  - standards set for documentation and style
  - new compilation procedure

## v0.2.0 - Quasiparticle tagging
  - Inclusion of "tagging" strategy to block quasiparticles 
  - Reenabling of the calculation of the divergence of the spin-orbit current density.
 
## v0.1.0 - Start of the documentation
  - Initial stages of the living documentation.
  - Started the automated testing on the master branch.
  - Moved the automated testing infrastructure from Gitlab CI/CD to Github actions.
  - Started the semantic versioning of code versions. 
