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

## v0.2.0 - Quasiparticle tagging
  - Inclusion of "tagging" strategy to block quasiparticles 
  - Reenabling of the calculation of the divergence of the spin-orbit current density.
 
## v0.1.0 - Start of the documentation
  - Initial stages of the living documentation.
  - Started the automated testing on the master branch.
  - Moved the automated testing infrastructure from Gitlab CI/CD to Github actions.
  - Started the semantic versioning of code versions. 
