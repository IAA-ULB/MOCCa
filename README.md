# MOCCa

[![License: AGPL v3](https://img.shields.io/badge/License-AGPL_v3-blue.svg)](https://www.gnu.org/licenses/agpl-3.0)

- describe licence
- ADD DOI INFO

The following is an extremely concise summary of the contents of this repository.

__MOCCa__ is a mean-field and linear response code using energy density functionals of the Skyrme type that represents the single-particle wavefunctions on a Lagrange mesh, written in FORTRAN.
__Hephaestos__ is a (collection of) Python modules/scripts that partially writes (parts of) the FORTRAN source code for Tantalus. 

Hephaestos, is in fact a code generator or preprocessor for Tantalus. For different choices regarding (a) types of functionals, (b) self-consistent symmetries and (c) type of many-body claculations, Hephaestos generates a version of Tantalus that is able to perform the calculations you wish to do. 

For anything more, please see the documentation. You can build it by typing 

    mkdocs serve

in your terminal and navigating to the adress it provides you with your favorite browser. Of course, you'll need a working installation of [mkdocs](https://www.mkdocs.org/) to do so. If you don't have that, you can also manually browse the Markdown documents in the `docs` directory.
