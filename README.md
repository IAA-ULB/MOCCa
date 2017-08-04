
--------------------------------------------------------------
# Tantalus & Hephaestos
--------------------------------------------------------------
Copyright W. Ryssens & M. Bender

--------------------------------------------------------------

## Tantalus

--------------------------------------------------------------

#### What is Tantalus?

Tantalus is the future of the mean-field codes using Lagrange mesh representation as developed originally in the Bonche-Flocard-Heenen collaboration. It is the successor to EV(B)8, CR8, EV4 and MOCCa, all rolled into one.

In short, it is a FORTRAN program that iteratively optimizes a given Skyrme-type energy density functional for a given number of protons and neutrons, possibly subject to constraints on the shape of the various nuclear densities. 

--------------------------------------------------------------

### Current status

* 
* 



### Known Problems 

1. Spin-orbit term currently crashes Hephaestos, even though it worked in the past.
2. Vector couplings: only one per term is admitted at this point in time. 
--------------------------------------------------------------
### Planned features

--------------------------------------------------------------

--------------------------------------------------------------

## Hephaestos

--------------------------------------------------------------
### What is Hephaestos?

Hephaestos is a (collection of) Python script(s) that partially writes the FORTRAN source code for Tantalus, depending on type of functional as well as intrinsic symmetries chosen. 

The generation of source code depends on a given list of terms that make up an energy density functional, Hephaestos writes the Tantalus source code that

* Calculates the mean-field densities needed (possibly contracted)
* Calculates the coupling constants associated with the various functional term
* Calculates and integrates the various energy densities
* Calculates the mean-fields associated with the densities
* Calculates the action of the single-particle Hamiltonian on the single-particle wave-functions.

plus some more decision making that depends on the type of functional considered, such as the number of derivatives of the single-particle wavefunctions to calculate. 

The second major function of Hephaestos is to adapt the Tantalus source code based on the intrinsic symmetry assumptions. 

--------------------------------------------------------------

### Current status

* Correct treatment for two-body Skyrme functionals (without pairing densities) up to N3LO level. This was tested self-consistently with MOCCa and WHISKY up to N2LO level and non-selfconsistently with WHISKY up to N3LO level.
* 



### Known Problems 

1. Spin-orbit term currently crashes Hephaestos, even though it worked in the past.
2. Vector couplings: only one per term is admitted at this point in time. 
--------------------------------------------------------------
### Planned features

* Determination of symmetries at compile-time
* 
--------------------------------------------------------------
### Other ideas

* Use the density ring structure to eliminate specific contractions of densities, specifically at N3LO level.
* 
--------------------------------------------------------------

