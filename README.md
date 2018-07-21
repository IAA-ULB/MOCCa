
--------------------------------------------------------------
# Tantalus & Hephaestos
Copyright W. Ryssens & M. Bender
Saturday, 21. July 2018 10:00PM 

--------------------------------------------------------------

## Tantalus
#### What is Tantalus?

Tantalus is the future of the mean-field codes using a 3D Lagrange mesh representation. Its
predecessors,  EV(B)8, CR8, EV4 and MOCCa, should all be replaced by Tantalus. 

### Currently supports

* Hartree-Fock calculations 
* EV8-like symmetries 
 
### Known Bugs
-  None that I know of at the moment.
 
## Hephaestos
---
### What is Hephaestos?

Hephaestos is a (collection of) Python modules that partially writes the FORTRAN source code for Tantalus, depending on type of functional as well as intrinsic symmetries chosen. 

The generation of source code depends on a given list of terms that make up an energy density functional, Hephaestos writes the Tantalus source code that

* Calculates the mean-field densities needed (possibly contracted)
*  Reads the input with respect to constants in the parameterization
* Calculates the coupling constants associated with the various functional term
* Calculates and integrates the various energy densities
* Calculates the mean-fields associated with the densities
* Calculates the action of the single-particle Hamiltonian on the single-particle wave-functions.

Plus a lot of decisions making that depends on the type of functional considered.

The second major function of Hephaestos is to adapt the Tantalus source code based on the intrinsic symmetry assumptions.  This is however not yet supported at all. 

### Current status

* Correct treatment for two-body Skyrme functionals (without pairing densities) up to N3LO level. 
* Intelligent input for functionals and parameterizations

### Known Problems (Hephaestos)

####Features
- [ ] **Vector products**: only one per term is admitted at this point in time. 

####Known Bugs
- [ ] **Densities identification:**  Hephaestos incorrectly determines the densities needed for these inputs   

    - E_C_I_NkSk_C_NqNq_NmSm  
    - E_C_I_NmSk_C_NqNq_NkSm    
  The code in that case decides to only calculate C_NqNq_NmSm, while the second term cannot be calculated with that contraction.
    The input in reversed order seems to work
    - E_C_I_NmSk_C_NqNq_NkSm
    - E_C_I_NkSk_C_NqNq_NmSm
   

