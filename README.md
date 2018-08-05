
----------
# Tantalus & Hephaestos
Copyright W. Ryssens, P.H. Heenen & M. Bender
Sunday, 05. August 2018 04:56PM 

### What is Tantalus?

Tantalus is the future of the mean-field codes using a 3D Lagrange mesh representation. Its predecessors,  EV(B)8, CR8, EV4 as well as MOCCa, should all be replaced by Tantalus. 

### What is Hephaestos?

Hephaestos is a (collection of) Python modules/scripts that partially writes the FORTRAN source code for Tantalus, depending on type of functional as well as intrinsic symmetries chosen. 

---------


### Known Problems 
#### Hephaestos 
- [ ] **Vector products**: only one per term is admitted at this point in time. 
- [ ] **Densities identification:**  Hephaestos incorrectly determines the densities needed for these inputs

		- E_C_I_NkSk_C_NqNq_NmSm  
    	- E_C_I_NmSk_C_NqNq_NkSm 
    
    The code in that case decides to only calculate C_NqNq_NmSm, while the second term cannot be calculated withthat contraction.
    The input in reversed order seems to work
    
		- E_C_I_NmSk_C_NqNq_NkSm
		- E_C_I_NkSk_C_NqNq_NmSm
	   

