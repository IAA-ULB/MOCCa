# Tantalus & Hephaestos
> Copyright W. Ryssens, P.H. Heenen & M. Bender

---- 
## Known bugs
---- 
### Tantalus/Fortran 

*  Canonical basis at finite temperature

    The canonical basis does not exist in general at finite temperature, as the transformation C in the Bogoliubov transformation no longer drops out.  
    In the language of Ring and Schuck, 

		U = D \bar{U} C, V  = D^* \bar{V} C.

    The density for instance becomes, 

		\rho = U f U^{\dagger} + V* (1-f) V^T
         =  D \bar{U}   C       f   C^{\dagger} \bar{U}^{\dagger} D^{\dagger}
         +  D \bar{V}^* C^*  (1-f)  C^{T} \bar{V}^{T} D^{\dagger} \, .
    
    When all the f are zero, all mentions of the matrix $C$ drop out and we 
    obtain

	    \rho = D \bar{V}^* \bar{V} D^{\dagger}\, ,

    meaning that rho is diagonal in the basis defined by D. However, if the f are different from zero (or rather, not all equal) then C doesn't drop out. The Bogoliubov transformation is still "at its core" a BCS transformation, but \rho and \kappa can not trivially be diagonalized as in the T=0 case. 

    I did not realize this when implementing finite T calculations in Tantalus, and currently the code naively diagonalizes rho. This is (I think at the moment of writing) not strictly wrong: Tantalus uses this diagonalization only to  efficiently sum the ordinary densities. For the pairing densities, Tantalus does the summation of the pairing densities in the HF basis. So, while the comments and intent was wrong, the final result and implementation should be right.
    
*  Angular momentum correction and HFB diagonalization
   The calculation of $$\langle J^2 \rangle$$ and the Belyaev moment of inertia  are unreliable, as figured out by G. Scamps and E. Olsen.   

   First, the openmp parallelization is buggy.

    >Dans la subroutine calcJ2andBelyaev_HFB, il faut remplacer :  
    >      !$OMP PARALLEL  
    >par  
    >      !$OMP PARALLEL shared(jx, jy, jz) private(i,ii,si,jj,j)  
    >pour que cela marche mieux.   



### Hephaestos/Python 

* Densities identification:
   Hephaestos incorrectly determines the densities needed for these inputs
	 	
		- E_C_I_NkSk_C_NqNq_NmSm  
		- E_C_I_NmSk_C_NqNq_NkSm 
		
   The code in that case decides to only calculate C_NqNq_NmSm, while the 
   second term cannot be calculated with that contraction. The input in reversed 
   order seems to work
	
		- E_C_I_NmSk_C_NqNq_NkSm
		- E_C_I_NkSk_C_NqNq_NmSm
		
* *assume_locality* - flag
  I need to investigate the behaviour of the assume_locality flag. Something is fishy, but I do not yet know enough. It definitely needs to be disabled by default though.

* 

### Overall 

---- 
## Wishlist/ToDos
---- 

### Tantalus/Fortran

 1. Streamline the error codes introduced by E. Olsen in the HFB solver. Extend to HF and BCS solvers.
 2. Implement lowest-qp style blocking also for T-broken calculations: the code does not distinguish between signatures yet.
 3. Test (and almost certainly correct) T-broken finite-temperature implementation
    * This requires for sure looking into P. Fanto's Pfaffian formula for the finite-T particle number projection.
 4. We should redo the Moments.f90 module to emphasize the role of the Lagrange multiplier more. After all, *all* constraints can be written as  
	  	R = E - \lambda \langle \hat{O} \rangle.
    Different constraints (quadratic, linear, augmented, .....) are simply different ways to calculate \lambda.


### Hephaestos/Python 
 1. The calculation of the action of the single-particle hamiltonian would be more efficient if derivatives of potentials were calculated.
	Example: right now the calculation of the part of the sp. hamiltonian associated with D_I_N proceeds as (modulo signs and factor i)
    	h psi \sim sum_m 0.5 * nabla_m ( F_I_Nm psi ) + 0.5 * F_I_Nm Nabla_m psi
    where the derivatives are calculated "on the fly" in the routine. However, it would be more efficient to precompute the derivative of the potential and use instead
		h psi \sim sum_m 0.5 * ( nabla_m F_I_Nm)  psi  + F_I_Nm Nabla_m psi

### Overall

 1. The implementation of the calculation of the pairing gaps is not optimal.
 	* This requires elimination of the delta_action routine.
 	* See email conversation with M. Bender on 07/11/2020.


