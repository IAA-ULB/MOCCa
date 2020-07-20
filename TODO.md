# Tantalus & Hephaestos
> Copyright W. Ryssens, P.H. Heenen & M. Bender


## Small to do's:

* We should redo the Moments.f90 module to emphasize the role of the Lagrange multiplier more. After all, *all* constraints can be written as  
  $$ R = E - \lambda \langle \hat{O} \rangle.$$
  Different constraints (quadratic, linear, augmented, .....) are simply 
  different ways to calculate $\lambda$.


## Known Bugs:

### Angular momentum corretion and HFB diagonalization
* The calculation of $$\langle J^2 \rangle$$ and the Belyaev moment of inertia  are unreliable, as figured out by G. Scamps and E. Olsen.   

  a) First, the openmp parallelization is buggy.
    
>Dans la subroutine calcJ2andBelyaev_HFB, il faut remplacer :  
>      !$OMP PARALLEL  
>par  
>      !$OMP PARALLEL shared(jx, jy, jz) private(i,ii,si,jj,j)  
>pour que cela marche mieux.   

 b)  I should replace the diagonalization routines by Lapack ones so that the diagonalization of the HFB problem is reproducible.
 
> J'ai creusé un peu plus cette histoire de résultats différents selon le  compilateur.  
> Pour résumer, j'ai des résultats différents si :
>-je mets qopenmp ou pas,
>-j'utilise gfortran ou ifort
>-je change les options d'optimisation de ifort

>Le point rassurant, c'est qu'ils convergent bien a terme vers la même valeur.
>-La première valeur du fichier d'output qui change c'est J^2.
>-Une déviation entre les calcul a lieu juste après la première >diagonalisation de la matrice HFB.
>-L'hamiltonien en entrée est identique, mais les fct d'ondes sont différentes, elle varient par exemple d'un signe, mais elle sont bien vecteur propre (je l'ai vérifié numériquement).

### Densities identification:
   Hephaestos incorrectly determines the densities needed for these inputs
	 	
		- E_C_I_NkSk_C_NqNq_NmSm  
		- E_C_I_NmSk_C_NqNq_NkSm 
		
   The code in that case decides to only calculate C_NqNq_NmSm, while the 
   second term cannot be calculated withthat contraction. The input in reversed 
   order seems to work
	
		- E_C_I_NmSk_C_NqNq_NkSm
		- E_C_I_NkSk_C_NqNq_NmSm
		

### *assume_locality* - flag

I need to investigate the behaviour of the assume_locality flag. Something is fishy, but I do not yet know enough.

### Canonical basis at finite temperature

The canonical basis does not exist in general at finite temperature, as the transformation C in the Bogoliubov transformation no longer drops out.  
In the language of Ring and Schuck, 
$$
U = D \bar{U} C, V  = D^* \bar{V} C.
$$

The density for instance becomes, 
$$
\rho = U f U^dagger + V* (1-f) V^T
     =  D \bar{U}   C       f   C^{\dagger} \bar{U}^{\dagger} D^{\dagger}
     +  D \bar{V}^* C^*  (1-f)  C^{T} \bar{V}^{T} D^{\dagger} \, .
$$

When all the $f$ are zero, all mentions of the matrix $C$ drop out and we 
obtain

$$
\rho = D \bar{V}^* \bar{V} D^{\dagger}\, ,
$$
%
meaning that rho is diagonal in the basis defined by D. However, if the 
$f$ are different from zero (or rather, not all equal) then C doesn't drop out. 
The Bogoliubov transformation is still "at its core" a BCS transformation, but 
\rho and \kappa can not trivially be diagonalized as in the T=0 case. 

I did not realize this when implementing finite T calculations in Tantalus, 
and currently the code naively diagonalizes rho. This is (I think at the moment
of writing) not strictly wrong: Tantalus uses this diagonalization only to 
efficiently sum the ordinary densities. For the pairing densities, Tantalus 
does the summation of the pairing densities in the HF basis. So, while the 
comments and intent was wrong, the final result and implementation should be 
right.
