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

I need to investigate the behaviour of the assume_locality flag. Something    if fishy, but I do not yet know enough.