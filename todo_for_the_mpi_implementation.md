# ToDo for the MPI implementation

## Goals
0) compiler nil8
1) run nil8 in parallel
2) compile full code

## Other/Various
[ ] remove temperature module(s)

## File-wise list of things todo

[ ] geninfo.f90 
	- symmetry block sizes for each rank
	- mapping between MPI rank and spwf number
	- spwf total numbers
[ ] timing.f90 
	- should work, except for the printing
[x] constants.f90
[ ] sphericalharmonics.f90
[ ] folding.f90
[ ] nil8
	- build wfs on separate processes
[x] derivatives.f90
[ ] wavefunctions.f90
	- array declaration
	- loops in all subroutines
	- orthonormalization procedure change
[x] pairingcutoffs.f90
[ ] parameterization.f90
 	- reading for iam == 0
[x] pairingstrengths.f90
[ ] basis_transform.f90
	- adapt loops
[ ] HF.f90
	- adapt loops
[ ] BCS.f90
	- adapt loops
[ ] HFB_gradient.f90, HFB_direct.f90 and HFB.f90
	- adapt loops
[x] pairing.f90
[ ] densities.f90
	- making an all_reduce over MPI ranks
	- Hephaestos matching
[ ] moments.f90
	- Input by just one rank
[x] coulomb.f90
[x] cranking.f90
[ ] momentsofinertia.f90	
	- change loops
	- allreduce (?)
[ ] transform.f90
	- change loops
[x] functional.f90
[ ] fission_MOI.f90
	- change loops
[ ] evolution.f90
	- change loops
[ ] scfiteration.f90
	- Input by just one rank
[ ] IO.f90
	- Input/output by just one rank
	- wavefunction file writing/reading per rank
[ ] convergence.f90
	- printing by just one rank
[ ] printing.f90
	- everything by just one rank
[ ] tantalus.f90
	- printing by just one rank
