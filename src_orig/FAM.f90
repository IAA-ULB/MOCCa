module fam

 !==============================================================================
 ! ________ _______  _        _ _________ _______  _                 _______
 !(  _____/(  ___  )( (      ) |\__   __/(  ___  )( \      |\     /|(  ____ \
 !| (      | (   ) ||  \    /  |   ) (   | (   ) || (      | )   ( || (    \/
 !| |___   | (___) ||   \  /   |   | |   | (___) || |      | |   | || (_____
 !|  ___)  |  ___  || (\ \/ /) |   | |   |  ___  || |      | |   | |(_____  )
 !| |      | (   ) || | \  / | |   | |   | (   ) || |      | |   | |      ) |
 !| |      | )   ( || )  \/  ( |   | |   | )   ( || (____/\| (___) |/\____) |
 !(_/      |/     \||/        \)   )_(   |/     \|(_______/(_______)\_______)
 !
 !  Copyright W. Ryssens & P. Demol
 !
 !------------------------------------------------------------------------------
 ! A FAM-QRPA implementation to complement MOCCa.
 !==============================================================================

  use densities
  use moments
  use fission_MOI

  implicit none

  !-----------------------------------------------------------------------------
  ! Define some FAM parameters
  real(KIND=dp) :: omega ! frequency of the perturbing field 
  real(KIND=dp) :: smear = 1.0_dp  ! complex smearing parameter, default 1.0 MeV
  real(KIND=dp) :: eta = 1.0e-3_dp ! small parameter entering derivatives, 
                                   ! default 10^-3
  !-----------------------------------------------------------------------------
  ! FAM amplitudes X, Y
  real(KIND=dp), allocatable :: X(:,:) ! forward amplitudes HF basis
  !                               | '-> sp index : hole
  !                               '-> sp index : particle
  real(KIND=dp), allocatable :: Y(:,:) ! backward amplitudes HF basis
  !                               | '-> sp index : hole
  !                               '-> sp index : particle
  !-----------------------------------------------------------------------------
  ! perturbed densities
  real(KIND=dp), allocatable :: drho(:,:)   ! preturbed normal density
  real(KIND=dp), allocatable :: dkappa(:,:) ! preturbed pairing density
  real(KIND=dp), allocatable :: dR(:,:)     ! preturbed generalised density
  !-----------------------------------------------------------------------------
  ! perturbed Hamiltonian
  real(KIND=dp), allocatable :: dH(:,:,:) ! perturbed Hamiltonian in HF basis
  !                                | | '-> 1: ph block, 2: hp block 
  !                                | '-> sp index : hole
  !                                '-> sp index : particle
  !-----------------------------------------------------------------------------
  ! external field
  real(KIND=dp), allocatable :: F(:,:,:)  ! perturbing external field in HF basis
  !                               | | '-> 1: ph block, 2: hp block 
  !                               | '-> sp index : hole
  !                               '-> sp index : particle
  integer :: l, m ! Principal and magnetic quantum number of the multipole moment

  ! Do we need more identifiers for electric vs mqgnetic and isovector 
  ! vs isoscalar


  contains

  subroutine inifam
    implicit none
    !---------------------------------------------------------------------------
    ! subroutine to initialise the FAM matrices, i.e.
    !---------------------------------------------------------------------------

    real(KIND=dp), allocatable :: SolidHarmHF(:,:)
    integer :: p, h
    real(KIND=dp) :: occ_h, occ_p, e_h, e_p
    logical :: ImPart

    print *, "Initialise FAM matrices" 

    allocate(drho(nwt,nwt))
    allocate(dkappa(nwt,nwt))
    allocate(dR(2*nwt,2*nwt))

    allocate(dH(nwt,nwt,2)) 
    allocate(F(nwt,nwt,2))


    allocate(X(nwt,nwt)) 
    allocate(Y(nwt,nwt))


    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Get the solid harmonics Q_lm(i,j) = < i | r^l Y_lm | j > expressed 
    ! in HF basis. 
    
    allocate(SolidHarmHF(nwt,nwt)) 

    ! Set external field to E2, hardcoded for now
    l = 2
    m = 0
    ImPart = .false. ! real (.false.) , imaginary (.true.) TBD later
   
    ! Calling a function in fission_MOI.f90
    SolidHarmHF = Qlm_spme(l, m, ImPart)

    ! TODO: write a general transfromation routine from the mesh to any 
    !       single-particle basis

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Define the external field F by selecting the particle-hole and 
    ! hole-particle subblocks of SolidHarmHF by multiplying by their 
    ! occupation, i.e. diagonal elements of rho in the canonical basis

    F = 0

    do h = 1, nwt
      occ_h = rho_can(h)
      if(occ_h < 1d-6) cycle
      do p = 1, nwt
        occ_p = 2.0 - rho_can(p) 
        ! degeneracy 2.0 must reduced if further symmetries are broken
        if(occ_p < 1d-6) cycle
        F(p,h,1) = occ_p * occ_h * SolidHarmHF(p,h) ! ph block F20(p,h)
        F(p,h,2) = occ_p * occ_h * SolidHarmHF(h,p) ! hp block F02(p,h)
      enddo
    enddo

    deallocate(SolidHarmHF)


    ! This can be improved by some element-wise products occ^T @ SolidHarmHF @ occ

    ! Note to future self: for QFAM this will be replaced by a transformation 
    ! to the qp basis. 

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! initialise the RPA amplitudes 

    X = -F(:,:,1)
    Y = -F(:,:,2)

    ! normalise with energy denominator
    do h = 1, nwt
      occ_h = rho_can(h)
      e_h = spenergies(h) 
      if(occ_h < 1d-6) cycle
      do p = 1, nwt
        occ_p = 2.0 - rho_can(p)
        e_p = spenergies(p) 
        if(occ_p < 1d-6) cycle
        X(p,h) = X(p,h) / (e_p - e_h - omega )
        Y(p,h) = Y(p,h) / (e_p - e_h + omega ) 
        ! print *, p, h, e_p, e_h, X(p,h), Y(p,h), F(p,h,1), F(p,h,2)
      enddo
    enddo

    ! TBD: X and Y are fully equivalent at this point. Is this expected?
    ! is this a consequency due to < i | Qlm | j > = < j | Qlm | i > 

    ! TODO: complex smearing

  end subroutine inifam

end module fam

program run_FAM

  use compilation
  use IO
  use Tantalus
  use fam

  implicit none

  ! integer :: ifail ! Future dev: required for HFB

  100 format &
     &  (/,8x,' _____________________________________________________________', &
     &   /,8x,'|                                                             |', &
     &   /,8x,'| MOCCa v2.0 =                                                |', &
     &   /,8x,'|                                                             |', &
     &   /,8x,'|    #####  ##   #     # #####   ##   #      #    #  ####     |', &
     &   /,8x,'|    #     #  #  ##   ##   #    #  #  #      #    # #         |', &
     &   /,8x,'|    #### #    # # # # #   #   #    # #      #    #  ####     |', &
     &   /,8x,'|    #    ###### #  #  #   #   ###### #      #    #      #    |', &
     &   /,8x,'|    #    #    # #     #   #   #    # #      #    # #    #    |', &
     &   /,8x,'|    #    #    # #     #   #   #    # ######  ####   ####     |', &
     &   /,8x,'|                                                             |', &
     &   /,8x,'|  Copyright  P.-H. Heenen, M. Bender, W. Ryssens & P. Demol  |', &
     &   /,8x,'|_____________________________________________________________|')

  print *
  print 100

  !------------------------------------------------------------------------------
  ! starting all timers
  ! 
  ! -> This is necessary since subroutines below make use of the timers
  ! 
  call initialize_all_timers

  !-----------------------------------------------------------------------------
  ! Read input from STDIN
  ! 
  ! For FAMQRPA, the code should read in addition:
  ! 
  ! -  the type of perturbing operator/external field: E1, E2, M1, M2, ...
  !    and more complicated stuff when targetting beta-decay
  !    Important note: we will need to distinguish
  ! -  the frequency \omega of the perturbing field
  ! -  the 'size' of the perturbation to perform the finite differencing
  ! -  a smearing parameter to avoid discontinuities at the poles of the 
  !    response function
  ! 
  call ReadInput()

  !-----------------------------------------------------------------------------
  ! Initalize the matrices for performing derivatives on the mesh
  call inilag()

  !------------------------------------------------------------------------------
  ! Read all information from a wf file
  call ReadWavefunction()

  !------------------------------------------------------------------------------
  ! Print all relevant input gleaned from STDIN and the wf file.
  call PrintInput()

  ! Provide memory for the derivatives of the spwfs
  call allocate_memory_derivatives(PairingType)


  ! Future dev: required for HFB
  ! ifail = 0
  ! call SolvePairing(pairingscheme, ifail)

  ! Derive all single-particle wavefunctions on the mesh
  if(store_derivatives) call deriveHF()

  ! Calculate the initial densities and the charge density (separately)
  ! call densit(SaveRho=.false.)

  call construct_canonical_basis(rho_pairing,kappa_pairing,rho_can,kappa_can)
  Density = densit(rho_can, kappa_pairing)

  call ConstructChargeDensity(Density) ! PD: necessary? 

  ! Adopt the relevant quantities to the centre-of-mass of the nucleus ! PD: necessary? 
  call adapt_com(Density)  

  call CalculateMoments(Density) ! Recalculate because the COM might have changed.

  ! initialise perturbed matrices
  call inifam()

  print *, "Reached the end successfully" 

  ! end of one FAM calculation;

end program run_FAM