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

  implicit none

  !------------------------------------------------------------------------------
  ! The FAM frequency of perturbing field
  real(KIND=dp) :: omega = 0.5
  !------------------------------------------------------------------------------
  ! Type of perturbing operator
  ! (0) : E0
  ! (1) : M0
  ! (2) : E1
  ! (3) : M1
  ! (4) : E2
  ! (5) : M2
  ! (6) : E3
  ! (7) : M3
  integer :: perturbationtype = 4
  !------------------------------------------------------------------------------

  contains

  subroutine inifam
    implicit none

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
  call densit(SaveRho=.false.)

  call ConstructChargeDensity(ChargeDensity) ! PD: necessary? 

  ! Adopt the relevant quantities to the centre-of-mass of the nucleus ! PD: necessary? 
  ! call adapt_com()

  call CalculateMoments()

  ! initialise perturbed matrices
  call inifam()

  print *, "Reached the end successfully" 

  ! end of one FAM calculation;

end program run_FAM