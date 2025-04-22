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
  ! Define some FAM parameters
  real(KIND=dp) :: omega ! frequency of the perturbing field 
  real(KIND=dp) :: smear = 1.0_dp ! complex smearing parameter, default 1.0 MeV
  real(KIND=dp) :: eta = 1.0e-3_dp ! ! small parameter entering derivatives, default 10^-3
  !------------------------------------------------------------------------------
  ! FAM amplitudes X, Y
  real(KIND=dp), allocatable :: X(:,:) ! forward amplitudes HF basis
  real(KIND=dp), allocatable :: Y(:,:) ! backward amplitudes HF basis
  !------------------------------------------------------------------------------
  ! perturbed densities
  real(KIND=dp), allocatable :: drho(:,:)   ! preturbed normal density
  real(KIND=dp), allocatable :: dkappa(:,:) ! preturbed pairing density
  real(KIND=dp), allocatable :: dR(:,:)     ! preturbed generalised density
  !------------------------------------------------------------------------------
  ! perturbed Hamiltonian
  real(KIND=dp), allocatable :: dH(:,:,:) ! perturbed Hamiltonian in HF basis
  !                               | '-'--> qp index 
  !                               -> 1: H^20, 2: H^02 
  !------------------------------------------------------------------------------
  ! external field
  real(KIND=dp), allocatable :: F(:,:,:)  ! perturbing external field in HF basis
  !                               | '-'--> qp index 
  !                               -> 1: F^20, 2: F^02 
  integer :: l, m ! Principal and magnetic quantum number of the multipole moment
  ! Do we need more identifiers for electric vs mqgnetic and isovector vs isoscalar


  contains

  subroutine inifam
    implicit none
    !------------------------------------------------------------------------------
    ! subroutine to initialise the FAM matrices, i.e.
    !------------------------------------------------------------------------------

    print *, "Initialise FAM matrices" 


    allocate(drho(nwt,nwt))
    allocate(dkappa(nwt,nwt))
    allocate(dR(2*nwt,2*nwt))

    allocate(dH(2,nwt,nwt)) 


    ! Set external field to E2, hardcoded for now
    l = 2
    m = 0
    allocate(F(2,nwt,nwt))


    allocate(SpherHarmMesh(nx,ny,nz,0:maxmoment,0:maxmoment,2))

    ! Generate the spherical harmonics Y^l_m(x,y,z) up to l=maxmoment (default:10)
    ! This could be reduced to just calling the necessary one
    call GenSphericalHarmonics(maxmoment,nx,ny,nz,meshx,meshy,meshz,           & 
    &                          SpherHarmMesh,quantisationaxis,secondaryaxis)

    allocate(drho(nx,ny,nz))

    ImPart = 0 ! 0 = Real, 1 = Im, => TBD later

    select case (perturbationtype)
      case (4)
        drho(:,:,:) = SpherHarmMesh(:,:,:,2,0,ImPart+1)
      case default
        print *, "only E2 (perturbationtype = 4) implemented so far"
    end select
    
    deallocate(spherharmmesh)
    

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