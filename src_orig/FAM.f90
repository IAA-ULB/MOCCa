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
  !                                | | '-> 1: dH^20, 2: dH^02 
  !                                '-'--> qp index 
  !------------------------------------------------------------------------------
  ! external field
  real(KIND=dp), allocatable :: F(:,:,:)  ! perturbing external field in HF basis
  !                               | | '-> 1: F^20, 2: F^02 
  !                               '-'--> qp index 
  integer :: l, m ! Principal and magnetic quantum number of the multipole moment
  ! Do we need more identifiers for electric vs mqgnetic and isovector vs isoscalar


  contains

  subroutine inifam
    implicit none
    !------------------------------------------------------------------------------
    ! subroutine to initialise the FAM matrices, i.e.
    !------------------------------------------------------------------------------

    real(KIND=dp), allocatable :: SpherHarmMesh(:,:,:,:,:,:)
    ! real(KIND=dp), allocatable :: harm_3D(:,:,:)
    real(KIND=dp) :: f_ph
    integer :: ImPart, p, h, i, j, k

    print *, "Initialise FAM matrices" 

    allocate(drho(nwt,nwt))
    allocate(dkappa(nwt,nwt))
    allocate(dR(2*nwt,2*nwt))

    allocate(dH(nwt,nwt,2)) 


    ! Set external field to E2, hardcoded for now
    l = 2
    m = 0
    allocate(F(nwt,nwt,2))


    allocate(SpherHarmMesh(nx,ny,nz,0:maxmoment,0:maxmoment,2))
    ! allocate(harm_3D(nx,ny,nz))

    ! Generate the spherical harmonics Y^l_m(x,y,z) up to l=maxmoment (default:10)
    ! This could be reduced to just calling the necessary one
    call GenSphericalHarmonics(maxmoment,nx,ny,nz,meshx,meshy,meshz,           & 
    &                          SpherHarmMesh,quantisationaxis,secondaryaxis)

    ! SpherHarmMesh(nx,ny,nz,l,m,2) = r^l Y_{lm}(x,y,z)
    !               -------      '-> (real, imaginary) 
    !                  '--> mesh coordinates (x, y, z)
    !
    ! /!\: double check with Wouter that the factor r^l is indeed included already
    !      such that they are in fact regular solid harmonics R_{lm} so that multipole
    !      moments are simply < Q_{lm} > = \int rho * R_{lm} .
    ! TBD: is it correct the imaginary part is zero when m=0 since 
    !          Y^l_m^dagger = (-1)^m * Y^l_-m 


    ImPart = 1 ! real (1) or imaginary (2) part => TBD later

    ! select case (perturbationtype)
    !   case (4)
    !     harm_3D(:,:,:) = SpherHarmMesh(:,:,:,l,m,ImPart) ! select l=2, m=0 component
    !   case default
    !     print *, "only E2 (perturbationtype = 4) implemented so far"
    ! end select

    ! transform SpherHarmMesh to HF basis
    ! is this basis transformation already implemented? Not in basis_transfrom.f90
    ! only some complicated (MPI) ones in wavefunctions.f90 
    ! if not, I think it should be something like this

    do p=1,nwt ! can be restricted to particles
      do h=1,nwt ! can be restricted to holes
        f_ph = 0.0_dp
        do k=1,nz 1 ! outerloop should be rightmost index
          do j=1,ny 
            do i=1,nx 
              f_ph = f_ph + dv * SpherHarmMesh(i,j,k,l,m,ImPart) * sum(HFpsi(i+(j-1)*nx+(k-1)*ny*nx,:,p) &
                & * HFpsi(i+(j-1)*nx+(k-1)*ny*nx,:,h)) 
            enddo
          enddo
        enddo
        F(1,p,h) = f_ph
        ! print *, p, h, f_ph
      enddo
    enddo

    ! can be restricted by running over allowed symmetry blocks


    do p=1,nwt ! can be restricted to particles
      do h=1,nwt ! can be restricted to holes
        f_ph = dv * sum(D_I_I(:,1) * sum(HFpsi(:,:,p) * HFpsi(:,:,h),2) ) 
              ! F(2,:,:) will host F02 in the future. 
        ! print *, p, h, f_ph
      enddo
    enddo

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