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
  use evolution

  implicit none

  !-----------------------------------------------------------------------------
  ! Define some FAM parameters
  !-----------------------------------------------------------------------------
  ! FAM energy frequencies
  real(KIND=dp) :: omega_fam  ! frequency of the perturbing field 
                              ! omega already defined as cranking frequency 
  ! A range of omega values can be passed by defining the min, max and stepsize
  ! i.e. omega = omega_min + k * omega_step < omega max for k=0,...
  real(KIND=dp) :: omega_min = 0.0_dp, omega_max = 0.0_dp
  real(KIND=dp) :: omega_step = 1.0_dp  ! Default stepsize of 1 MeV
  real(KIND=dp) :: smear = 1.0_dp  ! complex smearing parameter, default 0.5 MeV
  !    Note that the obtained strength is convoluted with a Lorentzian with FWHM 
  !    equal to Gamma = 2 * smear 
  real(KIND=dp) :: eta = 1.0e-3_dp ! small parameter entering derivatives, 
  !    Default currently set to 10-3. In the end, the strength should be 
  !    reasonably indepedent of the choice. 
  !      -> obsolete in explicit linearisation of the fields
  integer :: maxfamiter = 1000 ! maximal number of FAM iterations 
  !-----------------------------------------------------------------------------
  ! FAM amplitudes X, Y
  complex(KIND=dp), allocatable :: X(:,:) ! forward amplitudes HF basis
  !                                  | '-> sp index : hole
  !                                  '-> sp index : particle
  complex(KIND=dp), allocatable :: Y(:,:) ! backward amplitudes HF basis
  !                                  | '-> sp index : hole
  !                                  '-> sp index : particle
  !-----------------------------------------------------------------------------
  ! Perturbed densities
  ! /!\: contains the perturbation relative to the mean-field, e.g.
  !         rho_fam = rho_MF + drho
  complex(KIND=dp), allocatable :: drho(:,:)   ! perturbed normal density matrix
  complex(KIND=dp), allocatable :: dkappa(:,:) ! perturbed pairing density matrix
  ! complex(KIND=dp), allocatable :: dR(:,:)     ! perturbed generalised density matrix
  type(DensityVector)   :: dRs, dRa  ! perturbed densities in the mesh
  !                         |    '-> anti-symmetric part
  !                         '-> symmetric part
  type(PotentialVector) :: dFs, dFa  ! perturbed potentials on the mesh
  !                         |    '-> anti-symmetric part
  !                         '-> symmetric part
  ! type(DensityVector)   :: DensityPert  ! TO BE REMOVED
  ! type(PotentialVector) :: PotentialPert! TO BE REMOVED
  !-----------------------------------------------------------------------------
  ! unperturbed Hamiltonian and perturbed hamiltonian
  real(KIND=dp), allocatable :: H_unpert(:,:) ! unperturbed Hamiltonian in HF basis
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
  ! Do we need more identifiers for electric vs magnetic and isovector 
  ! vs isoscalar
  !-----------------------------------------------------------------------------
  ! convergence
  complex(KIND=dp), allocatable :: X_hist(:,:,:) ! history of X through FAM iters
  !                                       | | '-> sp index : hole
  !                                       | '-> sp index : particle
  !                                       '-> history index 
  complex(KIND=dp), allocatable :: Y_hist(:,:,:) ! history of Y through FAM iters
  !                                       | | '-> sp index : hole
  !                                       | '-> sp index : particle
  !                                       '-> history index 
  integer :: hist_max = 2 ! history size 
  integer :: hist_current_idx = 0 ! rolling index through the history
  ! notes: 
  !   Histories are implemented as circular buffers to mitigate copying data. 
  !   hist(hist_current_idx,:,:) contains the latest entry; the previous one can be
  !   accessed at idx = modulo(hist_current_idx - 2, hist_max) + 1). Rolling the
  !   index two steps back and then one forward is because mod gives values 
  !   0..hist_max-1 while fortran arrays use a 1-based index. 
  real(KIND=dp) :: tol_XY_conv = 1.0e-5_dp ! convergence tolerance for X and Y

  contains

  subroutine inifam(omega, PotentialsUnpert)
    implicit none
    !---------------------------------------------------------------------------
    ! Allocate the FAM objects, set the external field F and initialise the X
    ! and Y from first order, i.e. dH=0. 
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in) :: omega
    type(PotentialVector), intent(in) :: PotentialsUnpert
    real(KIND=dp), allocatable :: SolidHarmHF(:,:)
    integer :: p, h
    real(KIND=dp) :: occ_h, occ_p
    logical :: ImPart

    1 format(' S_',i1,i1,' (', f5.2, ') = ', es10.3)

    print *, "Initialise FAM matrices" 

    ! set omega frequency of perturbation
    omega_fam = omega

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! initialise the external field F

    if(.not.allocated(F)) then 
      allocate(F(nwt,nwt,2))

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      ! Get the solid harmonics Q_lm(i,j) = < i | r^l Y_lm | j > expressed 
      ! in HF basis. 
      
      allocate(SolidHarmHF(nwt,nwt)) 

      ! Set external field to E2, hardcoded for now
      ! l = 0
      ! m = 0
      ImPart = .false. ! real (.false.) , imaginary (.true.) 
      ! note: odd m and Im parts are not implemeted yet


      if(l==0) then
        SolidHarmHF = Rsq_spme()

        ! print *, 'Rsq'
        ! call print_spme_real(SolidHarmHF)

      else
        ! Calling a function in fission_MOI.f90
        SolidHarmHF = Qlm_spme(l, m, ImPart)

        ! Rescale, Qlm comes in units barn^(l/2)
        SolidHarmHF = SolidHarmHF * (100**(l/2.0)) 

        ! TODO: investigate sign change in the third row (column) wrt almost identical (row)
        ! SolidHarmHF(3,:) = -1.0 * SolidHarmHF(3,:)
        ! SolidHarmHF(:,3) = -1.0 * SolidHarmHF(:,3)
  
        ! print *, 'Q20'
        ! call print_spme_real(SolidHarmHF)

      endif
     

      ! note: 
      !   Stoitsov PRC 84 (2011) normalises the external field by a parameter
      !   alpha converting the units of the perturbation to MeV, and eventually 
      !   devides the obtained strength by alpha. 


      ! TODO: write a general transformation routine from the mesh to any 
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
          occ_p = 1.0 - rho_can(p) 
          if(occ_p < 1d-6) cycle
          F(p,h,1) = occ_p * occ_h * SolidHarmHF(p,h) ! ph block F20(p,h)
          F(p,h,2) = occ_p * occ_h * SolidHarmHF(h,p) ! hp block F02(p,h)
        enddo
      enddo
      ! This can be improved by some element-wise products occ^T @ SolidHarmHF @ occ
      
      ! Note to future self: for QFAM this will be replaced by a transformation 
      ! to the qp basis. 

      deallocate(SolidHarmHF)

    endif

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! set up the unperturbed Hamiltonian from the unperturbed potentials
    if(.not.allocated(H_unpert)) then 
      allocate(H_unpert(nwt,nwt))
      H_unpert = calc_sphamil(PotentialsUnpert, .false.)
    endif


    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! initialise perturbed Hamiltonian as 0
    if(.not.allocated(dH)) then 
      allocate(dH(nwt,nwt,2))
    endif

    dH = 0

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! initialise X and Y amplitudes and their history
    if(.not.allocated(X)) then
      allocate(X(nwt,nwt)) 
      allocate(Y(nwt,nwt))
    endif

    ! X and Y initialised from non-interacting response, i.e. setting dH = 0 
    ! in the FAM master
    call calculate_XY(dH)

    if(.not.allocated(X_hist)) then
      allocate(X_hist(hist_max,nwt,nwt)) 
      allocate(Y_hist(hist_max,nwt,nwt))
    endif

    X_hist = 0
    Y_hist = 0

    ! storing the initial x and Y in the history
    call store_XY_hist()

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! initialise the perturbed densities
    if(.not.allocated(drho)) then 
      allocate(drho(nwt,nwt))
      allocate(dkappa(nwt,nwt))
      ! allocate(dR(2*nwt,2*nwt))
    endif

    call build_perturbed_densities(X, Y, dRs, dRa)
    
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! initialise perturbed potentials as 0
    dFs = 0.0_dp * PotentialsUnpert
    dFa = 0.0_dp * PotentialsUnpert

    ! TO DO: replace by a better initialisation routine

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! print unperturbed strenght
    print 1, l, m, omega_fam, calc_strength()

  end subroutine inifam


  subroutine readfam(file_number)
    !---------------------------------------------------------------------------
    ! Read the namelist &fam/.
    !
    ! Input:
    !     file_number : channel number of opened file where to read from.
    !                   Optional. If not present, read from STDIN.
    !---------------------------------------------------------------------------
    integer(dp), intent(in),optional :: file_number
    real(KIND=dp) :: omega = -1.0_dp

    namelist /fam/      omega, omega_min, omega_max, omega_step,    &
    &                   smear, maxiter, l, m

    if(MPI_rank .eq. 0) then    
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      ! Reading the information on fam by the first MPI rank
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      print *, file_number
      if(present(file_number)) then
        read (unit=file_number, nml=fam)
      else
        read (unit=*, nml=fam)
      endif

      ! if a single fams frequency omega is passed, set min and max to omega
      if(omega .ne. -1.0_dp) then
        omega_min = omega
        omega_max = omega
      endif
    endif


  end subroutine

  subroutine printfam
    1 format ( 32('-'), ' FAM information ', 31('-'))
    2 format ( ' FAM frequency range:   ', /, &
    &          '    omega_min        = ', f10.3, /,  &
    &          '    omega_max        = ', f10.3, /,  &
    &          '    omega_step       = ', f10.3)
    21 format ('    complex smearing = ' ,f10.3)
    3 format ( ' Perturbing field:   ', /, &
    &          '    F = Q_', i1, i1)
    4 format ( ' Convergence:   ', /, &
    &          '    Maximal number of iterations: ',i8, /,  &
    &          '    ||X||, ||Y|| convergence  < ',es8.1)

    print 1
    print 2, omega_min, omega_max, omega_step
    print 21, smear
    print 3, l, m
    print 4, maxfamiter, tol_XY_conv

  end subroutine

  subroutine calculate_XY(dH)
    !---------------------------------------------------------------------------
    ! Compute the X and Y amplitudes from the FAM master equation
    !---------------------------------------------------------------------------

    implicit none

    real(KIND=dp), dimension(:,:,:), intent(in)  :: dH ! perturbed H in QP basis
    integer :: p, h
    real(KIND=dp) :: occ_h, occ_p, e_h, e_p
    ! complex(KIND=dp), allocatable :: denomX(:,:),  denomY(:,:)

    print *, "update X and Y"


    X = - (F(:,:,1) + dH(:,:,1))
    Y = - (F(:,:,2) + dH(:,:,2))

    ! normalise with energy denominator
    do h = 1, nwt
      occ_h = rho_can(h)
      e_h = spenergies(h) 
      if(occ_h < 1d-6) cycle
      do p = 1, nwt
        occ_p = 1.0 - rho_can(p) ! degeneracy is always 1 since T is broken
        e_p = spenergies(p) 
        if(occ_p < 1d-6) cycle
        X(p,h) = X(p,h) / (e_p - e_h - CMPLX(omega_fam,smear,KIND=dp) )
        Y(p,h) = Y(p,h) / (e_p - e_h + CMPLX(omega_fam,smear,KIND=dp) )
        ! print *, p, h, e_p, e_h, occ_h, occ_p, 
        ! print *, X(p,h), Y(p,h), F(p,h,1), F(p,h,2)
      enddo
    enddo

    ! print *, 'X'
    ! call print_spme_complex(X)

    ! print *, 'Y'
    ! call print_spme_complex(Y)

    ! print *, 'rho'
    ! call print_spme_real(rho_pairing)

    ! allocate(denomX(nwt,nwt))
    ! allocate(denomY(nwt,nwt))

    ! ! energy denominator
    ! do h = 1, nwt
    !   occ_h = rho_can(h)
    !   e_h = spenergies(h) 
    !   if(occ_h < 1d-6) cycle
    !   do p = 1, nwt
    !     occ_p = 1.0 - rho_can(p) ! degeneracy is always 1 since T is broken
    !     e_p = spenergies(p) 
    !     if(occ_p < 1d-6) cycle
    !     denomX(p,h) = (e_p - e_h - CMPLX(omega_fam,smear,KIND=dp) )
    !     denomY(p,h) = (e_p - e_h + CMPLX(omega_fam,smear,KIND=dp) )
    !     ! print *, p, h, e_p, e_h, occ_h, occ_p, 
    !     ! print *, X(p,h), Y(p,h), F(p,h,1), F(p,h,2)
    !   enddo
    ! enddo

    ! print *, 'denomX'
    ! call print_spme_complex(denomX)

    ! print *, 'denomY'
    ! call print_spme_complex(denomY)

    ! stop


    print *, '||F20|| = ', sqrt(sum( abs(F(:,:,1))**2) )
    print *, '||dH20|| = ', sqrt(sum( abs(dH(:,:,1))**2) )


  ! print *, '||X|| = ', sqrt(sum( abs(X(:,:))**2) )
  ! print *, '||Y|| = ', sqrt(sum( abs(Y(:,:))**2) )

  end subroutine calculate_XY

  subroutine print_all_fam_spmat()

    print *, 'X'
    call print_spme_complex(X)

    print *, 'Y'
    call print_spme_complex(Y)

    print *, 'drho'
    call print_spme_complex(X)

    print *, 'dH20'
    call print_spme_real(dH(:,:,1))

    print *, 'dH02'
    call print_spme_real(dH(:,:,2))
    
    print *, 'F20'
    call print_spme_real(F(:,:,1))

    print *, 'F02'
    call print_spme_real(F(:,:,2))

  end subroutine


  subroutine print_spme_real(A)
    implicit none
    real(kind=dp), intent(in) :: A(:,:)
    integer :: si, B, N, i

    si = 0
    do B=1,8
      N = HFBLocks(B)

      print *, 'BLOCK', B
      do i=si+1,si+N
        print '(99f10.5)',  A(i, si+1:si+N)
      enddo
      print *
      si = si + N
    enddo
    print *
    
  end subroutine print_spme_real

  subroutine print_spme_complex(A)
    implicit none
    complex(kind=dp), intent(in) :: A(:,:)
    integer :: si, B, N, i

    si = 0
    do B=1,8
      N = HFBLocks(B)

      print *, 'BLOCK', B
      do i=si+1,si+N
        print "(*('('sf8.5','sf8.5')':x))",  A(i, si+1:si+N)
      enddo
      print *
      si = si + N
    enddo
    print *
    
  end subroutine print_spme_complex

  subroutine store_XY_hist()
    !---------------------------------------------------------------------------
    ! Store the current X and Y into their histories. 
    !---------------------------------------------------------------------------

    ! roll the current index one step forward
    hist_current_idx = modulo(hist_current_idx, hist_max) + 1

    ! store X and Y in current spot
    X_hist(hist_current_idx, :, :) = X(:,:)
    Y_hist(hist_current_idx, :, :) = Y(:,:)

  end subroutine store_XY_hist

  subroutine mix_XY_linear(alpha)
    !---------------------------------------------------------------------------
    ! Simple linear mixing of the X and amplitudes, i.e. 
    !   X^[i] = alpha * X^[i] + (1-alpha) X^[i-1]
    ! No return. Changes are made to the current X and Y.
    !---------------------------------------------------------------------------

    real(KIND=dp), intent(in) :: alpha

    X = alpha * X + (1.0 - alpha) * X_hist(hist_current_idx, :, :) 
    Y = alpha * Y + (1.0 - alpha) * Y_hist(hist_current_idx, :, :) 

  end subroutine mix_XY_linear


  subroutine iniHFdensities()
    !---------------------------------------------------------------------------
    ! initialse the rho and kappa matrices in HF basis as (nwt, nwt) matrices
    ! these are coined as rho_pairing and kappa_pairing
    !---------------------------------------------------------------------------
    implicit none
    integer :: i
      
    allocate(rho_pairing(nwt,nwt))
    allocate(kappa_pairing(nwt,nwt))
  
    rho_pairing = 0
    kappa_pairing = 0
    
    do i=1,nwt
      rho_pairing(i,i) = rho_can(i)
    enddo
  
  end subroutine iniHFdensities


  subroutine build_perturbed_densities(X, Y, dRs, dRa)
    !---------------------------------------------------------------------------
    ! Build the perturbed mean-field densities.
    !---------------------------------------------------------------------------

    implicit none
    complex(KIND=dp), intent(in) :: X(:,:), Y(:,:)
    type(DensityVector), intent(out) :: dRs, dRa


    print *, "build perturbed densities"

    drho = X + transpose(Y) 
    dkappa = 0

    ! print *, 'drho'
    ! call print_spme_complex(drho)

    call densit_offdiag(drho, dkappa, dRs, dRa)

  end subroutine build_perturbed_densities

  subroutine build_dH_explicit(R, dRs, dRa)
    !---------------------------------------------------------------------------
    ! Build the perturbed single-particle Hamiltonian
    !---------------------------------------------------------------------------
    
    implicit none
    type(DensityVector), intent(in) :: R, dRs, dRa
    type(PotentialVector) :: dFsNew, dFaNew
    real(KIND=dp), allocatable :: HPert(:,:)
    integer :: i, h, p
    real(KIND=dp) :: occ_h, occ_p
    real(KIND=dp) :: alpha = 1.0d0 ! linear mixing coeff, disable : mixing of XY

    print *, "build perturbed hamiltonian using explicit linearisation"

    allocate(HPert(nwt,nwt))

    ! explicit linearisation of the fields
    call calc_perturbed_potentials(R, dRs, dRa, dFsNew, dFaNew)

    ! necessary? 
    call combine_potentials(dFsNew)

    ! linear mixing of the sym and antisym perturbed fields with previous iteration
    dFs = alpha * dFsNew  + (1.0 - alpha) * dFs
    dFa = alpha * dFaNew  + (1.0 - alpha) * dFa

    ! construct the sp hamiltonian
    HPert = calc_sphamil_me( HFpsi, HFdpsi, HFddpsi,dFs, dFa, .false.)

    ! print *, 'HPert'
    ! call print_spme_real(HPert)


    ! store ph and hp blocks in dH
    do h = 1, nwt
      occ_h = rho_can(h)
      if(occ_h < 1d-6) cycle
      do p = 1, nwt
        occ_p = 1.0 - rho_can(p) 
        if(occ_p < 1d-6) cycle
        dH(p,h,1) = occ_p * occ_h * HPert(p,h) ! ph block dH20(p,h)
        dH(p,h,2) = occ_p * occ_h * HPert(h,p) ! hp block dH02(p,h)
      enddo
    enddo

    ! print *, 'dH20'
    ! call print_spme_real(dH(:,:,1))

    ! print *, 'dH02'
    ! call print_spme_real(dH(:,:,2))


    ! print *,  sqrt(sum( abs(dH(:,:,1))**2) )

  end subroutine build_dH_explicit

  subroutine build_dH_findiff(Density, DensityPert)
    !---------------------------------------------------------------------------
    ! Build the perturbed single-particle Hamiltonian using finite difference
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! NOT OPERATIONAL
    !   -> use build_dH_explicit() instead.
    ! 
    ! Notes:
    !    This function is not operational at this point but is kept for potential 
    !    test in the future. In particular, it would allow to test 
    !    calc_perturbed_potentials.
    !---------------------------------------------------------------------------

    implicit none
    type(DensityVector), intent(in) :: Density, DensityPert
    type(DensityVector) :: DensityTot
    type(PotentialVector) :: PotentialTotNew
    real(KIND=dp), allocatable :: HPert(:,:)
    integer :: i, h, p
    real(KIND=dp) :: occ_h, occ_p

    print *, "build perturbed hamiltonian using finite difference"

    ! allocate(HPert(nwt,nwt))

    ! Calculate the total perturbed potentials, i.e. static mean-field + perturbation, 
    ! from the total perturbed densit, i.e. static mean-field + perturbation

    ! PotentialTotNew = calcPotentials(Density + eta * DensityPert)

    ! /!\ presently incorrect. This should be replaced by a function which treats
    ! symm and atisymm parts where the Rs = Density + eta * dRs and Ra = eta * dRa
    ! call calcPotentials(Density + eta * dRs, eta * dRa, Fs, Fa)

    ! call combine_potentials(PotentialPertNew)


    ! TODO mixing is still absent, it would require keeping track of the  
    ! unperturbed plus perturbed fields

    ! construct the sp hamiltonian
    ! HPert = calc_sphamil_me( HFpsi, HFdpsi, HFddpsi,PotentialTotNew, .false.)

    ! /!\ presently incorrect. Again, this should be replaced by a function which 
    ! treats symm and atisymm parts of the fields. 
    ! call calc_sphamil_me( HFpsi, HFdpsi, HFddpsi, Fs, Fa, .false.)


    ! ! compute dH by finite difference, i.e. subtract the unperturbed Hamiltonian
    ! HPert = HPert - H_unpert
    ! ! ... and devide by small parameter eta
    ! HPert = HPert / eta

    ! ! store ph and hp blocks in dH
    ! do h = 1, nwt
    !   occ_h = rho_can(h)
    !   if(occ_h < 1d-6) cycle
    !   do p = 1, nwt
    !     occ_p = 1.0 - rho_can(p) 
    !     if(occ_p < 1d-6) cycle
    !     dH(p,h,1) = occ_p * occ_h * HPert(p,h) ! ph block dH20(p,h)
    !     dH(p,h,2) = occ_p * occ_h * HPert(h,p) ! hp block dH02(p,h)
    !   enddo
    ! enddo

    ! print *,  sqrt(sum( abs(dH(:,:,1))**2) )

  end subroutine build_dH_findiff


  function calc_strength() result (S_out)
    !---------------------------------------------------------------------------
    ! Calculate the strength S(omega,F)
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! obtained from 
    !     S(omega,F) = - 1 /pi * Im Tr (F^dagger * drho)
    ! where 
    !    Tr (F^dagger * drho) = sum_ab (F^20_ab^* X_ab + F^02_ab^* Y_ab)
    ! 
    ! note: 
    !  - normalisation of external field may have to be taken into account
    !    S -> S/alpha
    !  - F is supposed to be real. If F is replaced by a complex field, the
    !    complex conjugation must be added
    !---------------------------------------------------------------------------

    complex(KIND=dp) :: S
    real(KIND=dp) S_out
    integer :: h, p
    real(KIND=dp) :: occ_h, occ_p

    S = 0
    do h = 1, nwt
      occ_h = rho_can(h)
      if(occ_h < 1d-6) cycle
      do p = 1, nwt
        occ_p = 1.0 - rho_can(p) 
        if(occ_p < 1d-6) cycle
        S = S + F(p,h,1) * X(p,h) + F(p,h,2) * Y(p,h)
      enddo
    enddo

    S_out = - S%im / pi

  end function calc_strength


  subroutine test_convergence(conv, div)
    !---------------------------------------------------------------------------
    ! Judge the convergence of the FAM iterations based on difference of X and Y
    ! with respect to previous iteration
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! The convergence measure corresponds to the Frobenius norm of the 
    ! difference fo the current X(i) (Y(i)) and the one of the previous 
    ! iteration X(i-1) (Y(i-1)) stored in X_hist and Y_hist, i.e.
    !      ||X(i) - X(i-1)|| / ||X(i)|| < tolerance
    ! 
    ! The Frobenius norm ||A|| is evaluated as sqrt(sum[abs(A(:,:))**2]) where
    ! the abs takes care of obtaining the modulus of the complex values.
    !---------------------------------------------------------------------------
    logical, intent(out) :: conv, div
    integer :: idx_prev
    real(KIND=dp) :: DX_norm, DY_norm, X_norm, Y_norm

    conv = .false.
    div = .false.

    X_norm = sqrt(sum( abs(X_hist(hist_current_idx,:,:))**2))
    Y_norm = sqrt(sum( abs(Y_hist(hist_current_idx,:,:))**2))


    print * , "||X|| = ", X_norm
    print * , "||Y|| = ", Y_norm

    if( (X_norm .ge. 1.0d3) .or. (Y_norm .ge. 1.0d3)) then
      div = .true.
    endif

    ! previous index in hist obtained by rolling back twice and adding one
    idx_prev = modulo(hist_current_idx - 2, hist_max) + 1

    DX_norm = sqrt( sum( abs(X_hist(hist_current_idx,:,:) - X_hist(idx_prev,:,:))**2) )
    DX_norm = DX_norm / X_norm

    DY_norm = sqrt( sum( abs(Y_hist(hist_current_idx,:,:) - Y_hist(idx_prev,:,:))**2) )
    DY_norm = DY_norm / Y_norm

    print * , "convergence: ||DX|| = ", DX_norm, "   ||DY|| = ", DY_norm

    if( (DX_norm<tol_XY_conv) .and. (DY_norm<tol_XY_conv)) then
      conv = .true.
    endif

  end subroutine test_convergence


  function Rsq_spme() result (Rsq)

    real(kind=dp) , allocatable :: Rsq(:,:)
    integer :: i, j, B, N, s
    
    allocate(Rsq(nwt, nwt))

    ! Initialize to zero
    Rsq = 0.0_dp

    s = 0
    do B = 1, 8
      N =  HFBlocks(B) ; if(N.eq.0) cycle
      ! Rsq is always block diagonal, for both P conserving and P broken states
      do i=1,N   
        do j=i,N
          ! me = Int d^3r Sum_sigma psi^*_i(r,sigma) psi_j(r,sigma) Rsq(r)
          Rsq(s+i,s+j) = sum(sum(HFpsi(:,:,s+i)*HFpsi(:,:,s+j),2) * sum(meshgrid**2,2)) * dv

          ! Rsq matrix elements are real symmetric
          Rsq(s+j,s+i) = Rsq(s+i,s+j)

        enddo
      enddo
      s = s + N
    enddo

  end function

end module fam
