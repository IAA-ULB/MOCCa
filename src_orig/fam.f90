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
  integer :: maxfamiter = 100 ! maximal number of FAM iterations 
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
  ! complex(KIND=dp), allocatable :: dR(:,:)   ! perturbed generalised density matrix
  type(DensityVector)   :: Runper    ! static mean-field densities on the mesh
  type(DensityVector)   :: dRs, dRa  ! perturbed densities on the mesh
  !                         |    '-> anti-symmetric part
  !                         '-> symmetric part
  type(PotentialVector) :: dFs, dFa  ! perturbed potentials on the mesh
  !                         |    '-> anti-symmetric part
  !                         '-> symmetric part
  !-----------------------------------------------------------------------------
  ! unperturbed Hamiltonian and perturbed hamiltonian
  real(KIND=dp), allocatable :: HUnper(:,:) ! unperturbed Hamiltonian in HF basis
  real(KIND=dp), allocatable :: dH(:,:,:)   ! ph and hp block of the perturbed
  !                                | | |      Hamiltonian in HF basis
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
  real(KIND=dp) :: XY_prec = 1.0e-10_dp ! convergence tolerance for X and Y
  !-----------------------------------------------------------------------------
  ! verbosity
  integer :: verbose = 0
  ! 0: very limited printing
  ! 1: printing some function calls
  ! 2: printing all sp matrices at each iteration


  contains
  

  subroutine inifam(omega, DensUnper, PotUnper)
    implicit none
    !---------------------------------------------------------------------------
    ! Allocate the FAM objects, set the external field F and initialise the X
    ! and Y from first order, i.e. dH=0. 
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in)          :: omega
    type(DensityVector), intent(in)    :: DensUnper
    type(PotentialVector), intent(in)  :: PotUnper
    real(KIND=dp), allocatable         :: SolidHarmHF(:,:)
    logical                            :: ImPart

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

        ! TODO: investigate signs in Q20 which seems suspicious in O16 nwt24 test case
        ! 3rd row/col in sym block 1 differs in sign wrt blocks 2, 5 and 6. 


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

      call get_ph_hp_blocks(SolidHarmHF, F(:,:,1), F(:,:,2))

      deallocate(SolidHarmHF)

    endif

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! set up the unperturbed Hamiltonian from the unperturbed potentials
    if(.not.allocated(Hunper)) then 
      allocate(Hunper(nwt,nwt))
      Hunper = calc_sphamil(PotUnper, .false.)
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
    ! store the unperturbed densities
    RUnper = DensUnper

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
    dFs = 0.0_dp * PotUnper
    dFa = 0.0_dp * PotUnper

    ! TODO: replace by a better initialisation routine
    ! This might require Hephaestos

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
    &                   smear, maxiter, l, m, XY_prec

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

      maxfamiter = maxiter

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
    print 4, maxfamiter, XY_prec

  end subroutine

  subroutine iterate_dHsp(dHsp_flat, dHspout_flat)
    !---------------------------------------------------------------------------
    ! Perform one FAM loop of the perturbed single-particle hamiltonian dH
    ! (in HF basis), which contain dh and ddelta (in the QFAM).  
    !---------------------------------------------------------------------------
    1 format('||dH_ph|| = ', es10.3, '     ||dH_hp|| = ', es10.3)
    2 format('||X|| = ', es10.3, '     ||Y|| = ', es10.3)
    3 format(' S_',i1,i1,' (', f5.2, ') = ', es10.3)

    implicit none
    real(KIND=dp), dimension(:), target, intent(in)   :: dHsp_flat
    real(KIND=dp), dimension(:), target, intent(out)  :: dHspout_flat

    real(KIND=dp), pointer :: dHsp(:,:), dHspout(:,:)

    integer       :: p, h
    real(KIND=dp) :: occ_h, occ_p

    ! pointer remapping for reshaping 1D flat arrays into 2D matrices
    dHsp(1:nwt,1:nwt) => dHsp_flat(:)
    dHspout(1:nwt,1:nwt) => dHspout_flat(:)

    ! get the ph and hp subblocks
    call get_ph_hp_blocks(dHsp, dH(:,:,1), dH(:,:,2))

    print 1, sqrt(sum( abs(dH(:,:,1))**2) ), sqrt(sum( abs(dH(:,:,2))**2) )

    ! calculate X and Y from the perturbed dH
    call calculate_XY(dH)
    print 2, sqrt(sum( abs(X(:,:))**2) ), sqrt(sum( abs(Y(:,:))**2) )

    print 3, l,m, omega_fam,  calc_strength()


    ! Apply simple linear mixing of X and Y. 
    ! call mix_XY_linear(lin_mix_coeff)
    ! -> this may be skipped when using GMRES

    call store_XY_hist()

    ! build the perturbed densities on the mesh dRs, dRa from X and Y
    call build_perturbed_densities(X, Y, dRs, dRa)

    ! explicit linearisation of the fields
    call calc_perturbed_potentials(RUnper, dRs, dRa, dFs, dFa)

    ! necessary? 
    call combine_potentials(dFs)

    ! construct the sp hamiltonian
    dHspout = calc_sphamil_me( HFpsi, HFdpsi, HFddpsi,dFs, dFa, .false.)

    if (verbose > 1) call print_all_fam_spmat()

  end subroutine iterate_dHsp

  subroutine calculate_XY(dH)
    !---------------------------------------------------------------------------
    ! Compute the X and Y amplitudes from the FAM master equation
    !---------------------------------------------------------------------------
    implicit none
    real(KIND=dp), intent(in)  :: dH(:,:,:) ! perturbed H in QP basis

    integer       :: p, h
    real(KIND=dp) :: occ_h, occ_p, e_h, e_p

    if (verbose > 0) print *, "update X and Y"


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
      enddo
    enddo

  end subroutine calculate_XY


  subroutine store_XY_hist()
    !---------------------------------------------------------------------------
    ! Store the current X and Y into their histories. 
    !---------------------------------------------------------------------------
    if (verbose > 0) print *, "store X and Y"

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

    if (verbose > 0) print *, "mix X and Y with alpha=", alpha

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
    complex(KIND=dp), intent(in)     :: X(:,:), Y(:,:)
    type(DensityVector), intent(out) :: dRs, dRa

    if (verbose > 0) print *, "build perturbed densities"

    drho = X + transpose(Y) 
    dkappa = 0

    call densit_offdiag(drho, dkappa, dRs, dRa)

  end subroutine build_perturbed_densities


  subroutine build_dH_explicit(R, dRs, dRa)
    !---------------------------------------------------------------------------
    ! Build the perturbed single-particle Hamiltonian
    !---------------------------------------------------------------------------
    1 format('||dH_ph|| = ', es10.3, '     ||dH_hp|| = ', es10.3)

    implicit none
    type(DensityVector), intent(in) :: R, dRs, dRa
    real(KIND=dp), allocatable :: dHsp(:,:)

    allocate(dHsp(nwt,nwt))

    ! explicit linearisation of the fields
    call calc_perturbed_potentials(R, dRs, dRa, dFs, dFa)

    ! necessary? 
    call combine_potentials(dFs)

    ! construct the sp hamiltonian
    dHsp = calc_sphamil_me( HFpsi, HFdpsi, HFddpsi,dFs, dFa, .false.)

   ! get the ph and hp subblocks
    call get_ph_hp_blocks(dHsp, dH(:,:,1), dH(:,:,2))

    print 1, sqrt(sum( abs(dH(:,:,1))**2) ), sqrt(sum( abs(dH(:,:,2))**2) )


  end subroutine build_dH_explicit


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

    1 format('||X|| = ', es10.3, '     ||Y|| = ', es10.3)
    2 format('Convergence: ', '||dX|| = ', es10.3, '     ||dY|| = ', es10.3)
    logical, intent(out) :: conv, div
    integer :: idx_prev
    real(KIND=dp) :: DX_norm, DY_norm, X_norm, Y_norm

    conv = .false.
    div = .false.

    X_norm = sqrt(sum( abs(X_hist(hist_current_idx,:,:))**2))
    Y_norm = sqrt(sum( abs(Y_hist(hist_current_idx,:,:))**2))


    ! print 1, X_norm, Y_norm

    if( (X_norm .ge. 1.0d3) .or. (Y_norm .ge. 1.0d3)) then
      div = .true.
    endif

    ! previous index in hist obtained by rolling back twice and adding one
    idx_prev = modulo(hist_current_idx - 2, hist_max) + 1

    DX_norm = sqrt( sum( abs(X_hist(hist_current_idx,:,:) - X_hist(idx_prev,:,:))**2) )
    DX_norm = DX_norm / X_norm

    DY_norm = sqrt( sum( abs(Y_hist(hist_current_idx,:,:) - Y_hist(idx_prev,:,:))**2) )
    DY_norm = DY_norm / Y_norm

    print 2, DX_norm, DY_norm

    if( (DX_norm < XY_prec) .and. (DY_norm < XY_prec)) then
      conv = .true.
    endif

  end subroutine test_convergence


  subroutine get_ph_hp_blocks(M, Mph, Mhp)
    !---------------------------------------------------------------------------
    ! Get the particle-hole and hole-particle subblocks of a one-body operator
    ! M. Occupation are obtained from the diagonal elements of rho_can. 
    ! 
    ! NOTE : 
    ! - the ordering of the sp labels is always the particle label first 
    ! and the hole label second, i.e. Mhp(p,h) and Mph(p,h). 
    ! - for simplicity, Mph and Mhp are of size (nwt,nwt). 
    !---------------------------------------------------------------------------

    implicit none
    real(KIND=dp), intent(in) :: M(:,:)
    real(KIND=dp), intent(out) :: Mph(:,:), Mhp(:,:)
    integer       :: p, h
    real(KIND=dp) :: occ_h, occ_p

    Mph = 0
    Mhp = 0

    do h = 1, nwt
      occ_h = rho_can(h)
      if(occ_h < 1d-6) cycle
      do p = 1, nwt
        occ_p = 1.0 - rho_can(p) 
        if(occ_p < 1d-6) cycle
        Mph(p,h) = occ_p * occ_h * M(p,h)
        Mhp(p,h) = occ_p * occ_h * M(h,p)
      enddo
    enddo

    ! This can be more efficient by using some mask and elementwise multiplication

    ! For QFAM this will have to be generalised to M20 and M02 obtained from a 
    ! Bogoliubov transformation to the qp basis. 

  end subroutine get_ph_hp_blocks



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

  end function Rsq_spme

  function norm_dH(dH) result(res)
    ! abstract template procedure dH -> real required for procedural argument to gmres
    ! to be updated to the objects of the dimensions of the perturbed
    ! sp hamiltonian dh and ddelta (in HF basis)
    real(KIND=dp), dimension(:), intent(in)  :: dH
    real(KIND=dp)                            :: res

    res = sqrt(sum(dH(:) * dH(:)))

  end function

  function ScProd_dH(dHl, dHr) result(res)
    ! abstract template procedure (dH,dH) -> complex required for procedural argument to gmres
    ! to be updated to the objects of the dimensions of the perturbed
    ! sp hamiltonian dh and ddelta (in HF basis)
    real(KIND=dp), dimension(:), intent(in)  :: dHl, dHr
    real(KIND=dp)                            :: res

    res = sum(dHl(:) *  dHr(:))

  end function

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

end module fam
