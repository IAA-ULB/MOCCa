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
  !------------------------------------------------------------------------------
  ! Hephaestos keywords
  ! 
  ! TR  : $TR
  ! NTR : $NTR
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
  real(KIND=dp) :: omega_min = 0.0_dp, omega_max = 30.0_dp
  real(KIND=dp) :: omega_step = 1.0_dp  ! Default stepsize of 1 MeV
  real(KIND=dp) :: smear = 1.0_dp  ! complex smearing parameter, default 0.5 MeV
  !    Note that the obtained strength is convoluted with a Lorentzian with FWHM 
  !    equal to Gamma = 2 * smear 
  !-----------------------------------------------------------------------------
  ! mixing strategy
  integer :: fam_mixingscheme = 0 ! 0 : GMRES (default)
  !                                 1 : linear mixing of dH
  integer :: fam_maxiter = 100 ! maximal number of FAM iterations 
  integer :: fam_maxhist = 30 ! maximal history size of GMRES 
  ! Coefficient for the linear mixing of FAM iterations
  real(KIND=dp) :: fam_lin_mix = 0.3_dp
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
  !         rho(omega) = rho_MF + drho(omega)
  complex(KIND=dp), allocatable :: drho(:,:)   ! perturbation to the normal density matrix
  complex(KIND=dp), allocatable :: dkappa(:,:) ! perturbation to the pairing density matrix
  ! complex(KIND=dp), allocatable :: dR(:,:)   ! perturbation to the generalised density matrix
  type(DensityVector)   :: Runper    ! static mean-field densities on the mesh
  type(DensityVector)   :: dRs, dRa  ! perturbation to the densities on the mesh
  !                         |    '-> anti-symmetric part
  !                         '-> symmetric part
  type(PotentialVector) :: dFs, dFa  ! perturbation to the potentials on the mesh
  !                         |    '-> anti-symmetric part
  !                         '-> symmetric part
  !-----------------------------------------------------------------------------
  ! unperturbed Hamiltonian and perturbed hamiltonian
  real(KIND=dp), allocatable :: HUnper(:,:) ! unperturbed Hamiltonian in HF basis
  ! -> currently not used, except for one routine in fam_testing.f90
  complex(KIND=dp), allocatable :: dH(:,:,:)   ! ph and hp block of the perturbing
  !                                   | | |      Hamiltonian in HF basis
  !                                   | | '-> 1: ph block, 2: hp block
  !                                   | '-> sp index : hole
  !                                   '-> sp index : particle
  complex(KIND=dp), allocatable :: dH_free_flat(:) ! free response of sp hamil
  !                                             '-> nwt x nwt
  !-----------------------------------------------------------------------------
  ! external field
  complex(KIND=dp), allocatable :: F(:,:,:)  ! perturbing external field in HF basis
  !                                  | | '-> 1: ph block, 2: hp block 
  !                                  | '-> sp index : hole
  !                                  '-> sp index : particle
  integer :: l, m ! anuglar momentum and projection quantum number of the multipole moment
  real(KIND=dp) :: eff_charge_n = 1.0_dp ! effective charges for neutrons in units of e
  real(KIND=dp) :: eff_charge_p = 1.0_dp ! effective charges for protons in units of e
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
  real(KIND=dp) :: fam_precision = 1.0e-5_dp ! convergence tolerance for X and Y
  !-----------------------------------------------------------------------------
  ! verbosity
  integer :: fam_verbose = 1
  ! 0: no printing. Used during GMRES as output would be confusing
  ! 1: limited printing. Used in the final FAM iteration once GMRES is converged
  !    (default)
  ! 2: printing all function calls. Useful for debugging. 
  ! 3: printing all sp matrices at each iteration. Useful for debugging. 

  interface get_ph_hp_blocks
    module procedure get_ph_hp_blocks_complex
    module procedure get_ph_hp_blocks_real
  end interface get_ph_hp_blocks


  contains
  

  subroutine inifam(omega, DensUnper, PotUnper)
    implicit none
    !---------------------------------------------------------------------------
    ! Allocate the FAM objects, set the external field F and initialise the X
    ! and Y from first order, i.e. dH=0. 
    !
    ! Input:
    !    omega      : frequency of the perturbing field
    !    DensUnper  : unperturbed densities on the mesh
    !    PotUnper   : unperturbed potentials on the mesh
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in)          :: omega
    type(DensityVector), intent(in)    :: DensUnper
    type(PotentialVector), intent(in)  :: PotUnper

    1 format(' S_',i1,i1,' (', f5.2, ') = ', es10.3)

    print *, "Initialise FAM matrices" 

    ! set omega frequency of perturbation
    omega_fam = omega

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! initialise the external field F

    if(.not.allocated(F)) then 
      F = get_f_LK(l, m, eff_charge_n, eff_charge_p)
    endif

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! set up the unperturbed Hamiltonian from the unperturbed potentials
    !if(.not.allocated(Hunper)) then
    !  allocate(Hunper(nwt,nwt))
    !  Hunper = calc_sphamil(PotUnper, .false.)
    ! endif


    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! initialise perturbed Hamiltonian as 0
    if(.not.allocated(dH)) then 
      allocate(dH(nwt,nwt,2))
    endif

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! initialise the free response of the sp hamiltonian
    if(.not.allocated(dH_free_flat)) then 
      allocate(dH_free_flat(nwt * nwt))
    endif

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! initialise X and Y amplitudes and their history
    if(.not.allocated(X)) then
      allocate(X(nwt,nwt)) 
      allocate(Y(nwt,nwt))
    endif

    if(.not.allocated(X_hist)) then
      allocate(X_hist(hist_max,nwt,nwt)) 
      allocate(Y_hist(hist_max,nwt,nwt))
    endif

    X_hist = 0
    Y_hist = 0

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
  

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! calculate free response by iterating FAM once starting from 0
    ! this also sets X, Y, drho, dkappa etc. to their respective free values

    dH_free_flat = 0
    call iterate_dHsp(dH_free_flat, dH_free_flat)


  end subroutine inifam


  subroutine readfam(file_number)
    !---------------------------------------------------------------------------
    ! Read the namelist &fam/.
    !
    ! Input:
    !     file_number : channel number of opened file where to read from.
    !                   Optional. If not present, read from STDIN.
    !---------------------------------------------------------------------------
    integer(dp), intent(in), optional :: file_number
    real(KIND=dp) :: omega = -1.0_dp
    integer       :: mixingscheme = 0
    integer       :: maxiter = 100, maxhist = 30

    namelist /fam/  omega, omega_min, omega_max, omega_step, smear, maxiter, &
    &               maxhist, l, m, fam_precision, mixingscheme, fam_lin_mix, &
    &               eff_charge_n, eff_charge_p


    if(MPI_rank .eq. 0) then
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      ! Reading the information on fam by the first MPI rank
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      if(present(file_number)) then
        read (unit=file_number, nml=fam)
      else
        read (unit=*, nml=fam)
      endif

      fam_maxiter = maxiter
      fam_maxhist = maxhist
      fam_mixingscheme = mixingscheme

      ! if a single fams frequency omega is passed, set min and max to omega
      if(omega .ne. -1.0_dp) then
        omega_min = omega
        omega_max = omega
      endif
    endif

  end subroutine readfam

  subroutine printfam
    1 format ( 32('-'), ' FAM information ', 31('-'))
    2 format ( ' FAM frequency range:   ', /, &
    &          '    omega_min        = ', f10.3, /,  &
    &          '    omega_max        = ', f10.3, /,  &
    &          '    omega_step       = ', f10.3, /,  &
    &          '    complex smearing = ' ,f10.3)
    3 format ( ' Perturbing field:   ', /, &
    &          '    F = Q_', i1, i1,/, &
    &          '    neutron eff charge = ', f10.3, ' e', /, &
    &          '    proton eff charge  = ', f10.3, ' e')
    41 format (' Convergence strategy: GMRES', /,  &
    &          '    max history size = ', i8, /,  &
    &          '    max # iterations = ', i8, /,  &
    &          '    res convergence  < ', es8.1)
    42 format (' Convergence strategy: linear mixing', /,  &
    &          '    mixing coef alpha = ', f10.3, /,  &
    &          '    max # iterations = ', i8, /,  &
    &          '    dh convergence   < ', es8.1)



    print 1
    print 2, omega_min, omega_max, omega_step, smear
    print 3, l, m, eff_charge_n, eff_charge_p
    if (fam_mixingscheme==0) print 41, fam_maxhist, fam_maxiter, fam_precision
    if (fam_mixingscheme==1) print 42, fam_lin_mix, fam_maxiter, fam_precision

  end subroutine

  subroutine iterate_dHsp(dHsp_flat, dHspout_flat)
    !---------------------------------------------------------------------------
    ! Perform one FAM loop of the perturbed single-particle hamiltonian dH
    ! (in HF basis), which contain dh and ddelta (in the QFAM).  
    !
    ! Input:
    !    dHsp_flat    : perturbed sp hamiltonian in HF basis as a flat array
    ! Output:
    !    dHspout_flat : updated perturbed sp hamiltonian in HF basis as a
    !                   flat array
    !---------------------------------------------------------------------------
    1 format('||dH_ph|| = ', es10.3, '     ||dH_hp|| = ', es10.3)
    2 format('||X|| = ', es10.3, '     ||Y|| = ', es10.3)
    3 format(' S_',i1,i1,' (', f5.2, ') = ', es10.3)

    implicit none
    complex(KIND=dp), dimension(:), target, intent(in)   :: dHsp_flat
    complex(KIND=dp), dimension(:), target, intent(out)  :: dHspout_flat

    complex(KIND=dp), pointer :: dHsp(:,:), dHspout(:,:)

    integer       :: p, h
    real(KIND=dp) :: occ_h, occ_p

    if (fam_verbose > 1) print *, "iterate_dH :: starting full FAM loop "

    ! pointer remapping for reshaping 1D flat arrays into 2D matrices
    dHsp(1:nwt,1:nwt) => dHsp_flat(:)
    dHspout(1:nwt,1:nwt) => dHspout_flat(:)

    ! get the ph and hp subblocks
    call get_ph_hp_blocks(dHsp, dH(:,:,1), dH(:,:,2))

    ! calculate X and Y from the perturbed dH
    call calculate_XY(dH)

    if (fam_verbose>0) then
      print 1, sqrt(sum( abs(dH(:,:,1))**2) ), sqrt(sum( abs(dH(:,:,2))**2) )
      print 2, sqrt(sum( abs(X(:,:))**2) ), sqrt(sum( abs(Y(:,:))**2) )
      print 3, l,m, omega_fam,  calc_strength()
    endif

    ! Apply simple linear mixing of X and Y. 
    ! call mix_XY_linear(lin_mix_coeff)
    ! -> this may be skipped when using GMRES

    call store_XY_hist()

    ! build the perturbed densities on the mesh dRs, dRa from X and Y
    call build_perturbed_densities(X, Y, dRs, dRa)

    ! explicit linearisation of the fields
    call calc_perturbed_potentials(RUnper, dRs, dRa, dFs, dFa)

    ! We add in all additional contributions to F_I_I that do not 
    !  result from the Skyrme functional.  
    call combine_potentials(dFs)
    call combine_potentials(dFa)

    ! construct the sp hamiltonian
    dHspout = calc_sphamil_me( HFpsi, HFdpsi, HFddpsi,dFs, dFa, .false.)

    if(fam_verbose > 2) call print_all_fam_spmat()

  end subroutine iterate_dHsp

  subroutine one_minus_T(dHsp_flat, dHspout_flat)
    !---------------------------------------------------------------------------
    ! The precedure iterate_dH constitutes an affine transformation 
    !    dH -> T(dH) + dH_free
    ! where T is a linear map. Fixed-point solutions of this affine problem are
    ! therefor also a solution of the standard linear problem
    !    (I - T) dH = dH_free. 
    ! Thus, (I-T) is the relevant linear operator to use in GMRES. One must 
    ! therefore compute
    !    (I - T) dH = dH - (T(dH) + dH_free) + dH_free 
    !               = dH - iterate_dH(dH) +  dH_free 
    !---------------------------------------------------------------------------

    implicit none
    complex(KIND=dp), dimension(:), intent(in)   :: dHsp_flat
    complex(KIND=dp), dimension(:), intent(out)  :: dHspout_flat

    if (fam_verbose > 1) print *, "compute (I-T) (dH)"

    call iterate_dHsp(dHsp_flat, dHspout_flat)

    dHspout_flat = dHsp_flat - dHspout_flat + dH_free_flat

    if(fam_verbose > 0) then
      print * , "||H_in||",   norm_dH(dHsp_flat)
      print * , "||H_out||",   norm_dH(dHspout_flat)
      print * , "||dH_free||",   norm_dH(dH_free_flat)
    endif

  end subroutine one_minus_T


  subroutine calculate_XY(dH)
    !---------------------------------------------------------------------------
    ! Compute the X and Y amplitudes from the FAM master equation
    !---------------------------------------------------------------------------
    implicit none
    complex(KIND=dp), intent(in)  :: dH(:,:,:) ! perturbed H in QP basis

    integer       :: p, h
    real(KIND=dp) :: occ_h, occ_p, e_h, e_p

    if (fam_verbose > 1) print *, "calculate_XY :: update X and Y"


    X = - (F(:,:,1) + dH(:,:,1))
    Y = - (F(:,:,2) + dH(:,:,2))

    ! normalise with energy denominator
    do h = 1, nwt
      occ_h = rho_can(h)
      e_h = spenergies(h) 
      if(occ_h < 1d-6) cycle
      do p = 1, nwt
$TR         occ_p = 2.0d0 - rho_can(p) ! degeneracy is 2 when T is conserved
$NTR        occ_p = 1.0d0 - rho_can(p) ! degeneracy is 1 when T is broken
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
    if (fam_verbose > 1) print *, "store_XY_hist :: store X and Y in hostory"

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

    if (fam_verbose > 1) print *, "mix_XY_linear :: linear mixing of X and Y with alpha=", alpha

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

    if (fam_verbose > 1) print *, "build_perturbed_densities :: "

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
    complex(KIND=dp), allocatable :: dHsp(:,:)

    if (fam_verbose > 1) print *, "build_dH_explicit :: "


    allocate(dHsp(nwt,nwt))

    ! explicit linearisation of the fields
    call calc_perturbed_potentials(R, dRs, dRa, dFs, dFa)

    call combine_potentials(dFs)
    call combine_potentials(dFa)

    ! construct the sp hamiltonian
    dHsp = calc_sphamil_me( HFpsi, HFdpsi, HFddpsi,dFs, dFa, .false.)
!
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
    !---------------------------------------------------------------------------

    complex(KIND=dp) :: S
    real(KIND=dp) S_out
    integer :: h, p
    real(KIND=dp) :: occ_h, occ_p

    if (fam_verbose > 1) print *, "calc_strength :: S_lm where l= ", l, "m=", m


    S = 0
    do h = 1, nwt
      occ_h = rho_can(h)
      if(occ_h < 1d-6) cycle
      do p = 1, nwt
$TR         occ_p = 2.0d0 - rho_can(p) 
$NTR        occ_p = 1.0d0 - rho_can(p) 
        if(occ_p < 1d-6) cycle
        S = S + conjg(F(p,h,1)) * X(p,h) + conjg(F(p,h,2)) * Y(p,h)
      enddo
    enddo

$TR S = 2 * S ! Time-reversal factor 2

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
    2 format('Convergence: ', '||FAM(X) - X||/||X|| = ', es10.3, '     ||FAM(Y) - Y||/||Y|| = ', es10.3)
    logical, intent(out) :: conv, div
    integer :: idx_prev
    real(KIND=dp) :: DX_norm, DY_norm, X_norm, Y_norm

    if (fam_verbose > 1) print *, "test_convergence :: "


    conv = .false.
    div = .false.

    X_norm = sqrt(sum( abs(X_hist(hist_current_idx,:,:))**2))
    Y_norm = sqrt(sum( abs(Y_hist(hist_current_idx,:,:))**2))


    if( (X_norm .ge. 1.0d3) .or. (Y_norm .ge. 1.0d3)) then
      div = .true.
    endif

    ! previous index in hist obtained by rolling back twice and adding one
    idx_prev = modulo(hist_current_idx - 2, hist_max) + 1

    DX_norm = sqrt( sum( abs(X_hist(hist_current_idx,:,:) - X_hist(idx_prev,:,:))**2) )
    DX_norm = DX_norm / X_norm

    DY_norm = sqrt( sum( abs(Y_hist(hist_current_idx,:,:) - Y_hist(idx_prev,:,:))**2) )
    DY_norm = DY_norm / Y_norm

    if (fam_verbose > 1) print 2, DX_norm, DY_norm

    if( (DX_norm < fam_precision) .and. (DY_norm < fam_precision)) then
      conv = .true.
    endif

  end subroutine test_convergence


  function get_f_LK(L, K, eff_e_n, eff_e_p) result (f_LK_ph_hp)
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Get the particle-hole and hole-particle matrix elements of the multipole
    ! transition operators f_LK where f_LK(i,j) = < i | r^L Y_LK | j > 
    ! while the monopole operator is Q_00(i,j) = < i | r^2 Y_00 | j > 
    ! Only operational for even L at this point. s
    !
    ! INPUT:
    !     L, K          : multipolarity of the perturbing operator
    !     eff_e_n  : effective charge of neutrons (in units of e)
    !     eff_e_p  : effective charge of protons (in units of e)
    ! 
    ! REMARKS:  
    !   - Note that the code works with Re(Y_LK) and Im(Y_LK) which are NOT normalised; 
    !     they integrate to 1/2 when K != 0. 
    !
    !   - We define f^+_LK = 1/sqrt(2) r^L ( Y_LK + Y_L-K) = sqrt(2) * r^L Re(Y_LK), 
    !     when K = 2n > 0, which are normalised such that |f^+_LK|^2 integrates to 1
    !     over the unit sphere. 
    !     The code gives back f^+_LK for now. Since f_LK and f_L-K would give identical strengths 
    !     for axial even-even nuclei when L is even, f^-=0.
    !
    !   - Note that if eff_e_n = eff_e_p, the operator is of isoscalar type, 
    !     if eff_e_n=-eff_e_p, the operator purely isovector. In certain 
    !     applications, e.g. isovector dipole excitation, one choses eff_e_p = N/A
    !     and eff_e_n = -Z/A such that one only has eff_e_n ~ - eff_e_p, 
    !     but still calls the operator isovector. 
    ! 
    !   - One might add a normalisation to the external field F -> F / alpha in order to have 
    !     dh_free of order 1. Due to linearity of all FAM steps, this then needs to be 
    !     compensated as X -> alpha X , Y -> alpha Y, dh -> alpha * dh, ..., and 
    !     S -> alpha^2 S
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


    implicit none
    complex(KIND=dp), allocatable :: f_LK_ph_hp(:,:,:)
    integer, intent(in) :: L, K
    real(KIND=dp), intent(in) :: eff_e_n, eff_e_p
    logical :: ImPart
    complex(KIND=dp), allocatable :: f_LK_spme(:,:)
      
    allocate(f_LK_spme(nwt,nwt)) 
    allocate(f_LK_ph_hp(nwt,nwt,2)) 

   

    if(l==0) then
      f_LK_spme = Rsq_spme()
    else 
      ! Calling a function in fission_MOI.f90, which returns <i|r^L Re(Y_LK)|j> 
      ! in strange fission units barn^(l/2) = (100 fm^2)^(l/2)
      ImPart = .false. ! real (.false.) , imaginary (.true.) 
      f_LK_spme = Qlm_spme(L, K, ImPart)
      
      ! Rescale f_LK_spme to express in units of fm^l
      f_LK_spme = f_LK_spme * (100**(l/2.0)) 

      ! Renormalise with sqrt(2) if K is not 0
      if(K.ne.0) f_LK_spme = f_LK_spme * sqrt(2.0)


      ! Multiply the operator by the effective charges 
      f_LK_spme(1:nwn,1:nwn) = eff_e_n * f_LK_spme(1:nwn,1:nwn)
      f_LK_spme(nwn+1:,nwn+1:) = eff_e_p * f_LK_spme(nwn+1:,nwn+1:)

      ! TODO: investigate signs in Q20 which seems suspicious in O16 nwt24 test case
      ! 3rd row/col in sym block 1 differs in sign wrt blocks 2, 5 and 6. 


      endif
     
      if(fam_verbose > 2) then
        print *, 'f^+_LK'
       call print_spme_complex(f_LK_spme)
     endif

      ! note: 
      !   Stoitsov PRC 84 (2011) normalises the external field by a parameter
      !   alpha converting the units of the perturbation to MeV, and eventually 
      !   devides the obtained strength by alpha. 


      ! TODO: write a general transformation routine from the mesh to any 
      !       single-particle basis

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      ! Define the external field F by selecting the particle-hole and 
      ! hole-particle subblocks of f_LK by multiplying by their 
      ! occupation, i.e. diagonal elements of rho in the canonical basis

    call get_ph_hp_blocks(f_LK_spme, f_LK_ph_hp(:,:,1), f_LK_ph_hp(:,:,2))

    deallocate(f_LK_spme)

  end function


  subroutine get_ph_hp_blocks_complex(M, Mph, Mhp)
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
    complex(KIND=dp), intent(in) :: M(:,:)
    complex(KIND=dp), intent(out) :: Mph(:,:), Mhp(:,:)
    integer       :: p, h
    real(KIND=dp) :: occ_h, occ_p

    Mph = 0
    Mhp = 0

    do h = 1, nwt
      occ_h = rho_can(h)
      if(occ_h < 1d-6) cycle
      do p = 1, nwt
$TR         occ_p = 2.0d0 - rho_can(p)
$NTR        occ_p = 1.0d0 - rho_can(p) 
        if(occ_p < 1d-6) cycle
        Mph(p,h) = M(p,h)
        Mhp(p,h) = M(h,p)
      enddo
    enddo

    ! This can be more efficient by using some mask and elementwise multiplication

    ! For QFAM this will have to be generalised to M20 and M02 obtained from a 
    ! Bogoliubov transformation to the qp basis. 

  end subroutine get_ph_hp_blocks_complex

  subroutine get_ph_hp_blocks_real(M, Mph, Mhp)
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
$TR         occ_p = 2.0d0 - rho_can(p)
$NTR        occ_p = 1.0d0 - rho_can(p) 
        if(occ_p < 1d-6) cycle
        Mph(p,h) = M(p,h)
        Mhp(p,h) = M(h,p)
      enddo
    enddo

    ! This can be more efficient by using some mask and elementwise multiplication

    ! For QFAM this will have to be generalised to M20 and M02 obtained from a
    ! Bogoliubov transformation to the qp basis.

  end subroutine get_ph_hp_blocks_real


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
    complex(KIND=dp), dimension(:), intent(in)  :: dH
    real(KIND=dp)                            :: res

    res = sqrt(sum(abs(dH(:))**2))

  end function

  function ScProd_dH(dHl, dHr) result(res)
    ! abstract template procedure (dH,dH) -> complex required for procedural argument to gmres
    ! to be updated to the objects of the dimensions of the perturbed
    ! sp hamiltonian dh and ddelta (in HF basis)
    complex(KIND=dp), dimension(:), intent(in)  :: dHl, dHr
    complex(KIND=dp)                            :: res

    ! res = sum(conjg(dHl(:)) * dHr(:))
    res = sum(dHl(:) * conjg(dHr(:)))

  end function

  subroutine print_all_fam_spmat()

    print *, 'X'
    call print_spme_complex(X)

    print *, 'Y'
    call print_spme_complex(Y)

    print *, 'drho'
    call print_spme_complex(X + transpose(Y))

!     print *, 'drho_sym'
!     call print_spme_complex(X + transpose(Y) + transpose(X) + Y)
!
!     print *, 'drho_antisym'
!     call print_spme_complex(X + transpose(Y) - transpose(X) - Y)


    print *, 'dH20'
    call print_spme_complex(dH(:,:,1))

    print *, 'dH02'
    call print_spme_complex(dH(:,:,2))

!     print *, 'F20'
!     call print_spme_real(F(:,:,1))
!
!     print *, 'F02'
!     call print_spme_real(F(:,:,2))

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
        print "(*( '(',g12.5,',',g12.5,')',:))",  A(i, si+1:si+N)
      enddo
      print *
      si = si + N
    enddo
    print *
    
  end subroutine print_spme_complex

end module fam
