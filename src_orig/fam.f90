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
  ! A FAM-(Q)RPA implementation to complement MOCCa.
  !------------------------------------------------------------------------------
  ! Hephaestos keywords
  ! 
  ! TR  : $TR
  ! NTR : $NTR
  ! TAUPRESENT : $TAUPRESENT
  !==============================================================================

  use densities
  use moments
  use fission_MOI
  use evolution
  use pairing

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
  ! FAM strength
  complex(KIND=dp) :: strength_complex = CMPLX(0.0_dp,0.0_dp,KIND=dp)
  !    the complex strength S(w,F) = Tr(F^dagger drho(w))
  real(KIND=dp) :: strength = 0.0_dp
  !    the transition strength (aka dB/dw) obtained as - 1/pi * Im(strength_complex)
  real(KIND=dp) :: ewsr = 0.0_dp ! energy weighted sum rule
  !-----------------------------------------------------------------------------
  ! mixing strategy
  integer :: fam_mixingscheme = 0 ! 0 : GMRES (default)
  !                                 1 : linear mixing of dH
  integer :: fam_maxiter = 100 ! maximal number of FAM iterations 
  integer :: fam_maxhist = 30 ! maximal history size of GMRES 
  ! Coefficient for the linear mixing of FAM iterations
  real(KIND=dp) :: fam_lin_mix = 0.3_dp
  !-----------------------------------------------------------------------------
  ! Attention: storage convention for FAM matrices 
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !  When z-signature is conserved, all (Q)FAM matrices come in either of the 
  !   following two forms 
  ! 
  !    A =  ( A_+    0 )     or     B = (  0         B_{+-} )
  !         ( 0    A_- )                ( -B^T_{+-}  0      )
  !
  ! where the division is made between the states with \eta = +-i.
  !
  ! Examples of A-shape are: \delta \rho, \delta H^{11}, ....
  !          of B-shape are: \delta H^{20}, X, Y, ....
  !
  ! When time-reversal is conserved, we do not store the complete matrices 
  ! A and B; rather we store half of them. Our convention is to store 
  !  A_+ and B_{+-}, i.e. the top-left of A and the top-right of B. 
  !
  ! Note that the (Q)FAM equations for the "reduced" objects A_+ and B_{+-}
  ! are not quite identical to the (Q)FAM equations for the complete matrices
  ! A and B; some signs and/or transposes appear in not-so-obvious places.
  !-----------------------------------------------------------------------------
  ! FAM amplitudes X, Y
  complex(KIND=dp), allocatable :: X(:,:) ! forward amplitudes in sp (FAM) or  
  !                                         qp (QFAM) basis, size (nwt,nwt)
  complex(KIND=dp), allocatable :: Y(:,:) ! backward amplitudes in sp (FAM) or 
  !                                         qp (QFAM) basis, size (nwt,nwt)
  ! - In absence of pairing : (FAM)
  !   X and Y contain ph ans hp elements respectively. They are ordered as X(p,h) 
  !   and Y(p,h) /!\ where p is a unoccupied sp index h is an occupied sp index.
  !   They are allocated ove the complete basis size (nwt,nwt)
  ! - In presence of pairing : (QFAM)
  !   X and Y are stored in the quasi-particle basis, their size is (nwt,nwt).
  !
  !-----------------------------------------------------------------------------
  ! Perturbed densities
  ! /!\: perturbations are always RELATIVE to the static mean-field, e.g.
  !         rho(omega) = rho_MF + drho(omega)
  complex(KIND=dp), allocatable :: drho(:,:)         ! perturbation to the normal density matrix in HF basis
  complex(KIND=dp), allocatable :: dkappa_plus(:,:)  ! perturbations to the pairing density matrix in HF basis
  complex(KIND=dp), allocatable :: dkappa_minus(:,:) 
  type(DensityVector)   :: Runper    ! static mean-field densities on the mesh
  type(DensityVector)   :: dRs, dRa  ! perturbation to the particle-hole densities on the mesh
  !                         |    '-> anti-symmetric part
  !                         '-> symmetric part
  type(PotentialVector) :: dFs, dFa  ! perturbation to the particle-hole potentials on the mesh
  !                         |    '-> anti-symmetric part
  !                         '-> symmetric part
  !
  type(DensityVector)   :: dR_pp_plus, dR_pp_minus  ! perturbation to the particle-particle densities on the mesh
  !                         |           '-> associated with kappa_minus
  !                         '-> associated with kappa^plus 
  type(PotentialVector) :: dF_pp_plus, dF_pp_minus  ! perturbation to the particle-particle potentials on the mesh
  !                         |           '-> associated with kappa_minus
  !                         '-> associated with kappa_plus
  !-----------------------------------------------------------------------------
  ! unperturbed Hamiltonian and perturbed hamiltonian
  real(KIND=dp), allocatable :: Hunper(:,:) ! unperturbed Hamiltonian in HF basis
  ! -> currently not used, except for one routine in fam_testing.f90
  complex(KIND=dp), allocatable :: dH(:,:,:) ! perturbed Hamiltonian in sp (FAM) or
  !                                   | | |    qp (QFAM) basis, size (nwt,nwt,2)
  !                                   | | '-> 1 : ph/20 or 2 : hp/02 component 
  !                                   | '-> sp/qp index
  !                                   '-> sp/qp index
  ! - FAM, no pairing: dH(:,:,1) and dH(:,:,2) 
  !     contain ph and hp elements of the perturbed single-particle hamiltonian in the HF basis
  !     They are allocated over the complete basis size, i.e. dimension (nwt,nwt). 
  ! - QFAM, pairing: dH(:,:,1) and dH(:,:,2) 
  !     contain the 20 and 02 elements of the perturbed HFB hamiltonian in the qp basis. 
  !     These remain matrices of dimension (nwt,nwt).
  !   
  complex(KIND=dp), allocatable :: dH_free_flat(:) ! free response of Hamiltonian in the HF basis
  ! The free response is obtained by performing one complete FAM loop starting from dH=0
  ! - In absence of pairing : (FAM)
  !   dH_free_flat contains the free perturbed sp hamiltonian dh(:,:) in HF basis. 
  !   This is a flat array of length (nwt x nwt).  
  ! - In presence of pairing : (QFAM)
  !   dH_free_flat contains three matrices in a stacked fashion in this order
  !      1.    d\Delta^+ 
  !      2.    dh 
  !      3.  - d\Delta^-,*
  !   The whole is a flat array of length (nwt x nwt x 3).
  !-----------------------------------------------------------------------------
  ! external field
  complex(KIND=dp), allocatable :: F(:,:,:)  ! perturbed external field in sp (FAM) or
  !                                   | | |    qp (QFAM) basis, size (nwt,nwt,2)
  ! - same remark as dH(:,:,:)        | | '-> 1 : ph/20 or 2 : hp/02 component 
  !                                   | '-> sp/qp index
  !                                   '-> sp/qp index
  ! Type of perturbing operator
  !   'multipole'       = multipole moment Q_{\ell m}
  !   'particle number' = particle number operator N
  character(len=20) :: operator_type = 'multipole'
  integer :: l = -1, m = -1 ! angular momentum and projection quantum number of the multipole moment
  real(KIND=dp) :: eff_charge_n = 1.0_dp ! effective charge for neutrons in units of e
  real(KIND=dp) :: eff_charge_p = 1.0_dp ! effective charge for protons in units of e
  !-----------------------------------------------------------------------------
  ! convergence
  complex(KIND=dp), allocatable :: X_hist(:,:,:) ! history of X through FAM iters
  !                                       | | '-> sp/qp index 
  !                                       | '-> sp/qp index
  !                                       '-> history index 
  complex(KIND=dp), allocatable :: Y_hist(:,:,:) ! history of Y through FAM iters
  !                                       | | '-> sp/qp index 
  !                                       | '-> sp/qp index 
  !                                       '-> history index 
  ! => REMARK: would it better to set the last index to be the history for memory contiguity
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
  ! XYtoF
  logical :: XYtoF = .false. ! compute F starting from X and Y. (default = .false.)
  !   This corresponds to the inverse problem of FAM and is mush easier to solve.
  !   It comes down to reading in X and Y and mutiplying with the QRPA matrix, 
  !   which can be achieved by ONE partial FAM iteration. 
  !-----------------------------------------------------------------------------
  ! verbosity
  integer :: fam_verbose = 1
  ! 0: no printing. Used during GMRES as output would be confusing
  ! 1: limited printing. Used in the final FAM iteration once GMRES is converged
  !    (default)
  ! 2: printing all function calls. Useful for debugging. 
  ! 3: printing all sp matrices at each iteration. Useful for debugging. 
  !-----------------------------------------------------------------------------
  ! Run the unit tests on start-up; this will not result in a FAM calculation!
  logical :: unit_test = .false.

  interface get_ph_hp_blocks
    module procedure get_ph_hp_blocks_complex
    module procedure get_ph_hp_blocks_real
  end interface get_ph_hp_blocks

contains

  subroutine inifam(omega, DensUnper, PotUnper, Finfile)
    !---------------------------------------------------------------------------
    ! Allocate the FAM objects and set the external field F. X, Y and perturbed 
    ! densities, fields and strength are computed from the free response, i.e. one 
    ! FAM loop starting from dH20 = dH02 = 0. 
    !
    ! Input:
    !    omega      : frequency of the perturbing field
    !    DensUnper  : unperturbed densities on the mesh
    !    PotUnper   : unperturbed potentials on the mesh
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in)          :: omega
    type(DensityVector), intent(in)    :: DensUnper
    type(PotentialVector), intent(in)  :: PotUnper
    character(len=*), intent(in)       :: Finfile



    1 format(' S_',i1,i1,' (', f5.2, ') = ', es10.3)

    print *, "Initialise FAM matrices" 

    ! set omega frequency of perturbation
    omega_fam = omega

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! initialise the external field F
    if(.not.allocated(F)) then
      if (Finfile .ne. '') then 
        F = read_f(Finfile)
      else
        select case(trim(to_lower(operator_type)))
         ! lower to make the selection case insensitive
         ! trim to not bother about string length and possible trailing spaces
        case('multipole')
          F = get_f_LK(l, m, eff_charge_n, eff_charge_p)
        case('particle number')
          F = get_N(eff_charge_n, eff_charge_p)
        case DEFAULT
          call stp('Unrecognized operator_type!')
        end select
      endif
    endif

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! initialise perturbed Hamiltonian as 0
    if(.not.allocated(dH)) then 
      allocate(dH(nwt,nwt,2)) ! stores dHph (dH20), dHhp (dH02) in HF(B) basis for (Q)FAM
    endif

    dH = 0

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! initialise the free response of the sp hamiltonian
    if(.not.allocated(dH_free_flat)) then 
      if(pairingtype==0) then ! FAM
        allocate(dH_free_flat(nwt * nwt)) ! stores dh in HF basis
      else ! QFAM
        allocate(dH_free_flat(3 * nwt * nwt)) ! stores dh, d\Delta^+, -d\Delta^{-,*} in HF basis
      endif
    endif

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! initialise X and Y amplitudes and their history
    if(.not.allocated(X)) then
      allocate(X(nwt,nwt)) ! stores Xph (X20) in HF(B) basis for (Q)FAM
      allocate(Y(nwt,nwt)) ! stores Yhp (Y02) in HF(B) basis for (Q)FAM
    endif

    X = 0
    Y = 0

    if(.not.allocated(X_hist)) then
      allocate(X_hist(hist_max,nwt,nwt)) 
      allocate(Y_hist(hist_max,nwt,nwt))
    endif

    X_hist = 0
    Y_hist = 0

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! store the unperturbed densities
    Runper = DensUnper

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! initialise the perturbed densities
    if(.not.allocated(drho)) then 
      allocate(drho(nwt,nwt))
      allocate(dkappa_plus(nwt,nwt), dkappa_minus(nwt,nwt))
      ! todo : do not allocate dkappa in absence of pairing
      !        this requires to modify densit_offdiag to optional arguments
      drho         = 0.0d0 
      dkappa_plus  = 0.0d0 
      dkappa_minus = 0.0d0
    endif
  
  end subroutine inifam


  subroutine readfam(file_number)
    !---------------------------------------------------------------------------
    ! Read the namelist &fam/.
    !
    ! Input:
    !     file_number : channel number of opened file where to read from.
    !                   Optional. If not present, read from STDIN.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Note that I have to be careful here with default values. in case the variable 
    ! name differs from the one in input data, it gets always overwritten
    !---------------------------------------------------------------------------
    integer(dp), intent(in), optional :: file_number
    real(KIND=dp) :: omega = -1.0_dp
    integer       :: mixingscheme = 0
    integer       :: maxiter = 100, maxhist = 30

    namelist /fam/  omega, omega_min, omega_max, omega_step, smear, maxiter, &
    &               maxhist, l, m, fam_precision, mixingscheme, fam_lin_mix, &
    &               eff_charge_n, eff_charge_p, XYtoF, unit_test, operator_type

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
    5 format (' Compute F from X and Y')

    print 1
    print 2, omega_min, omega_max, omega_step, smear
    if (XYtoF) then
      print 5
    else
      print 3, l, m, eff_charge_n, eff_charge_p
      if (fam_mixingscheme==0) print 41, fam_maxhist, fam_maxiter, fam_precision
      if (fam_mixingscheme==1) print 42, fam_lin_mix, fam_maxiter, fam_precision
    endif
  
  end subroutine printfam

  subroutine iterate_dHsp(dHsp_flat, dHspout_flat)
    !---------------------------------------------------------------------------
    ! Perform one (Q)FAM loop.
    ! 
    ! Input:
    !    dHsp_flat    : perturbed hamiltonian in HF basis as a flat array
    ! Output:
    !    dHspout_flat : iterated perturbed hamiltonian in HF basis as a flat array
    !
    ! Reminder: for QFAM, the arrays store -d\Delta^{-,*}, NOT d\Delta.
    !           This is because GMRES works for LINEAR problems; the QFAM 
    !           equations are a linear function of -d\Delta^{-,*} but NOT of 
    !           d\Delta^{-}.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! 
    ! One full FAM iteration consists of 6 steps : 
    ! 
    ! (1) transform dh, ddelta+/- to QP basis               => dH20, dH02
    ! (2) compute XY from linear response equation          => X   , Y 
    ! (3) transform XY back to sp basis                     => drho, dkappa+/-
    ! (4) calculate perturbed densities on the mesh         => dRs, dRa (DensityVector)
    ! (5) compute perturbed fields on the mesh              => dFs, dFa (PotentialVector)
    ! (6) compute perturbed sp hamiltonian and pairing gaps => dh, d\Delta^{+}, -d\Delta^{-,*}
    !
    ! These are executed by calling two larger routines 
    !   FAM_dh_to_XY => steps (1) to (2)
    !   FAM_XY_to_dh => steps (3) to (6)
    !---------------------------------------------------------------------------
    1 format('||dH20||² = ', es10.3, '     ||dH02||² = ', es10.3)
    12 format('||dh||² = ', es10.3, '     ||ddelta+||² = ', es10.3, '     ||ddelta-||² = ', es10.3)
    2 format('||X||² = ', es10.3, '     ||Y||² = ', es10.3)
    22 format('||drho||² = ', es10.3, '     ||dkappa+||² = ', es10.3, '     ||dkappa-||² = ', es10.3)
    3 format(' S_',i1,i1,' (', f5.2, ') = ', es18.8)

    complex(KIND=dp), dimension(:), target, intent(in)   :: dHsp_flat
    complex(KIND=dp), dimension(:), target, intent(out)  :: dHspout_flat

    complex(KIND=dp), pointer :: dHsp(:,:,:), dHspout(:,:,:)

    real(KIND=dp) :: strength

    integer :: si, i, B, N, N2, T
    complex :: gauge 

    if (fam_verbose > 1) print *, "iterate_dH :: starting full FAM loop "

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! (1) unpack the flat vector to dh, ddelta+/- and transform to QP basis dH20 dH02

    if (pairingtype==0) then ! FAM

      ! pointer remapping for reshaping 1D flat arrays into one 2D matrices
      ! in absence of pairing, dHsp contains only normal field dhsp of sp hamiltonian in HF basis
      dHsp(1:nwt,1:nwt,1:1)    => dHsp_flat(:)
      dHspout(1:nwt,1:nwt,1:1) => dHspout_flat(:)

      ! get the ph and hp subblocks of the perturbed sp hamiltonian
      call get_ph_hp_blocks(dHsp(:,:,1), dH(:,:,1), dH(:,:,2))

      if (fam_verbose>1) print 1, sum( abs(dH(:,:,1))**2) , sum( abs(dH(:,:,2))**2) 

    else ! QFAM

      ! pointer remapping for reshaping 1D flat arrays into three 2D matrices
      !    dH(:,:,1) =  d\Delta^{-,*} 
      !    dH(:,:,2) = dh             
      !    dH(:,:,3) = -d\Delta^{-,*}  
      dHsp   (1:nwt,1:nwt,1:3) => dHsp_flat(:)
      dHspout(1:nwt,1:nwt,1:3) => dHspout_flat(:)

      ! transform the perturbed hamiltonian to the qp basis only interested in dH20 and dH02 components
      call transform_sp_to_qp(Bogoliubov, OTRsp=dHsp(:,:,1), OTLsp=dHsp(:,:,2), OBLsp=dHsp(:,:,3), & ! input 
      &                                   OTRqp=dH(:,:,1),   OBLqp=dH(:,:,2))                        ! output
      ! Attention: 1. there is NO (-CONJG) operation for dHsp(:,:,3), because this 
      !               routine takes -d\Delta^{-,*} as input.
      !            2. the routine spits out the 'bottom left' block of the full matrix, 
      !               which is \delta H^{02, T}. We could apply a transpose, but this 
      !               matrix is antisymmetric; exchanging signs works in both T-conserved 
      !               and T-broken cases.
      dH(:,:,2) = -         dH(:,:,2)
    endif

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! (2) compute X and Y amplitudes from linear response equation

    call calculate_XY(dH,X,Y)
    ! .... and store them into history
    call store_XY_hist(X,Y)

    if (fam_verbose>1) then
      print *, 'Verify antisymmetry of X and Y'
      print * , '||X + X^T|| = ', sum(abs(X+transpose(X))**2)
      print * , '||Y + Y^T|| = ', sum(abs(Y+transpose(Y))**2)
    endif

    if (fam_verbose>0) then
      print 2, sum( abs(X(:,:))**2) , sum( abs(Y(:,:))**2) 
      strength =  calc_strength()
      print 3, l,m, omega_fam, strength
    endif
  
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! (3) Obtain perturbed (pairing) density matrices in HF basis
    if (pairingtype==0) then ! FAM
      drho = X  + transpose(Y)
      dkappa_plus  = 0  
      dkappa_minus = 0
    else                     ! QFAM
      ! Attention for the subtleties of time-reversal invariance. 
      !
      ! When time-reversal is not conserved, i.e. in the "complete" calculation: 
      ! the perturbed density matrix is
      !
      !     \delta \mathcal{R} =    (   0     X   )
      !                             (   Y^T   0   )
      !
      !    which is why the input to the routine has Y^T

$NTR  call transform_qp_to_sp(Bogoliubov, OTRqp=X, OBLqp=transpose(Y), &  
$NTR  &                       OTRsp=dkappa_plus, OTLsp=drho, OBLsp=dkappa_minus)
    
      ! When time-reversal is conserved: the perturbed density matrix is 
      !
      !       \delta \mathcal{R} =    (   0     X   )   = ( 0  0        0    X_r )
      !                               (   Y^T   0   )     ( 0  0       -X_r  0   )
      !                                                   ( 0  -Y^T_r   0    0   )
      !                                                   ( Y^T_r       0    0   )
      !
      !     because we store the top-right part Y_r of  Y = ( 0   Y_r ) 
      !                                                     (-Y_r 0   )
      !     and similar for X. 
      !

$TR   call transform_qp_to_sp(Bogoliubov, OTRqp=X, OBLqp=-Y, &  
$TR   &                       OTRsp=dkappa_plus, OTLsp=drho, OBLsp=dkappa_minus)

      ! The output of the qp -> sp transformation should be:
      ! 
      !      (  d\rho            d\kappa^+ )
      !      ( -d\kappa^{-,*}    -d\rho^{T}) 
      ! 
      ! Notice that the 'bottom left' corner of the resulting matrix is 
      !       -\delta \kappa^{-,*} 
      ! hence the -(CONJG) operation here, independent of whether or not
      ! time-reversal is conserved.
      dkappa_minus = - CONJG(dkappa_minus)
    endif
    
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! (4) Compute perturbed densities on the mesh
    call densit_offdiag(drho, dkappa_plus, dkappa_minus, dRs, dRa, dR_pp_plus, dR_pp_minus)

    if (fam_verbose > 0) then
      if(pairingtype==0) then
        print 22,  sum(abs(drho)**2), 0.0,  0.0
      else
        print 22,  sum(abs(drho)**2), sum(abs(dkappa_plus)**2),  sum(abs(dkappa_minus)**2)
      endif
    endif

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! (5) compute perturbed fields on the mesh

    ! explicit linearisation of the fields
    call calc_perturbed_potentials(RUnper, dRs, dRa, dR_pp_plus, dR_pp_minus, &
    &                                      dFs, dFa, dF_pp_plus, dF_pp_minus)
    ! We add in all additional contributions to F_I_I that do not 
    !  result from the Skyrme functional.  
    call combine_potentials(dFs)
    call combine_potentials(dFa)
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! (6) calculate perturbed hamiltonian and pairing in the HF basis
    if (pairingtype==0) then ! FAM
      
      ! construct the sp hamiltonian in HF basis
      dHspout(:,:,1) = calc_sphamil_me( HFpsi, HFdpsi, HFddpsi,dFs, dFa, .false.)

      if (fam_verbose > 0) print 12,  sum(abs(dHsp(:,:,1))**2), 0.0, 0.0

    else ! QFAM
      ! construct the sp hamiltonian + pairing fields in HF basis
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! \delta h, perturbation of the single-particle hamiltonian
      dHspout(:,:,2) = calc_sphamil_me( HFpsi, HFdpsi, HFddpsi, dFs, dFa, .false.)

      ! \delta \Delta^{+} -> stored 'as-is'
      dHspout(:,:,1) = calc_delta_me(   HFpsi, HFdpsi, HFddpsi, dF_pp_plus , .false.)

      ! \delta \Delta^{-}
      dHspout(:,:,3) = calc_delta_me(   HFpsi, HFdpsi, HFddpsi, dF_pp_minus, .false.)
      ! -> stored as (- \delta \Delta^{-,*} )
      dHspout(:,:,3) = - CONJG(dHspout(:,:,3))

    endif

    if(fam_verbose > 2) then
       call print_all_fam_spmat()
       print *, 'OUTPUT of dHSP_iterate'
       print *, 'dh'
       print *, '||dh||² = ', sum(abs(dHsp(:,:,2))**2)
       if (pairingtype==0) then
          call print_spme_complex(dHsp(:,:,2))
       else
          call print_spme_complex_superblock(dHsp(:,:,2))
       endif
       print *, 'dDelta+'
       print *, '||dDelta+||² = ', sum(abs(dHsp(:,:,1))**2)
       if (pairingtype==0) then
          call print_spme_complex(dHsp(:,:,1))
       else
          call print_spme_complex_superblock(dHsp(:,:,1))
       endif
       print *, 'dDelta-'
       print *, '||dDelta-||² = ', sum(abs(dHsp(:,:,3))**2)
       if (pairingtype==0) then
          call print_spme_complex(dHsp(:,:,3))
       else
          call print_spme_complex_superblock(dHsp(:,:,3))
       endif
    endif
  end subroutine iterate_dHsp

  subroutine FAM_dh_to_XY(dHsp_flat, X_local, Y_local)
    !---------------------------------------------------------------------------
    ! Obtain the X-Y (Q)FAM amplitudes starting from the induced mean-fields. 
    !
    ! Input:
    !    dHsp_flat        : perturbed hamiltonian in HF basis as a flat array
    ! Output:
    !    X_local, Y_local : (Q)FAM amplitudes
    !
    ! Reminder: for QFAM, the arrays store -d\Delta^{-,*}, NOT d\Delta.
    !           This is because GMRES works for LINEAR problems; the QFAM 
    !           equations are a linear function of -d\Delta^{-,*} but NOT of 
    !           d\Delta^{-}.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! 
    ! These are two of the steps from a full FAM iteration : 
    ! 
    ! (1) transform dh, ddelta+/- to QP basis               => dH20, dH02
    ! (2) compute XY from linear response equation          => X   , Y 
    !---------------------------------------------------------------------------
    1 format('||dH20||² = ', es10.3, '     ||dH02||² = ', es10.3)
    12 format('||dh||² = ', es10.3, '     ||ddelta+||² = ', es10.3, '     ||ddelta-||² = ', es10.3)
    2 format('||X||² = ', es10.3, '     ||Y||² = ', es10.3)
    22 format('||drho||² = ', es10.3, '     ||dkappa+||² = ', es10.3, '     ||dkappa-||² = ', es10.3)
    3 format(' S_',i1,i1,' (', f5.2, ') = ', es18.8)

    complex(KIND=dp), intent(out)          :: X_local(:,:), Y_local(:,:)
    complex(KIND=dp), target, intent(in)   :: dHsp_flat(:)
    complex(KIND=dp), pointer              :: dHsp(:,:,:)

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! (1) unpack the flat vector to dh, ddelta+/- and transform to QP basis dH20 dH02

    if (pairingtype==0) then ! FAM

      ! pointer remapping for reshaping 1D flat arrays into one 2D matrices
      ! in absence of pairing, dHsp contains only normal field dhsp of sp hamiltonian in HF basis
      dHsp(1:nwt,1:nwt,1:1)    => dHsp_flat(:)

      ! get the ph and hp subblocks of the perturbed sp hamiltonian
      call get_ph_hp_blocks(dHsp(:,:,1), dH(:,:,1), dH(:,:,2))

      if (fam_verbose>1) print 1, sum( abs(dH(:,:,1))**2) , sum( abs(dH(:,:,2))**2) 

    else ! QFAM

      ! pointer remapping for reshaping 1D flat arrays into three 2D matrices
      !    dH(:,:,1) =  d\Delta^{-,*} 
      !    dH(:,:,2) = dh             
      !    dH(:,:,3) = -d\Delta^{-,*}  
      dHsp   (1:nwt,1:nwt,1:3) => dHsp_flat(:)

      ! transform the perturbed hamiltonian to the qp basis only interested in dH20 and dH02 components
      call transform_sp_to_qp(Bogoliubov, OTRsp=dHsp(:,:,1), OTLsp=dHsp(:,:,2), OBLsp=dHsp(:,:,3), & ! input 
      &                                   OTRqp=dH(:,:,1),   OBLqp=dH(:,:,2))                        ! output
      ! Attention: 1. there is NO (-CONJG) operation for dHsp(:,:,3), because this 
      !               routine takes -d\Delta^{-,*} as input.
      !            2. the routine spits out the 'bottom left' block of the full matrix, 
      !               which is \delta H^{02, T}. We could apply a transpose, but this 
      !               matrix is antisymmetric; exchanging signs works in both T-conserved 
      !               and T-broken cases.
      dH(:,:,2) = -         dH(:,:,2)
    endif

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! (2) compute X and Y amplitudes from linear response equation

    call calculate_XY(dH,X,Y)

  end subroutine FAM_dh_to_XY
  !subroutine partial_FAM_XY_to_dH(X, Y, dHspout_flat)
  !  !---------------------------------------------------------------------------
  !  ! Perform a partial FAM loop, starting from X and Y get the induced perturbed
  !  ! Hamiltonian
  !  ! 
  !  ! Input:
  !  !    X, Y         :  X Y FAM amplitudes 
  !  ! Output:
  !  !    dHspout_flat : iterated perturbed hamiltonian in HF basis as a flat array
  !  ! 
  !  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  !  !
  !  ! Note: 
  !  !   In FAM dHsp_flat only contains dh, while in QFAM it stacks dh, ddelta+ and ddelta-
  !  ! 
  !  !   This function only is almost copy-paste of iterate_dHsp, skipping steps (1) and (2)
  !  ! 
  !  ! (3) transform XY to sp basis                    => drho, dkappa+/-
  !  ! (4) calculate perturbed densities on the mesh   => dRs, dRa (DensityVector)
  !  ! (5) compute perturbed fields on the mesh        => dFs, dFa (PotentialVector)
  !  ! (6) compute perturbed sp hamiltonian and paring => dh, ddelta+/-
  !  !
  !  !---------------------------------------------------------------------------
  !  
  !  1 format('||dH20||² = ', es10.3, '     ||dH02||² = ', es10.3)
  !  12 format('||dh||² = ', es10.3, '     ||ddelta+||² = ', es10.3, '     ||ddelta-||² = ', es10.3)
  !  2 format('||X||² = ', es10.3, '     ||Y||² = ', es10.3)
  !  22 format('||drho||² = ', es10.3, '     ||dkappa+||² = ', es10.3, '     ||dkappa-||² = ', es10.3)
  !  3 format(' S_',i1,i1,' (', f5.2, ') = ', es18.8)
!
!    complex(KIND=dp), dimension(:,:), intent(in) :: X, Y
!    complex(KIND=dp), dimension(:), target, intent(out)  :: dHspout_flat
!    complex(KIND=dp), pointer :: dHspout(:,:,:)!!
!!!

!    if (fam_verbose > 1) print *, "partial_FAM_XY_to_dH :: starting partial FAM loop from X and Y"!!
!
!    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!    ! Set up the pointer remap!
!
!    if (pairingtype==0) then ! FAM!
!
!      ! pointer remapping for reshaping 1D flat arrays into one 2D matrices
!      ! in absence of pairing, dHsp contains only normal field dhsp of sp hamiltonian in HF basis
 !     dHspout(1:nwt,1:nwt,1:1) => dHspout_flat(:)


!    else ! QFAM

!      ! pointer remapping for reshaping 1D flat arrays into three 2D matrices
!      ! dHsp contains sp hamiltonian in HF basis: normal field + two pairing fields [ddelta+, dh, ddelta-]
!      !    dH(:,:,1) = ddelta+ = dH20, dH(:,:,2) = dh = dH11, dH(:,:,3) = ddelta- = dH02 
!      dHspout(1:nwt,1:nwt,1:3) => dHspout_flat(:)!
!
!    endif

!    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!    ! (3) Obtain perturbed (pairing) density matrices in HF basis
!    if (pairingtype==0) then ! FAM
!      drho = X  + transpose(Y)
!      dkappa_plus  = 0  
!      dkappa_minus = 0
!    else ! QFAM
!      call transform_qp_to_sp(Bogoliubov, O20qp=X, O02qp=Y, O20sp=dkappa_plus, O11sp=drho, O02sp=dkappa_minus)
!    endif!!

!    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!    ! (4) Compute perturbed densities on the mesh
!    call densit_offdiag(drho, dkappa_plus, dkappa_minus, dRs, dRa, dR_pp_plus, dR_pp_minus)
!
!    if (fam_verbose > 0) then
!      if(pairingtype==0) then
!        print 22,  sum(abs(drho)**2), 0.0,  0.0
!      else
!        print 22,  sum(abs(drho)**2), sum(abs(dkappa_plus)**2),  sum(abs(dkappa_minus)**2)
!      endif
!    endif
!    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!    ! (5) compute perturbed fields on the mesh
!    ! explicit linearisation of the fields
!    call calc_perturbed_potentials(RUnper, dRs, dRa, dR_pp_plus, dR_pp_minus, dFs, dFa, dF_pp_plus, dF_pp_minus)!!

!    ! We add in all additional contributions to F_I_I that do not 
!    !  result from the Skyrme functional.  
!    call combine_potentials(dFs)
!    call combine_potentials(dFa)!

!    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!    ! (6) calculate perturbed hamiltonian and pairing in the HF basis!

!    if (pairingtype==0) then ! FAM
!      
!      ! construct the sp hamiltonian in HF basis
!      dHspout(:,:,1) = calc_sphamil_me( HFpsi, HFdpsi, HFddpsi,dFs, dFa, .false.)!!
!
!    else ! QFAM!!
!
!      ! construct the sp hamiltonian + pairing fields in HF basis
!      dHspout(:,:,1) = calc_delta_me(   HFpsi, HFdpsi, HFddpsi, dF_pp_plus, .false.)
!      dHspout(:,:,2) = calc_sphamil_me( HFpsi, HFdpsi, HFddpsi, dFs, dFa, .false.)
!      dHspout(:,:,3) = calc_delta_me(   HFpsi, HFdpsi, HFddpsi, dF_pp_minus, .false.)!!
!
!    endif!
!
!    if(fam_verbose > 2) call print_all_fam_spmat()
!
!  end subroutine partial_FAM_XY_to_dH


!  subroutine Multiply_XY_with_QRPAmat(X, Y, omega, F, dHsp_flat_in)
!    !---------------------------------------------------------------------------
!    ! Multiply X and Y by the QRPA matrix by performing one adjusted FAM loop. 
!    ! i.e.
!    !        (E - omega) * X + dH20(X, Y) = - F20
!    !        (E + omega) * Y + dH02(X, Y) = - F02
!    !  
!    ! Input:
!    !    X, Y     :  X Y input amplitudes 
!    !    omega    :  frequency used in the linear response
!    !    dHsp_flat_in (optional)  :  flat array of the perturbed hamiltonian in the HF basis 
!    ! Output:
!    !    F20, F02 :  induced external field 
!    ! 
!    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!    !
!    !  (1) Compute the induced perturbed hamiltonian dH in the HF basis by calling 
!    !      partial_FAM_XY_to_dH()            =>  dHsp_flat = {dh, ddelta+, ddelta-}
!    !      -> this step gets skipped if dHsp is passed to this routine 
!    !  (2) Transform dHsp to the QP basis    =>  dH20, dH02
!    !  (3) Compute the resulting field F by linear response, i.e. the Eq. above
!    !
!    !---------------------------------------------------------------------------
!    
!    1 format('||dH20||² = ', es10.3, '     ||dH02||² = ', es10.3)
!    12 format('||dh||² = ', es10.3, '     ||ddelta+||² = ', es10.3, '     ||ddelta-||² = ', es10.3)
!    2 format('||X||² = ', es10.3, '     ||Y||² = ', es10.3)
!    22 format('||drho||² = ', es10.3, '     ||dkappa+||² = ', es10.3, '     ||dkappa-||² = ', es10.3)
!    3 format(' S_',i1,i1,' (', f5.2, ') = ', es18.8)
!
!    implicit none
!    complex(KIND=dp), intent(in) :: X(:,:), Y(:,:)
!    complex(KIND=dp), intent(in) :: omega
!    complex(KIND=dp), intent(out) :: F(:,:,:)
!    complex(KIND=dp), optional, intent(in) :: dHsp_flat_in(:)
!    complex(KIND=dp), allocatable, target :: dHsp_flat(:)
!    complex(KIND=dp), pointer :: dHsp(:,:,:)!!
!
!    if (fam_verbose > 1) print *, "Multiply_with_QRPAmat :: compute the external field induced by XY"!
!
!    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!    ! (1) if not present, compute the induced perturbed hamiltonian dH in the HF basis
!
!    if (.not. present(dHsp_flat_in)) then
!
!      if (.not. allocated(dHsp_flat)) then
!        if(pairingtype==0) then
!          allocate(dHsp_flat(nwt * nwt))
!        else
!          allocate(dHsp_flat(3 * nwt * nwt))
!        endif
!      endif!
!
!      call partial_FAM_XY_to_dH(X, Y, dHsp_flat)
!
!    else
!      dHsp_flat = dHsp_flat_in
!    endif

!   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!   ! (2) unpack via pointer remap and transfrom dH to the qp basis
!
!    if (pairingtype==0) then ! FAM
!
!      ! pointer remapping for reshaping 1D flat arrays into one 2D matrices
!      ! in absence of pairing, dHsp contains only normal field dhsp of sp hamiltonian in HF basis
!      dHsp(1:nwt,1:nwt,1:1) => dHsp_flat(:)

!      ! get the ph and hp subblocks of the perturbed sp hamiltonian
!      call get_ph_hp_blocks(dHsp(:,:,1), dH(:,:,1), dH(:,:,2))!
!
!      if (fam_verbose>1) print 1, sum( abs(dH(:,:,1))**2) , sum( abs(dH(:,:,2))**2) 
!
!    else ! QFAM
!
!      ! pointer remapping for reshaping 1D flat arrays into three 2D matrices
!      ! dHsp contains sp hamiltonian in HF basis: normal field + two pairing fields [ddelta+, dh, ddelta-]
!      !    dH(:,:,1) = ddelta+ = dH20, dH(:,:,2) = dh = dH11, dH(:,:,3) = ddelta- = dH02 
!      dHsp(1:nwt,1:nwt,1:3) => dHsp_flat(:)
!
!      ! transform the perturbed hamiltonian to the qp basis only interested in dH20 and dH02 components
!      call transform_sp_to_qp(Bogoliubov, O20sp=dHsp(:,:,1), O11sp=dHsp(:,:,2), O02sp=dHsp(:,:,3), O20qp=dH(:,:,1), O02qp=dH(:,:,2))
!
!      if (fam_verbose>1) print 1, sum( abs(dH(:,:,1))**2) , sum( abs(dH(:,:,2))**2) 
!
!    endif

!    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!    ! (3) compute induced external field F20 F02
!    
!    call compute_F_from_XYdH(X, Y, dH, omega, F)
!
!
!  end subroutine Multiply_XY_with_QRPAmat

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
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Input:
    !    dHsp_flat    : perturbed sp/qp hamiltonian in HF basis as a flat array
    ! Output:
    !    dHspout_flat : (I - T) * dHsp_flat
    !---------------------------------------------------------------------------
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

  subroutine calculate_XY(dH,X, Y)
    !---------------------------------------------------------------------------
    ! Compute the X and Y amplitudes from the (Q)FAM master equation. 
    !
    ! Note: 
    ! - In the absence of pairing, X, Y, dH, F store the particle-hole subblocks 
    !    and the loops are only over particle-hole pairs. 
    ! - In presence of pairing, X, Y, dH and F store quasiparticle matrix 
    !    elements and the loops run over the complete quasiparticle basis. 
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Input:
    !   dH   : perturbed single-particle (FAM) or HFB (QFAM) Hamiltonian
    !          complex array of size (nwt,nwt,2)
    ! Output:
    !    X,Y : (Q) FAM amplitudes
    !          complex arrays of size (nwt,nwt)
    ! 
    !---------------------------------------------------------------------------
    complex(KIND=dp), intent(out) :: X(nwt,nwt), Y(nwt,nwt) 
    complex(KIND=dp), intent(in)  :: dH(nwt,nwt,2) 

    integer       :: i, j, si, si2, N, N2, B, T, degeneracy

    if (fam_verbose > 1) print *, "calculate_XY :: update X and Y"

    X = - (F(:,:,1) + dH(:,:,1))
    Y = - (F(:,:,2) + dH(:,:,2))

    ! normalise with energy denominator

    if(pairingtype==0) then ! FAM : difference of particle and hole energy
      
      $TR   degeneracy = 2 ! degeneracy of sp states in case of T conservation
      $NTR  degeneracy = 1 ! degeneracy of sp states in case of T conservation
    
      si = 0
      do B=1,8,2
        N  = HFblocks(B)    ; if(N.eq.0) cycle 
        N2 = HFblocks(B+1)
        T = N + N2
        ! run over particle-hole pairs. hole (j) as outer, particle (i) as inner loop
        do j = 1, T
          if(rho_can(si+j) < 1d-6) cycle  ! skip if j is not a hole state
          do i = 1, T
            if(abs(degeneracy - rho_can(si+i)) < 1d-6) cycle  ! skip if i is not a particle state
            X(si+i,si+j) = X(si+i,si+j) / (spenergies(si+i) - spenergies(si+j) - CMPLX(omega_fam,smear,KIND=dp) )
            Y(si+i,si+j) = Y(si+i,si+j) / (spenergies(si+i) - spenergies(si+j) + CMPLX(omega_fam,smear,KIND=dp) )
          enddo
        enddo
        si = si+T
      enddo 
    
    else ! QFAM : sum of two qp energy
    
      ! loop over 4 isospin-parity (IP) block (signature unresolved)
      ! We require two start indices
      ! si  determines the start of the block in qp-basis of dimension nwt   -> X, Y
      ! si2 determines the start of the block in qp-basis of dimension 2*nwt -> qpenergies (-Emax,..., -E1, E1,..., Emax)
      si = 0; si2 = 0
      do B=1,8,2
        N  = HFblocks(B)    ; if(N.eq.0) cycle 
        N2 = HFblocks(B+1)
        T = N + N2
        do j = 1, T
          do i = 1, T

            ! fetch qpenergies from second half (si2 + T), i.e. positive qp spectrum 

            X(si+i,si+j) = X(si+i,si+j) / (qpenergies(si2+T+i) + qpenergies(si2+T+j) - CMPLX(omega_fam,smear,KIND=dp) )
            Y(si+i,si+j) = Y(si+i,si+j) / (qpenergies(si2+T+i) + qpenergies(si2+T+j) + CMPLX(omega_fam,smear,KIND=dp) )
          enddo
        enddo
        si  = si  +   T ! move start index by size of IP block 
        si2 = si2 + 2*T ! move start index by twice the size of IP block
      enddo
    endif

    ! TODO: could be optimised by pre-storing the energy denominator and doing a simple 
    !       elementwise multiplication

    if(fam_verbose > 2) then
      print * , "||X||",   sum(abs(X(:,:)**2))
      print * , "||Y||",   sum(abs(Y(:,:)**2))
    endif

  end subroutine calculate_XY


  subroutine compute_F_from_XYdH(X, Y, dH, omega, F)
    !---------------------------------------------------------------------------
    ! Calculate the external field F induced by X, Y and dH at frequency omega
    ! from the linear response equation. 
    !        F20 = - dH20 -(E - omega) * X 
    !        F02 = - dH02 -(E + omega) * Y
    !---------------------------------------------------------------------------
    implicit none
    complex(KIND=dp), intent(in)  :: dH(:,:,:) ! perturbed H in QP basis
    complex(KIND=dp), intent(in)  :: X(:,:), Y(:,:) ! X, Y in QP basis
    complex(KIND=dp), intent(in)  :: omega ! complex frequency
    complex(KIND=dp), intent(out) :: F(:,:,:) ! induced external field F in QP basis

    integer       :: i, j, si, si2, N, N2, B, T, degeneracy

    if (fam_verbose > 1) print *, "compute_F_from_XYdH ::"

    F = - dH

    if(pairingtype==0) then ! FAM : difference of particle and hole energy

      $TR   degeneracy = 2 ! degeneracy of sp states in case of T conservation
      $NTR  degeneracy = 1 ! degeneracy of sp states in case of T conservation

      si = 0
      do B=1,8,2
        N  = HFblocks(B)    ; if(N.eq.0) cycle 
        N2 = HFblocks(B+1)
        T = N + N2
        ! run over particle-hole pairs. hole (j) as outer, particle (i) as inner loop
        do j = 1, T
          if(rho_can(si+j) < 1d-6) cycle  ! skip if j is not a hole state
          do i = 1, T
            if(abs(degeneracy - rho_can(si+i)) < 1d-6) cycle  ! skip if i is not a particle state

            F(si+i,si+j,1) = F(si+i,si+j,1) - X(si+i,si+j) * (spenergies(si+i) - spenergies(si+j) - omega )
            F(si+i,si+j,2) = F(si+i,si+j,2) - Y(si+i,si+j) * (spenergies(si+i) - spenergies(si+j) + omega )
          enddo
        enddo
        si = si+T
      enddo 

    
    else ! QFAM : sum of two qp energy
    
      ! loop over 4 isospin-parity (IP) block (signature unresolved)
      ! We require two start indices
      ! si  determines the start of the block in qp-basis of dimension nwt   -> X, Y
      ! si2 determines the start of the block in qp-basis of dimension 2*nwt -> qpenergies (-Emax,..., -E1, E1,..., Emax)
      si = 0; si2 = 0
      do B=1,8,2
        N  = HFblocks(B)    ; if(N.eq.0) cycle 
        N2 = HFblocks(B+1)
        T = N + N2
        do j = 1, T
          do i = 1, T

            ! fetch qpenergies from second half (si2 + T), i.e. positive qp spectrum 

            F(si+i,si+j,1) = F(si+i,si+j,1) - X(si+i,si+j) * (qpenergies(si2+T+i) + qpenergies(si2+T+j) - omega )
            F(si+i,si+j,2) = F(si+i,si+j,2) - Y(si+i,si+j) * (qpenergies(si2+T+i) + qpenergies(si2+T+j) + omega )
          enddo
        enddo
        si  = si  +   T ! move start index by size of IP block 
        si2 = si2 + 2*T ! move start index by twice the size of IP block
      enddo
    endif

    if(fam_verbose > 2) then
      print * , "||F20||^2 = ",   sum(abs(F(:,:,1)**2))
      print * , "||F02||^2 = ",   sum(abs(F(:,:,2)**2))
    endif

  end subroutine compute_F_from_XYdH

  subroutine store_XY_hist(X,Y)
    !---------------------------------------------------------------------------
    ! Store the current X and Y into their histories. 
    !
    ! Input:
    !  X, Y: (Q)FAM amplitudes to be stored
    !        complex arrays of size (nwt,nwt)
    !---------------------------------------------------------------------------
    complex(KIND=dp), intent(in) :: X(:,:), Y(:,:)

    if (fam_verbose > 1) print *, "store_XY_hist :: store X and Y in history"

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
      
    if(.not.allocated(rho_pairing)) allocate(rho_pairing(nwt,nwt))
    if(.not.allocated(kappa_pairing)) allocate(kappa_pairing(nwt,nwt))
  
    rho_pairing   = 0
    kappa_pairing = 0
    
    do i=1,nwt
      rho_pairing(i,i) = rho_hf(i)
    enddo
  
  end subroutine iniHFdensities

  function calc_strength() result (res)
    !---------------------------------------------------------------------------
    ! Calculate the strength S(omega,F) and store output in strength and 
    ! strength_complex and return strength
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! strength_complex is defined as 
    !     strength_complex = Tr (F^dagger * drho)
    !                      = sum_ab (F^20_ab^* X_ab + F^02_ab^* Y_ab)
    ! while the strength  
    !     strength = -1/pi * strength_complex
    ! 
    ! note: 
    !  - normalisation of external field may have to be taken into account
    !    S -> S/alpha
    !  - in case of FAM, F(:,:,1) contains the ph block and F(:,:,2) contains
    !    the hp block which differ is F if not Hermitian
    !  - in case of QFAM, F(:,:,1) contains the F20 block in qp basis and F(:,:,2)
    !    contains the F02 block which differ if F is not Hermitian
    !---------------------------------------------------------------------------

    complex(KIND=dp) :: S = 0
    real(KIND=dp) :: res
    integer :: i, j, si, B, N, N2, T

    if (fam_verbose > 1) print *, "calc_strength :: S_lm where l= ", l, "m=", m

    S = 0


    if(pairingtype==0) then ! FAM
      si = 0
      ! loop over 8 isospin-parity-signature (IPS) block 
      do B=1,8
        N  = HFblocks(B)    ; if(N.eq.0) cycle 
        ! run over particle-hole pairs. hole (j) as outer, particle (i) as inner loop
        do j = si+1, si+N
          if(rho_can(j) < 1d-6) cycle  ! skip if j is not a hole state
          do i = si+1, si+N
            S = S + conjg(F(i,j,1)) * X(i,j) + conjg(F(i,j,2)) * Y(i,j)
          enddo
        enddo
        si = si+N
      enddo
    
    else ! QFAM
     ! loop over 4 isospin-parity (IP) block (signature unresolved)
      si = 0
      do B=1,8,2
        N  = HFblocks(B)    ; if(N.eq.0) cycle 
        N2 = HFblocks(B+1)
        T = N + N2
        ! Loop over all possible pairs (i,j)
        ! Note: we do not leverage symmetry here - the representation
        !  of the matrices F,X,Y in memory depends on the the conservation
        !  or breaking of T
        do j = si+1, si+T
          do i = si+1, si+T
             S = S + conjg(F(i,j,1)) * X(i,j) + conjg(F(i,j,2)) * Y(i,j)
          enddo
        enddo
        si  = si + T 
      enddo

    endif
    S = 0.5 * S
    $TR S = 2 * S ! Time-reversal factor 2


    strength_complex = S 
    strength = - IMAG(strength_complex) / pi

    ! return the strength
    res = strength

  end function calc_strength


  subroutine calc_strength_decomp(S_complex, strength)
    !---------------------------------------------------------------------------
    ! Calculate the complex response and the strength decomposed into
    ! different symmetry channels. For now, this assumes that the perturbing
    ! operator must respect all symmetries, i.e. diagonal in tau,pi,z-sign. 
    ! In the future, applying the idea for a non-trivial perturbation operator
    ! would require to loop over the blocks in a (partial) off-diagonal way, 
    ! e.g. pi=-pi' when l is odd. 
    !---------------------------------------------------------------------------

    complex(KIND=dp), intent(out) :: S_complex(8) 
    real(KIND=dp), intent(out) :: strength(8)
    integer :: i, j, B, N, N2, si, T
    real(KIND=dp) :: occ_h, occ_p

    if (fam_verbose > 1) print *, "calc_strength_decomp :: S_lm where l= ", l, "m=", m

    if (mod(l,2) == 1 .or. mod(m,2) == 1) then
      print *, "NOT IMPLEMENTED :: calc_strength_decomp() not applicable when l or m is odd"
      ! print *, "calling calc_strength() instead"
      ! call calc_strength()
      return
    endif

    S_complex = 0
    strength = 0


   if(pairingtype==0) then ! FAM
      si = 0
      ! loop over 8 isospin-parity-signature (IPS) block 
      do B=1,8
        N  = HFblocks(B)    ; if(N.eq.0) cycle 
        ! run over particle-hole pairs. hole (j) as outer, particle (i) as inner loop
        do j = si+1, si+N
          if(rho_can(j) < 1d-6) cycle  ! skip if j is not a hole state
          do i = si+1, si+N
            S_complex(B) = S_complex(B) + conjg(F(i,j,1)) * X(i,j) + conjg(F(i,j,2)) * Y(i,j)
          enddo
        enddo
        si = si+N
      enddo

    
    else ! QFAM
    
      ! loop over 4 isospin-parity (IP) block (signature unresolved)
      si = 0
      do B=1,8,2
        N  = HFblocks(B)    ; if(N.eq.0) cycle 
        N2 = HFblocks(B+1)
        T = N + N2
        ! Loop over all possible combinations of (i,j)
        ! - - - - - - - - - - - - - - - - - - -
        ! Note: we do not leverage symmetry here - the representation
        !  of the matrices F,X,Y in memory depends on the the conservation
        !  or breaking of T
        do j = si+1,si+T
          do i = si+1,si+T
             S_complex(B) = S_complex(B) + conjg(F(i,j,1)) * X(i,j) &
                  &                      + conjg(F(i,j,2)) * Y(i,j)
         enddo
        enddo
        print *
        S_complex(B) = S_complex(B) / 2.0d0

        si  = si + T
      enddo
    endif

$TR    S_complex(:) = 2.0 * S_complex(:) ! Time-reversal factor 2
    strength(:) = - IMAG(S_complex(:)) / pi

    if (fam_verbose > 0) then
      print *, 'Decomposed strength : '
      print * , 'S_n+ : (', strength(1), ' , ', strength(2), ' )'
      print * , 'S_n- : (', strength(3), ' , ', strength(4), ' )'
      print * , 'S_p+ : (', strength(5), ' , ', strength(6), ' )'
      print * , 'S_p- : (', strength(7), ' , ', strength(8), ' )'
      print * , 'S_tot : ', sum(strength(:))
    endif
  end subroutine calc_strength_decomp

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

    if (fam_verbose > 0) print 2, DX_norm, DY_norm

    if( (DX_norm < fam_precision) .and. (DY_norm < fam_precision)) then
      conv = .true.
    endif

  end subroutine test_convergence

  function get_f_LK(L, K, eff_e_n, eff_e_p) result (f_LK_qpme)
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
    integer, intent(in)       :: L, K
    real(KIND=dp), intent(in) :: eff_e_n, eff_e_p
    logical :: ImPart
    complex(KIND=dp), allocatable :: f_LK_qpme(:,:,:)
    complex(KIND=dp), allocatable :: f_LK_spme(:,:)

    if (fam_verbose > 1) print *, "get_f_LK :: "
      
    allocate(f_LK_spme(nwt,nwt)) 
    allocate(f_LK_qpme(nwt,nwt,2)) 

   

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

      endif

    ! TODO: refactor this; selecting particle-hole or quasiparticle parts
    !       of the perturbing operator and (ii) multiplying by effective charges
    !       is independent on our particular choice of perturbing operator and can
    !       thus be made universal...


      ! Multiply the operator by the effective charges 
      f_LK_spme(1:nwn,1:nwn) = eff_e_n * f_LK_spme(1:nwn,1:nwn)
      f_LK_spme(nwn+1:,nwn+1:) = eff_e_p * f_LK_spme(nwn+1:,nwn+1:)

     
      if(fam_verbose > 2) then
        print *, 'f^+_LK'
       call print_spme_complex(f_LK_spme)
       print *, '||f||²', sum(abs(f_LK_spme)**2)
       call print_spme_complex_superblock(f_LK_spme)
     endif

      ! note: 
      !   Stoitsov PRC 84 (2011) normalises the external field by a parameter
      !   alpha converting the units of the perturbation to MeV, and eventually 
      !   devides the obtained strength by alpha. 
    if (pairingtype==0) then ! FAM
      ! Define the external field F by selecting the particle-hole and 
      ! hole-particle subblocks of f_LK by multiplying by their 
      ! occupation, i.e. diagonal elements of rho in the canonical basis
      call get_ph_hp_blocks(f_LK_spme, f_LK_qpme(:,:,1), f_LK_qpme(:,:,2))

      if(fam_verbose > 2) then
        print *, ' f_LK_ph'
        call print_spme_complex_superblock( f_LK_qpme(:,:,1))
        print *, ' f_LK_hp'
        call print_spme_complex_superblock( f_LK_qpme(:,:,2))
      endif

    else ! QFAM
      
      ! Define the external field F as the qpme obtained by performing a bogolibov 
      ! transformation and storing the F^20 anf F^02 components
      call transform_sp_to_qp(Bogoliubov, OTLsp=f_LK_spme, &                                ! Input 
      &                                   OTRqp=f_LK_qpme(:,:,1), OBLqp=f_LK_qpme(:,:,2))   ! Output
      ! TODO: update Attention: the output of this routine is F^{02,T}!
$NTR  f_LK_qpme(:,:,2) = TRANSPOSE(f_LK_qpme(:,:,2))
$TR   f_LK_qpme(:,:,2) = -         f_LK_qpme(:,:,2)

      if(fam_verbose > 2) then
        print *, ' f_LK_qpme(:,:,1)'
        call print_spme_complex_superblock( f_LK_qpme(:,:,1))
        print *, ' f_LK_qpme(:,:,2)'
        call print_spme_complex_superblock( f_LK_qpme(:,:,2))
      endif

    endif

    if(fam_verbose > 1) then

      print *, '||F(:,:,1)||²', sum(abs(f_LK_qpme(:,:,1))**2)
      print *, '||F(:,:,2)||²', sum(abs(f_LK_qpme(:,:,2))**2)

    endif

    deallocate(f_LK_spme)

  end function get_F_LK

  function get_N(eff_e_n, eff_e_p) result(Nqpme)
    !-----------------------------------------------------------------------------
    ! Get the particle-hole and hole-particle matrix elements of the particle 
    ! number operator N.
    ! 
    ! Input:
    !     eff_e_n  : effective charge of neutrons (in units of e)
    !     eff_e_p  : effective charge of protons (in units of e)
    !
    ! Output:
    !     Nqpme    : matrix elements of N, either particle-hole (FAM)
    !                                      or 2qp (QFAM)
    !-----------------------------------------------------------------------------
    real(KIND=dp), intent(in)     :: eff_e_n, eff_e_p
    complex(KIND=dp), allocatable :: Nqpme(:,:,:), Nspme(:,:)
    integer                       :: i

    allocate(Nspme(nwt,nwt))
    allocate(Nqpme(nwt,nwt,2))

    ! Build the single-particle matrix elements of N
    Nspme = 0.0d0
    do i=1,nwt
       Nspme(i,i) = 1.0d0
    enddo

    ! Multiply by effective charges
    Nspme(1:nwn,1:nwn)   = eff_e_n * Nspme(1:nwn,1:nwn)
    Nspme(nwn+1:,nwn+1:) = eff_e_p * Nspme(nwn+1:,nwn+1:)

    if (pairingtype==0) then ! FAM
      call get_ph_hp_blocks(Nspme, Nqpme(:,:,1), Nqpme(:,:,2))
    else ! QFAM
      call transform_sp_to_qp(Bogoliubov, OTLsp=Nspme, &
           &                  OTRqp=Nqpme(:,:,1), OBLqp=Nqpme(:,:,2))
$NTR  Nqpme(:,:,2) = TRANSPOSE(Nqpme(:,:,2))
$TR   Nqpme(:,:,2) = -         Nqpme(:,:,2)
    endif

    if(fam_verbose > 2) then
      print *, ' f_LK_qpme(:,:,1)'
      call print_spme_complex_superblock( Nqpme(:,:,1))
      print *, ' f_LK_qpme(:,:,2)'
      call print_spme_complex_superblock( Nqpme(:,:,2))
    endif


    print *, '||F(:,:,1)||²', sum(abs(Nqpme(:,:,1))**2)
    print *, '||F(:,:,2)||²', sum(abs(Nqpme(:,:,2))**2)

  end function get_N

  function calc_EWSR(R) result (ewsr)
    !---------------------------------------------------------------------------
    ! Compute the energy-weighted sum rule from a ground-state expectation value. 
    ! This value should equal the first-moment of the strength function, i.e.
    !     ewsr = m_1(F) = int_0^inf dE E S(E, F). 
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Notes:
    !   - Expressions are based on N. Hinohara, PRC 100, 024310 (2019)
    !   - The ewsr for EDFs in principle consists of three parts :
    !
    !         ewsr = m1_kin + kappa + G_LGSB
    !
    !     * m1_kin accounts for the kinetic part
    !     * kappa is an enhancement term 
    !          -> absent for isoscalar perturbations
    !     * G_LGSB accounts for local-gauge-symmetry breaking effects. 
    !      /!\ : G_LGSB is currenlty NOT implemented
    !            the EDF param is thus assumed to respect LGS
    !   - Presently only implemented for monopole (L=0,K=0) and quadruole 
    !     perturbations (L=2, K=0) and (L=2, K=2), in fact it is for 
    !     Q22+ = 1/sqrt(2) (Q_22 + Q_2,-2).  
    !   - Isoscalar(vector) character of the perturbation is dealt with via the 
    !     effective charges. While for isoscalar both are positive and 
    !     approximately equal, for isovector effective charges differ in sign but 
    !     are close in magnitude.
    !---------------------------------------------------------------------------
    type(DensityVector), intent(in) :: R
    real(KIND=dp) :: ewsr
    real(KIND=dp) :: m1kin=0, kappa=0, Ctau0=0, Ctau1=0
    type(Moment), pointer  :: moment_ptr, r2_ptr

    ewsr = 0

    if (fam_verbose > 1) print *, "calc_EWSR :: calculate the energy-weighted sum rule m1"

    ! - - - - - - - - - - - - - - - - - - - - - - - - -
    ! l = 0, monopole
    if(l == 0) then
      
      if (fam_verbose > 1) print *, "monopole"
      
      r2_ptr => FindMoment(-2,0,.false.) ! pointer to <r^2>
      m1kin = 4.0 * hbm(1) * (eff_charge_n**2 * r2_ptr%Value(1) + eff_charge_p**2 * r2_ptr%Value(2) )
      ! note that r2_ptr%Value contains a factor N (or Z)


    ! - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! l = 1  dipole (same for m=0 and m=1)
    else if(l == 1) then
      
      if (fam_verbose > 1) print *, "dipole"

      m1kin = (3.0 / (4.0 * pi)) * hbm(1) * ( Neutrons * eff_charge_n**2  + Protons *eff_charge_p**2 ) 
     
    ! - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! l = 2, m = 0, axial quadrupole
    else if(l == 2 .and. m==0) then
      
      if (fam_verbose > 1) print *, "quadrupole Q20"

      ! Ony implemented for axial nuclei, return 0 if beta22 > 0.1 
      moment_ptr => FindMoment(2,2,.false.) ! pointer to <Q_22>
      if (moment_ptr%beta(4) > 0.1) then
        print *, 'NOT IMPLEMENTED: quadrupole EWSR assumes axial shape '
        return
      endif

    
      r2_ptr => FindMoment(-2,0,.false.) ! pointer to <r^2>
      moment_ptr => FindMoment(2,0,.false.) ! pointer to <Q_20>

      m1kin = (5.0 / (2.0 * pi)) * hbm(1) * ( &
        &   eff_charge_n**2 * r2_ptr%Value(1) * (1. + sqrt(5./(4.*pi)) * moment_ptr%beta(1)) &
        & + eff_charge_p**2 * r2_ptr%Value(2) * (1. + sqrt(5./(4.*pi)) * moment_ptr%beta(2))  )
      ! note that r2_ptr%Value contains a factor N (or Z)

    ! - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! l = 2, m = 2, quadrupole
    else if(l == 2 .and. m==2) then
      
      if (fam_verbose > 1) print *, "quadrupole Q22+ = 1/sqrt(2) (Q_22 + Q_2,-2)"

      ! Ony implemented for axial nuclei, return 0 if beta22 > 0.1 
      moment_ptr => FindMoment(2,2,.false.) ! pointer to <Q_22>
      if (moment_ptr%beta(4) > 0.1) then
        print *, 'NOT IMPLEMENTED: quadrupole EWSR assumes axial shape '
        return
      endif

      r2_ptr => FindMoment(-2,0,.false.) ! pointer to <r^2>
      moment_ptr => FindMoment(2,0,.false.) ! pointer to <Q_20> 
      ! I know this seems suspicious but actually do need beta_20 to calculate the sumrule for Q22

      m1kin = (5.0 / (2.0 * pi)) * hbm(1) * ( &
        &   eff_charge_n**2 * r2_ptr%Value(1) * (1. - sqrt(5./(16.*pi)) * moment_ptr%beta(1)) &
        & + eff_charge_p**2 * r2_ptr%Value(2) * (1. - sqrt(5./(16.*pi)) * moment_ptr%beta(2))  )
      ! note that r2_ptr%Value contains a factor N (or Z)

    ! - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! any other l, m
    else 
      print *, 'NOT IMPLEMENTED: only monopole (l=0) and axial quadrupole (l=2, m=0) EWSR implemented for now'
      return
    endif


    !--------------------------------------------------------------
    ! include enhancement factor kappa for isovector pertubations
    !--------------------------------------------------------------

    ! skip if Ctau0 and Ctau1 are not present (e.g. in LO functionals)
#if( $TAUPRESENT )

    if(eff_charge_n * eff_charge_p < 0) then ! IV if product of eff charges is negative

        ! Coupling constant of E_D_I_I_D_Nm_Nm, ['0', '0'] (aka C^tau_0)
        Ctau0 = +1.0/2.0*Cc(t1,x1,+1,0,0)+1.0/2.0*Cc(t2,x2,-1,0,0) 

        ! Coupling constant of E_D_I_I_D_Nm_Nm, ['1', '1'] (aka C^tau_1)
        Ctau1 = +1.0/2.0*Cc(t1,x1,+1,0,1)+1.0/2.0*Cc(t2,x2,-1,0,1) 


        ! init kappa as the common prefactor
        kappa = (Ctau0 - Ctau1) * (eff_charge_n - eff_charge_p)**2

        ! - - - - - - - - - - - - - - - - - - - - - - - - -
        ! l = 0, monopole
        if(l == 0) then
          
          ! integral over the mesh of 4 (x^2 + y^2 + z^2) * rho_n * rho_p
          kappa = kappa * 4.0 * sum( (meshgrid(:,1)**2 + meshgrid(:,2)**2 + meshgrid(:,3)**2) &
            &                       * R%D_I_I(:,1) * R%D_I_I(:,2)) * dv

        ! - - - - - - - - - - - - - - - - - - - - - - - - - 
        ! l = 1  dipole (same for m=0 and m=1)
        else if(l == 1) then

          ! integral over the mesh of 3/(4pi) rho_n * rho_p
          kappa = kappa * (3.0 / (4.0 * pi)) * sum(R%D_I_I(:,1) * R%D_I_I(:,2)) * dv
        
        ! - - - - - - - - - - - - - - - - - - - - - - - - - 
        ! l = 2, m = 0, axial quadrupole Q20
        else if(l == 2 .and. m==0) then

          ! integral over the mesh of 5/(4pi) (x^2 + y^2 + 4*z^2) * rho_n * rho_p
          kappa = kappa * (5.0 / (4.0 * pi) ) &
            &     * sum( (meshgrid(:,1)**2 + meshgrid(:,2)**2 + 4.*meshgrid(:,3)**2) &
            &            * R%D_I_I(:,1) * R%D_I_I(:,2)) * dv

        ! - - - - - - - - - - - - - - - - - - - - - - - - - 
        ! l = 2, m = 2, quadrupole Q22+ = 1/sqrt(2) (Q_22 + Q_2,-2)
        else if(l == 2 .and. m==2) then

          ! integral over the mesh of 15/(4pi) (x^2 + y^2) * rho_n * rho_p
          kappa = kappa * (15.0 / (4.0 * pi) ) * sum( (meshgrid(:,1)**2 + meshgrid(:,2)**2 ) &
            &                                          * R%D_I_I(:,1) * R%D_I_I(:,2)) * dv
          
        endif
          
        if (fam_verbose > 1) print *, "enhancement factor kappa = ", kappa / m1kin
      
    endif
#endif 

    ewsr = m1kin + kappa

    print *, "Energy-weighted sum rule : m1 = ", ewsr

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
$TR         occ_p = 2.0d0 - rho_can(p)   ! WR: Is this not superfluous? I mean, occ_h and occ_p do not actually enter the result? 
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
$TR         occ_p = 2.0d0 - rho_can(p) ! WR: Is this not superfluous? I mean, occ_h and occ_p do not actually enter the result? 
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

  subroutine transform_sp_to_qp(Bogo, OTRsp, OTLsp, OBLsp, OTRqp, OTLqp, OBLqp)
    !---------------------------------------------------------------------------
    ! Transform a matrix representation of an operator from the single-particle
    ! to the quasiparticle basis. This routine assumes that the Bogoliubov 
    ! transformation is real and that the operator commutes with time-reversal, 
    ! parity, and z-signature.
    ! 
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !
    ! Conventions: 
    !  1. the matrix representation of an operator that is bilinear in
    !     fermionic annihilation/creation operators is
    !
    !        O = 1/2 ( c^dagger c ) ( o^TL  o^TR   ) ( c         )
    !                               ( o^BL -o^TL,T ) ( c^\dagger )
    !
    !     where the {c, c^{\dagger}} can be particle or quasiparticle operators.
    !     The nomenclature of the variables is the following 
    !       TL = top left 
    !       TR = top right
    !       BL = bottom left
    !
    !     This convention is ALMOST the same as the standard operator notation
    !       O^TL =  O^11
    !       O^TR =  O^20 
    !       O^BL = -O^02   <---- sneaky minus sign 
    !
    !     If O is hermitian, then we have that
    !       ->  o^{TL} = + o^{TL, \dagger}
    !         ->  o^{BL} = - o^{TR,*}
    !       but this is not true in general.
    ! 
    ! 2. when time-reversal is explicitly conserved, the represented submatrices
    !     are halved; in that case, the matrices take the following form:
    ! 
    !       O^{TL} = ( O^{TL} 0      )   
    !                ( 0      O^{TL} )                       
    !     
    !       O^{TR} = ( 0      O^{TR} ) 
    !                (-O^{TR} 0      )
    ! 
    !       O^{BL} = ( 0      O^{BL} )
    !                (-O^{BL} 0      )
    !
    !     This is the labelling this routine adopts, i.e. the inputs concern 
    !     the top-left part of O^{TL} and the top-right parts of both O^{TR} 
    !     O^{BL}.
    !
    ! If these conventions do not appear natural to you, realize they have 
    ! been adopted to make the useage of this routines straightforward; i.e. 
    ! the user should not worry about (possibly symmetry-dependent) signs in 
    ! the inputs of this routine.
    ! 
    ! With these conventions, the formulas are 
    ! 
    !   O^{TL}_qp =     U^{\dagger} O^{TL} V^*  +     U^\dagger   O^{TR}   V 
    !             + \xi V^{\dagger} O^{BL} U    -     V^{\dagger} O^{TL,T} V   
    !
    !   O^{TR}_qp = \xi U^{\dagger} O^{TL} V^*  +     U^{\dagger} O^{TR}   U^*
    !             +     V^{\dagger} O^{BL} V    -     V^{\dagger} O^{TL,T} U^*  
    !    
    !   O^{BL}_qp =     V^T         O^{TL} U    +     V^T         O^{TR}   V 
    !             +     U^T         O^{BL} U    - \xi U^T         O^{TL}   V
    !  where \xi = +1(-1) is time-reversal is broken(converved).     
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input: 
    ! - Bogo: Bogoliubov transformation matrix
    ! - OTRsp, OTLsp, OBLsp: complex 2D arrays 
    !                        matrix elements in single-particle space 
    ! Output:
    ! - OTRqp, OTLqp, OBLqp: complex 2D arrays 
    !                        matrix elements in quasiparticle space
    !---------------------------------------------------------------------------     
    real(KIND=dp), intent(in)                       :: Bogo(:, :)
    complex(KIND=dp), intent(in), optional, target  :: OTRsp(:, :), OTLsp(:, :), OBLsp(:, :)
    complex(KIND=dp), intent(out), optional, target :: OTRqp(:, :), OTLqp(:, :), OBLqp(:, :)

    complex(KIND=dp), pointer     :: OTR_sp_p(:, :), OTL_sp_p(:, :), OBL_sp_p(:, :)
    complex(KIND=dp), pointer     :: OTR_qp_p(:, :), OTL_qp_p(:, :), OBL_qp_p(:, :)

    real(KIND=dp), allocatable    :: Ub(:, :), Vb(:, :)
    integer                       :: B, N, N2, si, sb, T, i
    real(KIND=dp)                 :: Tphase

    if (present(OTRqp)) OTRqp = 0._dp
    if (present(OTLqp)) OTLqp = 0._dp
    if (present(OBLqp)) OBLqp = 0._dp

$NTR Tphase = +1.0_dp
$TR  Tphase = -1.0_dp

    ! si determines the start of the block in sp-basis of dimension nwt
    ! sb determines the start of the block in qp-basis of dimension 2*nwt 
    !   -> Bogo contains all HFB eigenvectors ordered with increasing QPE (-Emax,..., -E1, E1,..., Emax)   
    si = 0 ; sb = 0
    do B=1,8,2
        N  = HFblocks(B)    ; if(N.eq.0) cycle 
        N2 = HFblocks(B+1)
        T = N + N2
  
        ! Getting the U and V out to make the formulas explicit
        ! and the matrix multiplications memory-local
        Ub = Bogo(sb  +1:sb+  T,sb+T+1:sb+2*T)
        Vb = Bogo(sb+T+1:sb+2*T,sb+T+1:sb+2*T)

        ! Pointers to make the equations below more compact
        if(present(OTRsp)) OTR_sp_p => OTRsp(si+1:si+T,si+1:si+T)
        if(present(OTLsp)) OTL_sp_p => OTLsp(si+1:si+T,si+1:si+T)
        if(present(OBLsp)) OBL_sp_p => OBLsp(si+1:si+T,si+1:si+T)

        if(present(OTRqp)) OTR_qp_p => OTRqp(si+1:si+T,si+1:si+T)
        if(present(OTLqp)) OTL_qp_p => OTLqp(si+1:si+T,si+1:si+T)
        if(present(OBLqp)) OBL_qp_p => OBLqp(si+1:si+T,si+1:si+T)

        if(present(OTLqp)) then
          if(present(OTLsp)) then 
              OTL_qp_p = OTL_qp_p +          matmul(transpose(Ub),matmul(          OTL_sp_p , Ub)) ! + U^\dagger O^11   U
              OTL_qp_p = OTL_qp_p -          matmul(transpose(Vb),matmul(transpose(OTL_sp_p), Vb)) ! - V^\dagger O^11,T V
          endif 
          if(present(OTRsp)) then 
              OTL_qp_p = OTL_qp_p +          matmul(transpose(Ub),matmul(          OTR_sp_p , Vb)) ! + U^\dagger O^20   V 
          endif 
          if(present(OBLsp)) then 
              OTL_qp_p = OTL_qp_p + Tphase * matmul(transpose(Vb),matmul(          OBL_sp_p , Ub)) ! + V^\dagger O^02   U
          endif 
        endif 

        if(present(OTRqp)) then 
          if(present(OTLsp)) then 
              OTR_qp_p = OTR_qp_p + Tphase * matmul(transpose(Ub),matmul(          OTL_sp_p , Vb)) ! + U^\dagger O^11   V^*
              OTR_qp_p = OTR_qp_p -          matmul(transpose(Vb),matmul(transpose(OTL_sp_p), Ub)) ! - V^\dagger O^11,T U^*
          endif 
          if(present(OTRsp)) then 
              OTR_qp_p = OTR_qp_p +          matmul(transpose(Ub),matmul(          OTR_sp_p , Ub)) ! + U^\dagger O^20   U^*
          endif 
          if(present(OBLsp)) then 
              OTR_qp_p = OTR_qp_p +          matmul(transpose(Vb),matmul(          OBL_sp_p , Vb)) ! + V^\dagger O^02   V^*
          endif
        endif 

        if(present(OBLqp)) then
          if(present(OTLsp)) then 
              OBL_qp_p = OBL_qp_p +          matmul(transpose(Vb),matmul(          OTL_sp_p , Ub)) ! + V^T       O^11   U
              OBL_qp_p = OBL_qp_p - Tphase * matmul(transpose(Ub),matmul(transpose(OTL_sp_p), Vb)) ! - U^T       O^11,T U
          endif
          if(present(OTRsp)) then 
              OBL_qp_p = OBL_qp_p +          matmul(transpose(Vb),matmul(          OTR_sp_p , Vb)) ! + V^T       O^20   V
          endif 
          if(present(OBLsp)) then 
              OBL_qp_p = OBL_qp_p +          matmul(transpose(Ub),matmul(          OBL_sp_p , Ub)) ! + U^T       O^02   U
          endif
        endif 

        si = si +  T
        sb = sb +2*T      
    enddo 

   end subroutine transform_sp_to_qp

   subroutine transform_qp_to_sp(Bogo, OTRqp, OTLqp, OBLqp, OTRsp, OTLsp, OBLsp)
    !---------------------------------------------------------------------------
    ! Transform a matrix representation of an operator from the quasiparticle
    ! to the single basis. This routine assumes that the Bogoliubov 
    ! transformation is real and that the operator commutes with time-reversal, 
    ! parity, and z-signature.
    ! 
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !
    ! Conventions: 
    !  1. the matrix representation of an operator that is bilinear in
    !     fermionic annihilation/creation operators is
    !
    !        O = 1/2 ( c^dagger c ) ( o^TL  o^TR   ) ( c         )
    !                               ( o^BL -o^TL,T ) ( c^\dagger )
    !
    !     where the {c, c^{\dagger}} can be particle or quasiparticle operators.
    !     The nomenclature of the variables is the following 
    !       TL = top left 
    !       TR = top right
    !       BL = bottom left
    !
    !     This convention is ALMOST the same as the standard operator notation
    !       O^TL =  O^11
    !       O^TR =  O^20 
    !       O^BL = -O^02   <---- sneaky minus sign 
    !
    !     If O is hermitian, then we have that
    !       ->  o^{TL} = + o^{TL, \dagger}
    !         ->  o^{BL} = - o^{TR,*}
    !       but this is not true in general.
    ! 
    ! 2. when time-reversal is explicitly conserved, the represented submatrices
    !     are halved; in that case, the matrices take the following form:
    ! 
    !       O^{TL} = ( O^{TL} 0      )   
    !                ( 0      O^{TL} )                       
    !     
    !       O^{TR} = ( 0      O^{TR} ) 
    !                (-O^{TR} 0      )
    ! 
    !       O^{BL} = ( 0      O^{BL} )
    !                (-O^{BL} 0      )
    !
    !     This is the labelling this routine adopts, i.e. the inputs concern 
    !     the top-left part of O^{TL} and the top-right parts of both O^{TR} 
    !     O^{BL}.
    !
    ! If these conventions do not appear natural to you, realize they have 
    ! been adopted to make the useage of this routines straightforward; i.e. 
    ! the user should not worry about (possibly symmetry-dependent) signs in 
    ! the inputs of this routine.
    ! 
    ! With these conventions, the formulas are 
    ! 
    !   O^{TL} =     U   O^{TL}_qp U^{\dagger} + \xi U   O^{TR}_qp   U^{\dagger}
    !          +     V^* O^{BL}_qp U^{\dagger} -     V^* O^{TL,T}_qp V 
    !
    !   O^{TR} =     U   O^{TL}_qp U^T         +     U   O^{TR}_qp   U^T 
    !          +     V^* O^{BL}_qp U^\dagger   - \xi V^* O^{TL,T}_qp U^\dagger 
    !
    !   O^{BL} = \xi V   O^{TL}_qp U^{\dagger} +     V   O^{TR}_qp   V^T 
    !          +     U^* O^{BL}_qp U^{\dagger} -     U^* O^{TL,T}_qp U^T
    ! 
    !  where \xi = +1(-1) is time-reversal is broken(converved).     ! 
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input: 
    ! - Bogo: Bogoliubov transformation matrix
    ! - OTRqp, OTLqp, OBLqp: complex 2D arrays 
    !                        matrix elements in quasiparticle space 
    ! Output:
    ! - OTRsp, OTLsp, OBLsp: complex 2D arrays 
    !                        matrix elements in single-particle space
    !---------------------------------------------------------------------------   
    real(KIND=dp), intent(in)                       :: Bogo(:, :)
    complex(KIND=dp), intent(in ), optional, target :: OTRqp(:, :), OTLqp(:, :), OBLqp(:, :)
    complex(KIND=dp), intent(out), optional, target :: OTRsp(:, :), OTLsp(:, :), OBLsp(:, :)

    complex(KIND=dp), pointer     :: OTR_sp_p(:, :), OTL_sp_p(:, :), OBL_sp_p(:, :)
    complex(KIND=dp), pointer     :: OTR_qp_p(:, :), OTL_qp_p(:, :), OBL_qp_p(:, :)

    real(KIND=dp), allocatable    :: Ub(:, :), Vb(:, :)
    integer                       :: B, N, N2, si, sb, T, i
    real(KIND=dp)                 :: Tphase

   ! initialise the single-particle matrix elements to zero if they are present
    if(present(OTRsp)) OTRsp = 0._dp
    if(present(OTLsp)) OTLsp = 0._dp
    if(present(OBLsp)) OBLsp = 0._dp

$NTR  Tphase = +1.0_dp
$TR   Tphase = -1.0_dp

    si = 0 ; sb = 0
    do B=1,8,2
      N  = HFblocks(B)    ; if(N.eq.0) cycle 
      N2 = HFblocks(B+1)
      T = N + N2
  
      ! Getting the U and V out to make the formulas explicit
      ! and the matrix multiplications memory-local
      Ub = Bogo(sb  +1:sb+  T,sb+T+1:sb+2*T)
      Vb = Bogo(sb+T+1:sb+2*T,sb+T+1:sb+2*T)

      ! Pointers to make the equations below more compact
      if(present(OTRsp)) OTR_sp_p => OTRsp(si+1:si+T,si+1:si+T)
      if(present(OTLsp)) OTL_sp_p => OTLsp(si+1:si+T,si+1:si+T)
      if(present(OBLsp)) OBL_sp_p => OBLsp(si+1:si+T,si+1:si+T)

      if(present(OTRqp)) OTR_qp_p => OTRqp(si+1:si+T,si+1:si+T)
      if(present(OTLqp)) OTL_qp_p => OTLqp(si+1:si+T,si+1:si+T)
      if(present(OBLqp)) OBL_qp_p => OBLqp(si+1:si+T,si+1:si+T)

      if(present(OTLsp)) then 
        if(present(OTLqp)) then 
          OTL_sp_p = OTL_sp_p +          matmul( Ub, matmul(           OTL_qp_p , transpose(Ub))) ! + U   O^{11}   U^{\dagger}
          OTL_sp_p = OTL_sp_p -          matmul( Vb, matmul( transpose(OTL_qp_p), transpose(Vb))) ! - V^* O^{11},T V^{\dagger}
        endif 
        if(present(OTRqp)) then 
          OTL_sp_p = OTL_sp_p + Tphase * matmul( Ub, matmul(           OTR_qp_p , transpose(Vb))) ! + U   O^{20}   V^{\dagger}
        endif 
        if(present(OBLqp)) then 
          OTL_sp_p = OTL_sp_p +          matmul( Vb, matmul(           OBL_qp_p , transpose(Ub))) ! + V^* O^{02}   U^{\dagger}
        endif 
      endif 

      if(present(OTRsp)) then 
        if(present(OTLqp)) then 
          OTR_sp_p = OTR_sp_p +          matmul( Ub, matmul(           OTL_qp_p , transpose(Vb))) ! + U   O^{11}   V^{\dagger}
          OTR_sp_p = OTR_sp_p - Tphase * matmul( Vb, matmul( transpose(OTL_qp_p), transpose(Ub))) ! - V^* O^{11},T U^{\dagger} 
        endif 
        if(present(OTRqp)) then 
          OTR_sp_p = OTR_sp_p +          matmul( Ub, matmul(           OTR_qp_p , transpose(Ub))) ! + U   O^{20}   U^T
        endif 
        if(present(OBLqp)) then 
          OTR_sp_p = OTR_sp_p +          matmul( Vb, matmul(           OBL_qp_p , transpose(Vb))) ! + V^* O^{20}   V^{\dagger}
        endif 
      endif 

      if(present(OBLsp)) then 
        if(present(OTLqp)) then 
          OBL_sp_p = OBL_sp_p + Tphase * matmul( Vb, matmul(           OTL_qp_p , transpose(Ub))) ! + V   O^{11}   U^{\dagger} 
          OBL_sp_p = OBL_sp_p -          matmul( Ub, matmul( transpose(OTL_qp_p), transpose(Vb))) ! - U^* O^{11},T V^T
        endif 
        if(present(OTRqp)) then 
          OBL_sp_p = OBL_sp_p +          matmul( Vb, matmul(           OTR_qp_p , transpose(Vb))) ! + V   O^{20}   V^T
        endif
        if(present(OBLqp)) then 
          OBL_sp_p = OBL_sp_p +          matmul( Ub, matmul(           OBL_qp_p , transpose(Ub))) ! + U^* O^{20}   U^{\dagger}
        endif
      endif

      si = si +  T
      sb = sb +2*T
    enddo

  end subroutine transform_qp_to_sp

  subroutine transform_sp_to_qp_pd(Bogo, O20sp, O11sp, O02sp, O20qp, O11qp, O02qp)
    !---------------------------------------------------------------------------
    ! Performing quasi-particle transformation of a 1-body operator that 
    ! does not have to be 
    !  (i)   purely particle-hole or particle-particle / hole-hole 
    !  (ii)  hermitian 
    ! 
    ! Its generic form is schematically 
    !               O = o20sp + o11sp + o02sp. 
    ! 
    ! The function returns the matrix elements in of O in the operator in the 
    ! quasiparticle basis associated with the Bogoliubov transformation Bogo.
    !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Important notes:
    !  - inputs not specified are assumed to be zero. 
    !  - outputs not specified are not calculated, but do not necessarily vanish.
    !  - this routine works for REAL-valued Bogoliubov transformations...
    !  - ... but complex-valued operator matrix elements
    !  - this routine assumes that O is compatible with the symmetries of run!
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !    Bogo             : Bogoliubov transformation matrix W from sp to qp basis 
    !                       (2*nwt,2*nwt)
    !    O20sp (optional) : sp matrix elements of 20 operator component (nwt,nwt)
    !    O11sp (optional) : sp matrix elements of 11 operator component (nwt,nwt)
    !    O02sp (optional) : sp matrix elements of 02 operator component (nwt,nwt)
    !  
    !    Important convention:  TODO 
    ! 
    !
    !
    ! Output:
    !    O20qp (optional) : qp matrix elements of 20 operator component (nwt,nwt)
    !    O11qp (optional) : qp matrix elements of 11 operator component (nwt,nwt)
    !    O02qp (optional) : qp matrix elements of 02 operator component (nwt,nwt)
    !
    !    Important reminder:  TODO
    !   
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    !  
    ! Bogo contains the bogoliubov transformation W organised in block matrices
    ! where blocks have twice the size of HFblocks, i.e.
    ! 
    !              (  Wb         )                        (  Vb^*   Ub   )
    !    Bogo  =   (     Wb    : )                Wb  =   (              )
    !              (        ..Wb )                        (  Ub^*   Vb   )
    ! 
    ! Note that the block structure wrt Rz is non-trivial as it is antihermitian
    ! Hence matrices U and V have block structure in Rz
    ! 
    !              (  Ub(++)   0  )                         (   0    Vb(+-) )
    !       Ub  =  (              )                Vb   =   (               )
    !              (   0   Ub(--) )                         (  Vb(-+)   0   ) 
    !   
    ! The expressions for the QP matrix elements are 
    ! 
    !       O20qp = + T * U^{dagger} o11sp   V^* + U^{dagger} o20sp   U^* 
    !               - T * V^{dagger} o02sp^* V^* - V^{dagger} o11sp^T U^*
    ! 
    !       O11qp = + U^{dagger} o11sp   U   +     U^{dagger} o20sp   V
    !               - V^{dagger} o02sp^* U   -     V^{dagger} o11sp^T V
    ! 
    !       O02qp = - T * V^T    o11sp   U   - T * V^T        o20sp   V 
    !               +     U^T    o02sp^* U   +     U^T        o11sp^T V
    ! 
    ! If time-reversal is not leveraged in the code, then T = 1 and these 
    !  expressions can be interpreted straightforwardly in terms of the full 
    !  matrices.
    ! 
    ! If time-reversal is not leveraged, then T = -1 and these expressions 
    !  tackle the objects stored in the code - i.e. NOT the full matrices.
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in)                       :: Bogo(:,:)
    complex(KIND=dp), intent(in) , optional, target :: O20sp(:,:), O11sp(:,:), O02sp(:,:)
    complex(KIND=dp), intent(out), optional, target :: O20qp(:,:), O11qp(:,:), O02qp(:,:)

    real(KIND=dp), allocatable    :: Ub(:,:), Vb(:,:)
    complex(KIND=dp), pointer     :: O20b(:,:), O11b(:,:), O02b(:,:)
    complex(KIND=dp), pointer     :: O20b_qp(:,:), O11b_qp(:,:), O02b_qp(:,:)
    integer                       :: B, N, N2, si, sb, T, i
    real(KIND=dp)                 :: Tphase

    if(present(O20qp)) O20qp = 0._dp
    if(present(O11qp)) O11qp = 0._dp
    if(present(O02qp)) O02qp = 0._dp

$NTR Tphase = +1.0_dp
$TR  Tphase = -1.0_dp

    if (fam_verbose > 2) print *, "transform_sp_to_qp"

    ! si determines the start of the block in sp-basis of dimension nwt
    ! sb determines the start of the block in qp-basis of dimension 2*nwt 
    !   -> Bogo contains all HFB eigenvectors ordered with increasing QPE (-Emax,..., -E1, E1,..., Emax)   
    si = 0 ; sb = 0
    do B=1,8,2
      N  = HFblocks(B)    ; if(N.eq.0) cycle 
      N2 = HFblocks(B+1)
      T = N + N2
  
      ! Getting the U and V out to make the formulas explicit
      ! and the matrix multiplications memory-local
      Ub = Bogo(sb  +1:sb+  T,sb+T+1:sb+2*T)
      Vb = Bogo(sb+T+1:sb+2*T,sb+T+1:sb+2*T)

     !  print *, 'U', B, B+1
     !  do i=1,T
     !     print ('(99f10.3)'), Ub(i,1:T)
     !  enddo
     !  print *
     ! print *, 'V', B, B+1
     !  do i=1,T
     !     print ('(99f10.3)'), Vb(i,1:T)
     !  enddo
     !  print *

      ! Pointers to make the equations below more compact
      if(present(O20sp)) O20b => O20sp(si+1:si+T,si+1:si+T)
      if(present(O11sp)) O11b => O11sp(si+1:si+T,si+1:si+T)
      if(present(O02sp)) O02b => O02sp(si+1:si+T,si+1:si+T)

      if(present(O20qp)) O20b_qp => O20qp(si+1:si+T,si+1:si+T)
      if(present(O11qp)) O11b_qp => O11qp(si+1:si+T,si+1:si+T)
      if(present(O02qp)) O02b_qp => O02qp(si+1:si+T,si+1:si+T)

      if (fam_verbose > 2) print '(A, I3, I3, A, I5)', 'Blocks: ', B, B+1, ' with size', T
      if (fam_verbose > 2) print '(A, F10.2)',  '||U||^2 = ', sum(Ub(:,:) * Ub(:,:))
      if (fam_verbose > 2) print '(A, F10.2)',  '||V||^2 = ', sum(Vb(:,:) * Vb(:,:))

      !  O20qp = + T * U^{dagger} o11sp   V^* + U^{dagger} o20sp   U^* 
      !          - T * V^{dagger} o02sp^* V^* - V^{dagger} o11sp^T U^*
      if(present(O20qp)) then
        if(present(O11sp)) then
          O20b_qp = O20b_qp + Tphase * matmul(transpose(Ub),  matmul(          O11b , Vb))  ! + T * U^{dagger} o11sp   V^*
          O20b_qp = O20b_qp -          matmul(transpose(Vb),  matmul(transpose(O11b), Ub))  ! -     V^{dagger} o11sp^T U^*
        endif
        if(present(O20sp)) then
          O20b_qp = O20b_qp +          matmul(transpose(Ub),  matmul(          O20b , Ub)) ! +     U^{dagger} o20sp   U^*
        endif
        if(present(O02sp)) then
          O20b_qp = O20b_qp - Tphase * matmul(transpose(Vb),  matmul(          O02b , Vb)) ! - T * V^{dagger} o02sp^* V^*
        endif
      endif 

      !       O11qp = + U^{dagger} o11sp   U   + U^{dagger} o20sp   V
      !               - V^{dagger} o02sp^* U   - V^{dagger} o11sp^T V
      if(present(O11qp)) then
        if(present(O11sp)) then
          O11b_qp = O11b_qp + matmul(transpose(Ub),  matmul(          O11b , Ub)) ! + U^{dagger} o11sp   U
          O11b_qp = O11b_qp - matmul(transpose(Vb),  matmul(transpose(O11b), Vb)) ! - V^{dagger} o11sp^T V
        endif
        if(present(O20sp)) then
          O11b_qp = O11b_qp + matmul(transpose(Ub),  matmul(          O20b , Vb)) ! + U^{dagger} o20sp   V
        endif
        if(present(O02sp)) then
          O11b_qp = O11b_qp - matmul(transpose(Vb),  matmul(          O02b , Ub)) ! - V^{dagger} o02sp   U
        endif
      endif

      !       O02qp = - T * V^T        o11sp   U   - T * V^T        o20sp   V 
      !               + U^T        o02sp^* U   + U^T        o11sp^T V
      if(present(O02qp)) then
        if(present(O11sp)) then
          O02b_qp = O02b_qp - Tphase * matmul(transpose(Vb),  matmul(          O11b , Ub)) ! - T * V^T        o11sp   U
          O02b_qp = O02b_qp +          matmul(transpose(Ub),  matmul(transpose(O11b), Vb)) ! +     U^T        o11sp^T V
        endif
        if(present(O20sp)) then
          O02b_qp = O02b_qp - Tphase * matmul(transpose(Vb),  matmul(          O20b , Vb)) ! - T * V^T        o20sp   V 
        endif
        if(present(O02sp)) then
          O02b_qp = O02b_qp +          matmul(transpose(Ub),  matmul(          O02b , Ub)) ! +     U^T        o02sp^* U
        endif
      endif 

      si = si +  T
      sb = sb +2*T
    enddo

    if (fam_verbose > 2) then
      print *, 'Symmetry : '
      if(present(O20sp)) then
         print *, '    O20sp = + O20sp^T   : satisfied up to',  &
              & sum(abs(O20sp(:,:) - transpose(O20sp(:,:))))
      endif
      if(present(O02sp)) then
         print *, '    O02sp = + O02sp^T   : satisfied up to',  &
              & sum(abs(O02sp(:,:) - transpose(O02sp(:,:))))
      endif
      if(present(O20qp)) then
         print *, '    O20qp = + O20qp^T   : satisfied up to',  &
              & sum(abs(O20qp(:,:) - transpose(O20qp(:,:))))
      endif
      if(present(O02qp)) then
         print *, '    O02qp = + O02qp^T   : satisfied up to',  &
              & sum(abs(O02qp(:,:) - transpose(O02qp(:,:))))
      endif
      print *, 'Antisymmetry : '
      if(present(O20sp)) then
         print *, '    O20sp = - O20sp^T   : satisfied up to',  &
              & sum(abs(O20sp(:,:) + transpose(O20sp(:,:))))
      endif
      if(present(O02sp)) then
         print *, '    O02sp = - O02sp^T   : satisfied up to',  &
              & sum(abs(O02sp(:,:) + transpose(O02sp(:,:))))
      endif
      if(present(O20qp)) then
         print *, '    O20qp = - O20qp^T   : satisfied up to',  &
              & sum(abs(O20qp(:,:) + transpose(O20qp(:,:))))
      endif
      if(present(O02qp)) then
         print *, '    O02qp = - O02qp^T   : satisfied up to',  &
              & sum(abs(O02qp(:,:) + transpose(O02qp(:,:))))
      endif
      print *, 'Hermiticity : '
      if(present(O20sp) .and. present(O02sp)) then
         print *, '    O20sp = O02sp*   : satisfied up to', &
              & sum(abs(O20sp(:,:) - conjg(O02sp(:,:))))
      endif
      if(present(O11sp)) then
         print *, '    O11sp = O11sp^T^*   : satisfied up to', &
              & sum(abs(O11sp(:,:) - conjg(transpose(O11sp(:,:)))))
      endif
      if(present(O20qp) .and. present(O02qp)) then
         print *, '    O20qp = O02qp*   : satisfied up to', &
              & sum(abs(O20qp(:,:) - conjg(O02qp(:,:))))
      endif
      if(present(O11qp)) then
         print *, '    O11qp = O11qp^T^*   : satisfied up to', &
              & sum(abs(O11qp(:,:) - conjg(transpose(O11qp(:,:)))))
      endif
      print *, 'Anti-hermiticity : '
      if(present(O20sp) .and. present(O02sp)) then
         print *, '    O20sp = - O02sp*   : satisfied up to', &
              & sum(abs(O20sp(:,:) + conjg(O02sp(:,:))))
      endif
      if(present(O11sp)) then
         print *, '    O11sp = - O11sp^T^*   : satisfied up to', &
              & sum(abs(O11sp(:,:) + conjg(transpose(O11sp(:,:)))))
      endif
      if(present(O20qp) .and. present(O02qp)) then
         print *, '    O20qp = - O02qp*   : satisfied up to',  &
              & sum(abs(O20qp(:,:) + conjg(O02qp(:,:))))
      endif
      if(present(O11qp)) then
         print *, '    O11qp = - O11qp^T^*   : satisfied up to', &
              & sum(abs(O11qp(:,:) + conjg(transpose(O11qp(:,:)))))
      endif
    endif

  end subroutine transform_sp_to_qp_pd

  subroutine transform_qp_to_sp_pd(Bogo, O20qp, O11qp, O02qp, O20sp, O11sp, O02sp)
    !---------------------------------------------------------------------------
    ! Performing quasi-particle back transformation of a generic on 1-body operator 
    ! O = O20qp + O11qp + O02qp. The function returns the matrix elements in of 
    ! O in the operator in the sp basis.

    ! Input:
    !    Bogo             : Bogoliubov transformation matrix W from sp to qp basis 
    !                       (2*nwt,2*nwt)
    !    O20qp (optional) : qp matrix elements of 20 operator component (nwt,nwt)
    !    O11qp (optional) : qp matrix elements of 11 operator component (nwt,nwt)
    !    O02qp (optional) : qp matrix elements of 02 operator component (nwt,nwt)
    ! Output:
    !    O20sp (optional) : sp matrix elements of 20 operator component (nwt,nwt)
    !    O11sp (optional) : sp matrix elements of 11 operator component (nwt,nwt)
    !    O02sp (optional) : sp matrix elements of 02 operator component (nwt,nwt)
    ! 
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    !  
    ! Bogo contains the bogoliubov transformation W organised in block matrices
    ! where blocks have twice the size of HFblocks, i.e.
    ! 
    !              (  Wb         )                        (  Vb^*   Ub   )
    !    Bogo  =   (     Wb    : )                Wb  =   (              )
    !              (        ..Wb )                        (  Ub^*   Vb   )
    ! 
    ! Note that the block structure wrt Rz is non-trivial as it is antihermitian
    ! Hence matrices U and V have block structure in Rz
    ! 
    !              (  Ub(++)   0  )                         (   0    Vb(+-) )
    !       Ub  =  (              )                Vb   =   (               )
    !              (   0   Ub(--) )                         (  Vb(-+)   0   ) 
    ! 
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! REMARK : shouldnt this be equivalent to calling the subroutine 
    !   transform_sp_to_qp() for Bogo^dagger. 
    !---------------------------------------------------------------------------


    implicit none
    real(KIND=dp), intent(in)               :: Bogo(:,:)
    complex(KIND=dp), intent(in), optional  :: O20qp(:,:), O11qp(:,:), O02qp(:,:)
    complex(KIND=dp), intent(out), optional :: O20sp(:,:), O11sp(:,:), O02sp(:,:)

    real(KIND=dp), allocatable    :: Ub(:,:), Vb(:,:), rho(:,:), kappa(:,:)
    complex(KIND=dp), allocatable :: O20b(:,:), O11b(:,:), O02b(:,:)
    integer                       :: B, N, N2, si, sb, T, i
    real(KIND=dp)                 :: Tphase


    ! initialise Oijsp outputs to zero if they are present
    if(present(O20sp)) O20sp = 0._dp
    if(present(O11sp)) O11sp = 0._dp
    if(present(O02sp)) O02sp = 0._dp

$NTR Tphase = 1.0_dp
$TR  Tphase = -1.0_dp

    
    if (fam_verbose > 1) print *, "transform_qp_to_sp"


    ! si = O start index for O , sb = start index for bogo (increases twice as fast)

    si = 0 ; sb = 0
    do B=1,8,2
      N  = HFblocks(B)    ; if(N.eq.0) cycle 
      N2 = HFblocks(B+1)
      T = N + N2
  
      ! Getting the U and V out to make the formulas explicit
      ! and the matrix multiplications memory-local
      Ub = Bogo(sb  +1:sb+  T,sb+T+1:sb+2*T)
      Vb = Bogo(sb+T+1:sb+2*T,sb+T+1:sb+2*T)

      ! get the correct subblock if the present qp operator components
      if(present(O20qp)) O20b = O20qp(si+1:si+T,si+1:si+T)
      if(present(O11qp)) O11b = O11qp(si+1:si+T,si+1:si+T)
      if(present(O02qp)) O02b = O02qp(si+1:si+T,si+1:si+T)
  
      if (fam_verbose > 2) print '(A, I3, I3, A, I5)', 'Blocks: ', B, B+1, ' with size', T
      if (fam_verbose > 2) print '(A, F10.2)',  '||U||^2 = ', sum(Ub(:,:) * Ub(:,:))
      if (fam_verbose > 2) print '(A, F10.2)',  '||V||^2 = ', sum(Vb(:,:) * Vb(:,:))

      if(present(O20sp)) then
        if(present(O11qp)) then
          O20sp(si+1:si+T, si+1:si+T) = O20sp(si+1:si+T, si+1:si+T) + matmul(Ub,  matmul(          O11b , transpose(Vb)))
          O20sp(si+1:si+T, si+1:si+T) = O20sp(si+1:si+T, si+1:si+T) - Tphase * matmul(Vb,  matmul(transpose(O11b), transpose(Ub)))
        endif
        if(present(O20qp)) then
          O20sp(si+1:si+T, si+1:si+T) = O20sp(si+1:si+T, si+1:si+T) + matmul(Ub,  matmul(          O20b , transpose(Ub)))
        endif
        if(present(O02qp)) then
          O20sp(si+1:si+T, si+1:si+T) = O20sp(si+1:si+T, si+1:si+T) - Tphase * matmul(Vb,  matmul(          O02b , transpose(Vb)))
        endif
      endif 


      if(present(O11sp)) then
        if(present(O11qp)) then
          O11sp(si+1:si+T, si+1:si+T) = O11sp(si+1:si+T, si+1:si+T) + matmul(Ub,  matmul(          O11b , transpose(Ub)))
          O11sp(si+1:si+T, si+1:si+T) = O11sp(si+1:si+T, si+1:si+T) - matmul(Vb,  matmul(transpose(O11b), transpose(Vb)))
        endif
        if(present(O20qp)) then
          O11sp(si+1:si+T, si+1:si+T) = O11sp(si+1:si+T, si+1:si+T) + Tphase * matmul(Ub,  matmul(          O20b , transpose(Vb)))
        endif
        if(present(O02qp)) then
          O11sp(si+1:si+T, si+1:si+T) = O11sp(si+1:si+T, si+1:si+T) - Tphase * matmul(Vb,  matmul(          O02b , transpose(Ub)))
        endif
      endif 


      if(present(O02sp)) then
        if(present(O11qp)) then
          O02sp(si+1:si+T, si+1:si+T) = O02sp(si+1:si+T, si+1:si+T) - matmul(Vb,  matmul(          O11b , transpose(Ub)))
          O02sp(si+1:si+T, si+1:si+T) = O02sp(si+1:si+T, si+1:si+T) + Tphase * matmul(Ub,  matmul(transpose(O11b), transpose(Vb)))
        endif
        if(present(O20qp)) then
          O02sp(si+1:si+T, si+1:si+T) = O02sp(si+1:si+T, si+1:si+T) - Tphase * matmul(Vb,  matmul(          O20b , transpose(Vb)))
        endif
        if(present(O02qp)) then
          O02sp(si+1:si+T, si+1:si+T) = O02sp(si+1:si+T, si+1:si+T) + matmul(Ub,  matmul(          O02b , transpose(Ub)))
        endif
      endif

      si = si +  T
      sb = sb +2*T

    enddo
    

    if (fam_verbose > 2) then
      print *, 'Symmetry : '
      if(present(O20sp)) then
         print *, '    O20sp = + O20sp^T   : satisfied up to',&
              & sum(abs(O20sp(:,:) - transpose(O20sp(:,:))))
      endif
      if(present(O02sp)) then
         print *, '    O02sp = + O02sp^T   : satisfied up to',&
              &  sum(abs(O02sp(:,:) - transpose(O02sp(:,:))))
      endif
      if(present(O20qp)) then
         print *, '    O20qp = + O20qp^T   : satisfied up to',&
              &  sum(abs(O20qp(:,:) - transpose(O20qp(:,:))))
      endif
      if(present(O02qp)) then
         print *, '    O02qp = + O02qp^T   : satisfied up to',&
              &  sum(abs(O02qp(:,:) - transpose(O02qp(:,:))))
      endif
      print *, 'Antisymmetry : '
      if(present(O20sp)) then
         print *, '    O20sp = - O20sp^T   : satisfied up to',&
              &  sum(abs(O20sp(:,:) + transpose(O20sp(:,:))))
      endif
      if(present(O02sp)) then
         print *, '    O02sp = - O02sp^T   : satisfied up to',&
              &  sum(abs(O02sp(:,:) + transpose(O02sp(:,:))))
      endif
      if(present(O20qp)) then
         print *, '    O20qp = - O20qp^T   : satisfied up to',&
              &  sum(abs(O20qp(:,:) + transpose(O20qp(:,:))))
      endif
      if(present(O02qp)) then
         print *, '    O02qp = - O02qp^T   : satisfied up to',&
              &  sum(abs(O02qp(:,:) + transpose(O02qp(:,:))))
      endif
      print *, 'Hermiticity : '
      if(present(O20sp) .and. present(O02sp)) then
         print *, '    O20sp = O02sp*   : satisfied up to',&
              &  sum(abs(O20sp(:,:) - conjg(O02sp(:,:))))
      endif
      if(present(O11sp)) then
         print *, '    O11sp = O11sp^T^*   : satisfied up to',&
              &  sum(abs(O11sp(:,:) - conjg(transpose(O11sp(:,:)))))
      endif
      if(present(O20qp) .and. present(O02qp)) then
         print *, '    O20qp = O02qp*   : satisfied up to',&
              &  sum(abs(O20qp(:,:) - conjg(O02qp(:,:))))
      endif
      if(present(O11qp)) then
         print *, '    O11qp = O11qp^T^*   : satisfied up to',&
              &  sum(abs(O11qp(:,:) - conjg(transpose(O11qp(:,:)))))
      endif
      print *, 'Anti-hermiticity : '
      if(present(O20sp) .and. present(O02sp)) then
         print *, '    O20sp = - O02sp*   : satisfied up to',&
              &  sum(abs(O20sp(:,:) + conjg(O02sp(:,:))))
      endif
      if(present(O11sp)) then
         print *, '    O11sp = - O11sp^T^*   : satisfied up to',&
              &  sum(abs(O11sp(:,:) + conjg(transpose(O11sp(:,:)))))
      endif
      if(present(O20qp) .and. present(O02qp)) then
         print *, '    O20qp = - O02qp*   : satisfied up to',&
              &  sum(abs(O20qp(:,:) + conjg(O02qp(:,:))))
      endif
      if(present(O11qp)) then
         print *, '    O11qp = - O11qp^T^*   : satisfied up to',&
              &  sum(abs(O11qp(:,:) + conjg(transpose(O11qp(:,:)))))
      endif
    endif

  end subroutine transform_qp_to_sp_pd


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


  function read_f(Finfile) result(f_qpme)
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Read the quasi-particle matrix elements of the external field F20 and F02 
    ! from file named Finfile. 
    !
    ! INPUT:
    !     Finfile  : filename containing qpme of the external field F20 and F02
    !
    ! OUPUT: 
    !     f_qpme(:,:,:)  : complex qpme of F20 and F02 
    !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! note :
    !  -  The current version of the xy file written in IO.f90 contains the 
    !     perturbing filed F at the top of the file, followed by X, Y at several
    !     frequencies omega FAm was solved for. These different sections in the
    !     file are seperated by a line
    !     & omega =     [OMEGA]     [SMEAR]
    !  -  Also note that all the matrix elements larger than 1e-10 are stored in
    !     the file, antisymmetry is not exploited. I therefor explicitly
    !     anti-symmetrise the read in F to avoid noise wrt this symmetry. 
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


    implicit none
    character(len=*), intent(in) :: Finfile
    complex(KIND=dp), allocatable :: f_qpme(:,:,:)
    real(dp) :: F20_re, F20_im, F02_re, F02_im
    integer :: io, m, n, i
    character(len=256) :: line
    integer, parameter :: header_length = 14

    print *, "Reading F20 and F02 from ", Finfile
      
    allocate(f_qpme(nwt,nwt,2)) 
    f_qpme = 0


    ! Open the file
    open(unit=1, file=Finfile, status='old', iostat=io, action="read")
    if (io .ne. 0) then
      print *, 'Error opening file: ', Finfile
      return
    endif

    ! Read and print header lines 
    do i = 1, header_length
      read(1, '(A)') line
      if (fam_verbose > 1) print *, trim(line)
    enddo
    ! to do : could be a good idea to verify pqrqmeters in the header are compatible
    !         with the one read from stdin 

    do
      read(1, '(A)', iostat=io) line
      if (io < 0) exit  ! End of file
      if (line(1:1)=='&') exit ! End of F and beginning of XY section 

      ! Read sparse matrix data
      read(line, *) m, n, F20_re, F20_im, F02_re, F02_im
      ! print *, 'm = ', m, ' n = ', n, ' X = (', F20_re, ', ', F20_im, ') Y = (', F02_re, ', ', F02_im, ')'
      f_qpme(m,n,1) = dcmplx(F20_re, F20_im)
      f_qpme(m,n,2) = dcmplx(F02_re, F02_im)
    enddo

    close(1)

    if (pairingtype.ne.0) then ! QFAM

      ! explicitly antisymmetrise
      f_qpme(:,:,1) = 0.5 * (f_qpme(:,:,1) + transpose(f_qpme(:,:,1)))
      f_qpme(:,:,2) = 0.5 * (f_qpme(:,:,2) + transpose(f_qpme(:,:,2)))

    endif

    if(fam_verbose > 1) then

      print *, '||F20||² = ', sum(abs(f_qpme(:,:,1))**2)
      print *, '||F02||² = ', sum(abs(f_qpme(:,:,2))**2)

    endif

  end function read_f

  subroutine read_xy(XYinfile, X, Y)
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Read the quasi-particle matrix elements of the external field F20 and F02 
    ! from file named Finfile. 
    !
    ! INPUT:
    !     XYinfile  : filename containing qpme of X and Y
    !
    ! OUPUT: 
    !     X(:,:)    : complex qpme of X
    !     Y(:,:)    : complex qpme of Y
    !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! note :
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


    implicit none
    character(len=*), intent(in) :: XYinfile
    complex(KIND=dp), allocatable :: X(:,:),  Y(:,:)
    real(dp) :: X_re, X_im, Y_re, Y_im
    real(dp) :: omega_in, smear_in
    integer :: io, m, n, i
    character(len=256) :: line
    integer, parameter :: header_length = 14
    logical :: found_block = .false.

    print *, "Reading X and Y from ", XYinfile
      

    ! Open the file
    open(unit=1, file=XYinfile, status='old', iostat=io, action="read")
    if (io .ne. 0) then
      print *, 'Error opening file: ', XYinfile
      close(1)
      return
    endif

    ! Read and print header lines 
    do i = 1, header_length
      read(1, '(A)') line
      if (fam_verbose > 1) print *, trim(line)
    enddo
    ! to do : could be a good idea to verify pqrqmeters in the header are compatible
    !         with the one read from stdin 


    ! look for the block in the XY file starting with the separator line "& omega = [omega] [smear]"
    do
      read(1, '(A)', iostat=io) line
      
      ! reached end of file
      if (io < 0) then 
        print '(A, f8.3, A, f8.3)', 'End of file reached without finding good XY block for omega = ', omega_fam, ', smear = ', smear
        print *, 'Starting FAM solver with XY obtained from free response'
        close(1)
        return
      endif
      
      ! If the current line IS NOT a separator line (starting with '&'), then continue reading 
      if (line(1:1).ne.'&') then
       cycle 

      ! The current line IS a separator line
      else

        read(line, '(9X, F10.3, F10.3)') omega_in, smear_in

        ! Check if the separator line announces the XY block with correct frequency and smearing
        if ( (abs(omega_in - omega_fam) < 1e-6) .and. (abs(smear_in - smear) < 1e-6)) then 
          ! found the good block
          found_block = .true.
          exit
        else
          ! wrong block, continue reading
          cycle
        endif
      endif
    enddo


    ! if we found the correct block in the XYinfile, then fill the X and Y
    if (found_block) then
      X = 0
      Y = 0
      do
        read(1, '(A)', iostat=io) line
        if (io < 0) exit  ! End of file
        if (line(1:1)=='&') exit ! End of XY block 

        ! Read sparse matrix data
        read(line, *) m, n, X_re, X_im, Y_re, Y_im
        ! print *, 'm = ', m, ' n = ', n, ' X = (', X_re, ', ', X_im, ') Y = (', Y_re, ', ', Y_im, ')'
        X(m,n) = dcmplx(X_re, X_im)
        Y(m,n) = dcmplx(Y_re, Y_im)
      enddo

      close(1)
    endif

    if (pairingtype.ne.0) then ! QFAM

      ! explicitly antisymmetrise
      X = 0.5 * (X + transpose(X))
      Y = 0.5 * (Y + transpose(Y))

    
    endif

    if(fam_verbose > 1) then

      print *, '||X||² = ', sum(abs(X(:,:))**2)
      print *, '||Y||² = ', sum(abs(Y(:,:))**2)

    endif

  end subroutine read_xy

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
    print *, '||X||² = ', sum(abs(X)**2)
    if (pairingtype==0) then
      call print_spme_complex(X)
    else 
      call print_spme_complex_superblock(X)
    endif

    print *, 'Y'
    print *, '||Y||² = ', sum(abs(Y)**2)
    if (pairingtype==0) then
      call print_spme_complex(Y)
    else 
      call print_spme_complex_superblock(Y)
    endif

    print *, 'drho'
    print *, '||drho||² = ', sum(abs(drho)**2)
    call print_spme_complex(drho)

    if (pairingtype.ne.0) then
      print *, 'dkappa_plus'
      print *, '||dkappa_plus||² = ', sum(abs(dkappa_plus)**2)
      call print_spme_complex_superblock(dkappa_plus)

      print *, 'dkappa_minus'
      print *, '||dkappa_minus||² = ', sum(abs(dkappa_minus)**2)
      call print_spme_complex_superblock(dkappa_minus)
    endif

    print *, 'dh'
    print *, '||h||² = ', sum(abs(drho)**2)
    call print_spme_complex(drho)

    if (pairingtype.ne.0) then
      print *, 'dkappa_plus'
      print *, '||dkappa_plus||² = ', sum(abs(dkappa_plus)**2)
      call print_spme_complex_superblock(dkappa_plus)

      print *, 'dkappa_minus'
      print *, '||dkappa_minus||² = ', sum(abs(dkappa_minus)**2)
      call print_spme_complex_superblock(dkappa_minus)
    endif

    print *, 'dH20'
    print *, '||dH20||² = ', sum(abs(dH(:,:,1))**2)
    if (pairingtype==0) then
      call print_spme_complex(dH(:,:,1))
    else 
      call print_spme_complex_superblock(dH(:,:,1))
    endif

    print *, 'dH02'
    print *, '||dH02||² = ', sum(abs(dH(:,:,2))**2)
    if (pairingtype==0) then
      call print_spme_complex(dH(:,:,2))
    else 
      call print_spme_complex_superblock(dH(:,:,2))
    endif
  end subroutine print_all_fam_spmat
end module fam
