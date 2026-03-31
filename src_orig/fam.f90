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
  !   X and Y are stored in the quasi-particle basis, containing X20 and Y02 
  !   matrix elements. Only the positive qp spectrum is included such that size
  !   is still (nwt,nwt)
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
  ! - In absence of pairing : (FAM)
  !   dH(:,:,1) and dH(:,:,2) contain ph and hp elements of the perturbed hamiltonian in the HF basis
  !   They are allocated ove the complete basis size (nwt,nwt). 
  ! - In presence of pairing : (QFAM)
  !   dH(:,:,1) and dH(:,:,2) contain 20 and 02 elements of the perturbed hamiltonian in the HFG basis.
  !   Only the positive qp spectrum is included such that size is still (nwt,nwt)
  !   
  complex(KIND=dp), allocatable :: dH_free_flat(:) ! free response of Hamiltonian in the HF basis
  ! The free response is obtained by performing one complete FAM loop starting from dH=0
  ! - In absence of pairing : (FAM)
  !   dH_free_flat contains the free perturbed sp hamiltonian dh(:,:) in HF basis as a flat
  !   array of length (nwt x nwt).  
  ! - In presence of pairing : (QFAM)
  !   dH_free_flat stacks the free perturbed sp hamiltonian dh(:,:) and pairing fields ddelta_plus 
  !   and ddelta_minus in the HF basis as a flat array of length (nwt x nwt x 3)
  !-----------------------------------------------------------------------------
  ! external field
  complex(KIND=dp), allocatable :: F(:,:,:)  ! perturbed external field in sp (FAM) or
  !                                   | | |    qp (QFAM) basis, size (nwt,nwt,2)
  ! - same remark as dH(:,:,:)        | | '-> 1 : ph/20 or 2 : hp/02 component 
  !                                   | '-> sp/qp index
  !                                   '-> sp/qp index 
  integer :: l, m ! anuglar momentum and projection quantum number of the multipole moment
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

  subroutine inifam(omega, DensUnper, PotUnper, Finfile)
    implicit none
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
        F = get_f_LK(l, m, eff_charge_n, eff_charge_p)
      endif
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
      allocate(dH(nwt,nwt,2)) ! stores dHph (dH20), dHhp (dH02) in HF(B) basis for (Q)FAM
    endif

    dH = 0

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! initialise the free response of the sp hamiltonian
    if(.not.allocated(dH_free_flat)) then 
      if(pairingtype==0) then ! FAM
        allocate(dH_free_flat(nwt * nwt)) ! stores dh in HF basis
      else ! QFAM
        allocate(dH_free_flat(3 * nwt * nwt)) ! stores dh, ddelta+, ddelta- in HF basis
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
      ! todo : do not allocate dkappa in asbence of pairing
      !        this requires to modify densit_offdiag to optional arguments

    endif
  
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
    ! Perform one FAM loop of the perturbed single-particle hamiltonian dh
    ! (in HF basis), which contain dh and ddelta (in the QFAM). 
    ! 
    ! Input:
    !    dHsp_flat    : perturbed hamiltonian in HF basis as a flat array
    ! Output:
    !    dHspout_flat : iterated perturbed hamiltonian in HF basis as a flat array
    ! 
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    !
    ! Note: in FAM dHsp_flat only contains dh, while in QFAM it stacks dh, 
    !       ddelta+ and ddelta-
    ! 
    ! One full FAM iterations consists of 6 steps : 
    ! 
    ! (1) transform dh, ddelta+/- to QP basis         => dH20, dH02
    ! (2) compute XY from linear response equation    => X20, Y02
    ! (3) transform XY back to sp basis               => drho, dkappa+/-
    ! (4) calculate perturbed densities on the mesh   => dRs, dRa (DensityVector)
    ! (5) compute perturbed fields on the mesh        => dFs, dFa (PotentialVector)
    ! (6) compute perturbed sp hamiltonian and paring => dh, ddelta+/-
    !
    !---------------------------------------------------------------------------
    1 format('||dH20||² = ', es10.3, '     ||dH02||² = ', es10.3)
    12 format('||dh||² = ', es10.3, '     ||ddelta+||² = ', es10.3, '     ||ddelta-||² = ', es10.3)
    2 format('||X||² = ', es10.3, '     ||Y||² = ', es10.3)
    22 format('||drho||² = ', es10.3, '     ||dkappa+||² = ', es10.3, '     ||dkappa-||² = ', es10.3)
    3 format(' S_',i1,i1,' (', f5.2, ') = ', es18.8)

    implicit none
    complex(KIND=dp), dimension(:), target, intent(in)   :: dHsp_flat
    complex(KIND=dp), dimension(:), target, intent(out)  :: dHspout_flat

    complex(KIND=dp), pointer :: dHsp(:,:,:), dHspout(:,:,:)

    real(KIND=dp) :: strength

    if (fam_verbose > 1) print *, "iterate_dH :: starting full FAM loop "

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! (1) unpack the flat vector to dh, ddelta+/- and transform to QP basis dH20 dH02

    if (pairingtype==0) then ! FAM

      ! pointer remapping for reshaping 1D flat arrays into one 2D matrices
      ! in absence of pairing, dHsp contains only normal field dhsp of sp hamiltonian in HF basis
      dHsp(1:nwt,1:nwt,1:1) => dHsp_flat(:)
      dHspout(1:nwt,1:nwt,1:1) => dHspout_flat(:)

      ! get the ph and hp subblocks of the perturbed sp hamiltonian
      call get_ph_hp_blocks(dHsp(:,:,1), dH(:,:,1), dH(:,:,2))

      if (fam_verbose>1) print 1, sum( abs(dH(:,:,1))**2) , sum( abs(dH(:,:,2))**2) 

    else ! QFAM

      ! pointer remapping for reshaping 1D flat arrays into three 2D matrices
      ! dHsp contains sp hamiltonian in HF basis: normal field + two pairing fields [ddelta+, dh, ddelta-]
      !    dH(:,:,1) = ddelta+ = dH20, dH(:,:,2) = dh = dH11, dH(:,:,3) = ddelta- = dH02 
      dHsp(1:nwt,1:nwt,1:3) => dHsp_flat(:)
      dHspout(1:nwt,1:nwt,1:3) => dHspout_flat(:)

      ! transform the perturbed hamiltonian to the qp basis only interested in dH20 and dH02 components
      call transform_sp_to_qp(Bogoliubov, O20sp=dHsp(:,:,1), O11sp=dHsp(:,:,2), O02sp=dHsp(:,:,3), O20qp=dH(:,:,1), O02qp=dH(:,:,2))

      if (fam_verbose>1) print 1, sum( abs(dH(:,:,1))**2) , sum( abs(dH(:,:,2))**2) 

    endif

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! (2) compute X and Y amplitudes from linear response equation
   
    call calculate_XY(dH)

    call store_XY_hist()

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
    else ! QFAM
      call transform_qp_to_sp(Bogoliubov, O20qp=X, O02qp=Y, O20sp=dkappa_plus, O11sp=drho, O02sp=dkappa_minus)
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
    call calc_perturbed_potentials(RUnper, dRs, dRa, dR_pp_plus, dR_pp_minus, dFs, dFa, dF_pp_plus, dF_pp_minus)

    ! We add in all additional contributions to F_I_I that do not 
    !  result from the Skyrme functional.  
    call combine_potentials(dFs)
    call combine_potentials(dFa)

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! (6) calculate perturbed hamiltonian and pairing in the HF basis

    if (pairingtype==0) then ! FAM
      
      ! construct the sp hamiltonian in HF basis
      dHspout(:,:,1) = calc_sphamil_me( HFpsi, HFdpsi, HFddpsi,dFs, dFa, .false.)

      if (fam_verbose > 0) print 12,  sum(abs(dHspout(:,:,1))**2), 0.0, 0.0

    else ! QFAM

      ! construct the sp hamiltonian + pairing fields in HF basis
      dHspout(:,:,1) = calc_delta_me(   HFpsi, HFdpsi, HFddpsi, dF_pp_plus, .false.)
      dHspout(:,:,2) = calc_sphamil_me( HFpsi, HFdpsi, HFddpsi, dFs, dFa, .false.)
      dHspout(:,:,3) = calc_delta_me(   HFpsi, HFdpsi, HFddpsi, dF_pp_minus, .false.)

      if (fam_verbose > 0) print 12,  sum(abs(dHspout(:,:,2))**2), sum(abs(dHspout(:,:,1))**2),  sum(abs(dHspout(:,:,3))**2)

    endif

    if(fam_verbose > 2) call print_all_fam_spmat()

  end subroutine iterate_dHsp


  subroutine partial_FAM_XY_to_dH(X, Y, dHspout_flat)
    !---------------------------------------------------------------------------
    ! Perform a partial FAM loop, starting from X and Y get the induced perturbed
    ! Hamiltonian
    ! 
    ! Input:
    !    X, Y         :  X Y FAM amplitudes 
    ! Output:
    !    dHspout_flat : iterated perturbed hamiltonian in HF basis as a flat array
    ! 
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    !
    ! Note: 
    !   In FAM dHsp_flat only contains dh, while in QFAM it stacks dh, ddelta+ and ddelta-
    ! 
    !   This function only is almost copy-paste of iterate_dHsp, skipping steps (1) and (2)
    ! 
    ! (3) transform XY to sp basis                    => drho, dkappa+/-
    ! (4) calculate perturbed densities on the mesh   => dRs, dRa (DensityVector)
    ! (5) compute perturbed fields on the mesh        => dFs, dFa (PotentialVector)
    ! (6) compute perturbed sp hamiltonian and paring => dh, ddelta+/-
    !
    !---------------------------------------------------------------------------
    
    1 format('||dH20||² = ', es10.3, '     ||dH02||² = ', es10.3)
    12 format('||dh||² = ', es10.3, '     ||ddelta+||² = ', es10.3, '     ||ddelta-||² = ', es10.3)
    2 format('||X||² = ', es10.3, '     ||Y||² = ', es10.3)
    22 format('||drho||² = ', es10.3, '     ||dkappa+||² = ', es10.3, '     ||dkappa-||² = ', es10.3)
    3 format(' S_',i1,i1,' (', f5.2, ') = ', es18.8)

    implicit none
    complex(KIND=dp), dimension(:,:), intent(in) :: X, Y
    complex(KIND=dp), dimension(:), target, intent(out)  :: dHspout_flat
    complex(KIND=dp), pointer :: dHspout(:,:,:)


    if (fam_verbose > 1) print *, "partial_FAM_XY_to_dH :: starting partial FAM loop from X and Y"

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Set up the pointer remap

    if (pairingtype==0) then ! FAM

      ! pointer remapping for reshaping 1D flat arrays into one 2D matrices
      ! in absence of pairing, dHsp contains only normal field dhsp of sp hamiltonian in HF basis
      dHspout(1:nwt,1:nwt,1:1) => dHspout_flat(:)


    else ! QFAM

      ! pointer remapping for reshaping 1D flat arrays into three 2D matrices
      ! dHsp contains sp hamiltonian in HF basis: normal field + two pairing fields [ddelta+, dh, ddelta-]
      !    dH(:,:,1) = ddelta+ = dH20, dH(:,:,2) = dh = dH11, dH(:,:,3) = ddelta- = dH02 
      dHspout(1:nwt,1:nwt,1:3) => dHspout_flat(:)

    endif

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! (3) Obtain perturbed (pairing) density matrices in HF basis

    if (pairingtype==0) then ! FAM
      drho = X  + transpose(Y)
      dkappa_plus  = 0  
      dkappa_minus = 0
    else ! QFAM
      call transform_qp_to_sp(Bogoliubov, O20qp=X, O02qp=Y, O20sp=dkappa_plus, O11sp=drho, O02sp=dkappa_minus)
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
    call calc_perturbed_potentials(RUnper, dRs, dRa, dR_pp_plus, dR_pp_minus, dFs, dFa, dF_pp_plus, dF_pp_minus)

    ! We add in all additional contributions to F_I_I that do not 
    !  result from the Skyrme functional.  
    call combine_potentials(dFs)
    call combine_potentials(dFa)

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! (6) calculate perturbed hamiltonian and pairing in the HF basis

    if (pairingtype==0) then ! FAM
      
      ! construct the sp hamiltonian in HF basis
      dHspout(:,:,1) = calc_sphamil_me( HFpsi, HFdpsi, HFddpsi,dFs, dFa, .false.)

    else ! QFAM

      ! construct the sp hamiltonian + pairing fields in HF basis
      dHspout(:,:,1) = calc_delta_me(   HFpsi, HFdpsi, HFddpsi, dF_pp_plus, .false.)
      dHspout(:,:,2) = calc_sphamil_me( HFpsi, HFdpsi, HFddpsi, dFs, dFa, .false.)
      dHspout(:,:,3) = calc_delta_me(   HFpsi, HFdpsi, HFddpsi, dF_pp_minus, .false.)

    endif

    if(fam_verbose > 2) call print_all_fam_spmat()

  end subroutine partial_FAM_XY_to_dH


  subroutine Multiply_XY_with_QRPAmat(X, Y, omega, F)
    !---------------------------------------------------------------------------
    ! Multiply X and Y by the QRPA matrix by performing one adjusted FAM loop. 
    ! i.e.
    !        (E - omega) * X + dH20(X, Y) = - F20
    !        (E + omega) * Y + dH02(X, Y) = - F02
    !  
    ! Input:
    !    X, Y     :  X Y input amplitudes 
    !    omega    :  frequency used in the linear response
    ! Output:
    !    F20, F02 :  induced external field 
    ! 
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    !
    !  (1) Compute the induced perturbed hamiltonian dH in the HF basis by calling 
    !      partial_FAM_XY_to_dH()            =>  dHsp_flat = {dh, ddelta+, ddelta-}
    !  (2) Transform dHsp to the QP basis    =>  dH20, dH02
    !  (3) Compute the resulting field F by linear response, i.e. the Eq. above
    !
    !---------------------------------------------------------------------------
    
    1 format('||dH20||² = ', es10.3, '     ||dH02||² = ', es10.3)
    12 format('||dh||² = ', es10.3, '     ||ddelta+||² = ', es10.3, '     ||ddelta-||² = ', es10.3)
    2 format('||X||² = ', es10.3, '     ||Y||² = ', es10.3)
    22 format('||drho||² = ', es10.3, '     ||dkappa+||² = ', es10.3, '     ||dkappa-||² = ', es10.3)
    3 format(' S_',i1,i1,' (', f5.2, ') = ', es18.8)

    implicit none
    complex(KIND=dp), intent(in) :: X(:,:), Y(:,:)
    complex(KIND=dp), intent(in) :: omega
    complex(KIND=dp), intent(out) :: F(:,:,:)
    complex(KIND=dp), allocatable, target :: dHsp_flat(:)
    complex(KIND=dp), pointer :: dHsp(:,:,:)


    if (fam_verbose > 1) print *, "Multiply_with_QRPAmat :: compute the external field induced by XY"


    if(.not. allocated(dHsp_flat)) then
      if(pairingtype==0) then
        allocate(dHsp_flat(nwt * nwt))
      else
        allocate(dHsp_flat(3 * nwt * nwt))
      endif
    endif


    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! (1) Compute the induced perturbed hamiltonian dH in the HF basis
    
    call partial_FAM_XY_to_dH(X, Y, dHsp_flat)


    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! (2) unpack via pointer remap and transfrom dH to the qp basis

    if (pairingtype==0) then ! FAM

      ! pointer remapping for reshaping 1D flat arrays into one 2D matrices
      ! in absence of pairing, dHsp contains only normal field dhsp of sp hamiltonian in HF basis
      dHsp(1:nwt,1:nwt,1:1) => dHsp_flat(:)

      ! get the ph and hp subblocks of the perturbed sp hamiltonian
      call get_ph_hp_blocks(dHsp(:,:,1), dH(:,:,1), dH(:,:,2))

      if (fam_verbose>1) print 1, sum( abs(dH(:,:,1))**2) , sum( abs(dH(:,:,2))**2) 

    else ! QFAM

      ! pointer remapping for reshaping 1D flat arrays into three 2D matrices
      ! dHsp contains sp hamiltonian in HF basis: normal field + two pairing fields [ddelta+, dh, ddelta-]
      !    dH(:,:,1) = ddelta+ = dH20, dH(:,:,2) = dh = dH11, dH(:,:,3) = ddelta- = dH02 
      dHsp(1:nwt,1:nwt,1:3) => dHsp_flat(:)

      ! transform the perturbed hamiltonian to the qp basis only interested in dH20 and dH02 components
      call transform_sp_to_qp(Bogoliubov, O20sp=dHsp(:,:,1), O11sp=dHsp(:,:,2), O02sp=dHsp(:,:,3), O20qp=dH(:,:,1), O02qp=dH(:,:,2))

      if (fam_verbose>1) print 1, sum( abs(dH(:,:,1))**2) , sum( abs(dH(:,:,2))**2) 

    endif

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! (3) compute induced external field F20 F02
    
    call compute_F_from_XYdH(X, Y, dH, omega, F)


  end subroutine Multiply_XY_with_QRPAmat


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
    ! 
    ! Input:
    !    dHsp_flat    : perturbed hamiltonian in HF basis as a flat array
    ! Output:
    !    dHspout_flat : (I - T) * dHsp_flat
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
    ! Compute the X and Y amplitudes from the FAM master equation. 
    ! In absense of pairing, X, Y, dH, F store ph subblocks and loops are only 
    ! over ph pairs. In presence of pairing, X, Y, dH and F store qp matrix 
    ! elements and loops run over the complete qp basis.  
    !---------------------------------------------------------------------------
    implicit none
    complex(KIND=dp), intent(in)  :: dH(:,:,:) ! perturbed H in QP basis

    integer       :: i, j, si, si2, N, N2, B, T

    if (fam_verbose > 1) print *, "calculate_XY :: update X and Y"

    X = - (F(:,:,1) + dH(:,:,1))
    Y = - (F(:,:,2) + dH(:,:,2))


    ! normalise with energy denominator

    if(pairingtype==0) then ! FAM : difference of particle and hole energy
      si = 0
      do B=1,8,2
        N  = HFblocks(B)    ; if(N.eq.0) cycle 
        N2 = HFblocks(B+1)
        T = N + N2
        ! run over particle-hole pairs. hole (j) as outer, particle (i) as inner loop
        do j = 1, T
          if(rho_can(si+j) < 1d-6) cycle  ! skip if j is not a hole state
          do i = 1, T
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

    integer       :: i, j, si, si2, N, N2, B, T

    if (fam_verbose > 1) print *, "compute_F_from_XYdH ::"

    F = - dH

    if(pairingtype==0) then ! FAM : difference of particle and hole energy
      si = 0
      do B=1,8,2
        N  = HFblocks(B)    ; if(N.eq.0) cycle 
        N2 = HFblocks(B+1)
        T = N + N2
        ! run over particle-hole pairs. hole (j) as outer, particle (i) as inner loop
        do j = 1, T
          if(rho_can(si+j) < 1d-6) cycle  ! skip if j is not a hole state
          do i = 1, T
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


  subroutine store_XY_hist()
    !---------------------------------------------------------------------------
    ! Store the current X and Y into their histories. 
    !---------------------------------------------------------------------------
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
      
    allocate(rho_pairing(nwt,nwt))
    allocate(kappa_pairing(nwt,nwt))
  
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
        ! loop over unique qp pairs, i.e. j < i
        do j = si+1, si+T
          do i = j+1, si+T
             S = S + conjg(F(i,j,1)) * X(i,j) + conjg(F(i,j,2)) * Y(i,j)
          enddo
        enddo
        si  = si + T 
      enddo

    endif

    $TR S = 2 * S ! Time-reversal factor 2


    strength_complex = S 
    strength = - strength_complex%im / pi

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
        ! loop over unique qp pairs, i.e. j < i
        do j = si+1, si+T
          do i = j+1, si+T
             S_complex(B) = S_complex(B) + conjg(F(i,j,1)) * X(i,j) + conjg(F(i,j,2)) * Y(i,j)
          enddo
        enddo
        si  = si + T 
      enddo

    endif


$TR    S_complex(:) = 2.0 * S_complex(:) ! Time-reversal factor 2
    strength(:) = - S_complex(:)%im / pi

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


    implicit none
    integer, intent(in) :: L, K
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


      ! TODO: investigate signs in Q20 which seems suspicious in O16 nwt24 test case
      ! 3rd row/col in sym block 1 differs in sign wrt blocks 2, 5 and 6. 


      endif


      ! Multiply the operator by the effective charges 
      f_LK_spme(1:nwn,1:nwn) = eff_e_n * f_LK_spme(1:nwn,1:nwn)
      f_LK_spme(nwn+1:,nwn+1:) = eff_e_p * f_LK_spme(nwn+1:,nwn+1:)

     
      if(fam_verbose > 2) then
        print *, 'f^+_LK'
       call print_spme_complex(f_LK_spme)
       print *, '||f||²', sum(abs(f_LK_spme)**2)
     endif

      ! note: 
      !   Stoitsov PRC 84 (2011) normalises the external field by a parameter
      !   alpha converting the units of the perturbation to MeV, and eventually 
      !   devides the obtained strength by alpha. 


      ! TODO: write a general transformation routine from the mesh to any 
      !       single-particle basis


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
      ! transformation and storing the F^20 anf F^02 comnpnents
      call transform_sp_to_qp(Bogoliubov, O11sp=f_LK_spme, O20qp=f_LK_qpme(:,:,1), O02qp=f_LK_qpme(:,:,2))

      ! ! Assuming that F is Hermitian, then F20 = F02^*
      ! f_LK_qpme(:,:,2) = conjg(f_LK_qpme(:,:,1))

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

  end function


  function calc_EWSR() result (ewsr)
    !---------------------------------------------------------------------------
    ! Compute the energy weighted sum rule from a ground-state 
    ! expectation value. When Thouless' theorem is applicable, then this value 
    ! should equal the first-moment of the strength function, i.e.
    ! m_1(F) = int_0^inf dE E S(E, F). 
    ! Expressions are taken from N. Hinohara PRC 91, 044323 (2015)
    !---------------------------------------------------------------------------
    real(KIND=dp) :: ewsr
    type(Moment), pointer  :: moment_ptr

    ewsr = 0

    ! isovector perturbations
    if(eff_charge_n .ne. eff_charge_p) then
      print *, 'NOT IMPLEMENTED: only isoscalar pertubations are implemented for now'
      ! this would require an enhancement factor kappa
      return
    endif

    ! isoscalar monopole
    if(l == 0) then
      moment_ptr => FindMoment(-2,0,.false.) ! pointer to <r_ch^2>
      ewsr =  4.0 * eff_charge_p**2 * hbm(1) * sum(moment_ptr%Value)
      ! note that moment_ptr%Value contains a factor A 

    ! isoscalar quadrupole
    else if(l == 2) then
      moment_ptr => FindMoment(-2,0,.false.) ! pointer to <r_ch^2>
      ewsr = (5.0 / (2.0 * pi)) * eff_charge_p**2 * hbm(1) * sum(moment_ptr%Value)
      ! note that moment_ptr%Value contains a factor A 


      ! deformation correction, still to be worked out for more general shapes. 
      print *, 'INCOMPLETE: deformation correction for EWSR assumes axial shape '

      ! For axial nuclei, correction with mass quadruple deformation beta20
      moment_ptr => FindMoment(2,0,.false.) ! pointer to <Q_20>
      ewsr = ewsr * (1 + sqrt(5./(4.*pi)) * moment_ptr%beta(4))

    else 
      print *, 'NOT IMPLEMENTED: only monopole (l=0) and quadrupole (l=2) EWSR implemented for now'
      return
    endif


    print *, "Energy weighted sum rule : m1 = ", ewsr

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



  subroutine transform_sp_to_qp(Bogo, O20sp, O11sp, O02sp, O20qp, O11qp, O02qp)
    !---------------------------------------------------------------------------
    ! Performing quasi-particle transformation of a generic on 1-body operator 
    ! O = O20sp + O11sp + O02sp. The function returns the matrix elements in of 
    ! O in the operator in the QP-basis.

    ! Input:
    !    Bogo             : Bogoliubov transformation matrix W from sp to qp basis 
    !                       (2*nwt,2*nwt)
    !    O20sp (optional) : sp matrix elements of 20 operator component (nwt,nwt)
    !    O11sp (optional) : sp matrix elements of 11 operator component (nwt,nwt)
    !    O02sp (optional) : sp matrix elements of 02 operator component (nwt,nwt)
    ! Output:
    !    O20qp (optional) : qp matrix elements of 20 operator component (nwt,nwt)
    !    O11qp (optional) : qp matrix elements of 11 operator component (nwt,nwt)
    !    O02qp (optional) : qp matrix elements of 02 operator component (nwt,nwt)
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
    ! QP matrix elements are obtained from 
    ! 
    !       O20qp = + U^{dagger} o11sp   V^* + U^{dagger} o20sp   U^* 
    !               - V^{dagger} o02sp^* V^* - V^{dagger} o11sp^T U^*
    ! 
    !       O11qp = + U^{dagger} o11sp   U   + U^{dagger} o20sp   V
    !               - V^{dagger} o02sp^* U   - V^{dagger} o11sp^T V
    ! 
    !       O02qp = - V^T        o11sp   U   - V^T        o20sp   V 
    !               + U^T        o02sp^* U   + U^T        o11sp^T V
    ! 
    ! Note that for a hermitian operator O for which o20sp = o02sp^* and 
    ! o11sp^*=o11sp^T, the O02 component can be obtained trivially from O20
    ! by complex conjugation
    !      O20qp^* = + U^T  o11sp^* V + U^T  o20sp^*        U 
    !                - V^T  o02sp   V - V^T  o11sp^{dagger} U
    !              = O02qp
    ! 
    !---------------------------------------------------------------------------

    implicit none
    real(KIND=dp), intent(in)               :: Bogo(:,:)
    complex(KIND=dp), intent(in) , optional :: O20sp(:,:), O11sp(:,:), O02sp(:,:)
    complex(KIND=dp), intent(out), optional :: O20qp(:,:), O11qp(:,:), O02qp(:,:)

    real(KIND=dp), allocatable    :: Ub(:,:), Vb(:,:)
    complex(KIND=dp), allocatable :: O20b(:,:), O11b(:,:), O02b(:,:), temp(:,:)
    integer                       :: B, N, N2, si, sb, T, i
    real(KIND=dp)                 :: Tphase

    if(present(O20qp)) O20qp = 0._dp
    if(present(O11qp)) O11qp = 0._dp
    if(present(O02qp)) O02qp = 0._dp

$NTR Tphase = 1.0_dp
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

      if(present(O20sp)) O20b = O20sp(si+1:si+T,si+1:si+T)
      if(present(O11sp)) O11b = O11sp(si+1:si+T,si+1:si+T)
      if(present(O02sp)) O02b = O02sp(si+1:si+T,si+1:si+T)
  
      if (fam_verbose > 2) print '(A, I3, I3, A, I5)', 'Blocks: ', B, B+1, ' with size', T
      if (fam_verbose > 2) print '(A, F10.2)',  '||U||^2 = ', sum(Ub(:,:) * Ub(:,:))
      if (fam_verbose > 2) print '(A, F10.2)',  '||V||^2 = ', sum(Vb(:,:) * Vb(:,:))

      ! print * , 'U : '
      ! do i=1,T
      !   print "(99f10.5)",  Ub(i, 1:T)
      ! enddo
      
      ! print * , 'V : '
      ! do i=1,T
      !   print "(99f10.5)",  Vb(i, 1:T)
      ! enddo

      if(present(O20qp)) then
        if(present(O11sp)) then
          O20qp(si+1:si+T, si+1:si+T) = O20qp(si+1:si+T, si+1:si+T) + Tphase * matmul(transpose(Ub),  matmul(          O11b , Vb))
          O20qp(si+1:si+T, si+1:si+T) = O20qp(si+1:si+T, si+1:si+T) - matmul(transpose(Vb),  matmul(transpose(O11b), Ub))
          if (fam_verbose > 2) then
            temp = + Tphase * matmul(transpose(Ub),  matmul(          O11b , Vb))
            print * , ' O20 qp term 1 : ', temp(si+1,si+1)
            temp = - matmul(transpose(Vb),  matmul(transpose(O11b), Ub))
            print * , ' O20 qp term 2 : ', temp(si+1,si+1)
          endif
        endif
        if(present(O20sp)) then
          O20qp(si+1:si+T, si+1:si+T) = O20qp(si+1:si+T, si+1:si+T) + matmul(transpose(Ub),  matmul(          O20b , Ub)) 
          if (fam_verbose > 2) then
            temp = + matmul(transpose(Ub),  matmul(          O20b , Ub)) 
            print * , ' O20 qp term 3 : ', temp(si+1,si+1)
          endif
        endif
        if(present(O02sp)) then
          O20qp(si+1:si+T, si+1:si+T) = O20qp(si+1:si+T, si+1:si+T) - matmul(transpose(Vb),  matmul(          O02b , Vb))
          if (fam_verbose > 2) then
            temp = - matmul(transpose(Vb),  matmul(          O02b , Vb))
            print * , ' O20 qp term 4 : ', temp(si+1,si+1)
          endif
        endif
      endif 


      if(present(O11qp)) then
        if(present(O11sp)) then
          O11qp(si+1:si+T, si+1:si+T) = O11qp(si+1:si+T, si+1:si+T) + matmul(transpose(Ub),  matmul(          O11b , Ub))
          O11qp(si+1:si+T, si+1:si+T) = O11qp(si+1:si+T, si+1:si+T) - matmul(transpose(Vb),  matmul(transpose(O11b), Vb))
          if (fam_verbose > 2) then
            temp = + matmul(transpose(Ub),  matmul(          O11b , Ub))
            print * , ' O11 qp term 1 : ', temp(si+1,si+1)
            temp = - matmul(transpose(Vb),  matmul(transpose(O11b), Vb))
            print * , ' O11 qp term 2 : ', temp(si+1,si+1)
          endif
        endif
        if(present(O20sp)) then
          O11qp(si+1:si+T, si+1:si+T) = O11qp(si+1:si+T, si+1:si+T) + matmul(transpose(Ub),  matmul(          O20b , Vb))
          if (fam_verbose > 2) then
            temp = + matmul(transpose(Ub),  matmul(          O20b , Vb))
            print * , ' O11 qp term 3 : ', temp(si+1,si+1)
          endif
        endif
        if(present(O02sp)) then
          O11qp(si+1:si+T, si+1:si+T) = O11qp(si+1:si+T, si+1:si+T) - Tphase * matmul(transpose(Vb),  matmul(          O02b , Ub))
          if (fam_verbose > 2) then
            temp = - Tphase * matmul(transpose(Vb),  matmul(          O02b , Ub))
            print * , ' O11 qp term 4 : ', temp(si+1,si+1)
          endif
        endif
      endif


      if(present(O02qp)) then
        if(present(O11sp)) then
          O02qp(si+1:si+T, si+1:si+T) = O02qp(si+1:si+T, si+1:si+T) - Tphase * matmul(transpose(Vb),  matmul(          O11b , Ub))
          O02qp(si+1:si+T, si+1:si+T) = O02qp(si+1:si+T, si+1:si+T) + matmul(transpose(Ub),  matmul(transpose(O11b), Vb))
          if (fam_verbose > 2) then
            temp = - Tphase * matmul(transpose(Vb),  matmul(          O11b , Ub))
            print * , ' O02 qp term 1 : ', temp(si+1,si+1)
            temp = + matmul(transpose(Ub),  matmul(transpose(O11b), Vb))
            print * , ' O02 qp term 2 : ', temp(si+1,si+1)
          endif
        endif
        if(present(O20sp)) then
          O02qp(si+1:si+T, si+1:si+T) = O02qp(si+1:si+T, si+1:si+T) - Tphase * matmul(transpose(Vb),  matmul(          O20b , Vb))
          if (fam_verbose > 2) then
            temp = - Tphase * matmul(transpose(Vb),  matmul(          O20b , Vb))
            print * , ' O02 qp term 3 : ', temp(si+1,si+1)
          endif
        endif
        if(present(O02sp)) then
          O02qp(si+1:si+T, si+1:si+T) = O02qp(si+1:si+T, si+1:si+T) + matmul(transpose(Ub),  matmul(          O02b , Ub))
          if (fam_verbose > 2) then
            temp = + matmul(transpose(Ub),  matmul(          O02b , Ub))
            print * , ' O02 qp term 4 : ', temp(si+1,si+1)
          endif
        endif
      endif 

      si = si +  T
      sb = sb +2*T

    enddo

    if (fam_verbose > 2) then
      print *, 'Symmetry : '
      if(present(O20sp)) print *, '    O20sp = + O20sp^T   : satisfied up to',  sum(abs(O20sp(:,:) - transpose(O20sp(:,:))))
      if(present(O02sp)) print *, '    O02sp = + O02sp^T   : satisfied up to',  sum(abs(O02sp(:,:) - transpose(O02sp(:,:))))
      if(present(O20qp)) print *, '    O20qp = + O20qp^T   : satisfied up to',  sum(abs(O20qp(:,:) - transpose(O20qp(:,:))))
      if(present(O02qp)) print *, '    O02qp = + O02qp^T   : satisfied up to',  sum(abs(O02qp(:,:) - transpose(O02qp(:,:))))
      print *, 'Antisymmetry : '
      if(present(O20sp)) print *, '    O20sp = - O20sp^T   : satisfied up to',  sum(abs(O20sp(:,:) + transpose(O20sp(:,:))))
      if(present(O02sp)) print *, '    O02sp = - O02sp^T   : satisfied up to',  sum(abs(O02sp(:,:) + transpose(O02sp(:,:))))
      if(present(O20qp)) print *, '    O20qp = - O20qp^T   : satisfied up to',  sum(abs(O20qp(:,:) + transpose(O20qp(:,:))))
      if(present(O02qp)) print *, '    O02qp = - O02qp^T   : satisfied up to',  sum(abs(O02qp(:,:) + transpose(O02qp(:,:))))
      print *, 'Hermiticity : '
      if(present(O20sp) .and. present(O02sp)) print *, '    O20sp = O02sp*   : satisfied up to',  sum(abs(O20sp(:,:) - conjg(O02sp(:,:))))
      if(present(O11sp)) print *, '    O11sp = O11sp^T^*   : satisfied up to',  sum(abs(O11sp(:,:) - conjg(transpose(O11sp(:,:)))))
      if(present(O20qp) .and. present(O02qp)) print *, '    O20qp = O02qp*   : satisfied up to',  sum(abs(O20qp(:,:) - conjg(O02qp(:,:))))
      if(present(O11qp)) print *, '    O11qp = O11qp^T^*   : satisfied up to',  sum(abs(O11qp(:,:) - conjg(transpose(O11qp(:,:)))))
      print *, 'Anti-hermiticity : '
      if(present(O20sp) .and. present(O02sp)) print *, '    O20sp = - O02sp*   : satisfied up to',  sum(abs(O20sp(:,:) + conjg(O02sp(:,:))))
      if(present(O11sp)) print *, '    O11sp = - O11sp^T^*   : satisfied up to',  sum(abs(O11sp(:,:) + conjg(transpose(O11sp(:,:)))))
      if(present(O20qp) .and. present(O02qp)) print *, '    O20qp = - O02qp*   : satisfied up to',  sum(abs(O20qp(:,:) + conjg(O02qp(:,:))))
      if(present(O11qp)) print *, '    O11qp = - O11qp^T^*   : satisfied up to',  sum(abs(O11qp(:,:) + conjg(transpose(O11qp(:,:)))))
    endif
    
  end subroutine transform_sp_to_qp


  subroutine transform_qp_to_sp(Bogo, O20qp, O11qp, O02qp, O20sp, O11sp, O02sp)
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
      if(present(O20sp)) print *, '    O20sp = + O20sp^T   : satisfied up to',  sum(abs(O20sp(:,:) - transpose(O20sp(:,:))))
      if(present(O02sp)) print *, '    O02sp = + O02sp^T   : satisfied up to',  sum(abs(O02sp(:,:) - transpose(O02sp(:,:))))
      if(present(O20qp)) print *, '    O20qp = + O20qp^T   : satisfied up to',  sum(abs(O20qp(:,:) - transpose(O20qp(:,:))))
      if(present(O02qp)) print *, '    O02qp = + O02qp^T   : satisfied up to',  sum(abs(O02qp(:,:) - transpose(O02qp(:,:))))
      print *, 'Antisymmetry : '
      if(present(O20sp)) print *, '    O20sp = - O20sp^T   : satisfied up to',  sum(abs(O20sp(:,:) + transpose(O20sp(:,:))))
      if(present(O02sp)) print *, '    O02sp = - O02sp^T   : satisfied up to',  sum(abs(O02sp(:,:) + transpose(O02sp(:,:))))
      if(present(O20qp)) print *, '    O20qp = - O20qp^T   : satisfied up to',  sum(abs(O20qp(:,:) + transpose(O20qp(:,:))))
      if(present(O02qp)) print *, '    O02qp = - O02qp^T   : satisfied up to',  sum(abs(O02qp(:,:) + transpose(O02qp(:,:))))
      print *, 'Hermiticity : '
      if(present(O20sp) .and. present(O02sp)) print *, '    O20sp = O02sp*   : satisfied up to',  sum(abs(O20sp(:,:) - conjg(O02sp(:,:))))
      if(present(O11sp)) print *, '    O11sp = O11sp^T^*   : satisfied up to',  sum(abs(O11sp(:,:) - conjg(transpose(O11sp(:,:)))))
      if(present(O20qp) .and. present(O02qp)) print *, '    O20qp = O02qp*   : satisfied up to',  sum(abs(O20qp(:,:) - conjg(O02qp(:,:))))
      if(present(O11qp)) print *, '    O11qp = O11qp^T^*   : satisfied up to',  sum(abs(O11qp(:,:) - conjg(transpose(O11qp(:,:)))))
      print *, 'Anti-hermiticity : '
      if(present(O20sp) .and. present(O02sp)) print *, '    O20sp = - O02sp*   : satisfied up to',  sum(abs(O20sp(:,:) + conjg(O02sp(:,:))))
      if(present(O11sp)) print *, '    O11sp = - O11sp^T^*   : satisfied up to',  sum(abs(O11sp(:,:) + conjg(transpose(O11sp(:,:)))))
      if(present(O20qp) .and. present(O02qp)) print *, '    O20qp = - O02qp*   : satisfied up to',  sum(abs(O20qp(:,:) + conjg(O02qp(:,:))))
      if(present(O11qp)) print *, '    O11qp = - O11qp^T^*   : satisfied up to',  sum(abs(O11qp(:,:) + conjg(transpose(O11qp(:,:)))))
    endif

  end subroutine transform_qp_to_sp


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

    ! explicitly antisymmetrise
    f_qpme(:,:,1) = 0.5 * (f_qpme(:,:,1) + transpose(f_qpme(:,:,1)))
    f_qpme(:,:,2) = 0.5 * (f_qpme(:,:,2) + transpose(f_qpme(:,:,2)))

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
        print '(A, f8.3, A, f8.3)', 'End of file reached without finding good XY block for omega = ', omega_fam, ' ,smear = ', smear
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

      ! explicitly antisymmetrise
      X(:,:) = 0.5 * (X(:,:) + transpose(X(:,:)))
      Y(:,:) = 0.5 * (Y(:,:) + transpose(Y(:,:)))
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

  end subroutine


  subroutine print_spme_real(A)
    implicit none
    real(kind=dp), intent(in) :: A(:,:)
    integer :: si, B, N, i

    si = 0
    do B=1,8
      N = HFBlocks(B)

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
      N = HFBlocks(B)

      print *, 'BLOCK', B
      do i=si+1,si+N
        print "(*( '(',g12.5,',',g12.5,')',:))",  A(i, si+1:si+N)
      enddo
      print *
      si = si + N
    enddo
    print *
    
  end subroutine print_spme_complex

  subroutine print_spme_complex_superblock(A)
    implicit none
    complex(kind=dp), intent(in) :: A(:,:)
    integer :: si, B, N, i

    si = 0
    do B=1,8,2

      N  = HFblocks(B) + HFblocks(B+1) 

      print *, 'BLOCKS ', B ,' & ', B+1
      do i=si+1,si+N
        print "(*( '(',g12.5,',',g12.5,')',:))",  A(i, si+1:si+N)
      enddo
      print *
      si = si + N
    enddo
    print *
    
  end subroutine print_spme_complex_superblock

end module fam
