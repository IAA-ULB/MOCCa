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
  ! PBROKEN : $PBROKEN
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
  type(DensityVector), target   :: dRs, dRa  ! perturbation to the particle-hole densities on the mesh
  !                                 |    '-> anti-symmetric part
  !                                 '-> symmetric part
  type(PotentialVector) :: dFs, dFa  ! perturbation to the particle-hole potentials on the mesh
  !                         |    '-> anti-symmetric part
  !                         '-> symmetric part
  !
  type(DensityVector), target   :: dR_pp_plus, dR_pp_minus  ! perturbation to the particle-particle densities on the mesh
  !                                 |           '-> associated with kappa_minus
  !                                 '-> associated with kappa^plus 
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
  !   'Zcom'            = center-of mass z-coordinate
  !   'Zmomentum'       = center-of mass z momentum
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
        F = get_external_field(operator_type)
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


  subroutine printfam_init
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
  
  end subroutine printfam_init


  subroutine printfam_end(S_arr, num_iter, residual)
    real(KIND=dp), intent(in) :: S_arr(8)
    integer, intent(in) :: num_iter
    real(KIND=dp), intent(in) :: residual


    11 format(86('='))
    12 format(2x,74('-'))
    13 format(/,2x,32('-'), ' strength ', 32('-'),/)
    20 format ( "   total number of iterations : ", i4)
    21 format ( "                 fam residual :  ", es10.2)
    22 format ( "   omega = ", f8.3)
    23 format ( "   smear = ", f8.3)

    31 format ( '   Operator:   ', /, &
    &          '      F = Q_', i1, i1,/, &
    &          '      neutron eff charge = ', f10.3, ' e', /, &
    &          '      proton eff charge  = ', f10.3, ' e')
    32 format ( '   Operator:   ',30a) 
    33 format ( '   EWSR = ',  es16.6) 
    4 format (20x, '    neutron              proton                total')
    51 format ('    parity +  ', es20.6, es20.6, es20.6)
    52 format ('    parity -  ', es20.6, es20.6, es20.6)
    53 format ('    total     ', es20.6, es20.6)
    54 format ('    total strength :    ',30x, es20.6)



    print *
    print 11
    print 20, num_iter
    print 21, residual
    print *
    print 22, omega_fam
    print 23, smear
    
    if(operator_type=="multipole") then
      print 31, l, m, eff_charge_n, eff_charge_p
    else
      print 32, operator_type
    endif
    print 33, ewsr

    
    print 13
    print 4
    print 12
    print 51, sum(S_arr(1:2)), sum(S_arr(5:6)), sum(S_arr(1:2)) + sum(S_arr(5:6))
    print 52, sum(S_arr(3:4)), sum(S_arr(7:8)), sum(S_arr(3:4)) + sum(S_arr(7:8))
    print 53, sum(S_arr(1:4)), sum(S_arr(5:8))
    print 12
    print 54, sum(S_arr(1:8))
    print 12

    print 11

  end subroutine printfam_end


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
    1 format('||X||² = ', es10.3, '     ||Y||² = ', es10.3)

    complex(KIND=dp), dimension(:), target, intent(in)   :: dHsp_flat
    complex(KIND=dp), dimension(:), target, intent(out)  :: dHspout_flat

    complex(KIND=dp), pointer :: dHsp(:,:,:), dHspout(:,:,:)


    integer :: si, i, B, N, N2, T

    if (fam_verbose > 1) print *, "iterate_dH :: starting full FAM loop "

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Steps 1 & 2: from the perturbed hamiltonian in the HF basis to X,Y amplitudes
    call FAM_dh_to_XY(dHsp_flat,X,Y)
    ! .... but don't forget to store them into history for GMRES!
    call store_XY_hist(X,Y)

    if (fam_verbose>1) then
      print *, 'Verify antisymmetry of X and Y'
      print * , '||X + X^T|| = ', sum(abs(X+transpose(X))**2)
      print * , '||Y + Y^T|| = ', sum(abs(Y+transpose(Y))**2)
    endif

    if (fam_verbose>0) then
      print 1, sum( abs(X(:,:))**2) , sum( abs(Y(:,:))**2) 
      strength =  calc_strength()
    endif
  
    ! - - - - -  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Steps 3,4,5,6; from the X & Y amplitudes to the perturbed 
    !                 hamiltonian in the HF basis
    call  FAM_XY_to_dh(X,Y,dHspout_flat)

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
    call calculate_XY(dH,X_local,Y_local)

  end subroutine FAM_dh_to_XY


  subroutine FAM_XY_to_dh(X,Y,dHsp_flat)
    !-------------------------------------------------------------
    ! Compute the induced perturbation to the s.p./q.p. hamiltonian 
    ! from a set of X,Y (Q)FAM amplitudes. 
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input: 
    !  X,Y      : (Q)FAM amplitudes
    ! Output: 
    !  dHsp_flat: the perturbed hamiltonian in the HF-basis
    !             flattened, complex array
    !-------------------------------------------------------------

    12 format('||dh||² = ', es10.3, '     ||ddelta+||² = ', es10.3, '     ||ddelta-||² = ', es10.3)
    22 format('||drho||² = ', es10.3, '     ||dkappa+||² = ', es10.3, '     ||dkappa-||² = ', es10.3)

    complex(KIND=dp), intent(in)          :: X(:,:), Y(:,:)
    complex(KIND=dp), target, intent(out) :: dHsp_flat(:)
    complex(KIND=dp), pointer             :: dHsp(:,:,:) 

    if (pairingtype==0) then ! FAM
      ! pointer remapping for reshaping 1D flat arrays into one 2D matrices
      ! in absence of pairing, dHsp contains only normal field dhsp of sp hamiltonian in HF basis
      dHsp(1:nwt,1:nwt,1:1) => dHsp_flat(:)
    else ! QFAM
      ! pointer remapping for reshaping 1D flat arrays into three 2D matrices
      !    dH(:,:,1) =  d\Delta^{-,*} 
      !    dH(:,:,2) = dh             
      !    dH(:,:,3) = -d\Delta^{-,*}  
      dHsp(1:nwt,1:nwt,1:3) => dHsp_flat(:)
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

      $NTR call transform_qp_to_sp(Bogoliubov, OTRqp=X, OBLqp=transpose(Y), &  
      $NTR &                       OTRsp=dkappa_plus, OTLsp=drho, OBLsp=dkappa_minus)
    
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

      $TR  call transform_qp_to_sp(Bogoliubov, OTRqp=X, OBLqp=-Y, &  
      $TR  &                       OTRsp=dkappa_plus, OTLsp=drho, OBLsp=dkappa_minus)

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
      dHsp(:,:,1) = calc_sphamil_me( HFpsi, HFdpsi, HFddpsi,dFs, dFa, .false.)

      if (fam_verbose > 0) print 12,  sum(abs(dHsp(:,:,1))**2), 0.0, 0.0

    else ! QFAM
      ! construct the sp hamiltonian + pairing fields in HF basis
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! \delta h, perturbation of the single-particle hamiltonian
      dHsp(:,:,2) = calc_sphamil_me( HFpsi, HFdpsi, HFddpsi, dFs, dFa, .false.)

      ! \delta \Delta^{+} -> stored 'as-is'
      dHsp(:,:,1) = calc_delta_me(   HFpsi, HFdpsi, HFddpsi, dF_pp_plus , .false.)

      ! \delta \Delta^{-}
      dHsp(:,:,3) = calc_delta_me(   HFpsi, HFdpsi, HFddpsi, dF_pp_minus, .false.)
      ! -> stored as (- \delta \Delta^{-,*} )
      dHsp(:,:,3) = - CONJG(dHsp(:,:,3))

    endif

  end subroutine FAM_XY_to_dh


  subroutine Multiply_XY_with_QRPAmat(X, Y, omega, F, dHsp_flat_in)
    !---------------------------------------------------------------------------
    ! Multiply X and Y by the QRPA matrix by performing one adjusted FAM loop. 
    ! i.e.
    !        (E - omega) * X + dH20(X, Y) = - F20
    !        (E + omega) * Y + dH02(X, Y) = - F02
    !  
    ! Input:
    !    X, Y     :  X Y input amplitudes 
    !    omega    :  frequency used in the linear response
    !    dHsp_flat_in (optional)  :  flat array of the perturbed hamiltonian in the HF basis 
    ! Output:
    !    F20, F02 :  induced external field 
    ! 
    !---------------------------------------------------------------------------
    
    1 format('||dH20||² = ', es10.3, '     ||dH02||² = ', es10.3)

    complex(KIND=dp), intent(in) :: X(:,:), Y(:,:)
    complex(KIND=dp), intent(in) :: omega
    complex(KIND=dp), intent(out) :: F(:,:,:)
    complex(KIND=dp), optional, intent(in) :: dHsp_flat_in(:)
    complex(KIND=dp), allocatable, target :: dHsp_flat(:)
    complex(KIND=dp), pointer :: dHsp(:,:,:)!!

    if (fam_verbose > 1) print *, "Multiply_with_QRPAmat :: compute the external field induced by XY"!

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! (1) if not present, compute the induced perturbed hamiltonian dH in the HF basis
    if (.not. present(dHsp_flat_in)) then

      if (.not. allocated(dHsp_flat)) then
        if(pairingtype==0) then
          allocate(dHsp_flat(nwt * nwt))
        else
          allocate(dHsp_flat(3 * nwt * nwt))
        endif
      endif!

      call FAM_XY_to_dH(X, Y, dHsp_flat)

    else
      dHsp_flat = dHsp_flat_in
    endif

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   ! (2) unpack via pointer remap and transfrom dH to the qp basis
   if (pairingtype==0) then ! FAM

      ! pointer remapping for reshaping 1D flat arrays into one 2D matrices
      ! in absence of pairing, dHsp contains only normal field dhsp of sp hamiltonian in HF basis
      dHsp(1:nwt,1:nwt,1:1) => dHsp_flat(:)
      ! get the ph and hp subblocks of the perturbed sp hamiltonian
      call get_ph_hp_blocks(dHsp(:,:,1), dH(:,:,1), dH(:,:,2))!

      if (fam_verbose>1) print 1, sum( abs(dH(:,:,1))**2) , sum( abs(dH(:,:,2))**2) 

    else ! QFAM

      ! pointer remapping for reshaping 1D flat arrays into three 2D matrices
      ! dHsp contains sp hamiltonian in HF basis: normal field + two pairing fields [ddelta+, dh, ddelta-]
      !    dH(:,:,1) = ddelta+ = dH20, dH(:,:,2) = dh = dH11, dH(:,:,3) = ddelta- = dH02 
      dHsp(1:nwt,1:nwt,1:3) => dHsp_flat(:)
      
      ! Transform the perturbed hamiltonian to the qp basis only interested in dH20 and dH02 components
      call transform_sp_to_qp(Bogoliubov, OTRsp=dHsp(:,:,1), OTLsp=dHsp(:,:,2), OBLsp=dHsp(:,:,3), & ! input 
      &                                   OTRqp=dH(:,:,1),   OBLqp=dH(:,:,2))                        ! output
      ! Attention: 1. there is NO (-CONJG) operation for dHsp(:,:,3), because this 
      !               routine takes -d\Delta^{-,*} as input.
      !            2. the routine spits out the 'bottom left' block of the full matrix, 
      !               which is \delta H^{02, T}. We could apply a transpose, but this 
      !               matrix is antisymmetric; exchanging signs works in both T-conserved 
      !               and T-broken cases.
      dH(:,:,2) = -         dH(:,:,2)

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



  subroutine calc_strength_decomp(S_cmplx_arr, S_arr)
    !---------------------------------------------------------------------------
    ! Calculate the complex and the real strength S(omega,F) decomposed over 
    ! the symmetry channels. 
    !
    ! The complex strength is defined as 
    !     strength_complex = Tr (F^dagger * drho)
    !                      = 0.5 * sum_ab (F^20_ab^* X_ab + F^02_ab^* Y_ab)
    ! while the strength  
    !     strength = -1/pi * Im(strength_complex)
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !     -
    ! Output:
    !     S_cmplx_arr   : array of length 8 containing the strength_complex 
    !                     over (isospin, parity, z-sign) symmetry blocks 
    !     S_arr         : array of length 8 containing the real-valued strength 
    !                     over (isospin, parity, z-sign) symmetry blocks 
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Notes:
    !  - I follow the symmetry block conventions of the rest of the code. 
    !    isospin is always respected. If partity is broken, than all the strength
    !    is put in the positive parity block, etc.
    !---------------------------------------------------------------------------

    complex(KIND=dp), intent(out) :: S_cmplx_arr(8) 
    real(KIND=dp), intent(out) :: S_arr(8)

    integer :: i, j, B, N, N2, si, T
    real(KIND=dp) :: occ_h, occ_p

    if (fam_verbose > 1) print *, "calc_strength_decomp :: S_lm where l= ", l, "m=", m

    S_cmplx_arr = 0
    S_arr = 0


   if(pairingtype==0) then ! FAM
      si = 0
      ! loop over 8 isospin-parity-signature (IPS) block 
      do B=1,8
        N  = HFblocks(B)    ; if(N.eq.0) cycle 
        ! run over particle-hole pairs. hole (j) as outer, particle (i) as inner loop
        do j = si+1, si+N
          if(rho_can(j) < 1d-6) cycle  ! skip if j is not a hole state
          do i = si+1, si+N
            S_cmplx_arr(B) = S_cmplx_arr(B) + conjg(F(i,j,1)) * X(i,j) + conjg(F(i,j,2)) * Y(i,j)
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
             S_cmplx_arr(B) = S_cmplx_arr(B) + conjg(F(i,j,1)) * X(i,j) &
                  &                      + conjg(F(i,j,2)) * Y(i,j)
         enddo
        enddo
        S_cmplx_arr(B) = S_cmplx_arr(B) / 2.0d0 ! acount for double counting qp pairs

        si  = si + T
      enddo
    endif

    $TR S_cmplx_arr(:) = 2.0 * S_cmplx_arr(:) ! Time-reversal factor 2
    S_arr(:) = - IMAG(S_cmplx_arr(:)) / pi

    if (fam_verbose > 2) then
      print *, 'Decomposed strength : '
      print * , 'S_n+ : (', S_arr(1), ' , ', S_arr(2), ' )'
      print * , 'S_n- : (', S_arr(3), ' , ', S_arr(4), ' )'
      print * , 'S_p+ : (', S_arr(5), ' , ', S_arr(6), ' )'
      print * , 'S_p- : (', S_arr(7), ' , ', S_arr(8), ' )'
      print * , 'S_tot : ', sum(S_arr(:))
    endif

    ! set the global variables of the module
    strength_complex = sum(S_cmplx_arr(:))
    strength = sum(S_arr(:))

  end subroutine calc_strength_decomp


  function calc_strength() result (res)
    !---------------------------------------------------------------------------
    ! Alternative interface for calc_strength_decomp when only interested in the
    ! total strength. Calling calc_strength_decomp with unused dummy variables. 
    !---------------------------------------------------------------------------

    real(KIND=dp) :: res
    
    complex(KIND=dp) :: S_cmplx_dummy(8) = 0
    real(KIND=dp) :: S_dummy(8) = 0

    call calc_strength_decomp(S_cmplx_dummy, S_dummy)

    ! return the strength
    res = strength

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

    if (fam_verbose > 0) print 2, DX_norm, DY_norm

    if( (DX_norm < fam_precision) .and. (DY_norm < fam_precision)) then
      conv = .true.
    endif

  end subroutine test_convergence

  subroutine subtract_spurious_modes()
    !---------------------------------------------------------------------------
    ! Subtract the spurious modes from the X and Y amplitudes
    ! 
    ! Input:
    !    /
    ! Output:
    !    /
    ! 
    ! Remarks:
    !  - For now, only the subtraction of spurious translational mode is 
    !    subtracted, present when F is parity odd, i.e. L is odd and K = 0, 1
    !
    !  - This involves vacuum expectation values of commutators of 1B operators
    !    <[A,B]> which can be evaluated from their quasi-particle matrix
    !    elements as
    !
    !        <[A,B]> = 1/2 sum_ab(A20_ab B02_ab - B20_ab A02_ab)
    !
    !     -> if both operators are Hermitian such that A20 = A02*, then 
    !        <[A,B]> = Im (sum_ab(A20_ab B02_ab)) i
    !
    !  - Note the sign in the definition of the QRPA excitation operator 
    !    O^+ = X20 - Y02, such that 
    !        <[O^+,A]> = 1/2 sum_ab(X20_ab A02_ab + A20_ab Y02_ab)
    !    and causing several unexpected minus sign elsewhere
    !---------------------------------------------------------------------------

    complex(KIND=DP) :: Rz_qpme(nwt, nwt, 2)
    complex(KIND=DP) :: Pz_qpme(nwt, nwt, 2)
    complex(KIND=DP) :: comm_RP, comm_OR, comm_OP
    complex(KIND=DP) :: lambda_R, lambda_P
    real(KIND=DP) :: a


    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! skip when F=N or if F=Q_LK with even L or K>1
    if (operator_type == 'N' .or. (operator_type == 'multipole' .and. (mod(l,2)==0 .or. m>1)) ) then
      print *, 'No need to subtract translational spurious mode'
      return
    endif

    print *, 'Subtract translational spurious mode'

    a = calc_strength()
    print *, ' S prior =', strength_complex


    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! load the qpme of R and P
    Rz_qpme = get_external_field('Zcom')
    Pz_qpme = get_external_field('Zmomentum')

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Evaluation of commutator expectation value <[R,P]>
    comm_RP = 0.5 * (sum(Rz_qpme(:,:,1) * Pz_qpme(:,:,2)) - sum(Pz_qpme(:,:,1) * Rz_qpme(:,:,2)))
    $TR comm_RP = 2.0 * comm_RP ! account for absence of time-reversed states

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Evaluate the commutator <[O^+,R]>. 
    comm_OR = 0.5 * (sum(X(:,:) * Rz_qpme(:,:,2)) + sum(Rz_qpme(:,:,1) * Y(:,:)))
    $TR comm_OR = 2.0 * comm_OR ! account for absence of time-reversed states

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Evaluate the commutator <[O^+,P]>. 
    comm_OP = 0.5 * (sum(X(:,:) * Pz_qpme(:,:,2)) + sum(Pz_qpme(:,:,1) * Y(:,:)))
    $TR comm_OP = 2.0 * comm_OP ! account for absence of time-reversed states

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! compute lambda parameters 
    lambda_R = comm_OP / comm_RP
    lambda_P = - comm_OR / comm_RP

    print * , 'lambda_R = ', lambda_R
    print * , 'lambda_P = ', lambda_P


    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! subtract from X, Y
    X = X - lambda_R * Rz_qpme(:,:,1) - lambda_P * Pz_qpme(:,:,1)
    Y = Y + lambda_R * Rz_qpme(:,:,2) + lambda_P * Pz_qpme(:,:,2)

    a = calc_strength()
    print *, ' S after =', strength_complex


  end subroutine subtract_spurious_modes


  function get_external_field(op_type) result (f_qpme)
    !---------------------------------------------------------------------------
    ! Get quasi-particle matrix elements of the external field F based on 
    ! operator_type and possibly multipolarity l, m
    ! 
    ! Input:
    !    op_type : operator_type to be loaded. Current options are
    !               - 'multipole' :  Q_lm, the module variables m and l will be used 
    !               - 'particle number' : particle number operator 
    !               - 'zcom' : the z-coordinate operator  
    !               - 'zmomentum' : z-momentum operator
    ! Output:
    !    f_qpme  : quasi-particle matrix elements F20_mn and F02_mn of 
    !              the external field organised as a 3D complex array with 
    !              dimensions (nwt, nwt, 2)
    !                           |    |   '-> 1 : 20,  2 : 02 component 
    !                           |    '-> qp index
    !                           '-> qp index
    ! 
    ! Remarks:
    !  - in case of HF, qpme F20_mn and F02_mn reduce to Fph_ai and Fhp_ai, 
    !    the particle-hole and hole-particle subblocks of the spme, where 'a' is 
    !    an unoccupied sp index and 'i' is an occupied sp index
    !---------------------------------------------------------------------------
 
    character(len=*), intent(in) :: op_type
    complex(KIND=dp), allocatable :: f_qpme(:,:,:)
    complex(KIND=dp), allocatable :: f_spme(:,:)
    real(KIND=dp), allocatable :: nabla_spme(:,:,:,:)
    integer :: i, j

    if (fam_verbose > 1) print *, "get_external_field :: "
      
    allocate(f_qpme(nwt,nwt,2)) 
    allocate(f_spme(nwt,nwt)) 

    !----------------------------------------------------------------------------------
    ! 1) get the single particle matrux elements f_spme

    select case(trim(to_lower(op_type)))
     ! lower to make the selection case insensitive
     ! trim to not bother about string length and possible trailing spaces

     ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
     ! a) get the multipole operator Q_LK
      case('multipole')
        ! if l=0, then one needs r^2 rather then Q_00 for a monopole excitation
        if(l==0) then
          ! Get single-particle matrix elements of R in the HF basis
          f_spme = Rsq_spme()
        else 
          ! Get single-particle matrix elements of Q_lm in the HF basis
          !  -> calling a function in fission_MOI.f90, which returns <i|r^L Re(Y_LK)|j> 
          !     in strange fission units barn^(l/2) = (100 fm^2)^(l/2)
          f_spme = Qlm_spme(l, m, .false.)
          
          ! Convert f_spme to unit fm^l
          f_spme = f_spme * (100**(l/2.0)) 

          ! Normalise with sqrt(2) if K is not 0
          if(m.ne.0) f_spme = f_spme * sqrt(2.0)

          endif
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      ! b) get the particle number operator 
      case('particle number')
        ! Build the single-particle matrix elements of N
        f_spme = 0.0d0
        do i=1,nwt
           f_spme(i,i) = 1.0d0
        enddo

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      ! c) get the c.o.m. z-coordinate operator 
      case('zcom')
        ! Z_com = 1/A sum_i z_i is trivially related to Q_10 = sum_i sqrt(3/4pi) z_i

        ! Get single-particle matrix elements of Q_10 in the HF basis
        f_spme = Qlm_spme(1, 0, .false.)

        ! Convert f_spme to unit fm^l
        f_spme = f_spme * sqrt(100.0) 

        ! Cancel prefactor sqrt(3/4pi)
        f_spme = f_spme * sqrt( 4.0 * pi / 3.0) 

        ! multiply by 1/A
        f_spme = f_spme / (Neutrons + Protons)

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      ! d) get the c.o.m. z momentum = - i [hbar] nabla_z
      !    -> the factor hbar in the definition gets dropped
      case('zmomentum')

        allocate(nabla_spme(3,2,nwt,nwt))
        !                   | |  '---'-> sp indices 
        !                   | '-> real, imag
        !                   '-> x, y, z

        ! Get single-particle matrix elements of nabla in the HF basis
        nabla_spme = CompNablaMelements('HF')

        ! P_z = -i * nabla_z = IM(nabla_z) - Re(nabla_z) i
        f_spme = dcmplx(nabla_spme(3,2,:,:), -nabla_spme(3,1,:,:))



      case DEFAULT
        call stp('Unrecognized operator_type!')
    end select

    !----------------------------------------------------------------------------------
    ! 2) Multiply the single-particle matrix elements by the effective charges 
    f_spme(1:nwn,1:nwn) = eff_charge_n * f_spme(1:nwn,1:nwn)
    f_spme(nwn+1:,nwn+1:) = eff_charge_p * f_spme(nwn+1:,nwn+1:)

    !----------------------------------------------------------------------------------
    ! 3) convert spme to quasiparticle basis
    if (pairingtype==0) then ! FAM
      call get_ph_hp_blocks(f_spme, f_qpme(:,:,1), f_qpme(:,:,2))
    else ! QFAM
      call transform_sp_to_qp(Bogoliubov, OTLsp=f_spme, &
           &                  OTRqp=f_qpme(:,:,1), OBLqp=f_qpme(:,:,2))
      $NTR  f_qpme(:,:,2) = TRANSPOSE(f_qpme(:,:,2))
      $TR   f_qpme(:,:,2) = -         f_qpme(:,:,2)
    endif

    deallocate(f_spme)


    if(fam_verbose > 2) then
      print *, ' f_qpme(:,:,1)'
      call print_spme_complex_superblock( f_qpme(:,:,1))
      print *, ' f_qpme(:,:,2)'
      call print_spme_complex_superblock( f_qpme(:,:,2))
    endif


    print *, '||F(:,:,1)||²', sum(abs(f_qpme(:,:,1))**2)
    print *, '||F(:,:,2)||²', sum(abs(f_qpme(:,:,2))**2)


  end function get_external_field


  function calc_EWSR(R) result (res)
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
    real(KIND=dp) :: res
    real(KIND=dp) :: m1kin=0, kappa=0, Ctau0=0, Ctau1=0
    type(Moment), pointer  :: moment_ptr, r2_ptr


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

    ! set the global variable 
    ewsr = m1kin + kappa

    print *, "Energy-weighted sum rule : m1 = ", ewsr

    ! return ewsr
    res = ewsr

  end function calc_EWSR

  subroutine check_box_size()
    !---------------------------------------------------------------------------
    ! Check the values of drho and dkappa at the boundary of the box
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Notes:
    !   - We use a relative measure, drho(max(box)) / max(drho)
    !---------------------------------------------------------------------------

    integer :: it=1
    real(kind=DP) :: maximum, boundval
    complex(KIND=dp), pointer                  :: ddens3D(:,:,:,:)


    12 format(2x,78('-'))
    13 format(36('-'), ' box size check ', 36('-'))
    14 format(a45)

    2 format (4x,'Xmax = (nx+0.5)dx = ', f8.3, ' fm,  rel.max||d_dens(X=Xmax)|| = ', es10.3 )
    3 format (4x,'Ymax = (ny+0.5)dx = ', f8.3, ' fm,  rel.max||d_dens(Y=Ymax)|| = ', es10.3 )
    4 format (4x,'Zmax = (nz+0.5)dx = ', f8.3, ' fm,  rel.max||d_dens(Z=Zmax)|| = ', es10.3 )
    41 format (4x,'Zmin = -(nz+0.5)dx = ', f8.3, ' fm,  rel.max||d_dens(Z=Zmin)|| = ', es10.3 )
    

    print *
    print 13
    print 12

    print 14, 'drho_sym'
    print 12

    ddens3D(1:nx,1:ny,1:nz,1:2)   => dRs%D_I_I
  
    maximum = maxval(sum(abs(ddens3D(:,:,:,:)),4))

    print 2, meshX(nx) , maxval(sum(abs(ddens3D(nx,:,:,:)),3)) / maximum
    print 3, meshY(ny) , maxval(sum(abs(ddens3D(:,ny,:,:)),3)) / maximum
    print 4, meshZ(nz) , maxval(sum(abs(ddens3D(:,:,nz,:)),3)) / maximum
    $PBROKEN print 41, meshZ(1) , maxval(sum(abs(ddens3D(:,:,1,:)),3)) / maximum

    print 12
    print 14, 'drho_antisym'
    print 12

    ddens3D(1:nx,1:ny,1:nz,1:2)   => dRa%D_I_I
  
    maximum = maxval(sum(abs(ddens3D(:,:,:,:)),4))

    print 2, meshX(nx) , maxval(sum(abs(ddens3D(nx,:,:,:)),3)) / maximum
    print 3, meshY(ny) , maxval(sum(abs(ddens3D(:,ny,:,:)),3)) / maximum
    print 4, meshZ(nz) , maxval(sum(abs(ddens3D(:,:,nz,:)),3)) / maximum
    $PBROKEN print 41, meshZ(1) , maxval(sum(abs(ddens3D(:,:,1,:)),3)) / maximum

    print 12
    print 14, 'dkappa+'
    print 12

    ddens3D(1:nx,1:ny,1:nz,1:2)   => dR_pp_plus%DP_I_I
  
    maximum = maxval(sum(abs(ddens3D(:,:,:,:)),4))

    print 2, meshX(nx) , maxval(sum(abs(ddens3D(nx,:,:,:)),3)) / maximum
    print 3, meshY(ny) , maxval(sum(abs(ddens3D(:,ny,:,:)),3)) / maximum
    print 4, meshZ(nz) , maxval(sum(abs(ddens3D(:,:,nz,:)),3)) / maximum
    $PBROKEN print 41, meshZ(1) , maxval(sum(abs(ddens3D(:,:,1,:)),3)) / maximum

    print 12
    print 14, 'dkappa-'
    print 12

    ddens3D(1:nx,1:ny,1:nz,1:2)   => dR_pp_minus%DP_I_I
  
    maximum = maxval(sum(abs(ddens3D(:,:,:,:)),4))

    print 2, meshX(nx) , maxval(sum(abs(ddens3D(nx,:,:,:)),3)) / maximum
    print 3, meshY(ny) , maxval(sum(abs(ddens3D(:,ny,:,:)),3)) / maximum
    print 4, meshZ(nz) , maxval(sum(abs(ddens3D(:,:,nz,:)),3)) / maximum
    $PBROKEN print 41, meshZ(1) , maxval(sum(abs(ddens3D(:,:,1,:)),3)) / maximum
    print 12
    print *
    print *


  end subroutine check_box_size


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
        $TR   occ_p = 2.0d0 - rho_can(p)   ! WR: Is this not superfluous? I mean, occ_h and occ_p do not actually enter the result? 
        $NTR  occ_p = 1.0d0 - rho_can(p) 
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
        $TR  occ_p = 2.0d0 - rho_can(p) ! WR: Is this not superfluous? I mean, occ_h and occ_p do not actually enter the result? 
        $NTR occ_p = 1.0d0 - rho_can(p) 
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

    if (fam_verbose > 1) print *, "transform_sp_to_qp :: transform 1B operator from sp to qp basis"

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

    if (fam_verbose > 1) print *, "transform_qp_to_sp :: transform 1B operator from qp to sp basis"


    ! initialise the single-particle matrix elements to zero if they are present
    if(present(OTRsp)) OTRsp = 0._dp
    if(present(OTLsp)) OTLsp = 0._dp
    if(present(OBLsp)) OBLsp = 0._dp

    $NTR Tphase = +1.0_dp
    $TR  Tphase = -1.0_dp

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
