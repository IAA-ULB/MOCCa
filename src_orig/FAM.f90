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
  use FAM_testing, only : run_FAM_tests

  implicit none

  !-----------------------------------------------------------------------------
  ! Define some FAM parameters
  real(KIND=dp) :: omega_fam       ! frequency of the perturbing field 
  real(KIND=dp) :: smear = 1.0_dp  ! complex smearing parameter, default 0.5 MeV
  !    Note that the obtained strength is convoluted with a Lorentzian with FWHM 
  !    equal to double this complex shift
  real(KIND=dp) :: eta = 1.0e-3_dp ! small parameter entering derivatives, 
  !    Default currently set to 10-3. In the end, the strength should be 
  !    reasonably indepedent of the choice. 
  integer       :: maxfamiter = 10 ! maximal number of FAM iterations 
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

    ! set omega frequency of perturbation
    omega_fam = omega


    print *, "Initialise FAM matrices" 

    if(.not.allocated(drho)) then 
      allocate(drho(nwt,nwt))
      allocate(dkappa(nwt,nwt))
      ! allocate(dR(2*nwt,2*nwt))
    endif

    if(.not.allocated(X)) then
      allocate(X(nwt,nwt)) 
      allocate(Y(nwt,nwt))
    endif

    if(.not.allocated(X_hist)) then
      allocate(X_hist(hist_max,nwt,nwt)) 
      allocate(Y_hist(hist_max,nwt,nwt))
    endif

    X_hist=0
    Y_hist=0


    if(.not.allocated(F)) then 
      allocate(F(nwt,nwt,2))

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      ! Get the solid harmonics Q_lm(i,j) = < i | r^l Y_lm | j > expressed 
      ! in HF basis. 
      
      allocate(SolidHarmHF(nwt,nwt)) 

      ! Set external field to E2, hardcoded for now
      l = 2
      m = 0
      ImPart = .false. ! real (.false.) , imaginary (.true.) 
      ! note: odd m and Im parts are not implemeted yet


      if(l==0) then
        SolidHarmHF = Rsq_spme()

        print *, 'Rsq'
        call print_spme_real(SolidHarmHF)
        stop

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

      deallocate(SolidHarmHF)

    endif


    ! This can be improved by some element-wise products occ^T @ SolidHarmHF @ occ

    ! Note to future self: for QFAM this will be replaced by a transformation 
    ! to the qp basis. 

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
      dH(:,:,:) = 0
    endif

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! initialise perturbed potentials as 0
    dFs = 0.0_dp * PotentialsUnpert
    dFa = 0.0_dp * PotentialsUnpert

    ! TO DO: replace by a better initialisation routine


    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! initialise the RPA amplitudes 
    call calculate_XY()

    call store_XY_hist()

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! print unperturbed strenght
    print 1, l, m, omega_fam, calc_strength()

  end subroutine inifam


  subroutine calculate_XY()
    !---------------------------------------------------------------------------
    ! Compute the X and Y amplitudes from the FAM master equation
    !---------------------------------------------------------------------------

    implicit none
    integer :: p, h
    real(KIND=dp) :: occ_h, occ_p, e_h, e_p

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

    print *,  sqrt(sum( abs(dH(:,:,1))**2) )

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

    if( (X_norm .ge. 1.0d2) .or. (Y_norm .ge. 1.0d2)) then
      div = .true.
    endif

    ! previous index in hist obtained by rolling back twice and adding one
    idx_prev = modulo(hist_current_idx - 2, hist_max) + 1

    DX_norm = sqrt( sum( abs(X_hist(hist_current_idx,:,:) - X_hist(idx_prev,:,:))**2) )
    ! DX_norm = DX_norm / X_norm

    DY_norm = sqrt( sum( abs(Y_hist(hist_current_idx,:,:) - Y_hist(idx_prev,:,:))**2) )
    ! DY_norm = DY_norm / Y_norm

    print * , "convergence: ||DX|| = ", DX_norm, "   ||DY|| = ", DY_norm

    if( (DX_norm<tol_XY_conv) .and. (DY_norm<tol_XY_conv)) then
      conv = .true.
    endif

  end subroutine test_convergence


  function Rsq_spme() result (Rsq)

    real(kind=dp) , allocatable :: Rsq(:,:)
    integer :: i, j, k, Bi, Bj, Ni, Nj, si, sj
    
    allocate(Rsq(nwt, nwt))

    ! Initialize to zero
    Rsq = 0.0d0

    !--------------------------------------------------------------------------- 
    ! Loop over the neutron single-particle states
    !---------------------------------------------------------------------------
    si = 0
    do Bi = 1, 4
      Ni =  HFBlocks(Bi) ; if(Ni.eq.0) cycle
      do Bj = 1, 4
        Nj = HFBlocks(Bj) ; if(Nj.eq.0) cycle

        if(bj.eq.1) then
          sj = 0
        else
          sj = sum(HFblocks(1:Bj-1))
        endif

        ! Parity selection rule: only single-particle states of identical parity 
        !   (and signature) contribute    
! $PCONSERVED        if(Bi .ne. Bj)        cycle
        ! (THIS is only needed if parity conserved of course)
! $PBROKEN if( Bi .ne. Bj ) cycle

        if(Bi .ne. Bj)  cycle 

        do i=1,NI    
         do j=1,NJ ! PD: cant one restrict it to i, NJ
          ! me = Int d^3r Sum_sigma psi^*_i(r,sigma) psi_j(r,sigma) Rsq(r)
          Rsq(si+i,sj+j) = sum(sum(HFpsi(:,:,si+i)*HFpsi(:,:,sj+j),2) * sum(meshgrid**2,2)) * dv
          ! All these matrix elements are real if 
          ! (i)  we consider only real multipole moments
          ! (ii) time simplex is conserved
          ! 
          ! which means the matrix we store them in is symmetric
          Rsq(sj+j,si+i) = Rsq(si+i,sj+j)
         enddo
        enddo
      enddo
      si = si + NI
    enddo

    !--------------------------------------------------------------------------- 
    ! Loop over the proton single-particle states
    !---------------------------------------------------------------------------
    si = nwn
    do Bi = 5, 8
      Ni =  HFBlocks(Bi) ; if(Ni.eq.0) cycle
      do Bj = 5, 8
        Nj = HFBlocks(Bj) ; if(Nj.eq.0) cycle
        sj = sum(HFblocks(1:Bj-1))

        ! Parity selection rule: only single-particle states of identical parity 
        !   (and signature) contribute    
! $PCONSERVED        if(Bi .ne. Bj)        cycle
        ! (THIS is only needed if parity conserved of course)
! $PBROKEN if( Bi .ne. Bj ) cycle

        if(Bi .ne. Bj)  cycle 

        do i=1,NI    
         do j=1,NJ ! PD: cant one restrict it to i, NJ
          ! me = Int d^3r Sum_sigma psi^*_i(r,sigma) psi_j(r,sigma) Rsq(r)
          Rsq(si+i,sj+j) = sum(sum(HFpsi(:,:,si+i)*HFpsi(:,:,sj+j),2) * sum(meshgrid**2,2)) * dv
          ! All these matrix elements are real if 
          ! (i)  we consider only real multipole moments
          ! (ii) time simplex is conserved
          ! 
          ! which means the matrix we store them in is symmetric
          Rsq(sj+j,si+i) = Rsq(si+i,sj+j)
         enddo
        enddo
      enddo
      si = si + NI
    enddo

  end function Rsq_spme 


end module fam

program run_FAM

  use compilation
  use IO
  use Tantalus, only : print_header, initialize_all_timers
  use fam

  implicit none
  integer :: iteration
  logical :: is_converged, is_divergent
  real(kind=dp)  :: lin_mix_coeff=1.0d-2
  real(kind=dp) :: omega_curr, omega_min=20, omega_max=40, omega_step=+1.0
  integer :: omega_num, omega_index
  real(kind=dp), allocatable :: omega_arr(:), S_arr(:)
  integer, allocatable :: iter_arr(:)
  character(len=100) :: famfilename
  integer :: i, B, si,N

  ! integer :: ifail ! Future dev: required for HFB

  ! Print a nice header with all kinds of relevant info
  call print_header(.true.)

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
  ! -  the frequency \omega_fam of the perturbing field
  ! -  the 'size' of the perturbation to perform the finite differencing
  ! -  a smearing parameter to avoid discontinuities at the poles of the 
  !    response function
  ! 
  call ReadInput()

  !-----------------------------------------------------------------------------
  ! Initalize the matrices for performing derivatives on the mesh
  call inilag()

  !------------------------------------------------------------------------------
  ! Read all information from a .wf file
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

  !--------------------------------------------------------------------------------
  ! Step 0: build explicitly the matrix of the single-particle hamiltonian and
  !         diagonalize it within the subspace spanned by the spwfs read from file
  Density     = densit(rho_can, kappa_pairing)
  Potentials  = calcPotentials(Density)
  sphamil     = Calc_Sphamil(potentials, .true.)

  ! ! Testing printout - left for now
  ! print *, 'BEFORE'
  ! si = 0
  ! do B=1,8
  !   N = HFBLocks(B)

  !   print *, 'BLOCK', B, N
  !   do i=si+1,si+N
  !     print ('(99f10.3)'), sphamil(i, si+1:si+N)
  !   enddo
  !   print *
  !   si = si + N
  ! enddo

  call apply_subspace_rotation(sphamil, HFTransfo, spenergies)
  ! diagonalisation done; now recalculate other quantities
  if(store_derivatives) call deriveHF() ! and update derivatives
  Density     = densit(rho_can, kappa_pairing)
  Potentials  = calcPotentials(Density)
  sphamil     = Calc_Sphamil(potentials, .true.)

  ! ! Testing printout - left for now
  ! print *, 'AFTER'
  ! si = 0
  ! do B=1,8
  !   N = HFBLocks(B)

  !   print *, 'BLOCK', B, N
  !   do i=si+1,si+N
  !     print ('(99f10.3)'), sphamil(i, si+1:si+N)
  !   enddo
  !   print *
  !   si = si + N
  ! enddo
  ! stop

  ! print ('(99f10.3)'), rho_can(:)

  !
  ! Note: there is a silent assumption here that the HF-spectrum is sufficiently
  !       well-converged such that an explicit orthonormalisation will not change
  !       our mean-field state in any meaningful way. In the future, we might want
  !       to resolve the whole "pairing subproblem" again here and check that the
  !       structure does not vary too much.


  !---------------------------------------------------------------------------------
  ! construct the full HF densities rather than the merely the vector rho_can
  if (pairingtype .eq. 0) call iniHFdensities()


  !---------------------------------------------------------------------------------
  ! solving FAM for a range of omega frequencies

  omega_num = int((omega_max - omega_min) / omega_step) + 1

  allocate(omega_arr(omega_num))
  allocate(S_arr(omega_num))
  allocate(iter_arr(omega_num))
  iter_arr = 0

  omega_curr = omega_min

  do omega_index=1, omega_num

    !-------------------------------------------------------------------------------
    ! initialise FAM matrices end set perturbing external field
    call inifam(omega_curr, Potentials)

    if( calc_strength() .ge. 0.1) then
      lin_mix_coeff=1.0d-3
    else
      if( calc_strength() .ge. 0.001) then
        lin_mix_coeff=1.0d-2
      else
        lin_mix_coeff=1.0d-1
      endif
    endif

    ! Run all kinds of unit tests; should be made optional as this includes a stop statement
    ! call run_FAM_tests(X,Y)

    maxfamiter = 10000
    is_converged = .false.
    is_divergent = .false.

    ! Start of the iterations 
    do iteration=1, maxfamiter

      print *, "FAM iteration : ", iteration

      ! build the perturbed densities on the mesh dRs, dRa from X and Y
      call build_perturbed_densities(X, Y, dRs, dRa)

      ! build the perturbed hamiltonian using explicit linearisation of the field
      call build_dH_explicit(Density, dRs, dRa)

      ! calculate X and Y from the perturbed sp Hamil dH
      call calculate_XY()
      
      ! Apply simple linear mixing of X and Y. 
      ! call mix_XY_linear(lin_mix_coeff)
      ! To be replaced with something more fancy in the future

      call store_XY_hist()

      ! call print_all_fam_spmat()

      print *, " S(", omega_curr, ") = ",  calc_strength()

      ! Exit the loop if convergence is achieved.
      if (iteration > 1) then ! at least two iterations to be able to compare
       call test_convergence(is_converged, is_divergent)
        if(is_converged) then
          print *, "Hooray! FAM is converged! "
          ! omega_arr(omega_index) = omega_curr
          ! S_arr(omega_index) = calc_strength()
          iter_arr(omega_index) = iteration
          exit
        endif
        if(is_divergent) then
          print *, "FAM diverges, exiting"
          ! omega_arr(omega_index) = omega_curr
          ! S_arr(omega_index) = calc_strength()
          iter_arr(omega_index) = -iteration
          exit
        endif
      endif
      if (iteration == maxfamiter) then
        print *, "Reached maximal number of iterations, ", maxfamiter
        iter_arr(omega_index) = -maxfamiter
      endif
    enddo

    omega_arr(omega_index) = omega_curr
    S_arr(omega_index) = calc_strength()


    print *, " S(", omega_arr(omega_index), ") = ", S_arr(omega_index)


    omega_curr = omega_curr + omega_step

  enddo


  write (famfilename, fmt='(a2,2i1,a4)') "S_", l, m, ".fam"

  print *, omega_arr
  print *, S_arr
  print *, iter_arr
  call write_fam_strength(omega_arr, S_arr, iter_arr, l, m, famfilename)

  print *, "Reached the end successfully" 


  ! end of one FAM calculation;

end program run_FAM
