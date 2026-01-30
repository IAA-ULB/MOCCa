module fam_testing
  !--------------------------------------------------------------------------------
  ! A selection of tests to check the development of the FAM functionality and
  ! guard it against regressions.
  !--------------------------------------------------------------------------------
  use densities
  use moments
  use Coulombmod, only : solve_coulomb
  use pairing,    only :  rho_can, pairingtype, rho_pairing, kappa_pairing, FermiEnergy
  use fission_MOI
  use functional
  use evolution
  use fam
  use gmres

implicit none

contains
  subroutine run_FAM_tests(X,Y)
    !---------------------------------------------------------------------------------
    ! Catch-all routine to run all predefined unit tests for FAM routines.
    ! Input:
    !   X, Y : forward- and backward amplitudes to start testing perturbed quantities.
    ! Output:
    !   None
    !---------------------------------------------------------------------------------

    1 format ("Test = ", a20, " Success = ", i4)
    9 format (" ------- Start of the FAM testing routines -------")
    complex(KIND=dp), intent(in), allocatable :: X(:,:), Y(:,:)
    integer :: ifail


    print 9
    print *
    call test_HFMatrices_me(ifail)
    print 1, 'SPHAMIL_ME', ifail

    ! call test_potentials(X,Y,ifail)
    ! print 1, 'potentials', ifail

    !call test_potentials(X,Y,ifail)
    ! Attention: this testing routine has serious side effects on the state of the program.
    call test_densit_offdiag(ifail)
    print 1, 'DENSIT_OFFDIAG', ifail

    stop
  end subroutine run_FAM_tests
!
  subroutine test_potentials(X,Y,ifail)
    !------------------------------------------------------------------------
    ! TODO: describe
    !
    !
    !------------------------------------------------------------------------
    integer, intent(out)          :: ifail
    complex(KIND=dp), allocatable, intent(in)  :: X(:,:), Y(:,:)

    complex(KIND=dp), allocatable :: drho(:,:), dkappa(:,:), sphamil_me(:,:)
    real(KIND=dp)                 :: eta

    type(DensityVector)           :: R, dRa, dRs, dR_pp_plus, dR_pp_minus
    type(PotentialVector)         :: F, dFs, dFa, Fnew, dF_pp_minus, dF_pp_plus

    integer :: i


    allocate(drho(nwt,nwt))
    drho = 0.0d0
    drho = X + transpose(Y)

    R       = densit(rho_can, kappa_pairing)
    dRs     = densit_offdiag_ph_symmetric(drho)
    dRa     = densit_offdiag_ph_antisymmetric(drho)

    F       = calcpotentials(R)
    call calc_perturbed_potentials(R,dRs,dRa, dR_pp_plus, dR_pp_minus, dFs, dFa, dF_pp_plus, dF_pp_minus)

    eta = 0.000001
    Fnew    = calcpotentials(R + eta*dRs) + (-1.0d0) * F
    ! call print_deviations('F_I_I  -     symmetric', dFs%F_I_I       , Fnew%F_I_I/eta)

    ! call print_deviations('F_I_SX -     symmetric', dFs%F_I_S(:,1,:), Fnew%F_I_S(:,1,:)/eta)
    ! call print_deviations('F_I_SY -     symmetric', dFs%F_I_S(:,2,:), Fnew%F_I_S(:,2,:)/eta)
    ! call print_deviations('F_I_SZ -     symmetric', dFs%F_I_S(:,3,:), Fnew%F_I_S(:,3,:)/eta)

    ! call print_deviations('G_I_NX -     symmetric', dFs%G_I_N(:,1,:), Fnew%G_I_N(:,1,:)/eta)
    ! call print_deviations('G_I_NY -     symmetric', dFs%G_I_N(:,2,:), Fnew%G_I_N(:,2,:)/eta)
    ! call print_deviations('G_I_NZ -     symmetric', dFs%G_I_N(:,3,:), Fnew%G_I_N(:,3,:)/eta)


    Fnew    = calcpotentials(R + eta*dRa) + (-1.0d0) * F
    ! call print_deviations('F_I_I  - antisymmetric', dFa%F_I_I       , Fnew%F_I_I/eta)

    ! call print_deviations('F_I_SX - antisymmetric', dFa%F_I_S(:,1,:), Fnew%F_I_S(:,1,:)/eta)
    ! call print_deviations('F_I_SY - antisymmetric', dFa%F_I_S(:,2,:), Fnew%F_I_S(:,2,:)/eta)
    ! call print_deviations('F_I_SZ - antisymmetric', dFa%F_I_S(:,3,:), Fnew%F_I_S(:,3,:)/eta)

    ! call print_deviations('G_I_NX - antisymmetric', dFa%G_I_N(:,1,:), Fnew%G_I_N(:,1,:)/eta)
    ! call print_deviations('G_I_NY - antisymmetric', dFa%G_I_N(:,2,:), Fnew%G_I_N(:,2,:)/eta)
    ! call print_deviations('G_I_NZ - antisymmetric', dFa%G_I_N(:,3,:), Fnew%G_I_N(:,3,:)/eta)

    print *

  end subroutine test_potentials

  subroutine print_deviations(name, ref, findiff)
    !-----------------------------------------------------------------------
    ! TODO: describe
    !
    !-----------------------------------------------------------------------
    character(len=*), intent(in) :: name
    complex(KIND=dp), intent(in) :: ref(:,:), findiff(:,:)
    integer :: i

    print *, name, ' real part:'
    print *, '-------------------------'
    do i=1,nx
      print ('(i3, 2f10.3, es12.3)'), i, DBLE(ref(i,1)), DBLE(findiff(i,1)), DBLE(ref(i,1))- DBLE(findiff(i,1))
    enddo
    print *
    print *, name, ' imaginary part:'
    print *, '-------------------------'
    do i=1,nx
      print ('(i3, 2f10.3, es12.3)'), i, DBLE(ref(i,1)), DBLE(findiff(i,1)), DBLE(ref(i,1))- DBLE(findiff(i,1))
    enddo
    print *

  end subroutine print_deviations

  subroutine test_HFmatrices_me(ifail)
    !--------------------------------------------------------------------------------------
    ! Test whether the matrix elements of the single-particle Hamiltonian and the pairing 
    ! gaps in the HF basis when calculated in two different ways. For this purpose, we first
    ! calculate the densities and potentials as usual in a mean-field code, and then
    ! use 
    ! 
    ! (i)  apply_sphamil + calc_gaps: standard mean-field procedure
    !
    ! (ii) calc_sphamil_me + calc_delta_me : FAM-like calculation of matrix elements.
    !
    ! Although slightly wasteful in terms of CPU resources, this routine never assumes that
    ! any part of the matrix of the single-particle hamiltonian is hermitian/symmetric.
    !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Output:
    !  ifail : 0 if succesful, 1 if a deviation above 1e-10 has been detected
    !---------------------------------------------------------------------------------------
    integer, intent(out)       :: ifail
    complex(KIND=dp), allocatable :: sphamil_me(:,:), sphamil_orig(:,:), delta_orig(:,:), delta_me(:,:)
    complex(KIND=dp), allocatable :: hpsi(:,:), drho(:,:), dkappa(:,:), dkappa_plus(:,:), dkappa_minus(:,:)
    real(KIND=dp), allocatable    :: dev(:,:), dev_gaps(:,:)
    real(KIND=dp)                 :: stabfactor(2)
    integer                       :: si, B, N, i, it, j, N2, T
    type(PotentialVector)         :: dFs, dFa, F, dF_pp_minus, dF_pp_plus
    type(DensityVector)           :: R, dRs, dRa, R_pp_plus, R_pp_minus, dR_pp_plus, dR_pp_minus

    ifail = 0
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! (i) Ordinary mean-field-like calculation
    R  = densit(rho_can, kappa_pairing)
    F = calcpotentials(R)

    ! Important: apply_sphamil requires that the potentials in F are NOT combined!
    !            don't call combine_potentials(F) here!
    allocate(sphamil_orig(nwt,nwt)); sphamil_orig = 0.0d0
    si = 0
    do B=1,8
      N = HFBlocks(B)
      it = -1
      if(B .ge. 5) it = +1
      do j=si+1,si+N
        hpsi = apply_sphamil(HFPsi(:,:,j), HFdPsi(:,:,:,j), HFddPsi(:,:,:,j), &
        &                     sx(:,j), sy(:,j), sz(:,j), it ,.false. ,F)
        do i=si+1,si+N
          sphamil_orig(i,j) = sum(HFpsi(:,:,i) * hpsi)*dv
        enddo
      enddo
      si = si + N
    enddo

    ! Calculation of the pairing gaps
    stabfactor = 0
    call calcHFBgaps(FermiEnergy, stabfactor, F)
    delta_orig = HFBGaps

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! (ii) FAM-like calculation
    allocate(drho(nwt,nwt)) ; drho = 0.0d0
    drho = 0.0d0
    if(pairingtype.ne.0) then
      allocate(dkappa_plus(nwt,nwt))  ; dkappa_plus  = 0.0d0
      allocate(dkappa_minus(nwt,nwt)) ; dkappa_minus = 0.0d0
    endif
    ! Calculate all relevant densities and potentials WITHOUT perturbation (drho = 0), 
    ! mostly to allocate/intialize potentials that should be zero
    call densit_offdiag(drho, dkappa_plus, dkappa_minus, dRs, dRa, dR_pp_plus, dR_pp_minus)     
    call calc_perturbed_potentials(R, dRs, dRa, dR_pp_plus, dR_pp_minus, dFs, dFa, dF_pp_plus, dF_pp_minus)   

    call combine_potentials(F) ! calc_sphamil_me expects the potentials to be combined !
    ! .... and feed the result into the spwf sandwhiches
    sphamil_me = calc_sphamil_me( HFpsi, HFdpsi, HFddpsi, F, dFa, .false.)
    ! .... and add the matrix elements of the kinetic energy
    sphamil_me = sphamil_me + kinetic_me(HFpsi, HFdpsi, hfddpsi)
    ! Calculate the pairing gaps
    delta_me = calc_delta_me(HFpsi, HFdpsi, HFddpsi, F, .false.)
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! (iii)Check result and print output if needed
    si = 0
    do B = 1,8
      N = HFBlocks(B); if(N.eq.0) cycle
      ! The element-wise deviation
      dev       = abs(sphamil_orig(si+1:si+N, si+1:si+N) - sphamil_me(si+1:si+N, si+1:si+N))
      print *
      print *, 'BLOCK B=', B
      print *, "Maximal deviation of h = ", maxval(dev)
      if(maxval(dev)>1e-10) then
        ifail = 1
        print *, '------ Original calculation -------'
        do i=1,N
          print ('(99f10.3)'), sphamil_orig(si+i, si+1:si+N)
        enddo
        print *, '------ FAM-like calculation -------'
        do i=1,N
          print ('(99f10.3)'), sphamil_me(si+i, si+1:si+N)
        enddo
        print *
        print *, '------ difference           -------'
        do i=1,N
          print ('(99es10.2)'), dev(i, 1:N)
        enddo
        print *
      endif
      si = si + N
    enddo

    si = 0
    do B= 1,8,2
      N = HFBlocks(B); N2 = HFBlocks(B+1); T = N + N2; if(T.eq.0) cycle
      dev_gaps  = abs(delta_orig(si+1:si+T, si+1:si+T)   - delta_me(si+1:si+T, si+1:si+T))
      print *
      print *, 'BLOCKs B=', B, B+1
      print *, "Maximal deviation of Delta = ", maxval(dev_gaps)
      if(maxval(dev_gaps)>1e-10) then
        ifail = 1 
        print *, '------ Original calculation -------'
        do i=1,T
          print ('(99f10.3)'), DBLE(delta_orig(si+i, si+1:si+T))
        enddo
        print *, '------ FAM-like calculation -------'
        do i=1,T
          print ('(99f10.3)'), DBLE(delta_me(si+i, si+1:si+T))
        enddo
        print *
        print *, '------ difference           -------'
        do i=1,T
          print ('(99es10.2)'), DBLE(dev_gaps(i, 1:T))
        enddo
        print *
      endif
      si = si + T
    enddo
  end subroutine test_HFmatrices_me

  function kinetic_me(denpsi, dendpsi, denddpsi)
    !---------------------------------------------------------------------------------------
    ! Calculate the single-particle matrix elements (with 1-body COM included) for a
    ! set of single-particle wavefunctions on the mesh.
    !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !   denpsi  : a set of single-particle wavefunctions
    !   dendpsi : their first derivatives
    !   denddpsi: their second derivatives
    !
    ! Output:
    !   kinetic_me : the single-particle matrix elements of the kinetic energy with
    !                the effect of the one-body COM correction folded in
    !--------------------------------------------------------------------------------------
    real(KIND=dp), intent(in)         :: denpsi(:,:,:), dendpsi(:,:,:,:), denddpsi(:,:,:,:)
    complex(KIND=dp), allocatable     :: kinetic_me(:,:)

    integer                           :: it, B, si, N, wave_i, wave_j, i
    real(KIND=dp)                     :: reducedmass

    ! initialize
    allocate(kinetic_me(nwt,nwt)) ; kinetic_me = 0.0d0

    si = 0
    do B=1,8
      N = HFBlocks(B)

      !---------------------------------------------------------------------------
      ! Determine the isospin index
      if(B.ge.5) then
        it = 2
      else
        it = 1
      endif
      !---------------------------------------------------------------------------
      ! Reduced mass in case of self-consistent 1-body COM correction
      ! If doing pasta calculations, just skip.
      Reducedmass = 1.0_dp
      select case(COM1Body)
      case(0,1)
        Reducedmass = 1.0_dp
      case(2)
        Reducedmass = (1.0_dp-nucleonmass(it)/                                   &
        &                      (neutrons*nucleonmass(1)+protons*nucleonmass(2)))
      end select

      do wave_i=si+1,si+N
        do wave_j=si+1,si+N  ! Note: no assumption of hermeticity here!
           ! Action of the kinetic energy to the right
           kinetic_me(wave_i, wave_j) = kinetic_me(wave_i, wave_j) - hbm(it)* reducedmass &
           &                          * sum(denpsi(:,:,wave_j) &
           &                                *(   denddpsi(:,1,:,wave_i) &
           &                                   + denddpsi(:,4,:,wave_i) &
           &                                   + denddpsi(:,6,:,wave_i)))

          kinetic_me(wave_i, wave_j) = kinetic_me(wave_i, wave_j) * dv
        enddo
      enddo
      si = si + N
    enddo
  end function kinetic_me

  subroutine test_densit_offdiag(ifail)
    !-------------------------------------------------------------------------------
    ! This routine verifies that the density calculation of function densit
    ! and densit_offdiag result in the same densityvector; this is not a trivial
    ! observation as the former proceeds through sums in the canonical basis while
    ! the latter proceeds through a double sum in the Hartree-Fock basis.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Output:
    !  ifail : 0 if succesful, if a deviation above 1e-10 has been detected
    !-------------------------------------------------------------------------------
    integer, intent(out)          :: ifail
    type(DensityVector)           :: R_transformed, R_pp_transformed
    complex(KIND=dp), allocatable :: rho_test(:,:), kappa_test(:,:)
    real(KIND=dp), allocatable    :: transfo(:,:)
    real(KIND=dp)                 :: maxdev(3,2)
    integer                       :: i, it

    print *, '--------------- Test of densit_offdiag ---------------'
    if(pairingtype.eq.1) then
      allocate(rho_test(nwt,nwt))
      do i=1,nwt
            rho_test(i,i) = rho_can(i)
      enddo
    else
      rho_test   = rho_pairing
      kappa_test = kappa_pairing
    endif

    ! 1. Calculate the mean-field densities in the canonical basis
    Density     = densit(rho_can, kappa_pairing)
    ! 2. Generate a random unitary transformation with real coefficients
    transfo = gen_unitary_transform()
    ! 3. Transform the matrices rho and kappa, as well as the HFPsi array
    !    with this transformation
    !call mixup_rhokappa(rho_test, kappa_test, transfo)
    ! 4. Resum the densities with the offdiagonal routine
    R_transformed    = densit_offdiag_ph_symmetric(rho_test)
    R_pp_transformed = densit_offdiag_pp(kappa_test)

    ! 5. Compare
    ! TODO: make this more systematic!
    do it=1,2
        maxdev(1,it) = maxval(abs(Density%D_I_I(:,it)    - R_transformed%D_I_I  (:,it)))
        maxdev(2,it) = maxval(abs(Density%DP_I_I(:,it)   - R_pp_transformed%DP_I_I  (:,it)))

        !do i=1,nx
        !  print ('(i4,6es15.7)'), i, Density%D_I_I(i,it), R_transformed%D_I_I(i,it), &
        !  &                             Density%D_I_I(i,it) - R_transformed%D_I_I(i,it)
        !  print ('(i4,6es15.7)'), i, Density%DP_I_I(i,it), R_pp_transformed%DP_I_I(i,it), &
        !  &                             Density%DP_I_I(i,it) - R_pp_transformed%DP_I_I(i,it)
        !enddo
        !maxdev(2,it) = maxval(abs(Density%D_Nm_Nm(:,it) - R_transformed%D_Nm_Nm(:,it)))
        !maxdev(3,it) = maxval(abs(Density%C_I_Ns(:,:,:,it) - R_transformed%C_I_Ns(:,:,:,it)))
    enddo

    if(any(maxdev .gt. 1e-10)) then
        print *,'--------------------- Deviation in densities detected ---------------------'
        print *, 'Density    Neutrons      Protons'
        print ('(a4,6es15.7)'), 'RHO'       ,  maxdev(1,:), &
        &                                      DBLE(sum(Density%D_I_I(:,1))*dv), DBLE(sum(R_transformed%D_I_I(:,1)))*dv, &
        &                                      DBLE(sum(Density%D_I_I(:,2))*dv), DBLE(sum(R_transformed%D_I_I(:,2)))*dv
        print ('(a4,4es15.7)'), 'Tilde(Rho)',  maxdev(2,:)
        !print ('(a4,2es15.7)'), 'Jmn',  maxdev(3,:)
        print *,'---------------------------------------------------------------------------'
        ifail = 1
    else
        ifail = 0
        print *, 'Test passed successfully: no deviation detected.'
        print *, '---------------------------------------------------------------------------'
    endif

    ! 6. Restore the spwfs to their original condition in order to not mess with other tests
    !call mixup_rhokappa(rho_test, kappa_test, transfo)

  end subroutine test_densit_offdiag
!
!
!   subroutine build_dH_findiff(RUnper, dRs, dRa, eta)
!     !---------------------------------------------------------------------------
!     ! Build the perturbed single-particle Hamiltonian using finite difference
!     ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     ! NOT OPERATIONAL
!     !   -> use build_dH_explicit() instead.
!     !
!     ! Notes:
!     !    This function is not operational at this point but is kept for potential
!     !    test in the future. In particular, it would allow to test
!     !    calc_perturbed_potentials.
!     !---------------------------------------------------------------------------
!
!     1 format('||dH_ph|| = ', es10.3, '     ||dH_hp|| = ', es10.3)
!
!     implicit none
!     type(DensityVector), intent(in) :: RUnper, dRs, dRa
!     real(KIND=dp), intent(in)       :: eta     ! small finite diff. parameter
!     type(PotentialVector)           :: Fs, Fa
!     complex(KIND=dp), allocatable      :: HPert(:,:)
!
!     print *, "build perturbed Hamiltonian using finite difference"
!
!     allocate(HPert(nwt,nwt))
!
!     ! Calculate the total perturbed potentials, i.e. static mean-field + perturbation,
!     ! from the total perturbed densit, i.e. static mean-field + perturbation
!
!     ! call calcPotentials(RUnper + eta * dRs, eta * dRa, Fs, Fa)
!     ! => presently missing. /!\
!
!     call combine_potentials(Fs)
!
!     ! construct the sp hamiltonian
!     HPert = calc_sphamil_me( HFpsi, HFdpsi, HFddpsi, Fs,  Fa, .false.)
!
!     ! compute dH by finite difference, i.e. subtract the unperturbed Hamiltonian
!     HPert = HPert - Hunper
!     ! ... and devide by small parameter eta
!     HPert = HPert / eta
!
!     call get_ph_hp_blocks(HPert, dH(:,:,1), dH(:,:,2))
!
!     print 1, sqrt(sum( abs(dH(:,:,1))**2) ), sqrt(sum( abs(dH(:,:,2))**2) )
!
!   end subroutine build_dH_findiff



  subroutine test_gmres()
    !---------------------------------------------------------------------------
    ! Simply test for gmres: solve small linear problem
    !---------------------------------------------------------------------------

    implicit none
    integer, parameter :: n = 6
    complex(KIND=dp), dimension(n, n) :: A, Atmp
    complex(KIND=dp), dimension(n) :: b, x_explicit, x_choral
    integer :: i, info, iter, nbprod
    integer, dimension(n) :: ipiv
    real(KIND=dp) :: res

    print *, "test gmres"

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Define the matrix A
    A = reshape([ &
    dcmplx( 3.0_dp,  0.0_dp), dcmplx(-1.0_dp,  1.0_dp), dcmplx( 0.0_dp,  0.0_dp), &
    dcmplx( 0.0_dp,  0.0_dp), dcmplx( 0.0_dp,  0.0_dp), dcmplx( 0.0_dp,  0.0_dp), &
    dcmplx( 1.0_dp, -1.0_dp), dcmplx( 2.0_dp,  0.0_dp), dcmplx(-1.0_dp,  1.0_dp), &
    dcmplx( 0.0_dp,  0.0_dp), dcmplx( 0.0_dp,  0.0_dp), dcmplx( 0.0_dp,  0.0_dp), &
    dcmplx( 0.0_dp,  0.0_dp), dcmplx( 1.0_dp,  1.0_dp), dcmplx( 4.0_dp,  0.0_dp), &
    dcmplx(-1.0_dp,  1.0_dp), dcmplx( 0.0_dp,  0.0_dp), dcmplx( 0.0_dp,  0.0_dp), &
    dcmplx( 0.0_dp,  0.0_dp), dcmplx( 0.0_dp,  0.0_dp), dcmplx( 1.0_dp, -1.0_dp), &
    dcmplx( 5.0_dp,  0.0_dp), dcmplx(-1.0_dp,  1.0_dp), dcmplx( 0.0_dp,  0.0_dp), &
    dcmplx( 0.0_dp,  0.0_dp), dcmplx( 0.0_dp,  0.0_dp), dcmplx( 0.0_dp,  0.0_dp), &
    dcmplx( 1.0_dp,  1.0_dp), dcmplx( 6.0_dp,  0.0_dp), dcmplx(-1.0_dp,  1.0_dp), &
    dcmplx( 0.0_dp,  0.0_dp), dcmplx( 0.0_dp,  0.0_dp), dcmplx( 0.0_dp,  0.0_dp), &
    dcmplx( 0.0_dp,  0.0_dp), dcmplx( 1.0_dp, -1.0_dp), dcmplx( 7.0_dp,  0.0_dp)  &
    ], [n, n])

    Atmp = A

    b = [dcmplx(1.0_dp, 0.0_dp), dcmplx(1.0_dp, 0.0_dp), dcmplx(1.0_dp, 0.0_dp), &
         dcmplx(1.0_dp, 0.0_dp), dcmplx(1.0_dp, 0.0_dp), dcmplx(1.0_dp, 0.0_dp)]

    print * , "A : "
    do i=1,n
      print "(*('(', F8.5, ',', F8.5, ') ', :))",  A(:,i)
    enddo
    print * , "b : "
    do i=1,n
      print "(*('(', F8.5, ',', F8.5, ') ', :))",  b(i)
    enddo

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Compute the explicit solution x = A^{-1} * b

    x_explicit = b


    ! Call LAPACK routine ZGESV to solve the system
    call zgesv(n, 1, Atmp, n, ipiv, x_explicit, n, info)
    ! /!\ : this routine changes A 

    ! Print the results
    print *, "Explicit solution:"
    do i=1,n
      print "(*('(', F8.5, ',', F8.5, ') ', :))",  x_explicit(i)
    enddo
    print *, 'res : ', sqrt(sum(abs(b - matmul(A, x_explicit)) ** 2 ))

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Test my GMRES routine to solve Ax = b iteratively

    call alloc_gmres(multiply_by_A, b, 100, 6, 1e-6_dp, norm_2, ScalProd)

    call init_gmres(b)

    call iterate_gmres()
    call iterate_gmres()
    call iterate_gmres()
    call iterate_gmres()
    call iterate_gmres()
    call iterate_gmres()

    print *, "my GMRES solution:"
    do i=1,n
      print "(*('(', F8.5, ',', F8.5, ') ', :))",  x_gmres(i)
    enddo
    print *, 'res : ', gmres_res
    print *, 'res : ', sqrt(sum(abs(b - matmul(A, x_gmres)) ** 2 ))



  end subroutine test_gmres



  subroutine multiply_by_A(x_in, x_out)

    implicit none
    complex(KIND=dp), dimension(:), intent(in)   :: x_in
    complex(KIND=dp), dimension(:), intent(out)  :: x_out
    integer, parameter :: n = 6
    complex(KIND=dp), dimension(n, n) :: A
    integer :: i

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Define the matrix A
    A = reshape([ &
    dcmplx( 3.0_dp,  0.0_dp), dcmplx(-1.0_dp,  1.0_dp), dcmplx( 0.0_dp,  0.0_dp), &
    dcmplx( 0.0_dp,  0.0_dp), dcmplx( 0.0_dp,  0.0_dp), dcmplx( 0.0_dp,  0.0_dp), &
    dcmplx( 1.0_dp, -1.0_dp), dcmplx( 2.0_dp,  0.0_dp), dcmplx(-1.0_dp,  1.0_dp), &
    dcmplx( 0.0_dp,  0.0_dp), dcmplx( 0.0_dp,  0.0_dp), dcmplx( 0.0_dp,  0.0_dp), &
    dcmplx( 0.0_dp,  0.0_dp), dcmplx( 1.0_dp,  1.0_dp), dcmplx( 4.0_dp,  0.0_dp), &
    dcmplx(-1.0_dp,  1.0_dp), dcmplx( 0.0_dp,  0.0_dp), dcmplx( 0.0_dp,  0.0_dp), &
    dcmplx( 0.0_dp,  0.0_dp), dcmplx( 0.0_dp,  0.0_dp), dcmplx( 1.0_dp, -1.0_dp), &
    dcmplx( 5.0_dp,  0.0_dp), dcmplx(-1.0_dp,  1.0_dp), dcmplx( 0.0_dp,  0.0_dp), &
    dcmplx( 0.0_dp,  0.0_dp), dcmplx( 0.0_dp,  0.0_dp), dcmplx( 0.0_dp,  0.0_dp), &
    dcmplx( 1.0_dp,  1.0_dp), dcmplx( 6.0_dp,  0.0_dp), dcmplx(-1.0_dp,  1.0_dp), &
    dcmplx( 0.0_dp,  0.0_dp), dcmplx( 0.0_dp,  0.0_dp), dcmplx( 0.0_dp,  0.0_dp), &
    dcmplx( 0.0_dp,  0.0_dp), dcmplx( 1.0_dp, -1.0_dp), dcmplx( 7.0_dp,  0.0_dp)  &
    ], [n,n])

    x_out = matmul(A, x_in)

  end subroutine multiply_by_A


  subroutine test_gmres_affine()
    !---------------------------------------------------------------------------
    ! Simply test for gmres: look for a fixed-point of an affine transformation
    ! x -> Tx + x_free by writting it as a linear problem (1-T) x = x_free, 
    ! of the form Ax=b where A = 1-T and b = x_free. 
    !---------------------------------------------------------------------------

    implicit none
    integer, parameter :: n = 6
    complex(KIND=dp), dimension(n, n) :: T, Id, oneminusT, Atmp
    complex(KIND=dp), dimension(n) :: x_free, b, x_explicit, x_out
    integer :: i, info, iter, nbprod
    integer, dimension(n) :: ipiv
    real(KIND=dp) :: res

    print *, "test gmres"

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Define the matrix A
    T = reshape([ &
    dcmplx( 3.0_dp,  0.0_dp), dcmplx(-1.0_dp,  1.0_dp), dcmplx( 0.0_dp,  0.0_dp), &
    dcmplx( 0.0_dp,  0.0_dp), dcmplx( 0.0_dp,  0.0_dp), dcmplx( 0.0_dp,  0.0_dp), &
    dcmplx( 1.0_dp, -1.0_dp), dcmplx( 2.0_dp,  0.0_dp), dcmplx(-1.0_dp,  1.0_dp), &
    dcmplx( 0.0_dp,  0.0_dp), dcmplx( 0.0_dp,  0.0_dp), dcmplx( 0.0_dp,  0.0_dp), &
    dcmplx( 0.0_dp,  0.0_dp), dcmplx( 1.0_dp,  1.0_dp), dcmplx( 4.0_dp,  0.0_dp), &
    dcmplx(-1.0_dp,  1.0_dp), dcmplx( 0.0_dp,  0.0_dp), dcmplx( 0.0_dp,  0.0_dp), &
    dcmplx( 0.0_dp,  0.0_dp), dcmplx( 0.0_dp,  0.0_dp), dcmplx( 1.0_dp, -1.0_dp), &
    dcmplx( 5.0_dp,  0.0_dp), dcmplx(-1.0_dp,  1.0_dp), dcmplx( 0.0_dp,  0.0_dp), &
    dcmplx( 0.0_dp,  0.0_dp), dcmplx( 0.0_dp,  0.0_dp), dcmplx( 0.0_dp,  0.0_dp), &
    dcmplx( 1.0_dp,  1.0_dp), dcmplx( 6.0_dp,  0.0_dp), dcmplx(-1.0_dp,  1.0_dp), &
    dcmplx( 0.0_dp,  0.0_dp), dcmplx( 0.0_dp,  0.0_dp), dcmplx( 0.0_dp,  0.0_dp), &
    dcmplx( 0.0_dp,  0.0_dp), dcmplx( 1.0_dp, -1.0_dp), dcmplx( 7.0_dp,  0.0_dp)  &
    ], [n, n])

    print * , "T : "
    do i=1,n
      print "(*('(', F8.5, ',', F8.5, ') ', :))",  T(:,i)
    enddo

    Id = 0
    do i=1,n
      Id(i,i) = dcmplx( 1.0_dp, 0.0_dp)
    enddo

    oneminusT = Id - T

    print * , "A = 1-T : "
    do i=1,n
      print "(*('(', F8.5, ',', F8.5, ') ', :))",  oneminusT(:,i)
    enddo
  
    ! get the free part by calling the update on zero
    x_free = 0
    call affine_trafo(x_free, x_free) 

    print * , "x_free : "
    do i=1,n
      print "(*('(', F8.5, ',', F8.5, ') ', :))",  x_free(i)
    enddo

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Compute the explicit solution x = A^{-1} * b

    b = x_free

    Atmp = oneminusT

    ! Call LAPACK routine ZGESV to solve the system
    call zgesv(n, 1, Atmp, n, ipiv, b, n, info)
    ! /!\ : this routine changes A (-> triag matrix) and b (-> solution)
    x_explicit = b

    ! Print the results
    print *, "Explicit solution:"
    do i=1,n
      print "(*('(', F8.5, ',', F8.5, ') ', :))",  x_explicit(i)
    enddo

    ! test that the solution is a solution of Ax=b
    print *, '||b-Ax|| = || x_free - (1-T) x|| : ', sqrt(sum(abs(x_free - matmul(oneminusT, x_explicit)) ** 2 ))
    
    ! test that the solution is a fixed point of the affine transformation 
    call affine_trafo(x_explicit, x_out)

    print *, '||x - Tx + x_free|| : ', sqrt(sum(abs(x_explicit - x_out) ** 2 ))



    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Test my GMRES routine to solve Ax = b iteratively

    call alloc_gmres(affine_as_linear, x_free, 100, 6, 1e-6_dp, norm_2, ScalProd)

    call init_gmres(x_free)

    do i=1,6
      call iterate_gmres()
    enddo

    print *, "my GMRES solution:"
    do i=1,n
      print "(*('(', F8.5, ',', F8.5, ') ', :))",  x_gmres(i)
    enddo
    print *, 'GMRES res : ', gmres_res
    print *, '|| x_free - (1-T) x|| : ', sqrt(sum(abs(x_free - matmul(oneminusT, x_gmres)) ** 2 ))

    ! test that the solution is a fixed point of the affine transformation 
    call affine_trafo(x_gmres, x_out)

    print *, '||x - Tx + x_free|| : ', sqrt(sum(abs(x_gmres - x_out) ** 2 ))


    call dealloc_gmres()


  end subroutine test_gmres_affine

  subroutine affine_trafo(x_in, x_out)

    implicit none
    complex(KIND=dp), dimension(:), intent(in)   :: x_in
    complex(KIND=dp), dimension(:), intent(out)  :: x_out
    integer, parameter :: n = 6
    complex(KIND=dp), dimension(n, n) :: T
    complex(KIND=dp), dimension(n) :: x_free
    integer :: i

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Define the matrix A
    T = reshape([ &
    dcmplx( 3.0_dp,  0.0_dp), dcmplx(-1.0_dp,  1.0_dp), dcmplx( 0.0_dp,  0.0_dp), & 
    dcmplx( 0.0_dp,  0.0_dp), dcmplx( 0.0_dp,  0.0_dp), dcmplx( 0.0_dp,  0.0_dp), &
    dcmplx( 1.0_dp, -1.0_dp), dcmplx( 2.0_dp,  0.0_dp), dcmplx(-1.0_dp,  1.0_dp), & 
    dcmplx( 0.0_dp,  0.0_dp), dcmplx( 0.0_dp,  0.0_dp), dcmplx( 0.0_dp,  0.0_dp), &
    dcmplx( 0.0_dp,  0.0_dp), dcmplx( 1.0_dp,  1.0_dp), dcmplx( 4.0_dp,  0.0_dp), & 
    dcmplx(-1.0_dp,  1.0_dp), dcmplx( 0.0_dp,  0.0_dp), dcmplx( 0.0_dp,  0.0_dp), &
    dcmplx( 0.0_dp,  0.0_dp), dcmplx( 0.0_dp,  0.0_dp), dcmplx( 1.0_dp, -1.0_dp), & 
    dcmplx( 5.0_dp,  0.0_dp), dcmplx(-1.0_dp,  1.0_dp), dcmplx( 0.0_dp,  0.0_dp), &
    dcmplx( 0.0_dp,  0.0_dp), dcmplx( 0.0_dp,  0.0_dp), dcmplx( 0.0_dp,  0.0_dp), & 
    dcmplx( 1.0_dp,  1.0_dp), dcmplx( 6.0_dp,  0.0_dp), dcmplx(-1.0_dp,  1.0_dp), &
    dcmplx( 0.0_dp,  0.0_dp), dcmplx( 0.0_dp,  0.0_dp), dcmplx( 0.0_dp,  0.0_dp), & 
    dcmplx( 0.0_dp,  0.0_dp), dcmplx( 1.0_dp, -1.0_dp), dcmplx( 7.0_dp,  0.0_dp)  &
    ], [n,n])

    
    x_free = [dcmplx(1.0_dp, 0.0_dp), dcmplx(1.0_dp, 1.0_dp), dcmplx(2.0_dp, 0.0_dp), &
              dcmplx(6.0_dp, -1.0_dp), dcmplx(1.0_dp, 0.0_dp), dcmplx(3.0_dp, -2.0_dp)]

    x_out = matmul(T, x_in) + x_free

  end subroutine affine_trafo

  subroutine affine_as_linear(x_in, x_out)

    implicit none
    complex(KIND=dp), dimension(:), intent(in)   :: x_in
    complex(KIND=dp), dimension(:), intent(out)  :: x_out
    complex(KIND=dp), dimension(6) :: x_free

    ! Ax = (1-T)x = x - (Tx + x_free) + x_free = x - affine(x) + x_free

    x_free = [dcmplx(1.0_dp, 0.0_dp), dcmplx(1.0_dp, 1.0_dp), dcmplx(2.0_dp, 0.0_dp), &
              dcmplx(6.0_dp, -1.0_dp), dcmplx(1.0_dp, 0.0_dp), dcmplx(3.0_dp, -2.0_dp)]

    call affine_trafo(x_in, x_out)

    x_out = x_in - x_out + x_free

    print * , "||x_in||",   norm_dH(x_in)
    print * , "||x_out||",   norm_dH(x_out)
    print * , "||x_free||",   norm_dH(x_free)

  end subroutine affine_as_linear

 function norm_2(vec) result(norm)
    complex(KIND=dp), dimension(:), intent(in)  :: vec
    real(KIND=dp)                               :: norm

    norm = sqrt(sum(abs(vec(:)) ** 2))

  end function

  function ScalProd(vec_l, vec_r) result(res)
    ! Innner product on complex vector space.
    ! NOTE : we follow the maths convention where the inner product is linear 
    !        in the first component, while the typical physics convention assumes 
    !        linearity in the second component. 
    complex(KIND=dp), dimension(:), intent(in)  :: vec_l, vec_r
    complex(KIND=dp)                            :: res

    res = sum(vec_l(:) * conjg(vec_r(:)) )

  end function


  subroutine test_qptrafo()
    !---------------------------------------------------------------------------
    ! test the new qp trafo routines
    !---------------------------------------------------------------------------

   
    complex(KIND=dp), allocatable :: identity(:,:)
    complex(KIND=dp), allocatable :: Rsp(:,:,:), Rspback(:,:,:), Rqp(:,:,:),  Rph(:,:,:)

    real(KIND=dp), allocatable :: N20(:,:), N11(:,:)
    complex(KIND=dp), allocatable ::  N20new(:,:), N11new(:,:)

    integer :: i, si, T
      

    if(pairingtype==0) then
      print *, 'PairingType = 0 => HF :: Cannot perform QP trafo test. Printing ph matrix elements'
      
      allocate(Rsp(nwt,nwt,3))     ! contains the spme of a one-body operator F: f20_pq, f11_pq, f02_pq
  
      ! define Rsp as the spme of r^2
      Rsp(:,:,1) = 0
      Rsp(:,:,2) = Rsq_spme()
      Rsp(:,:,3) = 0

      allocate(Rph(nwt,nwt,2))     ! contains the qpme of a one-body operator F: F20_k1k2 and F11_k1k2
      Rph = 0

      call get_ph_hp_blocks(Rsp(:,:,2), Rph(:,:,1), Rph(:,:,2))
      
      print *, '||Fph (proton)||', sum(abs(Rph(nwn+1:nwt, nwn+1:nwt,1) * Rph(nwn+1:nwt, nwn+1:nwt,1)))
      ! call print_spme_complex( Rph(:,:,1))
      print *, 'Fph_ia : BLOCK 5 & 6'
      T = HFblocks(5)+HFblocks(6)
      do i=nwn+1,nwn+T
        print "(*( '(',g12.5,',',g12.5,')',:))",  Rph(i,nwn+1:nwn+T,1)
      enddo

      print '(99f10.5)', rho_can(nwn+1:nwn+T)

      return
    endif


    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! 1. some test for the particle number operator, comparing to existing routines

    print * , 'TEST 1 : comparing to exising routines for the particle number operator'


    ! construct the spme of the particle-number operator = identity 
    allocate(identity(nwt,nwt)) 
    identity = 0.0

    ! Set the diagonal elements to 1
    do i = 1, nwt
      identity(i, i) = 1.0
    end do


    ! transform to qpme N20
    allocate(N20new(nwt,nwt)) 
    allocate(N11new(nwt,nwt)) 
    call transform_sp_to_qp(Bogoliubov, O11sp=identity, O20qp=N20new, O11qp=N11new)

    ! get qpme of N20 from existing routine, for the neutron channel
    N20 = calcN20(Bogoliubov, HFblocks(1:4))

    ! get qpme of N11 from existing routine, for the neutron channel
    N11 = calcN11(Bogoliubov, HFblocks(1:4))

    print *, 'N20 (new routine) : BLOCK 1 & 2'
    T = HFblocks(1) + HFblocks(2) ! size of block1 + block2
    do  i=1,T
      print "(*( '(',f10.5,',',f10.5,')',:))",  N20new(i, 1:T)
    enddo


    print *, 'N20 (existing routine) : BLOCK 1 & 2'
    do i=1,T
      print "(99f10.5)",  N20(i, 1:T)
    enddo

    print *, 'diff '
    T = HFblocks(1) + HFblocks(2) ! size of block1 + block2
    do  i=1,T
      print "(*( '(',f10.5,',',f10.5,')',:))",  N20new(i, 1:T) - N20(i, 1:T)
    enddo

    
    print *, ' ||N20(existing routine) - N20(new)|| = ', &
    & sum(abs(N20(:,:) - N20new(1:nwn,1:nwn)))



    print *, 'N11 (new routine) : BLOCK 1 & 2'
    T = HFblocks(1) + HFblocks(2) ! size of block1 + block2
    do  i=1,T
      print "(*( '(',f10.5,',',f10.5,')',:))",  N11new(i, 1:T)
    enddo


    print *, 'N11 (existing routine) : BLOCK 1 & 2'
    do i=1,T
      print "(99f10.5)",  N11(i, 1:T)
    enddo

    print *, 'diff '
    T = HFblocks(1) + HFblocks(2) ! size of block1 + block2
    do  i=1,T
      print "(*( '(',f10.5,',',f10.5,')',:))",  N11new(i, 1:T) - N11(i, 1:T)
    enddo


    
    print *, ' ||N11(existing routine) - N11(new)|| = ', &
    & sum(abs(N11(:,:) - N11new(1:nwn,1:nwn)))



    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! 2. test unitarity by going from forth and back Bogo trafo: sp -> qp -> sp
    ! for the operator Rsq 

    print * , 'TEST 2 : perform forth and back Bogo trafo for Rsq: spme -> qpme -> spme'

    allocate(Rsp(nwt,nwt,3))     ! contains the spme of a one-body operator F: f20_pq, f11_pq, f02_pq
    allocate(Rqp(nwt,nwt,3))     ! contains the qpme of a one-body operator F: F20_k1k2 and F11_k1k2
    allocate(Rspback(nwt,nwt,3)) ! contains the spme of a one-body operator F: f20_pq, f11_pq, f02_pq

    ! define Rsp as the spme of r^2
    Rsp(:,:,1) = 0
    Rsp(:,:,2) = Rsq_spme()
    Rsp(:,:,3) = 0

    ! call transform_sp_to_qp(Bogoliubov, Rsp(:,:,1), Rsp(:,:,2), Rsp(:,:,3), Rqp(:,:,1),  Rqp(:,:,2),  Rqp(:,:,3))
    call transform_sp_to_qp(Bogoliubov, O11sp=Rsp(:,:,2), O20qp=Rqp(:,:,1), O11qp=Rqp(:,:,2), O02qp=Rqp(:,:,3))

    ! get qpme of F 
    ! no need to name the optional arguments if all are used 
    call transform_qp_to_sp(Bogoliubov, Rqp(:,:,1),  Rqp(:,:,2),  Rqp(:,:,3), Rspback(:,:,1), Rspback(:,:,2), Rspback(:,:,3))


    print *, 'F20 (sp) : BLOCK 1 & 2'
    T = HFblocks(1) + HFblocks(2) ! size of block1 + block2
    do  i=1,T
      print "(*( '(',f10.5,',',f10.5,')',:))",  Rsp(i,1:T, 1)
    enddo

    print *, 'F20 (qp) : BLOCK 1 & 2'
    T = HFblocks(1) + HFblocks(2) ! size of block1 + block2
    do  i=1,T
      print "(*( '(',f10.5,',',f10.5,')',:))",  Rqp(i, 1:T, 1)
    enddo

   
    print *, 'F20 (sp->qp->sp) : BLOCK 1 & 2'
    T = HFblocks(1) + HFblocks(2) ! size of block1 + block2
    do  i=1,T
      print "(*( '(',f10.5,',',f10.5,')',:))",  Rspback(i, 1:T, 1)
    enddo
    print *, ' ||F20 (sp) - F20 (sp->qp->sp)|| = ',  sum(abs(Rsp(:,:,1) - Rspback(:,:,1)))


    print *, 'F11 (sp) : BLOCK 1 & 2'
    T = HFblocks(1) + HFblocks(2) ! size of block1 + block2
    do  i=1,T
      print "(*( '(',f10.5,',',f10.5,')',:))",  Rsp(i,1:T, 2)
    enddo

    print *, 'F11 (qp) : BLOCK 1 & 2'
    T = HFblocks(1) + HFblocks(2) ! size of block1 + block2
    do  i=1,T
      print "(*( '(',f10.5,',',f10.5,')',:))",  Rqp(i, 1:T, 2)
    enddo

   
    print *, 'F11 (sp->qp->sp) : BLOCK 1 & 2'
    T = HFblocks(1) + HFblocks(2) ! size of block1 + block2
    do  i=1,T
      print "(*( '(',f10.5,',',f10.5,')',:))",  Rspback(i, 1:T, 2)
    enddo
    print *, ' ||F11 (sp) - F11 (sp->qp->sp)|| = ',  sum(abs(Rsp(:,:,2) - Rspback(:,:,2)))


    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! 3. test symmetry properties of the qpme

    print * , 'TEST 3 : perform symmetry checks on Rsq qpme'
    
    print *, 'Antisymmetry : '
    print *, '    R20_ab = -R20_ba   : satisfied up to',  sum(abs(Rqp(:,:,1) + transpose(Rqp(:,:,1))))
    print *, '    R02_ab = -R02_ba   : satisfied up to',  sum(abs(Rqp(:,:,3) + transpose(Rqp(:,:,3))))
    print *, 'Hermiticity : '
    print *, '    R11_ab =  R11_ba^* : satisfied up to',  sum(abs(Rqp(:,:,2) - conjg(transpose(Rqp(:,:,2)))))
    print *, '    R20_ab =  R02_ab^* : satisfied up to',  sum(abs(Rqp(:,:,1) - conjg(Rqp(:,:,3))))



    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! *. test whether F20 reduced to ph elements when the bogoliubov trafo is trivial
    ! Note that this is a bit tricky since the previous routines used in absence of pairing 
    ! get the occupations from rho_can which is diag(rho_hf) in that case but this in no longer
    ! true here since even in a trivial Bogo, indexing of the states can be altered

    print * , 'TEST * : check if zero-pairing limit of qp trafo: isolating to (ph,hp,pp,hh) blocks '

    print *, '||f11 (proton)||', sum(abs(Rsp(nwn+1:nwt, nwn+1:nwt,2) * Rsp(nwn+1:nwt, nwn+1:nwt,2)))

    print *, 'f11_pq : BLOCK 5 & 6'
    si =  nwn ! start index of block 5
    T = HFblocks(5) + HFblocks(6) ! size of block5 + block6

    do i=si+1,si+T
      print "(*( '(',g12.5,',',g12.5,')',:))",  Rsp(i,si+1:si+T,2)
    enddo


    
    print *, '||F20 (proton)||', sum(abs(Rqp(nwn+1:nwt, nwn+1:nwt,1) * Rqp(nwn+1:nwt, nwn+1:nwt,1)))

    print *, 'F20_k1k2 : BLOCK 5 & 6'
    do i=si+1,si+T
      print "(*( '(',g12.5,',',g12.5,')',:))",  Rqp(i,si+1:si+T,1)
    enddo


    allocate(Rph(nwt,nwt,2))     ! contains the qpme of a one-body operator F: F20_k1k2 and F11_k1k2
    Rph = 0


    ! assuming the bogo trafo is trivial in the proton block, the same result should be recovered from 
    ! the existing particle hole getters. 
    call get_ph_hp_blocks(Rsp(:,:,2), Rph(:,:,1), Rph(:,:,2))


    print *, '||Fph (proton)||', sum(abs(Rph(nwn+1:nwt, nwn+1:nwt,1) * Rph(nwn+1:nwt, nwn+1:nwt,1)))
    
    ! call print_spme_complex( Rph(:,:,1))
    
    print *, 'Fph_ia : BLOCK 5 & 6'
    do i=si+1,si+T
      print "(*( '(',g12.5,',',g12.5,')',:))",  Rph(i,si+1:si+T,1)
    enddo
    
    stop

  end subroutine test_qptrafo


  subroutine test_linearity_T()
    ! Test a complete FAM iteration is an affine transformation by calling iterate_dH on 
    ! a chosen linear combination a * dHa + b*dHb. One expects that
    ! FAM(a * dHa + b *dHb) = T(a * dHa + b * dHb) + dH_free
    !                       = a * T(dHa) + b * T(dHb) + dH_free
    ! such that
    ! FAM(a * dHa + b *dHb) - dH_free = a * (FAM(dHa) - dH_free) + b * (FAM(dHb) - dH_free)

    implicit none
    complex(KIND=dp), allocatable :: dHa(:), dHa_iter(:)
    complex(KIND=dp), allocatable :: dHb(:), dHb_iter(:)
    complex(KIND=dp), allocatable :: dHlincomb(:), dHlincomb_iter(:)
    real(KIND=dp),    allocatable :: rand_real(:), rand_imag(:)
    complex(KIND=dp), allocatable :: dH_free(:)
    complex(KIND=dp)              :: a, b
    integer                       :: i,j
    real                          :: diff_from_lin


    allocate(dHa(nwt*nwt))
    allocate(dHa_iter(nwt*nwt))
    allocate(dHb(nwt*nwt))
    allocate(dHb_iter(nwt*nwt))
    allocate(dHlincomb(nwt*nwt))
    allocate(dHlincomb_iter(nwt*nwt))
    allocate(rand_real(nwt*nwt))
    allocate(rand_imag(nwt*nwt))
    allocate(dH_free(nwt*nwt))


    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! get dH_free

    dH_free = 0
    call iterate_dHsp(dH_free, dH_free)

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! define dHa as random dcmplx matrix

    ! Seed the random number generator
    call random_seed()

    ! Generate random real and imaginary parts
    call random_number(rand_real)
    call random_number(rand_imag)

    ! Combine into complex matrix
    dHa = 0
    do i = 1, nwt
        do j = 1, nwt
            dHa( (j-1) * nwt + i ) = dcmplx(rand_real((j-1) * nwt + i), rand_imag((j-1) * nwt + i))
        end do
    end do

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! define dHb as random dcmplx matrix

    ! Generate random real and imaginary parts
    call random_number(rand_real)
    call random_number(rand_imag)

    ! Combine into complex matrix
    dHb = 0
    do i = 1, nwt
        do j = 1, nwt
            dHb( (j-1) * nwt + i ) = dcmplx(rand_real((j-1) * nwt + i), rand_imag((j-1) * nwt + i))
        end do
    end do


    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! dHlincomb as some linear combination of dHa and dHb
    
    a = dcmplx(6.0_dp, -1.0_dp) ! some random C scalar
    b = dcmplx(-2.0_dp, 0.5_dp) ! some random C scalar
    
    dHlincomb = a * dHa + b * dHb

    print * , "||dHa|| = ", norm_2(dHa)
    print * , "||dHb|| = ", norm_2(dHb)
    print * , "||dHlincomb|| = ", norm_2(dHlincomb)

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! calculate  FAM(dHa), FAM(dHb) and FAM(a * dHa + b *dHb)

    call iterate_dHsp(dHa, dHa_iter)
    call iterate_dHsp(dHb, dHb_iter)
    call iterate_dHsp(dHlincomb, dHlincomb_iter)

    print * , "||FAM(dHa)|| = ", norm_2(dHa_iter)
    print * , "||FAM(dHb)|| = ", norm_2(dHb_iter)
    print * , "||FAM(dHlincomb)|| = ", norm_2(dHlincomb_iter)

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! subtract the free response 

    dHa_iter = dHa_iter - dH_free
    dHb_iter = dHb_iter - dH_free
    dHlincomb_iter = dHlincomb_iter - dH_free

    print * , "||FAM(dHlincomb) - free|| = ", norm_2(dHlincomb_iter)
    print * , "||a(FAM(dHa) - free) + b(FAM(dHb) - free)|| = ", norm_2(a * dHa_iter + b * dHb_iter)


    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! verify linarity of T and hence affinity of FAM

    diff_from_lin = norm_2(dHlincomb_iter - (a * dHa_iter + b * dHb_iter))
    print * , "|| [FAM(dHlincomb) - free] - [a(FAM(dHa) - free) + b(FAM(dHb) - free)]|| = ", diff_from_lin

    print *, "the FAM iteration iterate_dHsp is an affine transformation upto a precision of ", diff_from_lin


  end subroutine

  subroutine test_linearity_FAM_coulomb()

    implicit none

    complex(KIND=dp), allocatable :: dHa(:), dHa_iter(:), dH_free(:)
    complex(KIND=dp), allocatable :: dHb(:), dHb_iter(:)
    complex(KIND=dp), allocatable :: dHlincomb(:), dHlincomb_iter(:)
    real(KIND=dp),    allocatable :: rand_real(:), rand_imag(:)
    type(potentialvector)         :: dFsym_a, dFsym_b, dFanti_a, dFanti_b, dFsym_lc, dFanti_lc, dFsym_free, dFanti_free
    type(densityvector)           :: dRsym_a, dRsym_b, dRanti_a, dRanti_b, dRsym_lc, dRanti_lc, dRsym_free, dRanti_free
    complex(KIND=dp)              :: a, b
    integer                       :: i,j
    real(KIND=dp)                 :: diff_from_lin_direct, diff_from_lin_exch, diff_from_lin


    allocate(dHa(nwt*nwt))
    allocate(dHa_iter(nwt*nwt))
    allocate(dHb(nwt*nwt))
    allocate(dHb_iter(nwt*nwt))
    allocate(dHlincomb(nwt*nwt))
    allocate(dHlincomb_iter(nwt*nwt))
    allocate(rand_real(nwt*nwt))
    allocate(rand_imag(nwt*nwt))
    allocate(dH_free(nwt*nwt))

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! get dH_free
    dH_free = 0
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! define dHa as random dcmplx matrix

    ! Seed the random number generator
    call random_seed()

    ! Generate random real and imaginary parts
    call random_number(rand_real)
    call random_number(rand_imag)

    ! Combine into complex matrix
    dHa = 0
    do i = 1, nwt
        do j = 1, nwt
            dHa( (j-1) * nwt + i ) = dcmplx(rand_real((j-1) * nwt + i), rand_imag((j-1) * nwt + i))
        end do
    end do

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! define dHb as random dcmplx matrix

    ! Generate random real and imaginary parts
    call random_number(rand_real)
    call random_number(rand_imag)

    ! Combine into complex matrix
    dHb = 0
    do i = 1, nwt
        do j = 1, nwt
            dHb( (j-1) * nwt + i ) = dcmplx(rand_real((j-1) * nwt + i), rand_imag((j-1) * nwt + i))
        end do
    end do

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! dHlincomb as some linear combination of dHa and dHb
    a = dcmplx(0.0_dp, 1.0_dp) ! some random C scalar
    b = dcmplx(1.0_dp, 0.0_dp) ! some random C scalar

    dHlincomb = a * dHa + b * dHb

    print * , "||dHa|| = ", norm_2(dHa)
    print * , "||dHb|| = ", norm_2(dHb)
    print * , "||dHlincomb|| = ", norm_2(dHlincomb)
    print * , "||CHECK|| = ", norm_2(dHlincomb - a * dHa - b * dHb)

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! calculate  FAM(dHa), FAM(dHb) and FAM(a * dHa + b *dHb)
    call iterate_partial_dHsp(dH_free, dH_free         , dRsym_free, dRanti_free, dFsym_free, dFanti_free)
    call iterate_partial_dHsp(dHa, dHa_iter            , dRsym_a , dRanti_a , dFsym_a , dFanti_a)
    call iterate_partial_dHsp(dHb, dHb_iter            , dRsym_b , dRanti_b , dFsym_b , dFanti_b)
    call iterate_partial_dHsp(dHlincomb, dHlincomb_iter, dRsym_lc, dRanti_lc, dFsym_lc, dFanti_lc)
!
!     print * , "||FAM(dHa)|| = ", norm_2(dHa_iter)
!     print * , "||FAM(dHb)|| = ", norm_2(dHb_iter)
!     print * , "||FAM(dHlincomb)|| = ", norm_2(dHlincomb_iter)

    dRsym_a  = dRsym_a    + (-1.0d0)* dRsym_free
    dRsym_b  = dRsym_b    + (-1.0d0)* dRsym_free
    dRsym_lc = dRsym_lc   + (-1.0d0)* dRsym_free
    
    dRanti_a  = dRanti_a  + (-1.0d0)* dRanti_free
    dRanti_b  = dRanti_b  + (-1.0d0)* dRanti_free
    dRanti_lc = dRanti_lc + (-1.0d0)* dRanti_free

    dFsym_a  = dFsym_a    + (-1.0d0)* dFsym_free
    dFsym_b  = dFsym_b    + (-1.0d0)* dFsym_free
    dFsym_lc = dFsym_lc   + (-1.0d0)* dFsym_free

    dFanti_a  = dFanti_a  + (-1.0d0)* dFanti_free
    dFanti_b  = dFanti_b  + (-1.0d0)* dFanti_free
    dFanti_lc = dFanti_lc + (-1.0d0)* dFanti_free

    
    ! symmetric parts
    print *, ' PROTON  DENSITIES'
    diff_from_lin = sum( ABS (a * dRsym_a%D_I_I(:,2) + b* dRsym_b%D_I_I(:,2) - dRsym_lc%D_I_I(:,2))**2 )
    print *, 'SYMMETRIC, REAL', diff_from_lin
    diff_from_lin = sum( ABS (a * dRanti_a%D_I_I(:,2) + b* dRanti_b%D_I_I(:,2) - dRanti_lc%D_I_I(:,2))**2 )
    print *, 'ANTISYMMETRIC, REAL', diff_from_lin

    print *, ' CHARGE DENSITIES'
    diff_from_lin = sum( DBLE (a * dRsym_a%chargedensity + b* dRsym_b%chargedensity - dRsym_lc%chargedensity)**2 )
    print *, 'SYMMETRIC, REAL', diff_from_lin
    diff_from_lin = sum( AIMAG (a * dRsym_a%chargedensity + b* dRsym_b%chargedensity - dRsym_lc%chargedensity)**2 )
    print *, 'SYMMETRIC, IMAG', diff_from_lin
    diff_from_lin = sum( abs (a * dRanti_a%chargedensity + b* dRanti_b%chargedensity - dRanti_lc%chargedensity)**2 )
    print *, 'ANTISYMMETRIC', diff_from_lin
    print *, 'POTENTIALS'
    diff_from_lin_direct = &
    & sum( abs (a * dFsym_a%Coulombpotential  + b* dFsym_b%Coulombpotential  - dFsym_lc%Coulombpotential)**2 )
    diff_from_lin_exch   = &
    & sum( abs (a * dFsym_a%Exchangepotential + b* dFsym_b%Exchangepotential - dFsym_lc%Exchangepotential)**2 )
     print *, 'SYMMETRIC, DIRECT ', diff_from_lin_direct
     print *, 'SYMMETRIC, EXCHANGE', diff_from_lin_exch
    diff_from_lin_direct = &
    & sum( abs (a * dFanti_a%Coulombpotential  + b* dFanti_b%Coulombpotential  - dFanti_lc%Coulombpotential)**2 )
    diff_from_lin_exch   = &
    & sum( abs (a * dFanti_a%Exchangepotential + b* dFanti_b%Exchangepotential - dFanti_lc%Exchangepotential)**2 )
     print *, 'ANTISYMMETRIC, DIRECT ', diff_from_lin_direct
     print *, 'ANTISYMMETRIC, EXCHANGE', diff_from_lin_exch

     ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     ! subtract the free response

     dHa_iter = dHa_iter - dH_free
     dHb_iter = dHb_iter - dH_free
     dHlincomb_iter = dHlincomb_iter - dH_free

     print * , "||FAM(dHlincomb) - free|| = ", norm_2(dHlincomb_iter)
     print * , "||a(FAM(dHa) - free) + b(FAM(dHb) - free)|| = ", norm_2(a * dHa_iter + b * dHb_iter)


    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! verify linarity of T and hence affinity of FAM

     diff_from_lin = norm_2(dHlincomb_iter - (a * dHa_iter + b * dHb_iter))
     print * , "|| [FAM(dHlincomb) - free] - [a(FAM(dHa) - free) + b(FAM(dHb) - free)]|| = ", diff_from_lin

     print *, "the FAM iteration iterate_dHsp is an affine transformation upto a precision of ", diff_from_lin
     stop

  end subroutine

  subroutine iterate_partial_dHsp(dHsp_flat, dHspout_flat, dRs, dRa, dFs, dFa)
    !---------------------------------------------------------------------------
    ! Perform one FAM loop of the perturbed single-particle hamiltonian dH
    ! (in HF basis), which contain dh and ddelta (in the QFAM).
    !---------------------------------------------------------------------------
    1 format('||dH_ph|| = ', es10.3, '     ||dH_hp|| = ', es10.3)
    2 format('||X|| = ', es10.3, '     ||Y|| = ', es10.3)
    3 format(' S_',i1,i1,' (', f5.2, ') = ', es10.3)

    implicit none
    complex(KIND=dp), dimension(:), target, intent(in)   :: dHsp_flat
    complex(KIND=dp), dimension(:), target, intent(out)  :: dHspout_flat

    complex(KIND=dp), pointer :: dHsp(:,:), dHspout(:,:)


    type(DensityVector), intent(out):: dRs, dRa
    type(DensityVector)             :: dR_pp_plus, dR_pp_minus ! temporary placeholders for this particular HF routine
    type(PotentialVector), intent(out):: dFs, dFa
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

    !call store_XY_hist()

    ! build the perturbed densities on the mesh dRs, dRa from X and Y
    call build_perturbed_densities(X, Y, dRs, dRa, dR_pp_plus, dR_pp_minus)

    ! explicit linearisation of the fields
    call calc_perturbed_potentials(RUnper, dRs, dRa, dR_pp_plus, dR_pp_minus, dFs, dFa, dF_pp_minus, dF_pp_plus)

    ! We add in all additional contributions to F_I_I that do not
    !  result from the Skyrme functional.
    call combine_potentials(dFs)
    call combine_potentials(dFa)

    ! construct the sp hamiltonian
    dHspout = calc_sphamil_me( HFpsi, HFdpsi, HFddpsi,dFs, dFa, .false.)

  end subroutine iterate_partial_dHsp


end module fam_testing
