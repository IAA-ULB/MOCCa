module fam_testing
  !--------------------------------------------------------------------------------
  ! A selection of tests to check the development of the FAM functionality and
  ! guard it against regressions.
  !--------------------------------------------------------------------------------
  use densities
  use moments
  use Coulombmod, only : SolveCoulomb
  use pairing,    only :  rho_can, pairingtype, rho_pairing, kappa_pairing
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

  !     call test_sphamil_me(ifail)
  !     print 1, 'SPHAMIL_ME', ifail

    call test_potentials(X,Y,ifail)
    print 1, 'potentials', ifail

    !call test_potentials(X,Y,ifail)
    ! Attention: this testing routine has serious side effects on the state of the program.
    !call test_densit_offdiag(ifail)
    !print 1, 'DENSIT_OFFDIAG', ifail

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

    type(DensityVector)           :: R, dRa, dRs
    type(PotentialVector)         :: F, dFs, dFa, Fnew

    integer :: i


    allocate(drho(nwt,nwt))
    drho = 0.0d0
    drho = X + transpose(Y)

    R       = densit(rho_can, kappa_pairing)
    dRs     = densit_offdiag_symmetric(drho, dkappa)
    dRa     = densit_offdiag_antisymmetric(drho,dkappa)

    F       = calcpotentials(R)
    call calc_perturbed_potentials(R,dRs,dRa, dFs, dFa)

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

  subroutine test_sphamil_me(ifail)
    !--------------------------------------------------------------------------------------
    ! Test whether the matrix elements of the single-particle Hamiltonian in the HF
    ! basis when calculated in two different ways.
    !
    ! (i)  densit + calc_potentials + apply_sphamil
    !      - - - - - - - - - - - - - - - - -
    !      the densities and potentials calculated as usual; with the latter
    !      applied to the spwfs as usual in the mean-field part of the code
    !
    ! (ii) densit_offdiag + calc_perturbed_potentials + calc_sphamil_me
    !      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !      a. the (offdiagonal) summation of densities with the mean-field density matrix
    !      b. the calculation of the potentials with calc_perturbed_potentials
    !         without perturbation
    !      c. the calculation of the matrix elements by explicit sandwiching
    !         of the potentials in calc_sphamil_me
    !      d. ... with the matrix elements of the kinetic energy added in manually!
    !
    ! Although slightly wasteful in terms of CPU resources, this routine never assumes that
    ! any part of the matrix of the single-particle hamiltonian is hermitian/symmetric.
    !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Output:
    !  ifail : 0 if succesful, 1 if a deviation above 1e-10 has been detected
    !---------------------------------------------------------------------------------------
    integer, intent(out)       :: ifail

    complex(KIND=dp), allocatable :: sphamil_me(:,:), sphamil_orig(:,:)
    complex(KIND=dp), allocatable :: hpsi(:,:), drho(:,:), dkappa(:,:)
    real(KIND=dp), allocatable    :: dev(:,:)
    integer                       :: si, B, N, i, it, j
    type(PotentialVector)         :: Fs, Fa, F
    type(DensityVector)           :: R, Rs, Ra

    ifail = 0
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! (i) Ordinary mean-field-like calculation
    R  = densit(rho_can, kappa_pairing)
    F = calcpotentials(R)

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

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! (ii) FAM-like calculation
    allocate(drho(nwt,nwt)) ; drho = 0.0d0
    ! Calculate the original densities through a non-diagonal summation
    do i=1,nwt
      drho(i,i) = rho_can(i)
    enddo
    call densit_offdiag(drho, dkappa, R , Ra)
    drho = 0.0d0
    call densit_offdiag(drho, dkappa, Rs, Ra)     ! No perturbation, drho = 0 in this call
    ! Calculate the potentials without perturbation
    call calc_perturbed_potentials(R, Rs, Ra, Fs, Fa)
    !- - - - - - - - - - - - - - - - -
    ! Convention for calc_sphamil_me !
    call combine_potentials(F)
    call combine_potentials(Fs)
    call combine_potentials(Fa)
    ! .... and feed the result into the spwf sandwhiches
    sphamil_me = calc_sphamil_me( HFpsi, HFdpsi, HFddpsi,F, Fa, .false.)
    ! .... and add the matrix elements of the kinetic energy
    sphamil_me = sphamil_me + kinetic_me(HFpsi, HFdpsi, hfddpsi)

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! (iii)Check result and print output if needed
    si = 0
    do B = 1,8
      N = HFBlocks(B)
      ! The element-wise deviation
      dev = abs(sphamil_orig(si+1:si+N, si+1:si+N) - sphamil_me(si+1:si+N, si+1:si+N))
      if(maxval(dev)>1e-10) then
        ifail = 1
        print *
        print *, 'BLOCK B=', B
        print *, "Maximal deviation = ", maxval(dev)
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
  end subroutine test_sphamil_me

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
    type(DensityVector)           :: R_transformed
    complex(KIND=dp), allocatable :: rho_test(:,:), kappa_test(:,:)
    real(KIND=dp), allocatable    :: transfo(:,:)
    real(KIND=dp)                 :: maxdev(3,2)
    integer                       :: i, it

    if(pairingtype.eq.1) then
      allocate(rho_test(nwt,nwt))
      do i=1,nwt
            rho_test(i,i) = rho_can(i)
      enddo
    else
      rho_test = rho_pairing
    endif

    ! 1. Calculate the mean-field densities in the canonical basis
    Density     = densit(rho_can, kappa_pairing)
    ! 2. Generate a random unitary transformation with real coefficients
    transfo = gen_unitary_transform()
    ! 3. Transform the matrices rho and kappa, as well as the HFPsi array
    !    with this transformation
    call mixup_rhokappa(rho_test, kappa_test, transfo)
    ! 4. Resum the densities with the offdiagonal routine
    R_transformed  = densit_offdiag_symmetric(rho_test, kappa_test)

    ! 5. Compare
    ! TODO: make this more systematic!
    do it=1,2
        maxdev(1,it) = maxval(abs(Density%D_I_I(:,it)   - R_transformed%D_I_I  (:,it)))
        !maxdev(2,it) = maxval(abs(Density%D_Nm_Nm(:,it) - R_transformed%D_Nm_Nm(:,it)))
!         maxdev(3,it) = maxval(abs(Density%C_I_Ns(:,:,:,it) - R_transformed%C_I_Ns(:,:,:,it)))
    enddo

    if(any(maxdev .gt. 1e-10)) then
        print *,'--------------------- Deviation in densities detected ---------------------'
        print *, 'Density    Neutrons      Protons'
        print ('(a4,2es15.7)'), 'RHO',  maxdev(1,:)
        print ('(a4,2es15.7)'), 'TAU',  maxdev(2,:)
        print ('(a4,2es15.7)'), 'Jmn',  maxdev(3,:)
        print *,'---------------------------------------------------------------------------'
        ifail = 1
    else
        ifail = 0
    endif

    ! 6. Restore the spwfs to their original condition in order to not mess with other tests
    call mixup_rhokappa(rho_test, kappa_test, transfo)

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
    dcmplx( 3.0_dp,  0.0_dp), dcmplx(-1.0_dp,  1.0_dp), dcmplx( 0.0_dp,  0.0_dp), dcmplx( 0.0_dp,  0.0_dp), dcmplx( 0.0_dp,  0.0_dp), dcmplx( 0.0_dp,  0.0_dp), &
    dcmplx( 1.0_dp, -1.0_dp), dcmplx( 2.0_dp,  0.0_dp), dcmplx(-1.0_dp,  1.0_dp), dcmplx( 0.0_dp,  0.0_dp), dcmplx( 0.0_dp,  0.0_dp), dcmplx( 0.0_dp,  0.0_dp), &
    dcmplx( 0.0_dp,  0.0_dp), dcmplx( 1.0_dp,  1.0_dp), dcmplx( 4.0_dp,  0.0_dp), dcmplx(-1.0_dp,  1.0_dp), dcmplx( 0.0_dp,  0.0_dp), dcmplx( 0.0_dp,  0.0_dp), &
    dcmplx( 0.0_dp,  0.0_dp), dcmplx( 0.0_dp,  0.0_dp), dcmplx( 1.0_dp, -1.0_dp), dcmplx( 5.0_dp,  0.0_dp), dcmplx(-1.0_dp,  1.0_dp), dcmplx( 0.0_dp,  0.0_dp), &
    dcmplx( 0.0_dp,  0.0_dp), dcmplx( 0.0_dp,  0.0_dp), dcmplx( 0.0_dp,  0.0_dp), dcmplx( 1.0_dp,  1.0_dp), dcmplx( 6.0_dp,  0.0_dp), dcmplx(-1.0_dp,  1.0_dp), &
    dcmplx( 0.0_dp,  0.0_dp), dcmplx( 0.0_dp,  0.0_dp), dcmplx( 0.0_dp,  0.0_dp), dcmplx( 0.0_dp,  0.0_dp), dcmplx( 1.0_dp, -1.0_dp), dcmplx( 7.0_dp,  0.0_dp)  &
    ], [n, n])

    ! A = reshape([ &
    ! dcmplx( 3.0_dp, 0.0_dp), dcmplx(-1.0_dp, 0.0_dp), dcmplx( 0.0_dp, 0.0_dp), dcmplx( 0.0_dp, 0.0_dp), dcmplx( 0.0_dp, 0.0_dp), dcmplx( 0.0_dp, 0.0_dp), &
    ! dcmplx( 1.0_dp, 0.0_dp), dcmplx( 2.0_dp, 0.0_dp), dcmplx(-1.0_dp, 0.0_dp), dcmplx( 0.0_dp, 0.0_dp), dcmplx( 0.0_dp, 0.0_dp), dcmplx( 0.0_dp, 0.0_dp), &
    ! dcmplx( 0.0_dp, 0.0_dp), dcmplx( 1.0_dp, 0.0_dp), dcmplx( 4.0_dp, 0.0_dp), dcmplx(-1.0_dp, 0.0_dp), dcmplx( 0.0_dp, 0.0_dp), dcmplx( 0.0_dp, 0.0_dp), &
    ! dcmplx( 0.0_dp, 0.0_dp), dcmplx( 0.0_dp, 0.0_dp), dcmplx( 1.0_dp, 0.0_dp), dcmplx( 5.0_dp, 0.0_dp), dcmplx(-1.0_dp, 0.0_dp), dcmplx( 0.0_dp, 0.0_dp), &
    ! dcmplx( 0.0_dp, 0.0_dp), dcmplx( 0.0_dp, 0.0_dp), dcmplx( 0.0_dp, 0.0_dp), dcmplx( 1.0_dp, 0.0_dp), dcmplx( 6.0_dp, 0.0_dp), dcmplx(-1.0_dp, 0.0_dp), &
    ! dcmplx( 0.0_dp, 0.0_dp), dcmplx( 0.0_dp, 0.0_dp), dcmplx( 0.0_dp, 0.0_dp), dcmplx( 0.0_dp, 0.0_dp), dcmplx( 1.0_dp, 0.0_dp), dcmplx( 7.0_dp, 0.0_dp)  &
    ! ], [n, n])

    ! A = reshape([ (3.0_dp, 0.0_dp), (-1.0_dp, 0.0_dp), (0.0_dp, 0.0_dp), &
    !             (-1.0_dp, 0.0_dp), (3.0_dp, 0.0_dp), (-1.0_dp, 0.0_dp), &
    !             (0.0_dp, 0.0_dp), (-1.0_dp, 0.0_dp), (3.0_dp, 0.0_dp) ], [n, n])

    Atmp = A

    b = [dcmplx(1.0_dp, 0.0_dp), dcmplx(1.0_dp, 0.0_dp), dcmplx(1.0_dp, 0.0_dp), dcmplx(1.0_dp, 0.0_dp), dcmplx(1.0_dp, 0.0_dp), dcmplx(1.0_dp, 0.0_dp)]

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
    ! Call the choral GMRES routine to solve Ax = b iteratively

    x_choral = b ! initial guess
    call do_gmres_choral(x_choral, iter, nbprod, res, b, multiply_by_A, norm_2, ScalProd, 1e-5_dp, 1, 5, 3)
    

    print *, "GMRES choral solution:"
    do i=1,n
      print "(*('(', F8.5, ',', F8.5, ') ', :))",  x_choral(i)
    enddo
    print *, 'res : ', sqrt(sum(abs(b - matmul(A, x_choral)) ** 2 ))

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Test my GMRES routine to solve Ax = b iteratively

    call alloc_gmres(multiply_by_A, b, 100, 6, 1e-6_dp, n, norm_2, ScalProd)

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
    complex(KIND=dp), dimension(:), target, intent(in)   :: x_in
    complex(KIND=dp), dimension(:), target, intent(out)  :: x_out
    integer, parameter :: n = 6
    complex(KIND=dp), dimension(n, n) :: A
    integer :: i

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Define the matrix A
    A = reshape([ &
    dcmplx( 3.0_dp,  0.0_dp), dcmplx(-1.0_dp,  1.0_dp), dcmplx( 0.0_dp,  0.0_dp), dcmplx( 0.0_dp,  0.0_dp), dcmplx( 0.0_dp,  0.0_dp), dcmplx( 0.0_dp,  0.0_dp), &
    dcmplx( 1.0_dp, -1.0_dp), dcmplx( 2.0_dp,  0.0_dp), dcmplx(-1.0_dp,  1.0_dp), dcmplx( 0.0_dp,  0.0_dp), dcmplx( 0.0_dp,  0.0_dp), dcmplx( 0.0_dp,  0.0_dp), &
    dcmplx( 0.0_dp,  0.0_dp), dcmplx( 1.0_dp,  1.0_dp), dcmplx( 4.0_dp,  0.0_dp), dcmplx(-1.0_dp,  1.0_dp), dcmplx( 0.0_dp,  0.0_dp), dcmplx( 0.0_dp,  0.0_dp), &
    dcmplx( 0.0_dp,  0.0_dp), dcmplx( 0.0_dp,  0.0_dp), dcmplx( 1.0_dp, -1.0_dp), dcmplx( 5.0_dp,  0.0_dp), dcmplx(-1.0_dp,  1.0_dp), dcmplx( 0.0_dp,  0.0_dp), &
    dcmplx( 0.0_dp,  0.0_dp), dcmplx( 0.0_dp,  0.0_dp), dcmplx( 0.0_dp,  0.0_dp), dcmplx( 1.0_dp,  1.0_dp), dcmplx( 6.0_dp,  0.0_dp), dcmplx(-1.0_dp,  1.0_dp), &
    dcmplx( 0.0_dp,  0.0_dp), dcmplx( 0.0_dp,  0.0_dp), dcmplx( 0.0_dp,  0.0_dp), dcmplx( 0.0_dp,  0.0_dp), dcmplx( 1.0_dp, -1.0_dp), dcmplx( 7.0_dp,  0.0_dp)  &
    ], [n,n])

    ! A = reshape([ &
    ! dcmplx( 3.0_dp, 0.0_dp), dcmplx(-1.0_dp, 0.0_dp), dcmplx( 0.0_dp, 0.0_dp), dcmplx( 0.0_dp, 0.0_dp), dcmplx( 0.0_dp, 0.0_dp), dcmplx( 0.0_dp, 0.0_dp), &
    ! dcmplx( 1.0_dp, 0.0_dp), dcmplx( 2.0_dp, 0.0_dp), dcmplx(-1.0_dp, 0.0_dp), dcmplx( 0.0_dp, 0.0_dp), dcmplx( 0.0_dp, 0.0_dp), dcmplx( 0.0_dp, 0.0_dp), &
    ! dcmplx( 0.0_dp, 0.0_dp), dcmplx( 1.0_dp, 0.0_dp), dcmplx( 4.0_dp, 0.0_dp), dcmplx(-1.0_dp, 0.0_dp), dcmplx( 0.0_dp, 0.0_dp), dcmplx( 0.0_dp, 0.0_dp), &
    ! dcmplx( 0.0_dp, 0.0_dp), dcmplx( 0.0_dp, 0.0_dp), dcmplx( 1.0_dp, 0.0_dp), dcmplx( 5.0_dp, 0.0_dp), dcmplx(-1.0_dp, 0.0_dp), dcmplx( 0.0_dp, 0.0_dp), &
    ! dcmplx( 0.0_dp, 0.0_dp), dcmplx( 0.0_dp, 0.0_dp), dcmplx( 0.0_dp, 0.0_dp), dcmplx( 1.0_dp, 0.0_dp), dcmplx( 6.0_dp, 0.0_dp), dcmplx(-1.0_dp, 0.0_dp), &
    ! dcmplx( 0.0_dp, 0.0_dp), dcmplx( 0.0_dp, 0.0_dp), dcmplx( 0.0_dp, 0.0_dp), dcmplx( 0.0_dp, 0.0_dp), dcmplx( 1.0_dp, 0.0_dp), dcmplx( 7.0_dp, 0.0_dp)  &
    ! ], [n, n])


    x_out = matmul(A, x_in)

  end subroutine multiply_by_A

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

end module fam_testing
