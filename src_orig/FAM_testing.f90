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

    complex(KIND=dp), intent(in) :: X(:,:), Y(:,:)
    integer :: ifail


    call test_sphamil_me(ifail)
    print 1, 'SPHAMIL_ME', ifail
    call test_eta_independence(X,Y,ifail)
    print 1, 'ETA_INDEPENDENCE', ifail

    ! Attention: this testing routine has serious side effects on the state of the program.
    !call test_densit_offdiag(ifail)
    !print 1, 'DENSIT_OFFDIAG', ifail

    stop
  end subroutine run_FAM_tests

  subroutine test_eta_independence(X,Y,ifail)
    !---------------------------------------------------------------------------
    ! TODO: document
    !
    !
    !---------------------------------------------------------------------------
    integer, intent(out)          :: ifail
    complex(KIND=dp), intent(in)  :: X(:,:), Y(:,:)

    integer                       :: ie, i, neta
    real(KIND=dp)                 :: eta
    type(DensityVector)           :: R, Rpert
    type(PotentialVector)         :: F, Fpert, Falt
    real(KIND=dp), allocatable    :: dev_F_Nm_Nm(:), dev_G_I_NS(:)
    complex(KIND=dp), allocatable :: drho(:,:), dkappa(:,:)
    complex(KIND=dp), allocatable :: perturbed_F_Nm_Nm(:,:,:)
    complex(KIND=dp), allocatable :: perturbed_G_I_NS(:,:,:,:,:)

    neta = 10
    allocate(drho(nwt,nwt))
    allocate(perturbed_F_Nm_Nm(mv,4,neta))
    allocate(perturbed_G_I_NS(mv,3,3,4,neta))
    allocate(dev_F_Nm_Nm(neta),dev_G_I_NS(neta))

    do ie = 1, neta
      eta = 0.01*ie

      !- - - - - - - - - - - - - - -- - - - - -
      ! Building the perturbed density matrix
      drho = 0.0d0
      do i=1,nwt
        drho(i,i) = rho_can(i)
      enddo
      drho = drho +  eta * (X + transpose(Y))

      ! Building the densities and the perturbed densities
      R       = densit(rho_can, kappa_pairing)
      RPert   = densit_offdiag(drho, dkappa)
      ! ... and the potentials and perturbed potentials
      F       = calcPotentials(R)
      FPert   = calcPotentials(RPert)

      perturbed_F_Nm_Nm(:,:,ie)    = (FPert%F_Nm_Nm - F%F_Nm_Nm)/eta
      perturbed_G_I_NS(:,:,:,:,ie) = (FPert%G_I_NS  - F%G_I_NS)/eta
    enddo

    do ie=1,neta
      dev_F_Nm_Nm(ie) = maxval(abs(perturbed_F_Nm_Nm(:,:,ie) - perturbed_F_Nm_Nm(:,:,1)))
      dev_G_I_NS(ie) = maxval(abs(perturbed_G_I_NS(:,:,:,:,ie) - perturbed_G_I_NS(:,:,:,:,1)))
    enddo

!     print *, 'DEVIATIONS'
!     print *, '     eta       F_Nm_Nm      G_I_NS'
!     do ie=1,neta
!           eta = 0.01*ie
!           print ('(99es12.3)'), eta, dev_F_Nm_Nm(ie), dev_G_I_NS(ie)
!     enddo

    ifail = 0
    if(any(dev_F_Nm_Nm.gt.1e-10)) ifail = 1
    if(any(dev_G_I_NS .gt.1e-10)) ifail = 1

  end subroutine test_eta_independence

!   subroutine test_linearity_response(X,Y,ifail)
!     !-------------------------------------------------------------------------------
!     ! Test whether the response of the individual fields to a small perturbation in
!     ! the densities is indeed linear w.r.t. to said densities.
!     ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     ! Output:
!     !  ifail : 0 if succesful, if a deviation above 1e-10 has been detected
!     !-------------------------------------------------------------------------------
!     integer, intent(out) :: ifail
!     real(KIND=dp)                 :: eta
!     type(DensityVector)           :: R, Rpert
!     type(PotentialVector)         :: F, Fpert, Falt
!     complex(KIND=dp), allocatable :: drho(:,:), dkappa(:,:), drho_sym(:,:), drho_asym(:,:)
!     complex(KIND=dp), intent(in)  :: X(:,:), Y(:,:)
!
!     integer :: i, si, B, N
!
!     ifail = 0
!
!     allocate(drho(nwt,nwt)); drho = 0.0d0
!     do i=1,nwt
!       drho(i,i) = rho_can(i)
!     enddo
!
!
!     eta =1
!     drho = drho +  eta * (X + transpose(Y))
!
!     drho_sym  = 0.5*(drho + transpose(drho))
!     drho_asym = 0.5*(drho - transpose(drho))
!
!     si = 0
!     do B=1,8
!       N = HFBlocks(B)
!
!       print *, 'B = ', B
!       print *, 'REAL PART X'
!       do i=1,N
!         print ('(99f10.3)'), DBLE(X(si+i, si+1:si+N))
!       enddo
!       print *, 'IMAG PART X'
!       do i=1,N
!         print ('(99f10.3)'), IMAG(X(si+i, si+1:si+N))
!       enddo
!       print *, 'REAL PART Y'
!       do i=1,N
!         print ('(99f10.3)'), DBLE(Y(si+i, si+1:si+N))
!       enddo
!       print *, 'IMAG PART Y'
!       do i=1,N
!         print ('(99f10.3)'), IMAG(Y(si+i, si+1:si+N))
!       enddo
!
!       print *, 'REAL PART DRHO'
!       do i=1,N
!         print ('(99f10.3)'), DBLE(drho(si+i, si+1:si+N))
!       enddo
!       print *, 'IMAG PART DRHO'
!       do i=1,N
!         print ('(99f10.3)'), IMAG(drho(si+i, si+1:si+N))
!       enddo
!       print *
!
!       print *, 'SYM PART DRHO'
!       do i=1,N
!         print ('(99f10.3)'), DBLE(drho_sym(si+i, si+1:si+N))
!       enddo
!       print *, 'ANTISYM PART DRHO'
!       do i=1,N
!         print ('(99f10.3)'), DBLE(drho_asym(si+i, si+1:si+N))
!       enddo
!       print *
!       print *
!
!       si = si + N
!     enddo
!
!
!     R       = densit(rho_can, kappa_pairing)
!     RPert   = densit_offdiag(drho, dkappa)
!
!     F       = calcPotentials(R)
!     FPert   = calcPotentials(RPert)
!
!     print *, 'ORIGINAL'
!     print *, '-------------------------------------------------------------------------------------------'
!     print *, 'ISOSCALAR D_Nm_Nm'
!     do i=1,nx
!       print ('(8f15.8)'),  R%D_I_I(i,3), F%F_Nm_Nm(i,3), DBLE(F%F_Nm_Nm(i,3))/DBLE(R%D_I_I(i,3)), &
!       &                                                  IMAG(F%F_Nm_Nm(i,3))/IMAG(R%D_I_I(i,3))
!     enddo
!     print *, 'ISOVECTOR D_Nm_Nm'
!     do i=1,nx
!       print ('(8f15.8)'),  R%D_I_I(i,4), F%F_Nm_Nm(i,4), DBLE(F%F_Nm_Nm(i,4))/DBLE(R%D_I_I(i,4)), &
!       &                                                  IMAG(F%F_Nm_Nm(i,4))/IMAG(R%D_I_I(i,4))
!     enddo
!
!     print *, 'PERTURBED'
!     print *, '-------------------------------------------------------------------------------------------'
!     print *, 'ISOSCALAR D_Nm_Nm'
!     do i=1,nx
!       print ('(8f15.8)'),  Rpert%D_I_I(i,3), FPert%F_Nm_Nm(i,3), DBLE(FPert%F_Nm_Nm(i,3))/DBLE(Rpert%D_I_I(i,3)), &
!       &                                                          IMAG(FPert%F_Nm_Nm(i,3))/IMAG(Rpert%D_I_I(i,3))
!     enddo
!     print *, 'ISOVECTOR D_Nm_Nm'
!     do i=1,nx
!       print ('(8f15.8)'),  Rpert%D_I_I(i,4), FPert%F_Nm_Nm(i,4), DBLE(FPert%F_Nm_Nm(i,4))/DBLE(Rpert%D_I_I(i,4)), &
!       &                                                          IMAG(FPert%F_Nm_Nm(i,4))/IMAG(Rpert%D_I_I(i,4))
!     enddo
!
!     print *, 'MANUAL DIFFERENCE'
!     print *, '------------------------------------------------------------------------------------------'
!     print *, 'ISOSCALAR D_Nm_Nm'
!     do i=1,nx
!       print ('(8f15.8)'),  Rpert%D_I_I(i,3)-R%D_I_I(i,3), FPert%F_Nm_Nm(i,3)-F%F_Nm_Nm(i,3), &
!       &     DBLE(FPert%F_Nm_Nm(i,3) - F%F_Nm_Nm(i,3))/DBLE(Rpert%D_I_I(i,3) - R%D_I_I(i,3)), &
!       &     IMAG(FPert%F_Nm_Nm(i,3) - F%F_Nm_Nm(i,3))/IMAG(Rpert%D_I_I(i,3) - R%D_I_I(i,3))
!     enddo
!
!     print *, 'OPERATOR DIFFERENCE'
!     print *, '--------------------------------------------------------------------------------------------'
!     RPert   = densit_offdiag(drho, dkappa)  + (-1.0d0)*R
!     Falt = (-1.0d0) * F
!     Fpert = Fpert + Falt
!     print *, 'ISOSCALAR D_Nm_Nm'
!     do i=1,nx
!       print ('(10f15.8)'),  Rpert%D_I_I(i,3), FPert%F_Nm_Nm(i,3), DBLE(FPert%F_Nm_Nm(i,3))/DBLE(Rpert%D_I_I(i,3)),&
!       &                                                           IMAG(FPert%F_Nm_Nm(i,3))/IMAG(Rpert%D_I_I(i,3))
!     enddo
!     print *, 'ISOVECTOR D_Nm_Nm'
!     do i=1,nx
!       print ('(8f15.8)'),  Rpert%D_I_I(i,4), FPert%F_Nm_Nm(i,4), DBLE(FPert%F_Nm_Nm(i,4))/DBLE(Rpert%D_I_I(i,4)), &
!       &                                                          IMAG(FPert%F_Nm_Nm(i,4))/IMAG(Rpert%D_I_I(i,4))
!     enddo
!
!
!   end subroutine test_linearity_response

  subroutine test_sphamil_me(ifail)
    !-------------------------------------------------------------------------------
    ! Test whether the matrix elements of the single-particle Hamiltonian in the HF
    ! basis are identical whether calculate through (i) apply_sphamil and
    ! (ii) calc_sphamil_me.
    !
    ! TODO: define failing case.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Output:
    !  ifail : 0 if succesful, if a deviation above 1e-10 has been detected
    !-------------------------------------------------------------------------------
    integer, intent(out)       :: ifail

    real(KIND=dp), allocatable :: sphamil_me(:,:), sphamil_recalc(:,:), hpsi(:,:)
    integer                    :: si, B, N, i, it, j
    type(PotentialVector)      :: F
    type(DensityVector)        :: R

    R  = densit(rho_can, kappa_pairing)
    call ConstructChargeDensity(R)
    F       = calcPotentials(R)

    ! Calculation with the new routine
    sphamil_me = calc_sphamil_me( HFpsi, HFdpsi, HFddpsi,F, .false.)

    allocate(sphamil_recalc(nwt,nwt)); sphamil_recalc = 0.0d0
    si = 0
    do B=1,8
      N = HFBlocks(B)
      it = -1
      if(B .ge. 5) it = +1
      do j=si+1,si+N
        hpsi = apply_sphamil(HFPsi(:,:,j), HFdPsi(:,:,:,j), HFddPsi(:,:,:,j), sx(:,j), sy(:,j), sz(:,j), it ,.false. ,F)
        do i=si+1,si+N
          sphamil_recalc(i,j) = sum(HFpsi(:,:,i) * hpsi)*dv
        enddo
      enddo
      si = si + N
    enddo

    ! print output
    si = 0
    do B = 1,8
      N = HFBLocks(B)
      print *, 'BLOCK B=', B
      print ('(99f12.5)'), spenergies(si+1:si+N)
      print *, 'SPH on file'
      do i=1,N
        print ('(99f12.5)'), sphamil_recalc(si+i, si+1:si+N)
      enddo
      print *, 'SPH from calc_sphamil_me'
      do i=1,N
        print ('(99f12.5)'), sphamil_me(si+i, si+1:si+N)
      enddo
      print *
      print *, 'Difference'
      do i=1,N
        print ('(99es12.2)'), abs(sphamil_recalc(si+i, si+1:si+N) - sphamil_me(si+i, si+1:si+N))
      enddo
      si = si + N
    enddo

    ifail = 0
  end subroutine test_sphamil_me

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
    R_transformed  = densit_offdiag(rho_test, kappa_test)

    ! 5. Compare
    ! TODO: make this more systematic!
    do it=1,2
        maxdev(1,it) = maxval(abs(Density%D_I_I(:,it)   - R_transformed%D_I_I  (:,it)))
        maxdev(2,it) = maxval(abs(Density%D_Nm_Nm(:,it) - R_transformed%D_Nm_Nm(:,it)))
        maxdev(3,it) = maxval(abs(Density%C_I_Ns(:,:,:,it) - R_transformed%C_I_Ns(:,:,:,it)))
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

end module fam_testing
