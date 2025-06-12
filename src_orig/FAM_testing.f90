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
!     call test_eta_independence(X,Y,ifail)
!     print 1, 'ETA_INDEPENDENCE', ifail

    !çcall test_potentials(X,Y,ifail)
    ! Attention: this testing routine has serious side effects on the state of the program.
    !call test_densit_offdiag(ifail)
    !print 1, 'DENSIT_OFFDIAG', ifail

    stop
  end subroutine run_FAM_tests
!
!   subroutine test_eta_independence(X,Y,ifail)
!     !---------------------------------------------------------------------------
!     ! TODO: document
!     !
!     !
!     !---------------------------------------------------------------------------
!     integer, intent(out)          :: ifail
!     complex(KIND=dp), intent(in)  :: X(:,:), Y(:,:)
!
!     integer                       :: ie, i, neta
!     real(KIND=dp)                 :: eta
!     type(DensityVector)           :: R, Rs, Ra
!     type(PotentialVector)         :: F, Fpert, Falt
!     real(KIND=dp), allocatable    :: dev_F_Nm_Nm(:), dev_G_I_NS(:)
!     complex(KIND=dp), allocatable :: drho(:,:), dkappa(:,:)
!     complex(KIND=dp), allocatable :: perturbed_F_Nm_Nm(:,:,:)
!     complex(KIND=dp), allocatable :: perturbed_G_I_NS(:,:,:,:,:)
!
!     neta = 10
!     allocate(drho(nwt,nwt))
!     allocate(perturbed_F_Nm_Nm(mv,4,neta))
!     allocate(perturbed_G_I_NS(mv,3,3,4,neta))
!     allocate(dev_F_Nm_Nm(neta),dev_G_I_NS(neta))
!
!     do ie = 1, neta
!       eta = 0.01*ie
!
!       !- - - - - - - - - - - - - - -- - - - - -
!       ! Building the perturbed density matrix
!       drho = 0.0d0
!       do i=1,nwt
!         drho(i,i) = rho_can(i)
!       enddo
!       drho = drho +  eta * (X + transpose(Y))
!
!       ! Building the densities and the perturbed densities
!       R       = densit(rho_can, kappa_pairing)
!       call densit_offdiag(drho, dkappa,Rs, Ra)
!       ! ... and the potentials and perturbed potentials
!       F       = calcPotentials(R)
!       FPert   = calcPotentials(RPert)
!
!       perturbed_F_Nm_Nm(:,:,ie)    = (FPert%F_Nm_Nm - F%F_Nm_Nm)/eta
!       perturbed_G_I_NS(:,:,:,:,ie) = (FPert%G_I_NS  - F%G_I_NS)/eta
!     enddo
!
!     do ie=1,neta
!       dev_F_Nm_Nm(ie) = maxval(abs(perturbed_F_Nm_Nm(:,:,ie) - perturbed_F_Nm_Nm(:,:,1)))
!       dev_G_I_NS(ie) = maxval(abs(perturbed_G_I_NS(:,:,:,:,ie) - perturbed_G_I_NS(:,:,:,:,1)))
!     enddo
!
! !     print *, 'DEVIATIONS'
! !     print *, '     eta       F_Nm_Nm      G_I_NS'
! !     do ie=1,neta
! !           eta = 0.01*ie
! !           print ('(99es12.3)'), eta, dev_F_Nm_Nm(ie), dev_G_I_NS(ie)
! !     enddo
!
!     ifail = 0
!     if(any(dev_F_Nm_Nm.gt.1e-10)) ifail = 1
!     if(any(dev_G_I_NS .gt.1e-10)) ifail = 1
!
!   end subroutine test_eta_independence
!
  subroutine test_potentials(X,Y,ifail)
    !------------------------------------------------------------------------
    !
    !
    !
    !------------------------------------------------------------------------
    integer, intent(out)          :: ifail
    complex(KIND=dp), intent(in)  :: X(:,:), Y(:,:)

    complex(KIND=dp), allocatable :: drho(:,:), dkappa(:,:), sphamil_me(:,:)

    type(DensityVector)           :: R, dRa, dRs
    type(PotentialVector)         :: F, dFs, dFa
    
    allocate(drho(nwt,nwt))
    drho = 0.0d0
    drho = X + transpose(Y)

    R       = densit(rho_can, kappa_pairing)
    dRs     = densit_offdiag_symmetric(drho, dkappa)
    dRa     = densit_offdiag_antisymmetric(drho,dkappa)
    
    F       = calcpotentials(R)
    call calc_perturbed_potentials(R,dRs,dRa, dFs, dFa)

    sphamil_me = calc_sphamil_me( HFpsi, HFdpsi, HFddpsi,dFs,dFa, .false.)


  end subroutine test_potentials

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

    complex(KIND=dp), allocatable :: sphamil_me(:,:), sphamil_recalc(:,:), hpsi(:,:), drho(:,:), dkappa(:,:)
    integer                    :: si, B, N, i, it, j
    type(PotentialVector)      :: Fs, Fa, F
    type(DensityVector)        :: R, Rs, Ra

    R  = densit(rho_can, kappa_pairing)
    call ConstructChargeDensity(R)

    allocate(drho(nwt,nwt)) ; drho = 0.0d0
    call densit_offdiag(drho, dkappa, Rs, Ra)
    F = calcpotentials(R)
    call calc_perturbed_potentials(R, Rs, Ra, Fs, Fa)

    ! Calculation with the new routine
    sphamil_me = calc_sphamil_me( HFpsi, HFdpsi, HFddpsi,F, Fa, .false.)
    sphamil_me = sphamil_me + kinetic_me(HFpsi, HFdpsi, hfddpsi)
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
      print ('(99f10.3)'), spenergies(si+1:si+N)
      print *, 'SPH on file'
      do i=1,N
        print ('(99f10.3)'), sphamil_recalc(si+i, si+1:si+N)
      enddo
      print *, 'SPH from calc_sphamil_me'
      do i=1,N
        print ('(99f10.3)'), sphamil_me(si+i, si+1:si+N)
      enddo
      print *
      print *, 'Difference'
      do i=1,N
        print ('(99es10.2)'), abs(sphamil_recalc(si+i, si+1:si+N) - sphamil_me(si+i, si+1:si+N))
      enddo
      si = si + N
    enddo

    ifail = 0
  end subroutine test_sphamil_me

  function kinetic_me(denpsi, dendpsi, denddpsi)
    !------------------------------------------------------------------------------------
    ! TODO: DOCUMENT
    !------------------------------------------------------------------------------------

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
#if(PASTA == 0)
      select case(COM1Body)
      case(0,1)
        Reducedmass = 1.0_dp
      case(2)
        Reducedmass = (1.0_dp-nucleonmass(it)/                                   &
        &                      (neutrons*nucleonmass(1)+protons*nucleonmass(2)))
      end select
#endif

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
