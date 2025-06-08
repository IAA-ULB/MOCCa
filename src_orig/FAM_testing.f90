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
  subroutine run_FAM_tests()
    !--------------------------------------------------------------------------------
    ! Catch-all routine to run all predefined unit tests for FAM routines.
    ! Input:
    !   None
    ! Output:
    !   None
    !--------------------------------------------------------------------------------
    integer :: ifail

    !call test_densit_offdiag(ifail)
    call test_linearity_response(ifail)
    call test_sphamil_me(ifail)

    stop
  end subroutine run_FAM_tests

  subroutine test_linearity_response(ifail)
    !-------------------------------------------------------------------------------
    ! Test whether the response of the individual fields to a small perturbation in
    ! the densities is indeed linear w.r.t. to said densities.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Output:
    !  ifail : 0 if succesful, if a deviation above 1e-10 has been detected
    !-------------------------------------------------------------------------------
    integer, intent(out) :: ifail

  end subroutine test_linearity_response

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

    real(KIND=dp), allocatable :: sphamil_me(:,:)
    integer                    :: si, B, N, i
    type(PotentialVector)      :: F
    type(DensityVector)        :: R

      do i=1,10
        print ('(99f10.3)'), sphamil(i, 1:10)
      enddo

    ! Unfortunately, the Coulomb potentials are not currently read from file
    ! hence we need to regenerate them explicitly.
    R  = densit(rho_can, kappa_pairing)
    call ConstructChargeDensity(R)
    F  = potentials_read
    call SolveCoulomb(R, F)

    ! Ensure that potentials match the conventions of sphamil_me with combine_potentials
    call combine_potentials(F)
    sphamil_me = calc_sphamil_me( HFpsi, HFdpsi, HFddpsi,F, .false.)

    ! print output
    si = 0
    do B = 1,8
      N = HFBLocks(B)
      print *, 'BLOCK B=', B
      print ('(99f12.5)'), spenergies(si+1:si+N)
      print *, 'SPH on file'
      do i=1,N
        print ('(99f12.5)'), sphamil(si+i, si+1:si+N)
      enddo
      print *, 'SPH from calc_sphamil_me'
      do i=1,N
        print ('(99f12.5)'), sphamil_me(si+i, si+1:si+N)
      enddo
      print *
      print *, 'Difference'
      do i=1,N
        print ('(99es12.2)'), abs(sphamil(si+i, si+1:si+N) - sphamil_me(si+i, si+1:si+N))
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
