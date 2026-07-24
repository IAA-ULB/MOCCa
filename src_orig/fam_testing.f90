!===============================================================================
!     __  __  ___   ____ ____
!    |  \/  |/ _ \ / ___/ ___|__ _
!    | |\/| | | | | |  | |   / _` |
!    | |  | | |_| | |__| |__| (_| |
!    |_|  |_|\___/ \____\____\__,_|
!
!    Copyright (C) 2026 W. Ryssens and M. Bender
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU Affero General Public License as published
!    by the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU Affero General Public License for more details.
!
!    You should have received a copy of the GNU Affero General Public License
!    along with this program.  If not, see <https://www.gnu.org/licenses/>.
!
!===============================================================================
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


  interface print_matrix_symmetries
    module procedure print_matrix_symmetries_complex
    module procedure print_matrix_symmetries_real
  end interface print_matrix_symmetries

contains

  subroutine run_FAM_tests()
    !---------------------------------------------------------------------------------
    ! Catch-all routine to run all predefined unit tests for FAM routines.
    !
    !  Notes
    !   - this routine assumes all of the setup to start FAM iterations has been done!
    !   - some of these tests dramatically change the state of the code, which is why
    !     the end of this routine features a call to exit
    !   - the code returns an exit status (0/1) to allow for this to be part of the
    !     MOCCa testing framework
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !   None
    ! Output:
    !   None
    !
    ! - Exit status = 0 => everything okay!
    !               = 1 => tests failed
    !---------------------------------------------------------------------------------

    1 format ("Test = ", a40, " Success = ", i4)
    9 format (" --------------------- Start of the FAM testing routines ---------------------")
   10 format (" --------------------- End of the FAM testing routines -----------------------")
   11 format (" Summary of results ")
   12 format (" All tests passed! ")
   13 format (" Some tests failed! ")

    integer :: ifail_sp_qp, ifail_HFme, ifail_potentials, ifail_linearity, ifail

    print 9
    print *

    print 10
    print 11

    if(PairingType .eq. 2) then
      call test_qptransfo(ifail_sp_qp)
      print 1, 'SP<->QP transformations'      , ifail_sp_qp
    endif
    call test_potentials(ifail_potentials)
    print 1, 'Linearisation of potentials'    , ifail_potentials
    call test_HFmatrices_me(ifail_HFme)
    print 1, 'Matrix elements of h and \Delta', ifail_HFme
    call test_linearity_T(ifail_linearity)
    print 1, 'Linearity of FAM iteration'     , ifail_linearity

  
    ifail = max(ifail_sp_qp, ifail_HFme, ifail_potentials, ifail_linearity)
    if(ifail.eq.0) print 12
    if(ifail.eq.1) print 13
    call exit(ifail)

  end subroutine run_FAM_tests

  subroutine test_potentials(ifail)
    !------------------------------------------------------------------------
    ! This subroutines tests the linearisation of the potentials
    !
    ! TODO: 
    !  - use Hephaestos to generalize this test to more potentials!
    !  - document this routine 
    ! 
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !  None 
    ! Output:
    !  ifail: 0 if all tests passed, 1 otherwise
    !------------------------------------------------------------------------
    integer, intent(out)          :: ifail
    complex(KIND=dp), allocatable :: drho(:,:), dkappa(:,:), dH(:,:,:)
    real(KIND=dp)                 :: eta
    logical                       :: check

    type(DensityVector)           :: R, dRa, dRs, dR_pp_plus, dR_pp_minus
    type(PotentialVector)         :: F, dFs, dFa, Fnew, dF_pp_minus, dF_pp_plus

    integer :: i

    ifail = 0

  1 format (100('-'))
 99 format (45('-'), ' Testing potentials ', 45('-') )

    print 99
    print *
    ! We have to be careful - a QFAM calculation does not by default construct the 
    !   canonical basis, but the mean-field routine densit requires it. 
    if(.not.allocated(canDpsi)) then
      allocate(candPsi  (nx*ny*nz, 3,4,nwt)) ! first order
      allocate(canddPsi (nx*ny*nz, 6,4,nwt)) ! full tensor second order
      call construct_canonical_basis(rho_pairing, kappa_pairing, rho_can, kappa_can)
    endif 

    ! Construct the perturbed density and anomalous density matrices 
    if(.not.allocated(drho))        allocate(drho(nwt,nwt))
    if(.not.allocated(dkappa_plus)) allocate(dkappa_plus(nwt,nwt), dkappa_minus(nwt,nwt))

    ! Set up free response X and Y amplitudes
    allocate(dH(  nwt,  nwt, 2)) 
    dH = 0.0d0
    call calculate_XY(dH,X,Y)

    ! Calculate the perturbed normal and anomalous density matrix
    if (pairingtype==0) then ! FAM
      drho = X  + transpose(Y)
      dkappa_plus  = 0  
      dkappa_minus = 0
    else ! QFAM
      call transform_qp_to_sp(Bogoliubov, OTRqp=X, OBLqp=Y,                                  & ! input 
      &                                   OTRsp=dkappa_plus, OTLsp=drho, OBLsp=dkappa_minus)   ! output
$TR   dkappa_minus = - dkappa_minus
    endif

    ! Construct the densities 
    R       = densit(rho_can, kappa_pairing)                                                 ! mean-field densities
    call densit_offdiag(drho, dkappa_plus, dkappa_minus, DRs, DRa, DR_pp_plus, DR_pp_minus)  ! perturbed densities

    ! Construct the potentials 
    F       = calcpotentials(R)                                               ! mean-field values
    call calc_perturbed_potentials(R, dRs, dRa, dR_pp_plus, dR_pp_minus, &    ! perturbed values
    &                                 dFs, dFa, dF_pp_plus, dF_pp_minus)

    ! small perturbation factor
    eta = 1d-8
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! 1. Test a perturbation of the symmetric part of the particle-hole perturbation
    print *
    Fnew    = calcpotentials(R + eta*dRs) + (-1.0d0) * F

    ! The perturbation of F_I_I is typically the hardest because of the density dependent term
    ifail = max(ifail, &
        & check_findiff_deviations('F_I_I  - symmetric', dFs%F_I_I       , Fnew%F_I_I/eta,        &
        & tol=1d-3))

    ifail = max(ifail, &
        & check_findiff_deviations('F_I_SX - symmetric', dFs%F_I_S(:,1,:), Fnew%F_I_S(:,1,:)/eta, &
        & tol=1d-9))
    ifail = max(ifail, &
        & check_findiff_deviations('F_I_SY - symmetric', dFs%F_I_S(:,2,:), Fnew%F_I_S(:,2,:)/eta, &
        & tol=1d-9))
    ifail = max(ifail, &
        & check_findiff_deviations('F_I_SZ - symmetric', dFs%F_I_S(:,3,:), Fnew%F_I_S(:,3,:)/eta, &
        & tol=1d-9))

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*
    ! 2. Test a perturbation of the antisymmetric part of the particle-hole perturbation
    Fnew  = calcpotentials(R + eta*dRa) + (-1.0d0) * F
    ifail = max(ifail, &
        & check_findiff_deviations('F_I_I - antisymmetric', dFa%F_I_I       , Fnew%F_I_I/eta, &
        & tol=1d-3))

    ifail = max(ifail, &
        & check_findiff_deviations('F_I_SX - antisymmetric', dFa%F_I_S(:,1,:), Fnew%F_I_S(:,1,:)/eta, &
        & tol=1d-9))
    ifail = max(ifail, &
        & check_findiff_deviations('F_I_SY - antisymmetric', dFa%F_I_S(:,2,:), Fnew%F_I_S(:,2,:)/eta, &
        & tol=1d-9))
    ifail = max(ifail, &
        & check_findiff_deviations('F_I_SZ - antisymmetric', dFa%F_I_S(:,3,:), Fnew%F_I_S(:,3,:)/eta, &
        & tol=1d-9))
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -* 
    ! 3. Test the perturbed pairing potentials
    Fnew  = calcpotentials(R + eta*DR_pp_plus) + (-1.0d0) * F
    ifail = max(ifail, &
        & check_findiff_deviations('FP_I_I  - plus' , DF_pp_plus%FP_I_I,  Fnew%FP_I_I/eta, &
        & tol=5d-7))

    Fnew  = calcpotentials(R + eta*DR_pp_minus) + (-1.0d0) * F
    ifail = max(ifail, &
        & check_findiff_deviations('FP_I_I  - minus', DF_pp_minus%FP_I_I, Fnew%FP_I_I/eta, &
        & tol=5d-7))

    print 1
    print *
  end subroutine test_potentials

  function check_findiff_deviations(name, ref, findiff, tol) result(ifail)
    !-----------------------------------------------------------------------
    ! Test that the potentials (ref, findiff) are equal, i.e.
    !
    !     | ref - f*findiff | < tol everywhere on the mesh
    !
    ! with f = 1 for a normal potential and f = 2 for a pairing potential.
    !
    ! The use-case of this routine is
    !       ref     =  potential calculated through FAM
    !       findiff =  potential calculated through the finite differencing
    !                  of the mean-field routine with explicitly perturbed
    !                  densities
    !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input :
    !    name   : string, used for printing
    !    ref    : first potential, probably a linear response potential
    !    findiff: potential obtained through finite differencing
    !    pairing: logical, whether it is a pairing potential or not
    !    tol    : real, tolerance to test
    ! Output :
    !    check  : logical, whether the test passed
    !-----------------------------------------------------------------------
    character(len=*), intent(in) :: name
    real(KIND=dp), intent(in)    :: tol
    complex(KIND=dp), intent(in) :: ref(:,:), findiff(:,:)
    logical                      :: check
    integer :: i, ifail

    check = .true.
    if( maxval(abs(ref - findiff)) > tol ) check = .false.

    if(.not. check) then
       print *, ' ------------------------------------------------------'
       print *, ' Comparison failed for ', name
       print *, ' Printing neutron values along the x-axis at j = k =1 '
       print *, name, ' real part - maxval = ', maxval(abs(DBLE(ref)))
       print *, '-------------------------'
       do i=1,nx
          print ('(i3, 2f10.3, es12.3)'), i, DBLE(ref(i,1)), DBLE(findiff(i,1)), &
               &                             DBLE(ref(i,1))- DBLE(findiff(i,1))
       enddo
       print *
       print *, name, ' imaginary part - maxval = ', maxval(abs(IMAG(ref)))
       print *, '-------------------------'
       do i=1,nx
          print ('(i3, 2f10.3, es12.3)'), i, IMAG(ref(i,1)), IMAG(findiff(i,1)), &
               &                             IMAG(ref(i,1)) -IMAG(findiff(i,1))
       enddo
       print *, ' ------------------------------------------------------'
       print *
    else
       print ('(a40,  1e15.3," < ", 1e15.3 )'), name,  maxval(abs(ref - findiff)), tol
    endif

    if (check) then
      ifail = 0
    else
      ifail = 1
    endif
  end function check_findiff_deviations

  subroutine test_HFmatrices_me(ifail)
    !--------------------------------------------------------------------------------------
    ! Test whether the matrix elements of the single-particle Hamiltonian and the pairing 
    ! gaps in the HF basis when calculated in two different ways. For this purpose, we first
    ! calculate the densities and potentials as usual in a mean-field code, and then use
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
    !  ifail : 0 if succesful, 1 if a deviation above 1d-10 has been detected
    !---------------------------------------------------------------------------------------
    integer, intent(out)       :: ifail
    complex(KIND=dp), allocatable :: sphamil_me(:,:), sphamil_orig(:,:), delta_orig(:,:), delta_me(:,:)
    complex(KIND=dp), allocatable :: hpsi(:,:), drho(:,:), dkappa(:,:), dkappa_plus(:,:), dkappa_minus(:,:)
    real(KIND=dp), allocatable    :: dev(:,:), dev_gaps(:,:)
    real(KIND=dp)                 :: stabfactor(2)
    integer                       :: si, B, N, i, it, j, N2, T
    type(PotentialVector)         :: dFs, dFa, F, dF_pp_minus, dF_pp_plus
    type(DensityVector)           :: R, dRs, dRa, R_pp_plus, R_pp_minus, dR_pp_plus, dR_pp_minus

  1 format (100('-'))
 99 format (25('-'), ' Matrix elements of h and Delta ', 25('-') )

    ifail = 0
    print 99
    print *
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! (i) Ordinary mean-field-like calculation

    ! We have to be careful - a QFAM calculation does not by default construct the 
    !   canonical basis, but the mean-field routine densit requires it.
    if(.not. allocated(candPsi)) then 
      allocate(candPsi  (nx*ny*nz, 3,4,nwt)) ! first order
      allocate(canddPsi (nx*ny*nz, 6,4,nwt)) ! full tensor second order
      call construct_canonical_basis(rho_pairing, kappa_pairing, rho_can, kappa_can)
    endif 

    R = densit(rho_can, kappa_pairing)
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
    if(pairingtype.eq.2) then
      stabfactor = 0
      call calcHFBgaps(FermiEnergy, stabfactor, F)
      delta_orig = HFBGaps
    endif 

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! (ii) FAM-like calculation
    allocate(drho(nwt,nwt)) ; drho = 0.0d0
    allocate(dkappa_plus(nwt,nwt))  ; dkappa_plus  = 0.0d0
    allocate(dkappa_minus(nwt,nwt)) ; dkappa_minus = 0.0d0
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
    ! (iii) Check result and print output if needed
    print *
    print *, 'Checking single-particle Hamiltonian matrix elements'
    print 1
    si = 0
    do B = 1,8
      N = HFBlocks(B); if(N.eq.0) cycle
      ! The element-wise deviation
      dev = abs(sphamil_orig(si+1:si+N, si+1:si+N) - sphamil_me(si+1:si+N, si+1:si+N))
      print ('(a15, i2, a20, es15.4)'), 'Block B = ', B, ' Max |h| = ', maxval(abs(sphamil_orig(si+1:si+N, si+1:si+N)))
      print ('(a15, i2, a20, es15.4)'), 'Block B = ', B, ' Max dev = ', maxval(dev)
      if(maxval(dev)>1d-10) then
        ifail = 1
        print *
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
      endif
      si = si + N
    enddo

    if(pairingtype.eq.2) then
      print *
      print *, 'Checking pairing gap matrix elements'
      print 1
      si = 0
      do B= 1,8,2
        N = HFBlocks(B); N2 = HFBlocks(B+1); T = N + N2; if(T.eq.0) cycle
        dev_gaps = abs(delta_orig(si+1:si+T, si+1:si+T) - delta_me(si+1:si+T, si+1:si+T))
        print ('(a15, i2, "-", i2, a20, es15.4)'), 'Blocks B = ', B, B+1, &
             & ' Max |Delta| = ', maxval(abs(delta_orig(si+1:si+T,si+1:si+T)))
        print ('(a15, i2, "-", i2, a20, es15.4)'), 'Blocks B = ', B, B+1, &
             &' Max dev = ', maxval(dev_gaps)
        if(maxval(dev_gaps)>1d-10) then
          ifail = 1
          print *
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
        endif
        si = si + T
      enddo
    endif
    print 1
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
    !  ifail : 0 if succesful, if a deviation above 1d-10 has been detected
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
    enddo

    if(any(maxdev .gt. 1d-10)) then
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

  subroutine test_L_Linv(X, Y, F, omega)
    !-------------------------------------------------------------------------------
    ! This routine tests the multiply_by_QRPAmat by checking that it is the inverse
    ! operator as solving FAM, i.e.
    !             *QRPAmat :     L : x -> L * x = f 
    !                 FAM  :  Linv : f -> Linv * f = x
    ! where x = (X, Y) and f = (F20,F02). Note that L and Linv also depends on omega. 
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! I dont want to bother to solve FAM here inside this function so please call
    ! this test when FAM is converged and pass the converged X, Y together with the
    ! external field F. 
    !-------------------------------------------------------------------------------
    complex(KIND=dp), intent(in)          :: X(:,:), Y(:,:)
    complex(KIND=dp), intent(in)          :: F(:,:,:)
    complex(KIND=dp), intent(in)          :: omega
    complex(KIND=dp), allocatable         :: Fout(:,:,:)

   ! allocate 
   allocate(Fout(nwt,nwt,2))


    print *, '--------------- Test L * Linv ---------------'

    print * , '||Fin||^2 = ', sum(abs(F))

    call Multiply_XY_with_QRPAmat(X, Y, omega, Fout)

    print * , '||Fout||^2 = ', sum(abs(Fout))
    print * , '||(Fin-Fout)||^2 = ',  sum(abs(F - Fout))

  end subroutine test_L_Linv

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

  subroutine test_qptransfo(ifail)
    !---------------------------------------------------------------------------
    ! Perform several tests for the SP->QP and QP-> SP transformation routines.
    !
    ! 1: Compare the transformation of N to existing routines
    ! 2: Compare the transformation of H to naive [W^{\dagger} H W]
    !     Only active when time-reversal is broken!
    ! 3: Unitarity - (SP->QP)(QP-SP) is identity (for H)
    ! 4: Symmetry properties of the resulting SPME and QPME
    ! 5: Unitarity - (SP->QP)(QP-SP) is identity for arbitrary (fermionic)
    !       operator matrix (with complex matrix elements)
    !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Potential upgrades
    !  - extend check 3 to an ARBITRARY one-body operator with random me
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !   None 
    ! Output:
    !   ifail: integer, 0 => success for any of the tests
    !                   1 => failure 
    !---------------------------------------------------------------------------
    integer, intent(out) :: ifail

    complex(KIND=dp), allocatable :: identity(:,:)
    complex(KIND=dp), allocatable :: Rsp(:,:,:), Rspback(:,:,:), Rqp(:,:,:),  Rph(:,:,:)

    real(KIND=dp), allocatable    :: N20(:,:), N11(:,:), H20(:,:), H11(:,:), H02(:,:)
    complex(KIND=dp), allocatable :: N20new(:,:), N11new(:,:), H20new(:,:), H11new(:,:), H02new(:,:)
    complex(KIND=dp), allocatable :: H20back(:,:), H11back(:,:), H02back(:,:)
$NTR  real(KIND=dp), allocatable    :: Hsp(:,:), Hqp_explicit(:,:)
    complex(KIND=dp), allocatable :: Hqp(:,:,:)

    ! Allocate random matrices and temporary arrays for random numbers
    complex(KIND=dp), allocatable :: O20sp(:,:), O11sp(:,:), O02sp(:,:)
    complex(KIND=dp), allocatable :: O20qp(:,:), O11qp(:,:), O02qp(:,:)
    complex(KIND=dp), allocatable :: O20sp_back(:,:), O11sp_back(:,:), O02sp_back(:,:)
    real(KIND=dp), allocatable :: rand_real(:,:), rand_imag(:,:)

    real(KIND=dp) :: ME_check(5)

    integer :: i, si, T, sb, N, N2, B
   
  1 format (100('-'))
  2 format ('Test ', i2, ' : ', a40)
 99 format (25('-'), ' Transformation routines (sp -> qp and qp -> sp) ', 25('-') )

    ifail  = 0

    print 99
    if(pairingtype==0) then
      print *, 'PairingType = 0 => HF :: Cannot perform QPtransfo test.'
    endif

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! 1. Compare the transformation of the particle number operator N to
    !        existing routines
    print 2, 1, ' Number operator N '
    print *

    ! construct the spme of the particle-number operator = identity
    allocate(identity(nwt,nwt)) 
    identity = 0.0

    ! Set the diagonal elements to 1
    do i = 1, nwt
      identity(i, i) = 1.0
    end do

    allocate(N20(nwt,nwt)   , N11(nwt,nwt))    ; N20    = 0.0d0 ; N11    = 0.0d0
    allocate(N20new(nwt,nwt), N11new(nwt,nwt)) ; N20new = 0.0d0 ; N11new = 0.0d0

    ! transform to qpme N20
    call transform_sp_to_qp(Bogoliubov, OTLsp=identity, OTRqp=N20new, OTLqp=N11new)

    ! get qpme of N20 from existing routine, for the neutron channel
    N20(    1:nwn,    1:nwn) = calcN20(Bogoliubov(      1:2*nwn,      1:2*nwn), HFblocks(1:4))
    N20(nwn+1:nwt,nwn+1:nwt) = calcN20(Bogoliubov(2*nwn+1:2*nwt,2*nwn+1:2*nwt), HFblocks(5:8))

    ! get qpme of N11 from existing routine, for the neutron channel
    N11(    1:nwn,    1:nwn) = calcN11(Bogoliubov(      1:2*nwn,      1:2*nwn), HFblocks(1:4))
    N11(nwn+1:nwt,nwn+1:nwt) = calcN11(Bogoliubov(2*nwn+1:2*nwt,2*nwn+1:2*nwt), HFblocks(5:8))

    print '(a50, 2es15.4)', ' ||N20(existing routine)||, || N20(new)|| = ', sum(abs(N20)), sum(abs(N20new))
    print '(a50, 2es15.4)', ' ||N11(existing routine)||, || N11(new)|| = ', sum(abs(N20)), sum(abs(N20new))
    print '(a50, es15.4)', ' ||N20(existing routine) - N20(new)|| = ', sum(abs(N20 - N20new))
    print '(a50, es15.4)', ' ||N11(existing routine) - N11(new)|| = ', sum(abs(N11 - N11new))

    if(sum(abs(N20 - N20new)) > 1d-10 .or. sum(abs(N11 - N11new)) > 1d-10) then
      print *, 'FAILURE!'
      ifail = 1
    else
      print *, 'SUCCESS!'
    endif
    print 1
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! 2. Compare the matrix multiplication of the HFB Hamiltonian with the
    print 2, 2, ' HFB Hamiltonian H'
$TR     print *, 'Test not valid when T is conserved'
    allocate(H20(nwt,nwt)    , H11(nwt,nwt))    ; H20    = 0.0d0 ; H11    = 0.0d0
    allocate(H20new(nwt,nwt) , H11new(nwt,nwt)) ; H20new = 0.0d0 ; H11new = 0.0d0
    allocate(H02(nwt,nwt)    , H02new(nwt,nwt)) ; H02new = 0.0d0
$NTR    allocate(Hsp(2*nwt,2*nwt), Hqp_explicit(2*nwt,2*nwt)) ; Hsp = 0.0d0

    ! transform to the quasiparticle representation of H
    call transform_sp_to_qp(Bogoliubov, OTLsp=dcmplx(sphamil), &
    &                                  OTRsp=dcmplx(HFBGaps), OBLsp=-dcmplx(HFBGaps), &
    &                                  OTRqp=H20new, OTLqp=H11new, OBLqp=H02new)

$NTR    ! Construct the HFB Hamiltonian in the sp space
$NTR    ! Note: we DO NOT use calcH20 and calcH11; these rely on a different memory layout
$NTR    Hsp = 0.0

$NTR    sb = 0; si = 0
$NTR    do B=1,8,2
$NTR      N  = HFBlocks(B)
$NTR      N2 = HFBlocks(B+1)
$NTR      T  = N + N2

$NTR      Hsp(sb  +1:sb+  T,sb  +1:sb+  T) =  sphamil(si+1:si+T,si+1:si+T)
$NTR      Hsp(sb+T+1:sb+2*T,sb+T+1:sb+2*T) = -sphamil(si+1:si+T,si+1:si+T)
$NTR      Hsp(sb+T+1:sb+2*T,sb  +1:sb+  T) = -HFBGaps(si+1:si+T,si+1:si+T)
$NTR      Hsp(sb  +1:sb+  T,sb+T+1:sb+2*T) = +HFBGaps(si+1:si+T,si+1:si+T)

$NTR      si = si +   T
$NTR      sb = sb + 2*T
$NTR    enddo
$NTR    ! Actual brute matrix multiplication : this is not correct when T is conserved
$NTR    ! as the array Bogoliubov does not contain the whole transformation
$NTR    Hqp_explicit = matmul(transpose(Bogoliubov), matmul(Hsp,Bogoliubov))
$NTR    sb = 0; si = 0
$NTR    do B=1,8,2
$NTR      N  = HFBlocks(B)
$NTR      N2 = HFBlocks(B+1)
$NTR      T  = N + N2

$NTR      H11(si+1:si+T, si+1:si+T) = Hqp_explicit(sb+1:sb+T,sb  +1:sb+  T)
$NTR      H20(si+1:si+T, si+1:si+T) = Hqp_explicit(sb+1:sb+T,sb+T+1:sb+2*T)
$NTR      H02(si+1:si+T, si+1:si+T) = Hqp_explicit(sb+T+1:sb+2*T,sb+1:sb+T)

$NTR      si = si +   T
$NTR      sb = sb + 2*T
$NTR    enddo

$NTR    print '(a50, 2es15.4)', ' ||H20(matmul)||   || H20(new)|| =', sum(abs(H20)), sum(abs(H20new))
$NTR    print '(a50, es15.4)',  ' ||H20(matmul)|| - || H20(new)|| =', sum(abs(H20)) - sum(abs(H20new))
$NTR    print '(a50, 2es15.4)', ' ||H11(matmul)||   || H11(new)|| =', sum(abs(H11)), sum(abs(H11new))
$NTR    print '(a50, es15.4)',  ' ||H11(matmul)|| - || H11(new)|| =', sum(abs(H11)) - sum(abs(H11new))
$NTR    print '(a50, 2es15.4)', ' ||H02(matmul)||   || H02(new)|| =', sum(abs(H02)), sum(abs(H02new))
$NTR    print '(a50, es15.4)',  ' ||H02(matmul)|| - || H02(new)|| =', sum(abs(H02)) - sum(abs(H02new))
$NTR    print *
$NTR    print *, 'Note : we do not check || Hmn - Hmnnew || because the QP reordering is not trivial.'
$NTR    print *, 'This means that this test is not sensitive to a global sign.'
$NTR    print *

    if(sum(abs(H20)) - sum(abs(H20new)) > 1d-10 .or. sum(abs(H11)) - sum(abs(H11new)) > 1d-10) then
      print *, 'FAILURE!'
      ifail = 1
    else
      print *, 'SUCCESS!'
    endif
    print 1
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! 3. Unitarity test for H
    print 2, 3, ' Unitarity of (SP->QP)(QP->SP) for H'
    print *
    allocate(H20back(nwt,nwt), H11back(nwt,nwt), H02back(nwt,nwt))

    call transform_qp_to_sp(Bogoliubov, OTLqp=H11new , OTRqp=H20new , OBLqp= H02new, &
    &                                   OTLsp=H11back, OTRsp=H20back, OBLsp= H02back)

    print *, ' ||   h   - h     (sp->qp->sp)|| = ',  sum(abs( sphamil - H11back))
    print *, ' || Delta - Delta (sp->qp->sp)|| = ',  sum(abs( HFBgaps - H20back))

    if(      sum(abs(sphamil - H11back))  > 1d-10 &
    &   .or. sum(abs(HFBgaps - H20back))  > 1d-10) then
      print *, 'FAILURE!'
      ifail = 1
    else
      print *, 'SUCCESS!'
    endif
    print 1
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! 5. test symmetry properties of the qpme
    print 2, 4, ' Symmetry properties of the SP and QP me'
    print *
    print * , '    SP ME    '
    print *,  ' ----------- '
    ! H20_sp and H02_sp should be related; H11 should be symmetric
    ME_check(1) = sum(abs(H11back - transpose(H11back)))
    ME_check(2) = sum(abs( H20back + H02back))

    print '(a50, 2es15.4)', ' || H11sp (sp->qp->sp) - H11sp^T (sp->qp->sp) || = ', ME_check(1)
    print '(a50, 2es15.4)', ' || H20sp (sp->qp->sp) + H02sp   (sp->qp->sp) || = ',  ME_check(2)

    print *
    print * , '    QP ME    '
    print *,  ' ----------- '

    ! H20_qp and H02_qp should be (anti)symmetric if T is conserved (broken)
    ME_check(3) = sum(abs(H11new - transpose(H11new)))
$NTR ME_check(4) = sum(abs( H20new + transpose(H20new)))
$TR  ME_check(4) = sum(abs( H20new - transpose(H20new)))
$NTR ME_check(5) = sum(abs( H02new + transpose(H02new)))
$TR  ME_check(5) = sum(abs( H02new - transpose(H02new)))

    print '(a50, 2es15.4)', '    H11_ab = +H11_ba   : satisfied up to', ME_check(3)
$NTR  print '(a50, 2es15.4)', '    H20_ab = -H20_ba   : satisfied up to',  ME_check(4)
$TR   print '(a50, 2es15.4)', '    H20_ab = +H20_ba   : satisfied up to',  ME_check(4)
$NTR  print '(a50, 2es15.4)', '    H02_ab = -H02_ba   : satisfied up to',  ME_check(5)
$TR   print '(a50, 2es15.4)', '    H02_ab = +H02_ba   : satisfied up to',  ME_check(5)

    if(any(ME_check > 1d-10)) then
      print *, 'FAILURE!'
      ifail = 1
    else
      print *, 'SUCCESS!'
    endif
    print 1

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! 5. Test that transform_sp_to_qp_wr and transform_qp_to_sp_wr are inverses
    print 2, 5, ' Inverse property: (QP->SP)(SP->QP) = identity'
    print *,    ' ... for (fermionic) random matrix'
    print *

    allocate(O20sp(nwt,nwt), O11sp(nwt,nwt), O02sp(nwt,nwt))
    allocate(O20qp(nwt,nwt), O11qp(nwt,nwt), O02qp(nwt,nwt))
    allocate(O20sp_back(nwt,nwt), O11sp_back(nwt,nwt), O02sp_back(nwt,nwt))

    ! Initialize to zero
    O20sp      = 0.0_dp; O11sp      = 0.0_dp; O02sp      = 0.0_dp
    O20sp_back = 0.0_dp; O11sp_back = 0.0_dp; O02sp_back = 0.0_dp

    ! Seed the random number generator
    call random_seed()

    ! Fill matrices with random complex values, block by block
    ! The transformation routines work block-by-block based on HFBlocks
    si = 0
    do B=1,8,2
      N = HFBlocks(B); if(N.eq.0) cycle
      N2 = HFBlocks(B+1)
      T = N + N2

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      ! Fill blocks of O11sp
      allocate(rand_real(N,N), rand_imag(N,N))
      call random_number(rand_real); call random_number(rand_imag)
      O11sp(si  +1:si+N,si+1:si+N) = dcmplx(rand_real, rand_imag)
      deallocate(rand_real, rand_imag)

      if(N2 .ne. 0) then
        allocate(rand_real(N2,N2), rand_imag(N2,N2))
        call random_number(rand_real); call random_number(rand_imag)
        O11sp(si+N+1:si+T,si+N+1:si+T) = dcmplx(rand_real, rand_imag)
        deallocate(rand_real, rand_imag)
      endif 
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      ! Fill blocks of O20sp - but respect fermionic antisymmetry
      allocate(rand_real(T,T), rand_imag(T,T))
      call random_number(rand_real); call random_number(rand_imag)
      if(N2 .ne. 0) then
        O20sp(si  +1:si+N,si+N+1:si+T) = dcmplx(rand_real(1:N,N+1:T), rand_imag(1:N,N+1:T)) 
        O20sp(si+N+1:si+T,si  +1:si+N) = - transpose(O20sp(si  +1:si+N,si+N+1:si+T))
      else
        O20sp(si  +1:si+N,si  +1:si+N) = dcmplx(rand_real, rand_imag) 
        O20sp(si  +1:si+N,si  +1:si+N) = O20sp(si  +1:si+N,si  +1:si+N) + transpose(O20sp(si  +1:si+N,si  +1:si+N))
      endif
      deallocate(rand_real, rand_imag)
      
      ! Fill blocks of O02sp - but respect fermionic antisymmetry
      allocate(rand_real(T,T), rand_imag(T,T))
      call random_number(rand_real);  call random_number(rand_imag)
      if(N2.ne.0) then
        O02sp(si  +1:si+N,si+N+1:si+T) = dcmplx(rand_real(1:N,N+1:T), rand_imag(1:N,N+1:T)) 
        O02sp(si+N+1:si+T,si  +1:si+N) = - transpose(O02sp(si+1:si+N,si+N+1:si+T))
      else 
        O02sp(si  +1:si+N,si  +1:si+N) = dcmplx(rand_real, rand_imag)
        O02sp(si  +1:si+N,si  +1:si+N) = O02sp(si  +1:si+N,si  +1:si+N) + transpose(O02sp(si  +1:si+N,si  +1:si+N))
      endif
      deallocate(rand_real, rand_imag)

      si = si + T
    enddo

    ! Transform to QP space
    call transform_sp_to_qp(Bogoliubov, OTRsp=O20sp, OTLsp=O11sp, OBLsp=O02sp, &
                           &            OTRqp=O20qp, OTLqp=O11qp, OBLqp=O02qp)

    !si = 0
    !do B=1,2,2
    !  N = HFBlocks(B); if(N.eq.0) cycle
    !  N2 = HFBlocks(B+1)
    !  T = N + N2!!

    !  print *, 'B = ', B , 'O20qp'
    !  do i=1,T
    !    print "(*( '(',g12.2,',',g12.2,')',:))", O20qp(si+i,si+1:si+T)
    !  enddo
    !  si = si + T
    !enddo 

    ! Transform back to SP space
    call transform_qp_to_sp(Bogoliubov, OTRqp=O20qp,      OTLqp=O11qp,      OBLqp=O02qp, &
                           &            OTRsp=O20sp_back, OTLsp=O11sp_back, OBLsp=O02sp_back)

    ! Check the differences
    print '("Maxval(|O20sp|)",es15.4)', maxval(abs(O20sp))
    print '("Maxval(|O11sp|)",es15.4)', maxval(abs(O11sp))
    print '("Maxval(|O02sp|)",es15.4)', maxval(abs(O02sp))
    print '(a50, es15.4)', ' || O20sp - (qp->sp->qp) O20sp || = ', sum(abs(O20sp - O20sp_back))
    print '(a50, es15.4)', ' || O11sp - (qp->sp->qp) O11sp || = ', sum(abs(O11sp - O11sp_back))
    print '(a50, es15.4)', ' || O02sp - (qp->sp->qp) O02sp || = ', sum(abs(O02sp - O02sp_back))

    !si = 0
    !do B=1,2,2
    !  N = HFBlocks(B); if(N.eq.0) cycle
    !  N2 = HFBlocks(B+1)
    !  T = N + N2!
!
!      print *, 'DEV', sum(abs(O20sp(si+1:si+T,si+1:si+T) - O20sp_back(si+1:si+T,si+1:si+T)))
!      print *, 'B = ', B , 'O20sp'
!      do i=1,T
!        print "(*( '(',g12.2,',',g12.2,')',:))", O20sp(si+i,si+1:si+T)
!      enddo
!      print *, 'B = ', B , 'O20sp_back'
!      do i=1,T
!        print "(*( '(',g12.2,',',g12.2,')',:))", O20sp_back(si+i,si+1:si+T)
!      enddo
!      print *, 'B = ', B , 'diff'
!      do i=1,T
!        print "(*( '(',g12.2,',',g12.2,')',:))", O20sp(si+i,si+1:si+T)- O20sp_back(si+i,si+1:si+T)
!      enddo
!      print *
!      si = si + T
!    enddo 

    if(    sum(abs(O20sp - O20sp_back)) > 1d-10 &
    & .or. sum(abs(O11sp - O11sp_back)) > 1d-10 &
    & .or. sum(abs(O02sp - O02sp_back)) > 1d-10) then
      print *, 'FAILURE!'
      ifail = 1
    else
      print *, 'SUCCESS!'
    endif
    print 1

    ! Clean up
    deallocate(O20sp, O11sp, O02sp, O20qp, O11qp, O02qp, O20sp_back, O11sp_back, O02sp_back)
  end subroutine test_qptransfo

  subroutine test_linearity_T(ifail)
    !--------------------------------------------------------------------------------------
    ! Test a complete FAM iteration is an affine transformation by calling iterate_dH on 
    ! a chosen linear combination a * dHa + b*dHb. One expects that
    ! FAM(a * dHa + b *dHb) = T(a * dHa + b * dHb) + dH_free
    !                       = a * T(dHa) + b * T(dHb) + dH_free
    ! such that
    ! FAM(a * dHa + b *dHb) - dH_free = a * (FAM(dHa) - dH_free) + b * (FAM(dHb) - dH_free)
    !
    !
    ! OUTPUT:
    !   ifail : integer
    !           0 if succesful, 1 if failed
    !--------------------------------------------------------------------------------------
    implicit none
    complex(KIND=dp), allocatable :: dHa(:), dHa_iter(:)
    complex(KIND=dp), allocatable :: dHb(:), dHb_iter(:)
    complex(KIND=dp), allocatable :: dHlincomb(:), dHlincomb_iter(:)
    real(KIND=dp),    allocatable :: rand_real(:), rand_imag(:)
    complex(KIND=dp), allocatable :: dH_free(:)
    complex(KIND=dp)              :: a, b
    integer                       :: i,j, allocsize
    integer, intent(out)          :: ifail 
    real                          :: diff_from_lin

  1 format (100('-'))
 99 format (25('-'), ' Linearity test of FAM iteration ', 25('-') )

    ifail = 0
    print 99
    print *

    if(pairingtype.ne.2) then
      allocsize = nwt*nwt
    else 
      allocsize = 3*nwt*nwt ! QRPA matrices are larger
    endif 
    allocate(dHa(allocsize))
    allocate(dHa_iter(allocsize))
    allocate(dHb(allocsize))
    allocate(dHb_iter(allocsize))
    allocate(dHlincomb(allocsize))
    allocate(dHlincomb_iter(allocsize))
    allocate(rand_real(allocsize))
    allocate(rand_imag(allocsize))
    allocate(dH_free(allocsize))
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
    dHa = DCMPLX(rand_real, rand_imag)
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! define dHb as random dcmplx matrix

    ! Generate random real and imaginary parts
    call random_number(rand_real)
    call random_number(rand_imag)

    ! Combine into complex matrix
    dHb = DCMPLX(rand_real, rand_imag)
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! dHlincomb as some linear combination of dHa and dHb
    
    a = dcmplx(6.0_dp, -1.0_dp) ! some random C scalar
    b = dcmplx(-2.0_dp, 0.5_dp) ! some random C scalar
    
    dHlincomb = a * dHa + b * dHb
    
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! calculate  FAM(dHa), FAM(dHb) and FAM(a * dHa + b *dHb)

    call iterate_dHsp(dHa, dHa_iter)
    call iterate_dHsp(dHb, dHb_iter)
    call iterate_dHsp(dHlincomb, dHlincomb_iter)

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! subtract the free response 

    dHa_iter = dHa_iter - dH_free
    dHb_iter = dHb_iter - dH_free
    dHlincomb_iter = dHlincomb_iter - dH_free

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! verify linarity of T and hence affinity of FAM
    diff_from_lin = norm_2(dHlincomb_iter - (a * dHa_iter + b * dHb_iter))
    print ('(a80, es15.4)'), '|| [FAM(dHlincomb) - free] - [a(FAM(dHa) - free) + b(FAM(dHb) - free)]|| = ', diff_from_lin
    print ('(a40, es15.4)'), 'FAM iteration is affine up to precision: ', diff_from_lin

    if(diff_from_lin > 1d-10) ifail = 1
    print 1

  end subroutine test_linearity_T


  subroutine test_RP_commutator()
    ! -----------------------------------------------------------------------
    ! In QP basis, the VEV of a commutator of two one-body operators is 
    !    <[A,B]> = <AB>_connected
    !            = 1/2 * sum_ab ( A^20_ab B^02_ab - A^20_ba B^02_ba )
    ! 
    ! For A=R_com and B=P_com, one expects 
    !    <[R,P]> = i
    ! (excluding the factor hbar inside P) 
    ! 

    complex(KIND=DP) :: Rz_qpme(nwt, nwt, 2)
    complex(KIND=DP) :: Pz_qpme(nwt, nwt, 2)
    complex(KIND=DP) :: Rz_spme(nwt, nwt)
    complex(KIND=DP) :: Pz_spme(nwt, nwt)
    complex(KIND=DP) :: RP_comm_spme(nwt, nwt), RP_comm_spme_exp(nwt,nwt)
    real(KIND=dp)    :: nabla_spme(3,2,nwt,nwt)
    complex(KIND=DP) :: commut0B
    complex(KIND=DP) :: term1RP, term1PR, term2RP, term2PR, term3RP, term3PR, term4RP, term4PR
    complex(KIND=DP) :: twobdens, tot
    integer          :: si, block, N, N2, i, j, k, T
    integer          :: a, b, c, d



    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! 1. evaluate [R, P] based on spme : 
    !
    ! <[R, P]> = sum_ab (sum_k R^11_ak P^11_kb - sum_k P^11_ak R^11_kb ) * rho_ba
    !          = Tr[  (R^11 P^11 - P^11 R^11) rho ]
    !
    ! Note that other terms in the two-body expectation value drop due to the commutator, 
    ! i.e. terms involving rho rho, kappa* kappa and rho rho* get canceled; terms appear 
    ! twice with opposite signs. For completeness, this is also verified.
    !

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! get the SPME for R

    ! get Q_10
    Rz_spme = Qlm_spme(1, 0, .false.)

    ! Convert Rz_spme to unit fm^l
    Rz_spme = Rz_spme * sqrt(100.0d0)

    ! Cancel prefactor sqrt(3/4pi)
    Rz_spme = Rz_spme * sqrt( 4.0d0 * pi / 3.0d0)

    ! multiply by 1/A
    Rz_spme = Rz_spme / (Neutrons + Protons)



    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! ... and then for P

    ! Get single-particle matrix elements of nabla in the HF basis
    nabla_spme = CompNablaMelements('HF')

    ! P_z = -i * nabla_z = IM(nabla_z) - Re(nabla_z) i
    Pz_spme = dcmplx(nabla_spme(3,2,:,:), -nabla_spme(3,1,:,:))


    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! print the symmetry properties of the matrices for debugging


    ! print * , '||Re(∇_z)||', sum(abs(nabla_spme(3,1,:,:)))
    ! print * , '||Im(∇_z)||', sum(abs(nabla_spme(3,2,:,:)))
    ! print * , '||∇_z||', sum(abs(dcmplx(nabla_spme(3,1,:,:), nabla_spme(3,2,:,:))))


    print * , 'R11 spme'
    call print_matrix_symmetries(Rz_spme)

    print * , '∇_z spme'
    call print_matrix_symmetries(dcmplx(nabla_spme(3,1,:,:), nabla_spme(3,2,:,:)))

    print * , 'P11 spme'
    call print_matrix_symmetries(Pz_spme)

    print * , 'rho'
    call print_matrix_symmetries(rho_pairing)

    print * , 'kappa'
    call print_matrix_symmetries(kappa_pairing)


    ! SP evaluation

    print * , 'Evaluation from spme : [R,P] = Tr[ (R^11 P^11 - P^11 R^11) rho ]' 

    term1RP = 0
    term1PR = 0
    
    do a = 1,nwt
      do b = 1,nwt
        do c = 1,nwt 
          term1RP = term1RP + Rz_spme(a,b) * Pz_spme(b,c) * rho_pairing(c, a)
          term1PR = term1PR - Pz_spme(a,b) * Rz_spme(b,c) * rho_pairing(c, a)
        enddo
      enddo
    enddo

    $TR term1RP = term1RP * 2.0d0
    $TR term1PR = term1PR * 2.0d0

    print * , '        R P rho :', term1RP 
    print * , '       -P R rho :', term1PR
    print * , '     naive    [R,P]  =', term1RP + term1PR

    print * 


    RP_comm_spme = matmul(Rz_spme, Pz_spme) - matmul(Pz_spme, Rz_spme)

    ! RP_comm_spme = matmul(RP_comm_spme, dcmplx(rho_pairing, 0_dp) )

    commut0B = 0
    do a = 1,nwt 
      do b = 1,nwt 
        commut0B = commut0B + RP_comm_spme(a,b) * rho_pairing(b, a)
      enddo
    enddo

    $TR commut0B = commut0B * 2.0d0


    print * , '     matmul   [R,P]  =', commut0B
    print * 


    print * , 'Evaluation from spme : [R,P] = Tr [(R^11 P^11 - P^11 R^11)*( rho rho  + kappa^* kappa + rho delta - rho rho^* ]'    


    term1RP = 0
    term1PR = 0
    term2RP = 0
    term2PR = 0
    term3RP = 0
    term3PR = 0
    term4RP = 0
    term4PR = 0
    
    do a = 1,nwt
      do b = 1,nwt
        do c = 1,nwt
          do d = 1,nwt
          
          term1RP = term1RP + Rz_spme(a,b) * Pz_spme(c,d) * rho_pairing(b,a) * rho_pairing(d,c)
          term1PR = term1PR - Pz_spme(a,b) * Rz_spme(c,d) * rho_pairing(b,a) * rho_pairing(d,c)

          term2RP = term2RP + Rz_spme(a,b) * Pz_spme(c,d) * kappa_pairing(c,a) * kappa_pairing(d,b)
          term2PR = term2PR - Pz_spme(a,b) * Rz_spme(c,d) * kappa_pairing(c,a) * kappa_pairing(d,b)
          
          term4RP = term4RP - Rz_spme(a,b) * Pz_spme(c,d) * rho_pairing(d,a) * rho_pairing(c,b)
          term4PR = term4PR + Pz_spme(a,b) * Rz_spme(c,d) * rho_pairing(d,a) * rho_pairing(c,b)
          
          enddo

          term3RP = term3RP + Rz_spme(a,c) * Pz_spme(c,b) * rho_pairing(b,a)
          term3PR = term3PR - Pz_spme(a,c) * Rz_spme(c,b) * rho_pairing(b,a)

        enddo
      enddo
    enddo

    $TR term1RP = term1RP * 2.0d0
    $TR term1PR = term1PR * 2.0d0
    $TR term2RP = term2RP * 2.0d0
    $TR term2PR = term2PR * 2.0d0
    $TR term3RP = term3RP * 2.0d0
    $TR term3PR = term3PR * 2.0d0
    $TR term4RP = term4RP * 2.0d0
    $TR term4PR = term4PR * 2.0d0


    print * , '+ R P rho rho       : ', term1RP 
    print * , '- P R rho rho       : ', term1PR
    print * , '+ R P kappa^* kappa : ', term2RP 
    print * , '- P R kappa^* kappa : ', term2PR
    print * , '+ R P rho delta     : ', term3RP 
    print * , '- P R rho delta     : ', term3PR
    print * , '+ R P rho rho^*     : ', term4RP 
    print * , '- P R rho rho^*     : ', term4PR
    print * 


    print * , '     [R,P]  =', term1RP + term1PR + term2RP + term2PR + term3RP + term3PR + term4RP + term4PR

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! 2. evaluate [R, P] based on qpme : 
    !
    ! <[R, P]> = 1/2 sum_ab ( R^20_ab P^02_ab - P^20_ab R^02_ab )
    !

    print *
    print *
    print * , 'Evaluation from qpme : [R,P] = 0.5 * sum_ab ( R^20_ab P^02_ab - P^20_ab R^02_ab )'    

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! get the QPME for R and P

    Rz_qpme = get_external_field('Zcom')

    Pz_qpme = get_external_field('Zmomentum')


    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! print some QPME symmetries

    print * , 'R20 qpme'
    call print_matrix_symmetries(Rz_qpme(:,:,1))

    print * , 'R02 qpme'
    call print_matrix_symmetries(Rz_qpme(:,:,2))

    print *, '||R20 - R02^*|| = ', sum(abs(Rz_qpme(:,:,1)-(conjg(Rz_qpme(:,:,2)))))

    print * , 'P20 qpme'
    call print_matrix_symmetries(Pz_qpme(:,:,1))

    print * , 'P02 qpme'
    call print_matrix_symmetries(Pz_qpme(:,:,2))

    print *, '||P20 - P02^*|| = ', sum(abs(Pz_qpme(:,:,1)-(conjg(Pz_qpme(:,:,2)))))

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! evaluate using naive loops

    term1RP = 0
    term1PR = 0

    do j = 1,nwt
      do i = 1,nwt

      term1RP = term1RP + Rz_qpme(i,j,1) * Pz_qpme(i,j,2)
      term1PR = term1PR - Pz_qpme(i,j,1) * Rz_qpme(i,j,2)

      enddo
    enddo

    term1RP = term1RP / 2.0d0
    term1PR = term1PR / 2.0d0

    $TR term1RP = term1RP * 2.0d0
    $TR term1PR = term1PR * 2.0d0

    print * 


    print * , '         0.5 R^20 P^02 :', term1RP 
    print * , '       - 0.5 P^20 R^02 :', term1PR
    print * , '     simple loops    [R,P]  =', term1RP + term1PR


    commut0B = dcmplx(0.0d0, imag(sum(Rz_qpme(:,:,1) * Pz_qpme(:,:,2))))
    $TR commut0B = commut0B * 2.0d0
    print * , ' imag(sum(R*P)) i    [R,P]  =', commut0B


    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! evaluate by looping over HFblocks
    !  -> THIS WILL NOT WORK WHEN PARITY IS CONSERVED SINCE P_z IS PARITY ODD
    commut0B = 0

    si = 0
    do block= 1,8,2
      N = HFBlocks(block); N2 = HFBlocks(block+1); T = N + N2; if(T.eq.0) cycle
      do j = si+1,si+T
        do i = si+1,si+T

        commut0B = commut0B + Rz_qpme(i,j,1) * Pz_qpme(i,j,2) - Pz_qpme(i,j,1) * Rz_qpme(i,j,2)

        enddo
      enddo
      si = si+T
    enddo

    commut0B = commut0B / 2.0d0

    $TR commut0B = commut0B * 2.0d0

    print * , '     via HFblocks    [R,P]  =', commut0B
    print * 

    call stp('end of test')

  end subroutine test_RP_commutator





  subroutine print_matrix_symmetries_complex(M)
    complex(kind=dp), intent(in) :: M (:,:)

    print *, 'This matrix is '

    if (sum(abs(M-conjg(M))) < 1d-10) print *, '    real'
    if (sum(abs(M+conjg(M))) < 1d-10) print *, '    imaginary'
    if (sum(abs(M-transpose(M))) < 1d-10) print *, '    symmetric'
    if (sum(abs(M+transpose(M))) < 1d-10) print *, '    antisymmetric'
    if (sum(abs(M-transpose(conjg(M)))) < 1d-10) print *, '    Hermitian'
    if (sum(abs(M+transpose(conjg(M)))) < 1d-10) print *, '    anti-Hermitian'

    print *

  end subroutine print_matrix_symmetries_complex

  subroutine print_matrix_symmetries_real(M)
    real(kind=dp), intent(in) :: M (:,:)

    call print_matrix_symmetries_complex(dcmplx(M ,0))

    print *

  end subroutine print_matrix_symmetries_real

end module fam_testing
