!===============================================================================
!     __  __  ___   ____ ____
!    |  \/  |/ _ \ / ___| ___|__ _
!    | |\/| | | | | |  | |_  / _` |
!    |  |\/| |_| | |__| |__| (_| |
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
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU Affero General Public License
!    along with this program.  If not, see <https://www.gnu.org/licenses/>.
!
!===============================================================================
module moments_testing
  !--------------------------------------------------------------------------------
  ! A selection of unit tests for the multipole-moment machinery of the
  ! moments module, in particular the symmetry-block matrix-element generator
  ! Qlm_spme_matrix and its dense expansion (design summary section 13).
  ! Styled after the fam_testing module.
  !
  ! The tests are armed with unit_test = .true. in the &mesh/ namelist: they
  ! act on the final mean-field state, and run_moments_tests ends the program
  ! afterwards - an armed flag will not result in a physics calculation.
  !
  ! Two tests for the multipole operators with l <= 4, all m and both the
  ! real and imaginary operator halves:
  !   (1) test_Qlm_dense: dump the dense representation of every Qlm matrix,
  !       via the print_dense helper of the block_matrix_module.
  !   (2) test_Qlm_trace: verify that taking the trace of the matrices with
  !       the density matrix - diagonal in the HF basis, with the canonical
  !       occupations rho_can on the diagonal (the convention of the cranking
  !       and momentsofinertia machinery) - reproduces the many-body
  !       multipole expectation values as calculated by CalculateMoments.
  !
  ! Both tests are exact identities up to floating-point reassociation: the
  ! matrix elements and the mesh moments are the same integrals, summed in
  ! a different order.
  !--------------------------------------------------------------------------------
  use compilation,        only: dp
  use moments,            only: Qlm_spme_matrix, FindMoment, MaxMoment
  use moments,            only: Moment
  use pairing,            only: rho_can, PairingType
  use wavefunctions,      only: nwn, nwt
  use block_matrix_module, only: sym_matrix

  implicit none

contains

  subroutine run_moments_tests()
    !------------------------------------------------------------------------------
    ! Catch-all routine to run all predefined unit tests for the moments module.
    !
    ! Notes
    !  - this routine assumes the mean-field run has been completed: it acts
    !    on the final HFPsi basis, occupations and moment values.
    !  - this routine contains an exit statement: arming the unit_test flag
    !    will not result in a physics calculation. The exit status allows
    !    this to be part of the MOCCa testing framework.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !   None
    ! Output:
    !   None
    !   Exit status = 0 => everything okay!
    !               = 1 => tests failed
    !------------------------------------------------------------------------------

    1 format ('Test = ', a40, ' Success = ', i4)
    9 format (' ---------------------- Start of the moments testing routines ---------------------')
    10 format (' ---------------------- End of the moments testing routines -----------------------')
    11 format (' Summary of results ')
    12 format (' All tests passed! ')
    13 format (' Some tests failed! ')

    integer :: ifail_dense, ifail_trace, ifail
    ! Maximum multipole order the tests act on: all l up to 4.
    integer, parameter :: lmax = 4

    print 9
    print *

    if (lmax .gt. MaxMoment) then
      print *, 'run_moments_tests: the build restricts moments to l =',      &
      &         MaxMoment, ', testing up to that order instead.'
    endif

    call test_Qlm_dense(min(lmax, MaxMoment), ifail_dense)
    print 1, 'Dense dump of the Qlm matrices'     , ifail_dense

    call test_Qlm_trace(min(lmax, MaxMoment), ifail_trace)
    print 1, 'Tr( rho_can Qlm ) vs moment value'   , ifail_trace

    print *
    print 10
    print 11

    ifail = max(ifail_dense, ifail_trace)
    if(ifail.eq.0) print 12
    if(ifail.eq.1) print 13
    call exit(ifail)

  end subroutine run_moments_tests

!===============================================================================
  subroutine test_Qlm_dense(lmax, ifail)
    !------------------------------------------------------------------------------
    ! Dump the dense representation of the single-particle multipole matrices
    ! for every multipole order l <= lmax, all magnetic quantum numbers
    ! 0 <= m <= l and both operator halves:
    !   Im = .false. : r^l Re Y_lm (stored numbers = Re < i | Q_lm | j >)
    !   Im = .true.  : r^l Im Y_lm (stored numbers = Im < i | Q_lm | j >)
    ! both in fm^l, through the print_dense helper of the block_matrix_module,
    ! which reconstructs the partner blocks and decodes the symmetry word.
    !
    ! This test is a visual-development tool: it prints the matrices and
    ! succeeds unless Qlm_spme_matrix fails loudly (stp) along the way.
    ! Intended for small bases only - print_dense materialises the full
    ! dense matrix.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !   lmax    : maximum multipole order to test
    ! Output:
    !   ifail   : 0 if all matrices could be built, 1 otherwise
    !------------------------------------------------------------------------------

    1 format (' ---------- Dense dump of the multipole matrices up to l =', i3, &
    &         ' ----------')

    integer, intent(in)  :: lmax
    integer, intent(out) :: ifail

    type(sym_matrix)      :: me
    character(len=80)     :: title
    logical               :: Im
    integer               :: l, m, ihalf

    print 1, lmax
    print *
    ifail = 0

    do l = 0, lmax
      do m = 0, l
        do ihalf = 1, 2
          Im = (ihalf .eq. 2)

          call Qlm_spme_matrix(l, m, Im, me)

          write(title, 11) l, m, merge('Im', 'Re', Im), l
          call me%print_dense(6, trim(title))
          print *
        enddo
      enddo
    enddo

    11 format ('Q_', i1, '_', i1, ' (', a2, ' Y half): stored half of the ',  &
    &           'matrix elements, fm^', i1)

  end subroutine test_Qlm_dense

!===============================================================================
  subroutine test_Qlm_trace(lmax, ifail)
    !------------------------------------------------------------------------------
    ! Verify that the trace of the multipole matrices with the density matrix
    ! of the many-body state reproduces the many-body multipole moment values
    ! as calculated in CalculateMoments:
    !
    !   Tr( rho Q_lm ) = sum_i rho_can(i) < i | Q_lm | i >   =   < Q_lm >
    !
    ! per isospin. The density matrix is diagonal in the HF basis, with the
    ! canonical occupations rho_can on the diagonal - the convention the
    ! cranking and momentsofinertia machinery use - and the isospin split
    ! follows the stored basis order: the first nwn spwfs are neutrons.
    !
    ! A symmetry-forbidden combination is either absent from the linked list
    ! of moments or carried in it with a vanishing value; absent moments are
    ! compared against the zero their selection rules impose on the matrix.
    !
    ! Both sides are the same mesh integrals in a different summation order
    ! (dense block products vs the density mesh sums), so the identity holds
    ! to floating-point reassociation; the tolerance is generous and the
    ! table is printed for inspection.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !   lmax    : maximum multipole order to test
    ! Output:
    !   ifail   : 0 if all moments traced correctly, 1 otherwise
    !------------------------------------------------------------------------------

    1 format (' --------- Trace of the multipole matrices against the occupation',          &
    &         ' density up to l =', i3, ' ---------')
    2 format (4x, ' l', 2x, ' m', 5x, 'half', 12x, 'Trace(n)', 14x, 'Value(n)',            &
    &         14x, 'Trace(p)', 14x, 'Value(p)')
    3 format (i6, i4, 5x, a2, 4es20.9)
    4 format (' largest deviation of the trace check : ', es12.3)
    5 format (' Test only valid for a density diagonal in the HF basis: skipped for HFB.')
    6 format (' rho_can is not available: has the mean-field run been completed?')

    integer, intent(in)  :: lmax
    integer, intent(out) :: ifail

    type(sym_matrix)       :: me
    type(Moment), pointer  :: pcur
    real(KIND=dp), allocatable :: D(:, :)
    real(KIND=dp)          :: tr(2), val(2), dev, tol, maxdev
    logical                :: Im
    integer                :: l, m, ihalf, i

    ifail = 0
    maxdev = 0.0d0

    ! The HFB density matrix is diagonal in the canonical basis, not in the
    !  HF basis these matrices live in: the trace identity is not tested.
    if (PairingType .eq. 2) then
      print 5
      return
    endif
    if (.not. allocated(rho_can)) then
      print 6
      ifail = 1
      return
    endif

    print 1, lmax
    print 2

    do l = 0, lmax
      do m = 0, l
        do ihalf = 1, 2
          Im = (ihalf .eq. 2)

          call Qlm_spme_matrix(l, m, Im, me)
          D = me%dense()

          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          ! The trace with the diagonal density: rho_can counts the spwf
          !  occupations in the stored HF-basis order, neutrons first.
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          tr = 0.0d0
          do i = 1, me%n
            if (i .le. nwn) then
              tr(1) = tr(1) + rho_can(i)*D(i, i)
            else
              tr(2) = tr(2) + rho_can(i)*D(i, i)
            endif
          enddo

          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          ! The many-body expectation value as calculated in moments.f90.
          !  A symmetry-forbidden moment is not in the linked list: its
          !  expectation value is consistently zero.
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          nullify(pcur)
          val = 0.0d0
          pcur => FindMoment(l, m, Im)
          if (associated(pcur)) val = pcur%Value
          nullify(pcur)

          dev = max( abs(tr(1) - val(1)), abs(tr(2) - val(2)) )
          tol = 1.0d-6*( 1.0d0 + abs(val(1)) + abs(val(2))                  &
          &                     + abs(tr(1)) + abs(tr(2)) )
          if (dev .gt. tol) ifail = 1
          maxdev = max(maxdev, dev/(1.0d0 + tol))

          print 3, l, m, merge('Im', 'Re', Im), tr(1), val(1), tr(2), val(2)
        enddo
      enddo
    enddo

    print 4, maxdev

  end subroutine test_Qlm_trace

end module
