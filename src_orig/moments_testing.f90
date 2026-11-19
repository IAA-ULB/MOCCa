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
  !
  ! One oracle, armed alongside the tests:
  !  (3) test_Qlm_fullbox: the full-box reconstruction oracle of the amended
  !      step 4 of design-summary section 10 - the reference path for the
  !      imaginary operator half. It recomputes selected matrix elements of
  !      both halves by independent mesh quadrature over the chamber-
  !      extended spinors, bypassing the chamber-weight table and the
  !      reflection-character estimation of the production fill. All three
  !      identities below are judged against physics, never against zeros.
  !--------------------------------------------------------------------------------
  use compilation,        only: dp
  use geninfo,            only: stp, nx, ny, nz, mv, dv, meshx_shifted,        &
  &                              meshy_shifted, meshz_shifted
  use geninfo,            only: reduX, reduY, reduZ
  use sphericalharmonics, only: generate_spherical_harmonics
  use moments,            only: Qlm_spme_matrix, FindMoment, MaxMoment
  use moments,            only: Moment
  use moments,            only: QuantisationAxis, SecondaryAxis
  use pairing,            only: rho_can, PairingType
  use wavefunctions,      only: nwn, nwt, HFPsi, sx, sy, sz
  use block_matrix_module, only: sym_matrix, sym_block

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

    integer :: ifail_dense, ifail_trace, ifail_fullbox, ifail
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

    call test_Qlm_fullbox(min(lmax, MaxMoment), ifail_fullbox)
    print 1, 'Full-box oracle vs stored Qlm'        , ifail_fullbox

    print *
    print 10
    print 11

    ifail = max(ifail_dense, ifail_trace, ifail_fullbox)
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

!===============================================================================
  subroutine test_Qlm_fullbox(lmax, ifail)
    !------------------------------------------------------------------------------
    ! Full-box reconstruction oracle for the single-particle multipole matrices
    ! (the amended step 4 of design-summary section 10): the reference path
    ! for the imaginary operator half, active for both halves.
    !
    ! For every multipole order l <= lmax, all m and both operator halves,
    ! selected matrix elements are recomputed INDEPENDENTLY of the Qlm_spme
    ! fill machinery: the spwfs are extended from the stored fundamental
    ! chamber to the full box with the per-component reflection signs
    ! sx/sy/sz (the folding convention of the derivative machinery), the
    ! operator mesh function is regenerated at the chamber coordinates, and
    ! the full complex matrix element < i | Q_lm | j > is obtained by plain
    ! mesh quadrature with the plain cell volume. The chamber-weight table
    ! and the numeric reflection-character estimation of the production
    ! fill are bypassed by construction; only the trusted raw inputs are
    ! shared.
    !
    ! Judged against physics, never against zeros, per element:
    !  (1) identity: the stored half agrees with the corresponding part
    !      (real/imaginary) of the full complex matrix element, at both the
    !      stored position and its partner position - the latter through the
    !      dense() partner reconstruction, so the storage-kind sign is
    !      validated along the way;
    !  (2) STy characteristic: the opposite part of the complex matrix
    !      element vanishes - purely real matrix elements for r^l Re Y_lm
    !      (STY_EVEN), purely imaginary ones for r^l Im Y_lm (STY_ODD);
    !  (3) Hermiticity: the reversed element is the complex conjugate.
    !
    ! Per stored block, the elements checked are the corner (1,1) - zero by
    ! construction for the antisymmetric halves, a genuine check through
    ! identity (1) - one element of the first row and of the first column,
    ! one hashed interior element and the largest-magnitude element, each
    ! in both orderings.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !   lmax    : maximum multipole order to test
    ! Output:
    !   ifail   : 0 if every checked element survived, 1 otherwise
    !------------------------------------------------------------------------------

    1 format (' --------- Full-box reconstruction oracle for the multipole',   &
    &         ' matrices up to l =', i3, ' ---------')
    2 format (4x, ' l', 2x, ' m', 5x, 'half', 3x, '#elem', 8x, 'identity',    &
    &         10x, 'STy half', 9x, 'Hermiticity')
    3 format (i6, i4, 5x, a2, i9, 4x, 3es14.5)
    4 format (' largest deviations of the full-box oracle -   identity:',     &
    &         es12.3, '   STy half:', es12.3, '   Hermiticity:', es12.3)
    5 format ('  fullbox oracle: ', a, ' mismatch for l =', i2, ' m =', i2,    &
    &         'half ', a2, ' block', i2, ' element (', i5, ',', i5, ')')
    6 format (i6, i4, 5x, a2, '     structurally empty matrix - skipped')

    integer, intent(in)  :: lmax
    integer, intent(out) :: ifail

    type(sym_matrix)          :: me
    real(KIND=dp), allocatable :: D(:, :)
    real(KIND=dp), allocatable :: Qc(:, :)
    integer, allocatable       :: flips(:, :)
    integer :: nch, l, m, ihalf, b, pk, npick, npairs, npairs_mom
    integer :: gi, gj
    integer, parameter :: PICK_MAX = 6
    integer :: ipick(2, PICK_MAX)
    real(KIND=dp) :: va_re, va_im, vb_re, vb_im
    real(KIND=dp) :: dev, tol, matmax
    real(KIND=dp) :: max_dev, max_sty, max_her, mom_dev, mom_sty, mom_her
    logical :: Im

    ifail  = 0
    npairs = 0
    max_dev = 0.0d0
    max_sty = 0.0d0
    max_her = 0.0d0

    print 1, lmax
    print 2

    do l = 0, lmax
      do m = 0, l
        do ihalf = 1, 2
          Im = (ihalf .eq. 2)

          call Qlm_spme_matrix(l, m, Im, me)
          if (me%nblocks .eq. 0) then
            print 6, l, m, merge('Im', 'Re', Im)
            cycle
          endif

          D  = me%dense()
          call fullbox_chambers(l, m, ihalf, nch, flips, Qc)

          matmax  = maxval(abs(D))
          dev     = 0.0d0
          npairs_mom = 0
          mom_dev = 0.0d0
          mom_sty = 0.0d0
          mom_her = 0.0d0

          do b = 1, me%nblocks
            call fullbox_block_picks(me%blocks(b), b, PICK_MAX, ipick, npick)
            do pk = 1, npick
              gi = me%blocks(b)%is + ipick(1, pk) - 1
              gj = me%blocks(b)%js + ipick(2, pk) - 1

              call fullbox_me(Qc, flips, gi, gj, va_re, va_im)
              call fullbox_me(Qc, flips, gj, gi, vb_re, vb_im)
              npairs = npairs + 1
              npairs_mom = npairs_mom + 1

              ! Tolerance: relative to the magnitudes at play, with a floor
              !  on the matrix scale so that structural zeros (diagonals of
              !  the antisymmetric halves, zero element picks) are judged
              !  against the noise of the quadrature rather than against
              !  exact zero.
              tol = 1.0d-9*(abs(va_re) + abs(va_im) + abs(vb_re) + abs(vb_im))&
              &   + 1.0d-10*matmax

              ! (1) identity against the stored half, both orderings - the
              !  partner side runs through the dense() reconstruction.
              ! (2) the STy-selected other part of the matrix element
              !  vanishes for both orderings.
              if (ihalf .eq. 1) then
                dev = abs(va_re - D(gi, gj)) + abs(vb_re - D(gj, gi))
                max_dev = max(max_dev, dev)
                mom_dev = max(mom_dev, dev)
                if (dev .gt. tol) then
                  print 5, 'stored-half', l, m, 'Re', b, gi, gj
                  ifail = 1
                endif
                dev = abs(va_im) + abs(vb_im)
                max_sty = max(max_sty, dev)
                mom_sty = max(mom_sty, dev)
                if (dev .gt. tol) then
                  print 5, 'STy-half  ', l, m, 'Re', b, gi, gj
                  ifail = 1
                endif
              else
                dev = abs(va_im - D(gi, gj)) + abs(vb_im - D(gj, gi))
                max_dev = max(max_dev, dev)
                mom_dev = max(mom_dev, dev)
                if (dev .gt. tol) then
                  print 5, 'stored-half', l, m, 'Im', b, gi, gj
                  ifail = 1
                endif
                dev = abs(va_re) + abs(vb_re)
                max_sty = max(max_sty, dev)
                mom_sty = max(mom_sty, dev)
                if (dev .gt. tol) then
                  print 5, 'STy-half  ', l, m, 'Im', b, gi, gj
                  ifail = 1
                endif
              endif

              ! (3) Hermiticity of the full complex matrix element:
              !      < j | Q | i > = conj( < i | Q | j > )
              dev = abs(va_re - vb_re) + abs(va_im + vb_im)
              max_her = max(max_her, dev)
              mom_her = max(mom_her, dev)
              if (dev .gt. tol) then
                print 5, 'Hermitian ', l, m, merge('Im', 'Re', Im), b, gi, gj
                ifail = 1
              endif
            enddo
          enddo

          print 3, l, m, merge('Im', 'Re', Im), npairs_mom, mom_dev,          &
          &        mom_sty, mom_her
        enddo
      enddo
    enddo

    print 4, max_dev, max_sty, max_her
    print *, ' number of matrix elements checked:', npairs

  end subroutine test_Qlm_fullbox

!===============================================================================
  subroutine fullbox_chambers(l, m, ihalf, nch, flips, Qc)
    !------------------------------------------------------------------------------
    ! Chamber setup of the full-box oracle: the operator mesh function
    ! evaluated freshly at the coordinates of EVERY chamber of the full box,
    ! with the per-chamber axis-flip masks and the plain cell volume folded
    ! in.
    !
    ! The stored mesh is the fundamental chamber: 2^nch folds with
    ! nch = reduX + reduY + reduZ, and dv carries the multiplicity, so the
    ! full-box cell volume is dv / 2^nch. For every chamber the spherical
    ! harmonics are regenerated with the mirrored coordinate arrays, so the
    ! reflection behaviour of the operator enters through the mesh values
    ! only.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !   l, m, ihalf : the operator slice r^l Re Y_lm (ihalf = 1) or
    !                 r^l Im Y_lm (ihalf = 2), conventions as Qlm_spme_matrix
    ! Output:
    !   nch         : number of chambers of the full box
    !   flips       : flips(1:3, 1:nch) - per-axis flip of each chamber
    !                 (0: stored mesh, 1: mirrored axis)
    !   Qc          : Qc(:, c) = operator values of chamber c, flattened over
    !                 the stored mesh points, in fm^l and multiplied by
    !                 dv / nch
    !------------------------------------------------------------------------------

    integer, intent(in)   :: l, m, ihalf
    integer, intent(out)  :: nch
    integer, allocatable, intent(out)  :: flips(:, :)
    real(KIND=dp), allocatable, intent(out) :: Qc(:, :)

    real(KIND=dp), allocatable :: SpherHarmMesh(:,:,:,:,:,:)
    integer :: c, ic, fx, fy, fz

    nch = 2**(reduX + reduY + reduZ)
    allocate(flips(3, nch), Qc(mv, nch))
    allocate(SpherHarmMesh(nx, ny, nz, 0:l, 0:l, 2))

    do c = 1, nch
      ! Chamber index decomposition over the folded axes only: axis bits of
      !  unfolded builds do not exist in the full box.
      ic = c - 1
      fx = 0
      fy = 0
      fz = 0
      if (reduX .eq. 1) then
        fx = mod(ic, 2)
        ic = ic/2
      endif
      if (reduY .eq. 1) then
        fy = mod(ic, 2)
        ic = ic/2
      endif
      if (reduZ .eq. 1) then
        fz = mod(ic, 2)
      endif
      flips(:, c) = (/ fx, fy, fz /)

      ! The operator values at the chamber coordinates: the COM-shifted mesh
      !  mirrors with the flipped axes only.
      call generate_spherical_harmonics(l, nx, ny, nz,                        &
      &       merge(-meshx_shifted, meshx_shifted, fx .eq. 1),                &
      &       merge(-meshy_shifted, meshy_shifted, fy .eq. 1),                &
      &       merge(-meshz_shifted, meshz_shifted, fz .eq. 1),                &
      &       SpherHarmMesh, QuantisationAxis, SecondaryAxis)
      Qc(:, c) = reshape(SpherHarmMesh(:,:,:,l,m,ihalf), (/ mv /))            &
      &        * (dv/nch)
    enddo

    deallocate(SpherHarmMesh)

  end subroutine fullbox_chambers

!===============================================================================
  subroutine fullbox_me(Qc, flips, gi, gj, me_re, me_im)
    !------------------------------------------------------------------------------
    ! One matrix element of the full-box oracle: the plain mesh quadrature of
    ! the full complex matrix element
    !
    !     < i | Q_lm | j > = sum_spin Int_fullbox d^3r  chi_i,s*(r) Q(r) chi_j,s(r)
    !
    ! over the chamber-extended spinors. At the full-box points of chamber c
    ! every spinor component k of state i takes the value of its mirror image
    ! in the stored chamber, multiplied by the product of the per-component
    ! reflection signs sx/sy/sz over the chamber's flipped axes (the folding
    ! convention of the derivative machinery).
    !
    ! With the real component storage of HFPsi - (2s-1) = Re and (2s) = Im of
    ! spin projection s - the summand reads
    !   Q (Re_i Re_j + Im_i Im_j) + i Q (Re_i Im_j - Im_i Re_j)
    ! per spin projection, the second bracket being the imaginary part that
    ! the STY_ODD operator half stores.
    !
    ! Qc carries the operator values per chamber with the full-box cell
    !  volume folded in, so the chamber sum below is the entire integral.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !   Qc, flips : the chamber state, as built by fullbox_chambers
    !   gi, gj    : global HF-basis indices of the row/column state
    ! Output:
    !   me_re, me_im : Re / Im of < i | Q_lm | j >
    !------------------------------------------------------------------------------

    real(KIND=dp), intent(in)  :: Qc(:, :)
    integer, intent(in)        :: flips(:, :)
    integer, intent(in)        :: gi, gj
    real(KIND=dp), intent(out) :: me_re, me_im

    integer :: c, k, s
    real(KIND=dp) :: si(4), sj(4)
    real(KIND=dp), allocatable :: Rei(:), Imi(:), Rej(:), Imj(:)

    if (gi .lt. 1 .or. gj .lt. 1 .or. gi .gt. nwt .or. gj .gt. nwt) then
      print *, 'fullbox_me: state index outside the basis: ', gi, gj
      call stp('State index outside the basis in fullbox_me.')
    endif

    allocate(Rei(mv), Imi(mv), Rej(mv), Imj(mv))
    me_re = 0.0d0
    me_im = 0.0d0

    do c = 1, size(Qc, 2)
      ! Extension signs of both states' four components for this chamber.
      si = 1.0d0
      sj = 1.0d0
      do k = 1, 4
        if (flips(1, c) .eq. 1) then
          si(k) = si(k)*sx(k, gi)
          sj(k) = sj(k)*sx(k, gj)
        endif
        if (flips(2, c) .eq. 1) then
          si(k) = si(k)*sy(k, gi)
          sj(k) = sj(k)*sy(k, gj)
        endif
        if (flips(3, c) .eq. 1) then
          si(k) = si(k)*sz(k, gi)
          sj(k) = sj(k)*sz(k, gj)
        endif
      enddo

      do s = 1, 2
        Rei = si(2*s-1)*HFPsi(:, 2*s-1, gi)
        Imi = si(2*s  )*HFPsi(:, 2*s  , gi)
        Rej = sj(2*s-1)*HFPsi(:, 2*s-1, gj)
        Imj = sj(2*s  )*HFPsi(:, 2*s  , gj)
        me_re = me_re + dot_product(Qc(:, c), Rei*Rej + Imi*Imj)
        me_im = me_im + dot_product(Qc(:, c), Rei*Imj - Imi*Rej)
      enddo
    enddo

    deallocate(Rei, Imi, Rej, Imj)

  end subroutine fullbox_me

!===============================================================================
  subroutine fullbox_block_picks(blk, b, pick_max, ipick, npick)
    !------------------------------------------------------------------------------
    ! Deterministic element selection of one stored block for the full-box
    !  oracle: the corner (1,1) - zero by construction for the antisymmetric
    !  halves - one element of the first row and of the first column, one
    !  hashed interior element and the largest-magnitude element.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !   blk       : one stored block of the matrix under test
    !   b         : its position in the block list, for the hashed interior
    !   pick_max  : size of the pick list handed in
    ! Output:
    !   ipick     : local (row, column) picks of the block, 1-based
    !   npick     : their number, <= pick_max
    !------------------------------------------------------------------------------

    type(sym_block), intent(in)  :: blk
    integer, intent(in)          :: b, pick_max
    integer, intent(out)        :: ipick(2, pick_max)
    integer, intent(out)        :: npick

    integer :: m1, n1, ij
    integer :: iloc(2)

    npick = 0
    m1 = size(blk%A, 1)
    n1 = size(blk%A, 2)

    npick = npick + 1
    ipick(1, npick) = 1
    ipick(2, npick) = 1
    if (n1 .ge. 2) then
      npick = npick + 1
      ipick(1, npick) = 1
      ipick(2, npick) = 2
    endif
    if (m1 .ge. 2) then
      npick = npick + 1
      ipick(1, npick) = 2
      ipick(2, npick) = 1
    endif
    npick = npick + 1
    ipick(1, npick) = 1 + mod(7*b + 3, m1)
    ipick(2, npick) = 1 + mod(11*b + 5, n1)

    ! Largest-magnitude element; maxloc returns the first maximum, so the
    !  pick is reproducible run to run.
    iloc = maxloc(abs(blk%A))
    ij = iloc(1)
    npick = npick + 1
    ipick(1, npick) = 1 + mod(ij - 1, m1)
    ipick(2, npick) = 1 + (ij - 1)/m1
    if (npick .gt. pick_max)                                                   &
    &  call stp('Pick list overflow in fullbox_block_picks.')

  end subroutine fullbox_block_picks

end module
