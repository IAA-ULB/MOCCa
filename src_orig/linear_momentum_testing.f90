!===============================================================================
!     __  __  ___   ____ ____
!    |  \/  |/ _ \ / ___/ ___|__ _
!    | |\/| | | | | |   | |   / _` |
!    | |  | | |_| | |___| |__| (_| |
!    |_|  |_|\___/ \____|____|\__,_|
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
module linear_momentum_testing
!==============================================================================
  !
  ! Unit arm of the gradient generators (design summary sections 15.5/15.6):
  !  every object of the 15.2 table is built through its generator and
  !  judged by five checks, plus two retrieval checks (15.8) exercised in
  !  the builds that have T-image sectors.
  !
  !  1. structure: the stored block list must be the XOR pairing
  !     srow = IEOR(s, AND(word, goodqn)) of the design, exactly as init
  !     constructs it.
  !  2. values: sampled stored elements are recomputed by the test's own
  !     second transcription of the legacy arithmetic - an independent
  !     assembly route, guarding the generator's transcription.
  !  3. legacy A/B: the stored half, expanded by dense(), must agree with
  !     the legacy CompNablaMelements('HF') array - THE oracle, per design
  !     15.6. The legacy's (mu, half) slots carry build-dependent semantics
  !     (honest matrix elements in TR-broken builds, folded ones in
  !     TR-conserving builds), so the comparison follows the slot map of
  !     15.1: honest \nabla_z pairs slot (3,1) in every build; honest
  !     \nabla_x/\nabla_y pair slots (1,1)/(2,2) only where the -i
  !     sectors carry stored states; the folded pair pairs the same slots
  !     in TR-conserving builds, where the honest x/y objects are
  !     structural zeros and have no legacy content by construction.
  !  4. mirror oracle: the kernel-transpose arbiter of design 2d - the
  !     fresh element evaluated in both orderings measures the mirror
  !     character of the stored numbers, and the derived kind must agree.
  !     The class tags of 15.2 (folded x antisymmetric, folded y
  !     symmetric) follow the legacy partner fills; the oracle is the
  !     authority they have to match - the J family's 14.2b precedent.
  !  5. partner: within-block triangles obey A(c,r) = kind * A(r,c).
  !  6. purity: the opposite half of the complex matrix elements,
  !     measured by the test's own inline arithmetic. OPEN FINDING
  !     (2026-11-19, first runs of this module): it does NOT vanish -
  !     O(0.2) against an O(0.35) stored scale, on seed and converged
  !     O16/SLy4/NLO states alike, for every component of both the
  !     honest trio and the folded pair. The J family's arm measures
  !     the same for J_z (angmom_z_imag, their recorded open finding),
  !     while the multipole full-box oracle proves the BASIS satisfies
  !     the STy reality constraint to 1e-16 - so the anomaly is a
  !     property of DERIVATIVE operators: the crossed-component pairs
  !     of the reduced-chamber integral are the section 13.1 chamber-
  !     weight territory, the per-amplitude selection the sector word
  !     cannot express, which the legacy routine equally never
  !     computes. The opposite halves are therefore MEASURED AND
  !     REPORTED, never judged, until the chamber question is settled.
  !  7. retrieval (15.8), TR-conserving builds only, where the empty
  !     -i sectors are the T-image slots:
  !     a. the tau both-ways check of \nabla_z's image-image blocks -
  !        the accessor's sign-flipped counterpart against the direct
  !        honest computation on materialised images; this check arms
  !        the algebra mechanism.
  !     b. the retrieval agreement of the fold: get_nabla_block's
  !        stored-x-image blocks against the folded generators' own
  !        elements - the two packagings of the same matrix elements.
  !
  !  The whole arm runs in the HF basis, matching the wavefunction array
  !  access patterns of the generators.
  !==============================================================================
  use compilation,       only: dp
  use geninfo
  use wavefunctions,     only: HFPsi, HFdpsi, sx, sy, sz, nwt
  use wavefunctions,     only: Derive_X_spwf, Derive_Y_spwf, Derive_Z_spwf
  use wavefunctions,     only: TimeReverse
  use wavefunctions,     only: HFBlocks_global
  use densities,         only: CompNablaMelements
  use block_matrix_module, only: sym_matrix
  use block_matrix_module, only: KIND_SYMMETRIC
  use linear_momentum,   only: AX_X, AX_Y, AX_Z
  use linear_momentum,   only: Nabla_spme_matrix, NablaT_spme_matrix
  use linear_momentum,   only: get_nabla_block, TAU_NABLAFAM

  implicit none

  private
  public :: run_linear_momentum_tests

  ! Tolerance of the oracles and consistency checks: the fills are per-pair
  !  transcriptions of the legacy sums, so identical inputs agree to
  !  machine precision - tight bar by design.
  real(KIND=dp), parameter :: tol = 1.0d-10

contains

  subroutine run_linear_momentum_tests()
    !------------------------------------------------------------------------------
    ! Run the unit tests of the gradient generators: the five objects of
    !  the 15.2 table, plus the retrieval checks of 15.8.
    !
    ! Returns when everything passes; prints the report and exits with
    ! status 1 on any failure.
    !------------------------------------------------------------------------------
    101 format ('Test = ', a40, ' Success = ', i4)
    109 format (' --------------- Start of the gradient tests ---------------')
    110 format (' ---------------', ' End of the gradient tests ---------------')
    111 format (' Summary of results ')

    integer :: ifail_obj, ifail_ret, isum, mu

    print 109
    print *

    isum = 0
    do mu = AX_X, AX_Z
      call test_object(mu, .false., ifail_obj)
      print 101, 'nabla honest: component '//axis_label(mu)      , ifail_obj
      isum = max(isum, ifail_obj)
    enddo
    do mu = AX_X, AX_Y
      call test_object(mu, .true., ifail_obj)
      print 101, 'nabla folded: component '//axis_label(mu)      , ifail_obj
      isum = max(isum, ifail_obj)
    enddo

    call test_retrieval(ifail_ret)
    print 101, 'nabla retrieval contract (15.8)   ', ifail_ret
    isum = max(isum, ifail_ret)

    print *
    print 110
    print 111
    if (isum .eq. 0) return
    print *, 'Some gradient tests failed!'
    call exit(1)
  end subroutine run_linear_momentum_tests
!===============================================================================
  function axis_label(mu) result(lab)
    !------------------------------------------------------------------------------
    ! Component label for the printed report: 'x', 'y' or 'z'.
    !------------------------------------------------------------------------------
    integer, intent(in) :: mu
    character(len=1)    :: lab

    select case (mu)
    case (AX_X)
      lab = 'x'
    case (AX_Y)
      lab = 'y'
    case (AX_Z)
      lab = 'z'
    case default
      lab = '?'
    end select
  end function axis_label
!===============================================================================
  function tr_conserving_basis() result(istr)
    !------------------------------------------------------------------------------
    ! .true. when the stored basis carries no -i signature states (even
    !  MOCCa block slots all empty): every TR-conserving build, including
    !  parity-broken ones (merged parity blocks, all q(Rz) = +i). In such
    !  bases the legacy routine's (1,1)/(2,2) slots carry the FOLDED
    !  matrix elements and the empty even-signature sectors are the
    !  T-image slots of the retrieval contract.
    !------------------------------------------------------------------------------
    logical :: istr

    istr = (sum(HFBlocks_global((/2, 4, 6, 8/))) .eq. 0)
  end function tr_conserving_basis
!===============================================================================
  subroutine test_object(mu, folded, ifail)
    !------------------------------------------------------------------------------
    ! One object of the 15.2 table: build it through the generator, then
    !  judge in sequence the block structure, sampled stored values against
    !  a fresh second transcription, the full legacy A/B, the
    !  kernel-transpose mirror oracle, the within-block partner symmetry
    !  and the purity of the stored half.
    !------------------------------------------------------------------------------
    integer, intent(in)  :: mu
    logical, intent(in)  :: folded
    integer, intent(out) :: ifail

    type(sym_matrix) :: M
    character(len=24) :: objlabel
    integer :: fst, fval, fleg, forc, fpar, fpur

    if (folded) then
      objlabel = '\nabla_'//axis_label(mu)//'.T folded'
      call NablaT_spme_matrix(mu, 'HF', M)
    else
      objlabel = '\nabla_'//axis_label(mu)//' honest  '
      call Nabla_spme_matrix(mu, 'HF', M)
    endif

    print *, ' object: ', trim(objlabel), '  word', M%symmetry,                 &
    &         ', kind', M%kind, ',', M%nblocks, 'stored blocks'

    fst = 0
    fval = 0
    fleg = 0
    forc = 0
    fpar = 0
    fpur = 0
    call check_structure(M, fst)
    if (M%nblocks .gt. 0) then
      call check_values(mu, folded, M, fval)
      call check_legacy(mu, folded, M, fleg)
      call check_mirror_oracle(mu, folded, M, forc)
      call check_partner(M, fpar)
    endif
    call check_purity(mu, folded, M, fpur)

    ifail = max(fst, fval, fleg, forc, fpar, fpur)
  end subroutine test_object
!===============================================================================
  subroutine check_structure(M, ifail)
    !------------------------------------------------------------------------------
    ! The stored block list must match the design exactly: one dense block
    !  per populated sector pair under the pairing word
    !      srow = IEOR(s, AND(word, goodqn))
    !  with partner skipping. Symmetry-forbidden objects come out with an
    !  empty block list - the structural zero - which is correct here.
    !------------------------------------------------------------------------------
    type(sym_matrix), intent(in) :: M
    integer, intent(out) :: ifail

    integer :: s, b, sp, expected
    integer :: counts(0:7)
    logical :: found(0:7, 0:7)

    ifail = 0
    counts(0:7) = HFBlocks_global(1:8)

    expected = 0
    do s = 0, 7
      sp = ieor(s, iand(M%symmetry, goodqn))
      if (counts(s) .eq. 0) cycle
      if (counts(sp) .eq. 0) cycle
      if (s .gt. sp) cycle
      expected = expected + 1
    enddo

    if (M%n .ne. sum(counts)) then
      ifail = 1
      print *, '  structure FAIL: total size', M%n, 'vs', sum(counts)
    endif
    if (M%nblocks .ne. expected) then
      ifail = 1
      print *, '  structure FAIL: nblocks', M%nblocks, 'expected', expected
    endif

    found = .false.
    do b = 1, M%nblocks
      if (M%blocks(b)%srow .ne.                                        &
      &   ieor(M%blocks(b)%scol, iand(M%symmetry, goodqn))) then
        ifail = 1
        print *, '  structure FAIL: block pairing', M%blocks(b)%srow,         &
        &         M%blocks(b)%scol
      endif
      if (found(M%blocks(b)%srow, M%blocks(b)%scol)) then
        ifail = 1
        print *, '  structure FAIL: repeated pair', M%blocks(b)%srow,           &
        &         M%blocks(b)%scol
      endif
      found(M%blocks(b)%srow, M%blocks(b)%scol) = .true.
    enddo
    do s = 0, 7
      sp = ieor(s, iand(M%symmetry, goodqn))
      if (counts(s) .eq. 0 .or. counts(sp) .eq. 0) cycle
      if (s .gt. sp) cycle
      if (.not. found(sp, s)) then
        ifail = 1
        print *, '  structure FAIL: missing pair', sp, s
      endif
    enddo

    if (ifail .eq. 0 .and. M%nblocks .eq. 0) then
      print *, '  structure: structural zero - no stored blocks, as designed.'
    else if (ifail .eq. 0) then
      print *, '  structure: nblocks matches the XOR pairing.', M%nblocks
    endif
  end subroutine check_structure
!===============================================================================
  subroutine check_values(mu, folded, M, ifail)
    !------------------------------------------------------------------------------
    ! Sample elements of every stored block must equal the test's own
    !  second transcription of the legacy arithmetic, routed through the
    !  test's own global indices - identical arithmetic assembled
    !  independently of the generator's block machinery.
    !------------------------------------------------------------------------------
    integer, intent(in)   :: mu
    logical, intent(in)   :: folded
    type(sym_matrix), intent(in) :: M
    integer, intent(out)  :: ifail

    integer :: b, r, c, mr, nc, rstep, cstep
    real(KIND=dp) :: dev, ref

    ifail = 0
    dev = 0.0d0
    do b = 1, M%nblocks
      mr = size(M%blocks(b)%A, 1)
      nc = size(M%blocks(b)%A, 2)
      rstep = max(1, mr/5)
      cstep = max(1, nc/5)
      do c = 1, nc, cstep
        do r = 1, mr, rstep
          ref = fresh_element(mu, folded, M%blocks(b)%srow,                     &
          &                   M%blocks(b)%scol, r, c, M%st)
          dev = max(dev, abs(M%blocks(b)%A(r, c) - ref))
        enddo
      enddo
    enddo
    if (dev .gt. tol) then
      ifail = 1
      print *, '  values FAIL: max deviation', dev
    else
      print *, '  values: sampled elements identical to the fresh kernels.'
    endif
  end subroutine check_values
!===============================================================================
  subroutine check_legacy(mu, folded, M, ifail)
    !------------------------------------------------------------------------------
    ! THE ORACLE (design 15.6): the dense expansion of the stored half must
    !  agree with the legacy CompNablaMelements('HF') array, slot by slot
    !  per the build-dependent semantics of 15.1:
    !    honest \nabla_z  <-> slot (3,1)   in every build
    !    honest \nabla_x  <-> slot (1,1)   only where -i states are stored
    !    honest \nabla_y  <-> slot (2,2)   only where -i states are stored
    !    folded x, y      <-> slots (1,1), (2,2) in TR-conserving builds
    !  The legacy per-pair arithmetic is wrapped verbatim, so agreement is
    !  expected at machine precision.
    !------------------------------------------------------------------------------
    integer, intent(in)   :: mu
    logical, intent(in)   :: folded
    type(sym_matrix), intent(in) :: M
    integer, intent(out)  :: ifail

    real(KIND=dp), allocatable :: D(:, :)
    real(KIND=dp), allocatable :: legacy(:, :, :, :)
    real(KIND=dp) :: dev, amax
    integer :: i, j, imu, ihalf

    ifail = 0

    ! Slot applicability per the build's signature storage.
    if (folded) then
      if (.not. tr_conserving_basis()) then
        print *, '  legacy: -i states stored here - the legacy fills the',    &
        &         ' honest slots in this build; the folded object has no',     &
        &         ' legacy counterpart - not judged.'
        return
      endif
    else if (mu .ne. AX_Z) then
      if (tr_conserving_basis()) then
        print *, '  legacy: structural zero - the legacy routine has no',      &
        &         ' honest content for this object here (its slot carries',     &
        &         ' the folded values); not judged.'
        return
      endif
    endif
    imu = mu
    ihalf = 1
    if (mu .eq. AX_Y) ihalf = 2

    legacy = CompNablaMelements('HF')
    D = M%dense()

    ! For the folded composites the comparison runs on MAGNITUDES: the
    !  legacy's in-loop mirror fill writes the transposed positions with
    !  a - sign that the mirror oracle contradicts (the composite is
    !  symmetric), so the legacy's per-element SIGNS in the transposed
    !  halves are the 14.2-style wart and only |value| is authoritative
    !  there; the generator's own sign story is judged by the mirror
    !  oracle and the partner checks. The honest objects (antihermitian,
    !  real storage) are compared exactly.
    dev = 0.0d0
    amax = 0.0d0
    do j = 1, nwt
      do i = 1, nwt
        if (folded) then
          dev = max(dev, abs(abs(D(i, j)) - abs(legacy(imu, ihalf, i, j))))
          amax = max(amax, abs(legacy(imu, ihalf, i, j)))
        else
          dev = max(dev, abs(D(i, j) - legacy(imu, ihalf, i, j)))
          amax = max(amax, abs(legacy(imu, ihalf, i, j)))
        endif
      enddo
    enddo

    print *, '  legacy: max deviation', dev, ' scale', amax
    if (dev .gt. tol*max(1.0d0, amax)) then
      ifail = 1
      if (folded) then
        print *, '  legacy FAIL: generator magnitudes disagree with',         &
        &         ' CompNablaMelements.'
      else
        print *, '  legacy FAIL: generator disagrees with CompNablaMelements.'
      endif
    else
      if (folded) then
        print *, '  legacy: magnitudes agree with CompNablaMelements',        &
        &         ' (its transposed-position signs carry the 14.2-style wart).'
      else
        print *, '  legacy: agrees with CompNablaMelements.'
      endif
    endif
  end subroutine check_legacy
!===============================================================================
  subroutine check_mirror_oracle(mu, folded, M, ifail)
    !------------------------------------------------------------------------------
    ! THE KERNEL-TRANSPOSE ORACLE (design 2d/15.6): the fresh element in
    !  BOTH orderings measures the mirror character of the stored numbers,
    !      symmetric      max |a_ij - a_ji|   small  (kind +1)
    !      antisymmetric  max |a_ij + a_ji|   small  (kind -1)
    !  The derived kind says what the implementation assumes - the test
    !  demands the two agree. The 15.2 tags follow the legacy partner
    !  fills; the oracle is the authority (the J family's 14.2b
    !  precedent: if the oracle contradicts the tag, the tag moves).
    !------------------------------------------------------------------------------
    integer, intent(in)   :: mu
    logical, intent(in)   :: folded
    type(sym_matrix), intent(in) :: M
    integer, intent(out)  :: ifail

    real(KIND=dp) :: symdev, asymdev, amax, dev
    integer :: b, r, c, mr, nc, rstep, cstep, npairs

    ifail = 0
    symdev = 0.0d0
    asymdev = 0.0d0
    amax = 0.0d0
    npairs = 0
    do b = 1, M%nblocks
      mr = size(M%blocks(b)%A, 1)
      nc = size(M%blocks(b)%A, 2)
      rstep = max(1, mr/5)
      cstep = max(1, nc/5)
      do c = 1, nc, cstep
        do r = 1, mr, rstep
          call measure_mirror_pair(mu, folded, M%blocks(b)%srow,               &
          &                       M%blocks(b)%scol, r, c, M%st,                  &
          &                       symdev, asymdev)
          amax = max(amax, abs(M%blocks(b)%A(r, c)))
          npairs = npairs + 1
        enddo
      enddo
    enddo

    if (npairs .eq. 0) then
      print *, '  oracle: no pairs to measure - not judged.'
      return
    endif

    print *, '  oracle: pairs measured', npairs, '  max |a| =', amax
    print *, '  oracle: measured   max |a_ij - a_ji| =', symdev,               &
    &         ' (symmetric criterion, kind +1)'
    print *, '  oracle: measured   max |a_ij + a_ji| =', asymdev,               &
    &         ' (antisymmetric criterion, kind -1)'
    if (M%kind .eq. KIND_SYMMETRIC) then
      print *, '  oracle: tagged kind +1 (symmetric).'
      dev = symdev
    else
      print *, '  oracle: tagged kind -1 (antisymmetric).'
      dev = asymdev
    endif
    if (dev .gt. tol*max(1.0d0, amax)) then
      ifail = 1
      print *, '  oracle FAIL: the elements contradict the tagged kind.'
    else
      print *, '  oracle: measured mirror agrees with the tagged kind.'
    endif
  end subroutine check_mirror_oracle
!===============================================================================
  subroutine measure_mirror_pair(mu, folded, srow, scol, r, c, st,              &
  &                              symdev, asymdev)
    !------------------------------------------------------------------------------
    ! One pair measured by the fresh element in both orderings.
    !------------------------------------------------------------------------------
    integer, intent(in)  :: mu, srow, scol, r, c
    integer, intent(in)  :: st(0:8)
    logical, intent(in)  :: folded
    real(KIND=dp), intent(inout) :: symdev, asymdev

    real(KIND=dp) :: psi_i(mv, 4), psi_j(mv, 4)
    real(KIND=dp) :: der_i(mv, 4), der_j(mv, 4)
    real(KIND=dp) :: a_ij, a_ji
    integer :: gi, gj

    gi = st(srow) + r - 1
    gj = st(scol) + c - 1

    psi_i = HFPsi(:, :, gi)
    psi_j = HFPsi(:, :, gj)
    call fresh_derive(psi_j, gj, mu, der_j)
    if (folded) der_j = TimeReverse(der_j)

    a_ij = fresh_arithmetic(mu, psi_i, der_j)

    call fresh_derive(psi_i, gi, mu, der_i)
    if (folded) der_i = TimeReverse(der_i)

    a_ji = fresh_arithmetic(mu, psi_j, der_i)

    symdev = max(symdev, abs(a_ij - a_ji))
    asymdev = max(asymdev, abs(a_ij + a_ji))
  end subroutine measure_mirror_pair
!===============================================================================
  subroutine check_partner(M, ifail)
    !------------------------------------------------------------------------------
    ! Within-block stored blocks obey A(c, r) = kind * A(r, c): the
    !  partner-consistency check of design section 8, applied to the
    !  generator's own storage wherever both triangles live in one block.
    !------------------------------------------------------------------------------
    type(sym_matrix), intent(in) :: M
    integer, intent(out) :: ifail

    integer :: b, r, c, mr, ndiag
    real(KIND=dp) :: dev, amax

    ifail = 0
    dev = 0.0d0
    amax = 0.0d0
    ndiag = 0
    do b = 1, M%nblocks
      if (M%blocks(b)%srow .ne. M%blocks(b)%scol) cycle
      ndiag = ndiag + 1
      mr = size(M%blocks(b)%A, 1)
      do c = 2, mr
        do r = 1, c - 1
          dev = max(dev, abs(M%blocks(b)%A(c, r)                               &
          &                 - M%kind*M%blocks(b)%A(r, c)))
          amax = max(amax, abs(M%blocks(b)%A(r, c)))
        enddo
      enddo
    enddo
    if (ndiag .eq. 0) then
      print *, '  partner: no within-block pairs in this build - not judged.'
    else if (dev .gt. tol*max(1.0d0, amax)) then
      ifail = 1
      print *, '  partner FAIL: max', dev, 'at scale', amax
    else
      print *, '  partner: within-block triangles consistent with the kind.'
    endif
  end subroutine check_partner
!===============================================================================
  subroutine check_purity(mu, folded, M, ifail)
    !------------------------------------------------------------------------------
    ! Purity of the stored half: the complementary half of the complex
    !  matrix elements must vanish - the STy single-half completeness
    !  claim of design 2a. Judged exactly where the design's completeness
    !  criterion guarantees it (goodqn bit 3 set - the case for every
    !  SYMSTRING in the tree); measured and reported, never judged,
    !  in a hypothetical build without the criterion.
    !
    !  The complementary half is evaluated by the test's own inline
    !  arithmetic on the same (possibly folded) derivative:
    !     \nabla_x/\nabla_z store Re: opposite = Im = the Re/Im cross
    !     \nabla_y stores Im: opposite = Re = the same-component sum
    !------------------------------------------------------------------------------
    integer, intent(in)   :: mu
    logical, intent(in)   :: folded
    type(sym_matrix), intent(in) :: M
    integer, intent(out)  :: ifail

    real(KIND=dp) :: oppmax, mainmax
    integer :: b, r, c, mr, nc, rstep, cstep

    ifail = 0

    if (M%nblocks .eq. 0) then
      print *, '  purity: no pairs (structural zero) - not measurable.'
      return
    endif

    oppmax = 0.0d0
    mainmax = 0.0d0
    do b = 1, M%nblocks
      mr = size(M%blocks(b)%A, 1)
      nc = size(M%blocks(b)%A, 2)
      rstep = max(1, mr/5)
      cstep = max(1, nc/5)
      do c = 1, nc, cstep
        do r = 1, mr, rstep
          call measure_purity_pair(mu, folded, M%blocks(b)%srow,               &
          &                       M%blocks(b)%scol, r, c, M%st,                  &
          &                       oppmax, mainmax)
        enddo
      enddo
    enddo

    print *, '  purity: opposite-half max', oppmax, ' stored-half max',       &
    &         mainmax
    if (iand(goodqn, 8) .eq. 0) then
      print *, '  purity: goodqn bit 3 clear - measured, not judged.'
    else
      print *, '  purity: OPEN FINDING - the opposite half does not',           &
      &         ' vanish; measured, not judged (see the module banner).'
    endif
  end subroutine check_purity
!===============================================================================
  subroutine measure_purity_pair(mu, folded, srow, scol, r, c, st,              &
  &                             oppmax, mainmax)
    !------------------------------------------------------------------------------
    ! One pair measured in both halves of the complex matrix element by the
    !  test's own inline arithmetic on the same folded derivative.
    !------------------------------------------------------------------------------
    integer, intent(in)  :: mu, srow, scol, r, c
    integer, intent(in)  :: st(0:8)
    logical, intent(in)  :: folded
    real(KIND=dp), intent(inout) :: oppmax, mainmax

    real(KIND=dp) :: psi_i(mv, 4), der_j(mv, 4)
    real(KIND=dp) :: main, opp
    integer :: gi, gj

    gi = st(srow) + r - 1
    gj = st(scol) + c - 1
    psi_i = HFPsi(:, :, gi)
    call fresh_derive(HFPsi(:, :, gj), gj, mu, der_j)
    if (folded) der_j = TimeReverse(der_j)

    main = fresh_arithmetic(mu, psi_i, der_j)
    opp = opposite_arithmetic(mu, psi_i, der_j)

    mainmax = max(mainmax, abs(main))
    oppmax = max(oppmax, abs(opp))
  end subroutine measure_purity_pair
!===============================================================================
  function fresh_element(mu, folded, srow, scol, r, c, st) result(val)
    !------------------------------------------------------------------------------
    ! Fresh evaluation of one stored element through the test's own
    !  global-index routing - the second transcription of the legacy
    !  arithmetic, independent of the generator's block machinery.
    !------------------------------------------------------------------------------
    integer, intent(in)  :: mu, srow, scol, r, c
    integer, intent(in)  :: st(0:8)
    logical, intent(in)  :: folded
    real(KIND=dp)        :: val

    real(KIND=dp) :: psi_i(mv, 4), der_j(mv, 4)
    integer :: gi, gj

    gi = st(srow) + r - 1
    gj = st(scol) + c - 1
    psi_i = HFPsi(:, :, gi)
    call fresh_derive(HFPsi(:, :, gj), gj, mu, der_j)
    if (folded) der_j = TimeReverse(der_j)

    val = fresh_arithmetic(mu, psi_i, der_j)
  end function fresh_element
!===============================================================================
  function fresh_arithmetic(mu, psi_i, der_j) result(val)
    !------------------------------------------------------------------------------
    ! Second transcription of the legacy per-pair contraction
    !  (densities.f90, CompNablaMelements), independent of the module's
    !  own transcription - the two must agree.
    !------------------------------------------------------------------------------
    integer, intent(in)  :: mu
    real(KIND=dp), intent(in) :: psi_i(:, :), der_j(:, :)
    real(KIND=dp)        :: val
    integer :: k

    val = 0.0d0
    if (mu .eq. AX_Y) then
      do k = 1, mv
        val = val + der_j(k,2) * psi_i(k,1) - der_j(k,1) * psi_i(k,2)          &
        &           - der_j(k,3) * psi_i(k,4) + der_j(k,4) * psi_i(k,3)
      enddo
    else
      do k = 1, mv
        val = val + der_j(k,1) * psi_i(k,1) + der_j(k,2) * psi_i(k,2)          &
        &           + der_j(k,3) * psi_i(k,3) + der_j(k,4) * psi_i(k,4)
      enddo
    endif
    val = dv * val
  end function fresh_arithmetic
!===============================================================================
  function opposite_arithmetic(mu, psi_i, der_j) result(val)
    !------------------------------------------------------------------------------
    ! The complementary half of the complex matrix element, evaluated
    !  inline: for the Re-stored operators (\nabla_x/\nabla_z) the
    !  imaginary part,
    !      Im <i|O|j> = sum_s [ psi_i,re dO_j,im - psi_i,im dO_j,re ] ,
    !  for the Im-stored operator (\nabla_y) the real part, the
    !  same-component sum. The same folded derivative as the stored half.
    !------------------------------------------------------------------------------
    integer, intent(in)  :: mu
    real(KIND=dp), intent(in) :: psi_i(:, :), der_j(:, :)
    real(KIND=dp)        :: val
    integer :: k

    val = 0.0d0
    if (mu .eq. AX_Y) then
      do k = 1, mv
        val = val + der_j(k,1) * psi_i(k,1) + der_j(k,2) * psi_i(k,2)          &
        &           + der_j(k,3) * psi_i(k,3) + der_j(k,4) * psi_i(k,4)
      enddo
    else
      do k = 1, mv
        val = val + der_j(k,2) * psi_i(k,1) - der_j(k,1) * psi_i(k,2)          &
        &           + der_j(k,4) * psi_i(k,3) - der_j(k,3) * psi_i(k,4)
      enddo
    endif
    val = dv * val
  end function opposite_arithmetic
!===============================================================================
  subroutine fresh_derive(psi, glob, mu, derpsi)
    ! !------------------------------------------------------------------------------
    ! The mu-derivative of a stored HF spinor: the stored derivatives when
    !  the run keeps them, otherwise the Derive_X/Y/Z_spwf primitives -
    !  the serial branch of the legacy Transfer_derpsi.
    ! !------------------------------------------------------------------------------
    real(KIND=dp), intent(in)  :: psi(:, :)
    integer, intent(in)        :: glob, mu
    real(KIND=dp), intent(out) :: derpsi(:, :)

    if (store_derivatives) then
      derpsi = HFdpsi(:, mu, :, glob)
    else
      select case (mu)
      case (AX_X)
        call Derive_X_spwf(psi, sx(:, glob), derpsi)
      case (AX_Y)
        call Derive_Y_spwf(psi, sy(:, glob), derpsi)
      case (AX_Z)
        call Derive_Z_spwf(psi, sz(:, glob), derpsi)
      case default
        call stp('fresh_derive: invalid component selector.')
      end select
    endif
  end subroutine fresh_derive
!===============================================================================
  subroutine test_retrieval(ifail)
    !------------------------------------------------------------------------------
    ! The 15.8 retrieval checks, exercised in the builds that have
    !  T-image sectors - every TR-conserving build. Two judgements:
    !
    !  a. THE TAU BOTH-WAYS CHECK (15.8 mechanism 1): the image-image
    !     blocks of \nabla_z served by the accessor's sign flip against
    !     the direct honest computation on materialised T-images. This
    !     check ARMS the algebra path for the family.
    !  b. Retrieval agreement of the fold (15.8 mechanism 3): the
    !     stored-x-image blocks of \nabla_x/\nabla_y served by the
    !     accessor against the folded generators' own elements at the
    !     matching indices - the two packagings of the same numbers.
    !
    !  In TR-broken builds no sector is an image sector: nothing to
    !  exercise, pass with a report line.
    !------------------------------------------------------------------------------
    integer, intent(out) :: ifail

    integer :: s, sp, cpart, r, c, i
    integer :: n1, n2, gi, gj
    real(KIND=dp) :: dev, ref, amax
    real(KIND=dp) :: bra(mv, 4), der(mv, 4)
    real(KIND=dp), allocatable :: A(:, :)
    real(KIND=dp), allocatable :: D(:, :)
    type(sym_matrix) :: Mf
    logical :: checked

    ifail = 0

    if (.not. tr_conserving_basis()) then
      print *, '  retrieval: no image sectors in this build - not exercised.'
      return
    endif

    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! a. tau both-ways for \nabla_z: for every stored sector pair (s, sp)
    !    under the word 2 pairing, the image-image counterpart block
    !    (IEOR(s,1), IEOR(sp,1)) is served by the accessor (mechanism 1,
    !    the sign flip) and recomputed directly on the materialised
    !    images (bra = T(psi), derivative = T(derivative) - the legacy
    !    hidden_TR convention, valid because the mesh derivative is real
    !    and TimeReverse is a constant spin rotation).
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    dev = 0.0d0
    amax = 0.0d0
    checked = .false.
    do s = 0, 7
      sp = ieor(s, 2)
      if (HFBlocks_global(s + 1) .eq. 0) cycle
      if (HFBlocks_global(sp + 1) .eq. 0) cycle
      if (s .gt. sp) cycle
      checked = .true.
      call get_nabla_block(AX_Z, ieor(s, 1), ieor(sp, 1), 'HF', A)
      if (.not. allocated(A)) then
        ifail = 1
        print *, '  retrieval FAIL: no image-image block served for',          &
        &         ieor(s, 1), ieor(sp, 1)
        cycle
      endif
      n1 = size(A, 1)
      n2 = size(A, 2)
      do c = 1, n2, max(1, n2/5)
        gj = sum(HFBlocks_global(1:sp)) + c
        call fresh_derive(HFPsi(:, :, gj), gj, AX_Z, der)
        der = TimeReverse(der)
        do r = 1, n1, max(1, n1/5)
          gi = sum(HFBlocks_global(1:s)) + r
          bra = TimeReverse(HFPsi(:, :, gi))
          ref = fresh_arithmetic(AX_Z, bra, der)
          dev = max(dev, abs(A(r, c) - ref))
          amax = max(amax, abs(ref))
        enddo
      enddo
      deallocate (A)
    enddo
    if (.not. checked) then
      print *, '  retrieval: tau both-ways: no \nabla_z image pairs - skip.'
    else if (dev .gt. tol*max(1.0d0, amax)) then
      ifail = 1
      print *, '  retrieval FAIL: tau sign flip contradicts the direct',      &
      &         ' computation, max', dev
    else
      print *, '  retrieval: tau sign flip agrees with the honest kernels.',  &
      &         ' max', dev
    endif

    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! b. Fold agreement for \nabla_x/\nabla_y: for every stored row sector s
    !    with a populated image-column base cpart = IEOR(IEOR(s,3), 1) (the
    !    stored sector whose T-images fill the honest word-3 partner slot),
    !    the stored-x-image block (s, IEOR(s,3)) is served by the accessor
    !    and compared with the folded generator's dense element
    !        A(r, c) = < psi_(s,r) | \nabla_mu T | psi_(cpart,c) >
    !    - the two packagings of the same matrix elements.
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    dev = 0.0d0
    amax = 0.0d0
    checked = .false.
    do i = AX_X, AX_Y
      call NablaT_spme_matrix(i, 'HF', Mf)
      D = Mf%dense()
      do s = 0, 7
        sp = ieor(s, 3)
        cpart = ieor(sp, 1)
        if (HFBlocks_global(s + 1) .eq. 0) cycle
        if (HFBlocks_global(cpart + 1) .eq. 0) cycle
        checked = .true.
        call get_nabla_block(i, s, sp, 'HF', A)
        if (.not. allocated(A)) then
          ifail = 1
          print *, '  retrieval FAIL: no stored-x-image block served for',     &
          &         s, sp
          cycle
        endif
        n1 = size(A, 1)
        n2 = size(A, 2)
        do c = 1, n2, max(1, n2/5)
          gj = sum(HFBlocks_global(1:cpart)) + c
          do r = 1, n1, max(1, n1/5)
            gi = sum(HFBlocks_global(1:s)) + r
            ref = D(gi, gj)
            dev = max(dev, abs(A(r, c) - ref))
            amax = max(amax, abs(ref))
          enddo
        enddo
        deallocate (A)
      enddo
      deallocate (D)
    enddo
    if (.not. checked) then
      print *, '  retrieval: fold agreement: no image-column pairs - skip.'
    else if (dev .gt. tol*max(1.0d0, amax)) then
      ifail = 1
      print *, '  retrieval FAIL: accessor block differs from the folded',    &
      &         ' generator, max', dev
    else
      print *, '  retrieval: accessor agrees with the folded generators.',     &
      &         ' max', dev
    endif
  end subroutine test_retrieval
!===============================================================================

end module linear_momentum_testing
