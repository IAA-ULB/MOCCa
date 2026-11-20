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
module angular_momentum_testing
!==============================================================================
  !
  ! Unit tests of the angular-momentum generators of the angular_momentum
  !  module (design summary section 14), armed with unit_test = .true. in
  !  the &mesh/ namelist exactly like the moments suite, and acting on the
  !  final mean-field state.
  !
  ! The centerpiece is the KERNEL-TRANSPOSE ORACLE of design section 14.3:
  !  the angmom kernels are pure functions of two spinors and one gradient,
  !  so evaluating K(wf_i, wf_j) and K(wf_j, wf_i) on real states measures
  !  the mirror character of each object's stored numbers - symmetric or
  !  antisymmetric - independent of every convention and consumer. The
  !  generators' class tags (and the storage kinds init derives from them)
  !  are pinned to what this oracle measures: the unit run PRINTS measured
  !  facts against the tags, so a class-table disagreement cannot pass
  !  silently. This is the verification path that settles the folded
  !  pair's kinds, whose class arithmetic does not transfer from the
  !  linear-operator derivations of design section 2a.
  !
  ! Checks per generated object (the honest trio J_x, J_y, J_z and the
  !  folded pair J_x.T, J_y.T), on every build:
  !   (1) structure: the stored block list matches the XOR pairing under
  !       the object's word and the goodqn mask; symmetry-forbidden
  !       objects come out as structural zeros and are judged here.
  !   (2) values: sample elements of every stored block equal a fresh
  !       kernel evaluation routed through the test's own global indices.
  !   (3) kernel-transpose oracle: sampled pairs per stored block are
  !       evaluated by the kernels in both orderings; the measured mirror
  !       character is printed against the derived kind - the oracle
  !       versus the sym_matrix implementation.
  !   (4) partner symmetry: within-block stored blocks obey
  !       A(c, r) = kind * A(r, c).
  !   (5) purity report of the stored half (honest trio only): the
  !       opposite-half kernel magnitudes on the same pairs are printed -
  !       single-chamber residues, not matrix elements, judged by (6).
  !       No opposite-half t-kernels exist, so the folded pair is not
  !       measurable.
  !   (6) full-box reconstruction (honest trio only): the deciding test
  !       of design section 14.2b. The complex matrix elements are
  !       recomputed over the chamber-extended spinors - the design 13.1
  !       step-4 pattern: the main half must equal the kernel evaluation
  !       on the stored chamber (chamber uniformity), the opposite half
  !       must vanish (the design 2a STy constraint, a statement about
  !       the chamber sum that no reduced-mesh evaluation can express).
  !       Green here adjudicates the O(1) opposite-half purity prints as
  !       the reduced-mesh residues of 14.2b, never the reverse.
  !
  ! Unlike run_moments_tests, this suite RETURNS on success (it runs
  !  first, so the moments tests still terminate the program afterwards);
  !  on failure it prints the report and exits with status 1.
  !
  ! Serial semantics, like the moments net (design sections 14.3 and 13.4).
  !==============================================================================
  use compilation,       only: dp
  use geninfo
  use wavefunctions,     only: HFPsi, HFdpsi, sx, sy, sz
  use wavefunctions,     only: Derive_X_spwf, Derive_Y_spwf, Derive_Z_spwf
  use wavefunctions,     only: angmom_x_real, angmom_x_imag
  use wavefunctions,     only: angmom_y_real, angmom_y_imag
  use wavefunctions,     only: angmom_z_real, angmom_z_imag
  use wavefunctions,     only: angmom_xt_real
  use wavefunctions,     only: angmom_yt_imag
  use wavefunctions,     only: HFBlocks_global
  use wavefunctions,     only: TimeReverse
  use block_matrix_module, only: sym_matrix, KIND_SYMMETRIC,                  &
  &                              KIND_ANTISYMMETRIC
  use block_matrix_module, only: HERMITIAN, ANTIHERMITIAN
  use block_matrix_module, only: STY_EVEN, STY_ODD
  use angular_momentum,  only: AX_X, AX_Y, AX_Z
  use angular_momentum,  only: Jmu_spme_matrix, get_J_block
  use angular_momentum,  only: TAU_ANGFAM

  implicit none

  private
  public :: run_angular_momentum_tests

  ! Tolerance of the oracle and of the consistency checks: the kernels
  !  accumulate over the full mesh, but identical inputs evaluated through
  !  the same arithmetic order agree bit for bit - tight bar by design.
  real(KIND=dp), parameter :: tol = 1.0d-10

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Arm-owned twin machinery (D-7, design 15.9 row 2): since the eager
  !  twin entry point retired behind the accessor, this test owns the
  !  folded reference arithmetic - the t-kernels through Jm_fresh_element
  !  - and builds the folded objects itself. The class labels below are
  !  the folded table of design 14.1: what the kernel-transpose oracle
  !  measured, never re-derived. Structure, value, mirror, partner and
  !  purity checks on these objects are unchanged; the accessor
  !  packaging is asserted against them in test_retrieval.
  type :: ang_twin_params
    integer :: st(0:8)     ! sector bounds of the twin under construction
    integer :: mu          ! component: AX_X / AX_Y
  end type ang_twin_params

contains

  subroutine run_angular_momentum_tests()
    !------------------------------------------------------------------------------
    ! Run the unit tests of the angular-momentum family: the honest trio
    !  through the module generator, the folded twins through the arm's
    !  own construction (D-7, 15.9 row 2 - the retired eager entry's
    !  arithmetic, the class labels the oracle asserted), each judged by
    !  the structure, value, mirror-oracle, partner-symmetry and purity
    !  checks, plus the retrieval checks of design 15.8.
    !
    ! Returns when everything passes; prints the report and exits with
    !  status 1 on any failure.
    !------------------------------------------------------------------------------
    101 format ('Test = ', a40, ' Success = ', i4)
    109 format (' -------------- Start of the angular-momentum tests --------------')
    110 format (' -------------- End of the angular-momentum tests ----------------')
    111 format (' Summary of results ')

    integer :: ifail_obj, isum, mu

    print 109
    print *

    isum = 0
    do mu = AX_X, AX_Z
      call test_object(mu, .false., ifail_obj)
      print 101, 'angmom honest: component '//axis_label(mu)      , ifail_obj
      isum = max(isum, ifail_obj)
    enddo
    do mu = AX_X, AX_Y
      call test_object(mu, .true., ifail_obj)
      print 101, 'angmom folded: component '//axis_label(mu)      , ifail_obj
      isum = max(isum, ifail_obj)
    enddo

    call test_retrieval(ifail_obj)
    print 101, 'angmom retrieval contract (15.8)   ', ifail_obj
    isum = max(isum, ifail_obj)

    print *
    print 110
    print 111
    if (isum .eq. 0) return
    print *, 'Some angular-momentum tests failed!'
    call exit(1)
  end subroutine run_angular_momentum_tests
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
  subroutine test_object(mu, folded, ifail)
    !------------------------------------------------------------------------------
    ! One object of the 14.1 table: build it through the generator, then
    !  judge in sequence the block structure, sampled stored values
    !  against fresh kernel evaluations, the kernel-transpose oracle
    !  against the derived kind, the within-block partner symmetry and -
    !  honest trio only - the single-chamber purity report plus the
    !  full-box reconstruction judgment. The printed
    !  detail lines identify the object under test.
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !   mu     : component selector, AX_X / AX_Y / AX_Z
    !   folded : .true. for the J_mu.T composite
    ! Output:
    !   ifail  : 0 all checks passed, 1 at least one failed
    !------------------------------------------------------------------------------
    integer, intent(in)  :: mu
    logical, intent(in)  :: folded
    integer, intent(out) :: ifail

    type(sym_matrix) :: M
    character(len=24) :: objlabel
    integer :: fst, fval, forc, fpar, fpur, fful

    if (folded) then
      objlabel = 'J_'//axis_label(mu)//'.T folded'
      call arm_build_twin(mu, M)
    else
      objlabel = 'J_'//axis_label(mu)//' honest  '
      call Jmu_spme_matrix(mu, M)
    endif

    print *, ' object: ', trim(objlabel), '  word', M%symmetry,                 &
    &         ', kind', M%kind, ',', M%nblocks, 'stored blocks'

    fst = 0
    fval = 0
    forc = 0
    fpar = 0
    fpur = 0
    fful = 0
    call check_structure(M, fst)
    if (M%nblocks .gt. 0) then
      call check_values(mu, folded, M, fval)
      call check_mirror_oracle(mu, folded, M, forc)
      call check_partner(M, fpar)
    endif
    call check_purity(mu, folded, M, fpur)
    if (.not. folded) call check_fullbox(mu, M, fful)

    ifail = max(fst, fval, forc, fpar, fpur, fful)
  end subroutine test_object
!===============================================================================
  subroutine check_structure(M, ifail)
    !------------------------------------------------------------------------------
    ! The stored block list must match the design exactly: one dense block
    !  per populated sector pair under the pairing word
    !      srow = IEOR(s, AND(word, goodqn))
    !  with partner skipping (one block per pair, scol <= srow).
    !  Symmetry-forbidden objects come out with an empty block list - the
    !  structural zero - which is the correct result here.
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

    ! Every stored block on an allowed pair, covering them all.
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
    ! Sample elements of every stored block must equal a fresh kernel
    !  evaluation routed through the test's own global indices - identical
    !  arithmetic assembled independently of the generator's block
    !  machinery.
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
          ref = Jm_fresh_element(mu, folded, M%blocks(b)%srow,                 &
          &                     M%blocks(b)%scol, r, c, M%st)
          dev = max(dev, abs(M%blocks(b)%A(r, c) - ref))
        enddo
      enddo
    enddo
    if (dev .gt. 0.0d0) then
      ifail = 1
      print *, '  values FAIL: max deviation', dev
    else
      print *, '  values: sampled elements identical to fresh kernels.'
    endif
  end subroutine check_values
!===============================================================================
  subroutine check_mirror_oracle(mu, folded, M, ifail)
    !------------------------------------------------------------------------------
    ! THE ORACLE: for sampled pairs of every stored block, the kernels are
    !  evaluated in BOTH orderings,
    !      a_ij = K(psi_i, psi_j, derpsi_j),   a_ji = K(psi_j, psi_i, derpsi_i)
    !  and the mirror character of the stored numbers is measured:
    !      symmetric      max |a_ij - a_ji|   small  (kind +1)
    !      antisymmetric  max |a_ij + a_ji|   small  (kind -1)
    !  The kernels are the authority (design 14.3): the printout says what
    !  they store, the derived kind says what the sym_matrix implementation
    !  assumes - and the test demands the two agree.
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
      print *, '  oracle FAIL: the kernels contradict the tagged kind.'
    else
      print *, '  oracle: measured mirror agrees with the tagged kind.'
    endif
  end subroutine check_mirror_oracle
!===============================================================================
  subroutine measure_mirror_pair(mu, folded, srow, scol, r, c, st,              &
  &                             symdev, asymdev)
    !------------------------------------------------------------------------------
    ! One pair measured by the kernels in both orderings.
    !------------------------------------------------------------------------------
    integer, intent(in)  :: mu, srow, scol, r, c
    integer, intent(in)  :: st(0:8)
    logical, intent(in)  :: folded
    real(KIND=dp), intent(inout) :: symdev, asymdev

    real(KIND=dp) :: psi_i(mv, 4), psi_j(mv, 4)
    real(KIND=dp) :: derpsi_i(mv, 3, 4), derpsi_j(mv, 3, 4)
    real(KIND=dp) :: a_ij, a_ji
    integer :: gi, gj

    gi = st(srow) + r - 1
    gj = st(scol) + c - 1

    psi_i = HFPsi(:, :, gi)
    psi_j = HFPsi(:, :, gj)
    call Jm_fresh_derive(psi_i, gi, derpsi_i)
    call Jm_fresh_derive(psi_j, gj, derpsi_j)

    a_ij = Jm_oracle_kernel(mu, folded, psi_i, psi_j, derpsi_j)
    a_ji = Jm_oracle_kernel(mu, folded, psi_j, psi_i, derpsi_i)

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
    ! Purity report of the stored half (honest trio): the opposite-half
    !  kernel measures the complementary part of the complex matrix
    !  elements - on the single stored chamber, with the dv multiplicity,
    !  which is valid only for chamber-uniform products. CAVEAT: the
    !  opposite-half kernels angmom_{x,z}_imag / angmom_y_real are
    !  production-unreferenced (this arm is their first caller in the
    !  tree), and their cross-character products are chamber-CANCELLING:
    !  the magnitudes printed here are O(1) chamber residues, NOT
    !  full-box matrix elements. They are measured, never judged by this
    !  report - the full-box reconstruction of check_fullbox below is
    !  the arbiter. No opposite-half t-kernels exist, so the folded pair
    !  is not measurable at all.
    !------------------------------------------------------------------------------
    integer, intent(in)   :: mu
    logical, intent(in)   :: folded
    type(sym_matrix), intent(in) :: M
    integer, intent(out)  :: ifail

    real(KIND=dp) :: oppmax, mainmax
    integer :: b, r, c, mr, nc, rstep, cstep

    ifail = 0

    if (folded) then
      print *, '  purity: no opposite-half t-kernel exists - not measurable.'
      return
    endif
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
          call measure_purity_pair(mu, M%blocks(b)%srow,                        &
          & M%blocks(b)%scol, r, c, M%st, oppmax, mainmax)
        enddo
      enddo
    enddo

    print *, '  purity: opposite-half max', oppmax, ' stored-half max',       &
    &         mainmax, ' (single-chamber residues: judged by check_fullbox)'
  end subroutine check_purity
!===============================================================================
  subroutine measure_purity_pair(mu, srow, scol, r, c, st, oppmax, mainmax)
    !------------------------------------------------------------------------------
    ! One pair measured by BOTH half kernels of the honest component mu:
    !  the stored half and the opposite half.
    !------------------------------------------------------------------------------
    integer, intent(in)  :: mu, srow, scol, r, c
    integer, intent(in)  :: st(0:8)
    real(KIND=dp), intent(inout) :: oppmax, mainmax

    real(KIND=dp) :: psi_i(mv, 4), psi_j(mv, 4)
    real(KIND=dp) :: derpsi_j(mv, 3, 4)
    integer :: gi, gj

    gi = st(srow) + r - 1
    gj = st(scol) + c - 1
    psi_i = HFPsi(:, :, gi)
    psi_j = HFPsi(:, :, gj)
    call Jm_fresh_derive(psi_j, gj, derpsi_j)

    select case (mu)
    case (AX_X)
      mainmax = max(mainmax, abs(angmom_x_real(psi_i, psi_j, derpsi_j)))
      oppmax = max(oppmax, abs(angmom_x_imag(psi_i, psi_j, derpsi_j)))
    case (AX_Y)
      mainmax = max(mainmax, abs(angmom_y_imag(psi_i, psi_j, derpsi_j)))
      oppmax = max(oppmax, abs(angmom_y_real(psi_i, psi_j, derpsi_j)))
    case (AX_Z)
      mainmax = max(mainmax, abs(angmom_z_real(psi_i, psi_j, derpsi_j)))
      oppmax = max(oppmax, abs(angmom_z_imag(psi_i, psi_j, derpsi_j)))
    end select
  end subroutine measure_purity_pair
!===============================================================================
  subroutine check_fullbox(mu, M, ifail)
    !------------------------------------------------------------------------------
    ! DECIDING TEST of design section 14.2b: the full-box reconstruction of
    !  the complex matrix elements
    !
    !      < i | J_mu | j >
    !
    !  over the chamber-extended spinors, for sampled pairs of every stored
    !  block - the design 13.1 step-4 pattern applied to the angmom family.
    !  Two facts are JUDGED per pair, closing what the single-chamber
    !  purity print of check_purity can only measure:
    !
    !   (a) identity: the main half is chamber-uniform - the full-box
    !       quadrature of fullbox_angmom_me (spinors extended across the
    !       fold planes with the per-component signs sx/sy/sz, gradients
    !       with the additional chain-rule sign, operator values fresh at
    !       every chamber's coordinates, plain cell volume dv/2^nch) must
    !       equal the plain kernel evaluation on the single stored chamber.
    !   (b) STy zero: the opposite half is chamber-cancelling - its
    !       full-box value must vanish: the design 2a characteristic
    !       constraint, a statement about the chamber sum that no
    !       reduced-mesh evaluation can express.
    !
    ! Green here adjudicates the O(1) opposite-half purity prints as the
    !  reduced-mesh residues predicted by 14.2b - the premise
    !  "jz_ij is always a real number" of momentsofinertia.f90 is judged,
    !  never assumed.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !   mu : honest component selector, AX_X / AX_Y / AX_Z
    !   M  : the honest object's sym_matrix; its stored blocks supply the
    !        sampled pairs, exactly as in check_purity
    ! Output:
    !   ifail : 0 every sampled pair passes both judgments,
    !           1 at least one failure occurred
    !------------------------------------------------------------------------------
    integer, intent(in)          :: mu
    type(sym_matrix), intent(in) :: M
    integer, intent(out)         :: ifail

    integer :: b, r, c, mr, nc, rstep, cstep, gi, gj, npair
    real(KIND=dp) :: fb_re, fb_im, fb_main, fb_opp, kern_main, dev, bar
    real(KIND=dp) :: maxid, maxop

    ifail = 0
    if (M%nblocks .eq. 0) then
      print *, '  fullbox: no pairs (structural zero) - not measurable.'
      return
    endif

    maxid = 0.0d0
    maxop = 0.0d0
    npair = 0
    do b = 1, M%nblocks
      mr = size(M%blocks(b)%A, 1)
      nc = size(M%blocks(b)%A, 2)
      rstep = max(1, mr/5)
      cstep = max(1, nc/5)
      do c = 1, nc, cstep
        do r = 1, mr, rstep
          gi = M%st(M%blocks(b)%srow) + r - 1
          gj = M%st(M%blocks(b)%scol) + c - 1
          kern_main = Jm_fresh_element(mu, .false., M%blocks(b)%srow,         &
          &                          M%blocks(b)%scol, r, c, M%st)
          call fullbox_angmom_me(mu, gi, gj, fb_re, fb_im)
          if (mu .eq. AX_Y) then
            fb_main = fb_im
            fb_opp = fb_re
          else
            fb_main = fb_re
            fb_opp = fb_im
          endif
          npair = npair + 1
          bar = tol*max(1.0d0, abs(kern_main))
          dev = abs(fb_main - kern_main)
          if (dev .gt. bar) ifail = 1
          maxid = max(maxid, dev)
          dev = abs(fb_opp)
          if (dev .gt. bar) ifail = 1
          maxop = max(maxop, dev)
        enddo
      enddo
    enddo

    print *, '  fullbox: identity max', maxid, ' opposite-half max', maxop,    &
    &        ' over', npair, 'pairs (judged)'
    if (ifail .eq. 0) then
      print *, '  fullbox: PASS - main half chamber-uniform, opposite half'   &
      &        //' vanishes over the full box.'
    else
      print *, '  fullbox: FAIL - a deviation above tolerance occurred.'
    endif
  end subroutine check_fullbox
!===============================================================================
  subroutine fullbox_angmom_me(mu, gi, gj, me_re, me_im)
    !------------------------------------------------------------------------------
    ! One full-box matrix element < psi_i | J_mu | psi_j > by plain mesh
    !  quadrature over ALL chambers of the folded mesh - the angmom
    !  counterpart of the full-box oracle of moments_testing (design
    !  13.1 step 4). It shares only trusted raw inputs with production:
    !  the stored spinors, their fresh gradients and the sx/sy/sz
    !  per-component reflection signs.
    !
    ! The stored mesh is the fundamental chamber; the full box is 2^nch
    !  copies with flipped axes (nch = reduX + reduY + reduZ, the dv
    !  multiplicity). At the full-box points of chamber c every spinor
    !  component k of a state takes its stored mirror value multiplied by
    !  the product of sx/sy/sz(k) over the chamber's flipped axes; every
    !  gradient field of component k acquires the additional chain-rule
    !  sign (-1)**(flip of its derivative axis); the operator coordinates
    !  are the mirrored COM-shifted mesh values. The full-box cell volume
    !  is dv/2^nch, so the chamber sum below is the entire integral.
    !
    ! With the real component storage - (2s-1)/(2s) = Re/Im of spin
    !  projection s - and a = Rei + i Imi, b = Rej + i Imj per spin, the
    !  per-spin summand of < a | (i V_mu + sigma_mu/2) | b > is
    !
    !     V_z = i(y d_x - x d_y)  V_x = i(z d_y - y d_z)
    !     V_y = i(x d_z - z d_x)
    !
    !     orbital Re: -Rei Vim + Imi Ver     orbital Im: Rei Ver + Imi Vim
    !
    !  plus the spin terms with the pairings sigma_mu dictates (z: same
    !  spin, x/y: crossed spins, y: with the (0,-i;i,0) i-factors).
    !------------------------------------------------------------------------------
    ! Input:
    !   mu : component selector, AX_X / AX_Y / AX_Z (honest components only)
    !   gi, gj : global HF-basis indices of the row/column state
    ! Output:
    !   me_re, me_im : Re / Im of < psi_i | J_mu | psi_j > over the full box
    !------------------------------------------------------------------------------
    integer, intent(in)  :: mu, gi, gj
    real(KIND=dp), intent(out) :: me_re, me_im

    real(KIND=dp) :: psi_i(mv, 4), psi_j(mv, 4), derpsi_j(mv, 3, 4)
    real(KIND=dp) :: Rei(mv), Imi(mv), Rej(mv), Imj(mv)
    real(KIND=dp) :: vxr(mv), vxi(mv), vyr(mv), vyi(mv), vzr(mv), vzi(mv)
    real(KIND=dp) :: vor(mv), voi(mv), xc(mv), yc(mv), zc(mv)
    real(KIND=dp) :: sri(4), srj(4), sgx, sgy, sgz, sgn, dvcell
    integer :: c, ic, fx, fy, fz, k, s, t, nch

    psi_i = HFPsi(:, :, gi)
    psi_j = HFPsi(:, :, gj)
    call Jm_fresh_derive(psi_j, gj, derpsi_j)

    nch = 2**(reduX + reduY + reduZ)
    dvcell = dv/nch
    me_re = 0.0d0
    me_im = 0.0d0

    do c = 1, nch
      ! Chamber decomposition over the folded axes only: axis bits of
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

      ! Per-component extension signs of both states and the chamber
      !  coordinates (mirrored COM-shifted mesh values).
      sri = 1.0d0
      srj = 1.0d0
      do k = 1, 4
        if (fx .eq. 1) then
          sri(k) = sri(k)*sx(k, gi)
          srj(k) = srj(k)*sx(k, gj)
        endif
        if (fy .eq. 1) then
          sri(k) = sri(k)*sy(k, gi)
          srj(k) = srj(k)*sy(k, gj)
        endif
        if (fz .eq. 1) then
          sri(k) = sri(k)*sz(k, gi)
          srj(k) = srj(k)*sz(k, gj)
        endif
      enddo
      sgx = merge(-1.0d0, 1.0d0, fx .eq. 1)
      sgy = merge(-1.0d0, 1.0d0, fy .eq. 1)
      sgz = merge(-1.0d0, 1.0d0, fz .eq. 1)
      xc = merge(-meshgrid_shifted(:, 1), meshgrid_shifted(:, 1), fx .eq. 1)
      yc = merge(-meshgrid_shifted(:, 2), meshgrid_shifted(:, 2), fy .eq. 1)
      zc = merge(-meshgrid_shifted(:, 3), meshgrid_shifted(:, 3), fz .eq. 1)

      do s = 1, 2
        Rei = sri(2*s-1)*psi_i(:, 2*s-1)
        Imi = sri(2*s  )*psi_i(:, 2*s  )
        ! Extended gradient fields of the ket spin s: the stored
        !  derivatives with the extension and chain-rule signs.
        vxr = srj(2*s-1)*sgx*derpsi_j(:, 1, 2*s-1)
        vxi = srj(2*s  )*sgx*derpsi_j(:, 1, 2*s  )
        vyr = srj(2*s-1)*sgy*derpsi_j(:, 2, 2*s-1)
        vyi = srj(2*s  )*sgy*derpsi_j(:, 2, 2*s  )
        vzr = srj(2*s-1)*sgz*derpsi_j(:, 3, 2*s-1)
        vzi = srj(2*s  )*sgz*derpsi_j(:, 3, 2*s  )

        ! Orbital part: i V_mu on the ket spin s, sandwiched by the bra.
        select case (mu)
        case (AX_Z)
          vor = yc*vxr - xc*vyr
          voi = yc*vxi - xc*vyi
        case (AX_X)
          vor = zc*vyr - yc*vzr
          voi = zc*vyi - yc*vzi
        case (AX_Y)
          vor = xc*vzr - zc*vxr
          voi = xc*vzi - zc*vxi
        case default
          print *, 'fullbox_angmom_me: components are AX_X, AX_Y and AX_Z.'
          call stp('Invalid component in fullbox_angmom_me.')
        end select
        me_re = me_re + dvcell*(dot_product(Imi, vor)                        &
        &                               - dot_product(Rei, voi))
        me_im = me_im + dvcell*(dot_product(Rei, vor)                        &
        &                               + dot_product(Imi, voi))

        ! Spin part: sigma_mu couples the bra spin s with ket spin t -
        !  the same spin for sigma_z, crossed for sigma_x / sigma_y.
        if (mu .eq. AX_Z) then
          t = s
        else
          t = 3 - s
        endif
        Rej = srj(2*t-1)*psi_j(:, 2*t-1)
        Imj = srj(2*t  )*psi_j(:, 2*t  )
        sgn = merge(1.0d0, -1.0d0, s .eq. 1)
        select case (mu)
        case (AX_Z)
          me_re = me_re + dvcell*0.5d0*sgn                                   &
          &                    *(dot_product(Rei, Rej) + dot_product(Imi, Imj))
          me_im = me_im + dvcell*0.5d0*sgn                                   &
          &                    *(dot_product(Rei, Imj) - dot_product(Imi, Rej))
        case (AX_X)
          me_re = me_re + dvcell*0.5d0                                       &
          &                    *(dot_product(Rei, Rej) + dot_product(Imi, Imj))
          me_im = me_im + dvcell*0.5d0                                       &
          &                    *(dot_product(Rei, Imj) - dot_product(Imi, Rej))
        case (AX_Y)
          me_re = me_re + dvcell*0.5d0*sgn                                   &
          &                    *(dot_product(Rei, Imj) - dot_product(Imi, Rej))
          me_im = me_im - dvcell*0.5d0*sgn                                   &
          &                    *(dot_product(Rei, Rej) + dot_product(Imi, Imj))
        end select
      enddo
    enddo
  end subroutine fullbox_angmom_me
!===============================================================================
  function Jm_fresh_element(mu, folded, srow, scol, r, c, st) result(val)
    !------------------------------------------------------------------------------
    ! Kernel evaluation of one element through the test's own global-index
    !  routing: spinors and the column gradient fetched here, dispatch to
    !  the kernel the object selects.
    !------------------------------------------------------------------------------
    integer, intent(in)  :: mu, srow, scol, r, c
    integer, intent(in)  :: st(0:8)
    logical, intent(in)  :: folded
    real(KIND=dp)        :: val

    real(KIND=dp) :: psi_i(mv, 4), psi_j(mv, 4)
    real(KIND=dp) :: derpsi_j(mv, 3, 4)
    integer :: gi, gj

    gi = st(srow) + r - 1
    gj = st(scol) + c - 1
    psi_j = HFPsi(:, :, gj)
    if (gi .eq. gj) then
      psi_i = psi_j
    else
      psi_i = HFPsi(:, :, gi)
    endif
    call Jm_fresh_derive(psi_j, gj, derpsi_j)
    val = Jm_oracle_kernel(mu, folded, psi_i, psi_j, derpsi_j)
  end function Jm_fresh_element
!===============================================================================
  function Jm_oracle_kernel(mu, folded, psi_i, psi_j, derpsi_j) result(val)
    !------------------------------------------------------------------------------
    ! Direct dispatch to the angmom kernel of the selected object - shared
    !  by the value, mirror-oracle and purity paths of this module.
    !------------------------------------------------------------------------------
    integer, intent(in)  :: mu
    logical, intent(in)  :: folded
    real(KIND=dp), intent(in) :: psi_i(:, :), psi_j(:, :), derpsi_j(:, :, :)
    real(KIND=dp)        :: val

    if (folded) then
      select case (mu)
      case (AX_X)
        val = angmom_xt_real(psi_i, psi_j, derpsi_j)
      case (AX_Y)
        val = angmom_yt_imag(psi_i, psi_j, derpsi_j)
      case default
        print *, 'Jm_oracle_kernel: folded components are AX_X and AX_Y.'
        call stp('Invalid folded component in the oracle.')
        val = 0.0d0
      end select
    else
      select case (mu)
      case (AX_X)
        val = angmom_x_real(psi_i, psi_j, derpsi_j)
      case (AX_Y)
        val = angmom_y_imag(psi_i, psi_j, derpsi_j)
      case (AX_Z)
        val = angmom_z_real(psi_i, psi_j, derpsi_j)
      case default
        print *, 'Jm_oracle_kernel: components are AX_X, AX_Y and AX_Z.'
        call stp('Invalid component in the oracle.')
        val = 0.0d0
      end select
    endif
  end function Jm_oracle_kernel
!===============================================================================
  subroutine Jm_fresh_derive(psi, glob, derpsi)
    !------------------------------------------------------------------------------
    ! Complete gradient of a stored spinor, the serial branch of
    !  Transfer_derpsi_complete: the stored derivatives when the run keeps
    !  them, otherwise the Derive_X/Y/Z_spwf primitives with the
    !  per-component reflection signs.
    !------------------------------------------------------------------------------
    real(KIND=dp), intent(in)  :: psi(:, :)
    integer, intent(in)        :: glob
    real(KIND=dp), intent(out) :: derpsi(:, :, :)

    if (store_derivatives) then
      derpsi = HFdpsi(:, :, :, glob)
    else
      call Derive_X_spwf(psi, sx(:, glob), derpsi(:, 1, :))
      call Derive_Y_spwf(psi, sy(:, glob), derpsi(:, 2, :))
      call Derive_Z_spwf(psi, sz(:, glob), derpsi(:, 3, :))
    endif
  end subroutine Jm_fresh_derive
!===============================================================================
  subroutine ang_twin_fill(srow, scol, A, params)
    !------------------------------------------------------------------------------
    ! block_filler callback of the arm-owned twin construction: the
    !  t-kernel fold through the test's own fresh transcription, one
    !  element per stored pair.
    !------------------------------------------------------------------------------
    integer, intent(in)   :: srow, scol
    real(KIND=dp), intent(out)  :: A(:, :)
    class(*), intent(in)  :: params

    integer :: r, c, m, n

    select type (params)
    type is (ang_twin_params)
      m = params%st(srow + 1) - params%st(srow)
      n = params%st(scol + 1) - params%st(scol)
      if (m .ne. size(A, 1) .or. n .ne. size(A, 2)) then
        print *, 'ang_twin_fill: params%st inconsistent with the block.'
        call stp('Sector bounds inconsistent in ang_twin_fill.')
      endif
      do c = 1, n
        do r = 1, m
          A(r, c) = Jm_fresh_element(params%mu, .true., srow, scol, r, c,        &
          &                      params%st)
        enddo
      enddo
    end select
  end subroutine ang_twin_fill
!===============================================================================
  subroutine arm_build_twin(mu, M)
    !------------------------------------------------------------------------------
    ! Construct the folded twin through the PUBLIC block machinery
    !  (init + fill_blocks + the arm's own t-kernel arithmetic): the
    !  retired eager entry point's object reconstruction. Class labels:
    !  the folded table of design 14.1 - both twins word 0; K_x
    !  HERMITIAN x STY_EVEN, K_y ANTIHERMITIAN x STY_ODD (herm the
    !  kind-completing label, D-1); kinds symmetric as the
    !  kernel-transpose oracle measures; the composite inherits the
    !  family tau.
    !------------------------------------------------------------------------------
    integer, intent(in)   :: mu
    type(sym_matrix), intent(out) :: M

    type(ang_twin_params) :: params
    integer :: herm_class, sty_class

    herm_class = HERMITIAN
    sty_class = STY_EVEN
    if (mu .eq. AX_Y) then
      herm_class = ANTIHERMITIAN
      sty_class = STY_ODD
    endif
    call M%init(HFBlocks_global, 0, herm_class, sty_class, TAU_ANGFAM,          &
    &           folded=.true., goodqn=goodqn)
    if (M%nblocks .eq. 0) return

    params%st = M%st
    params%mu = mu
    call M%fill_blocks(ang_twin_fill, params)
  end subroutine arm_build_twin
!===============================================================================
  function tr_conserving_basis() result(istr)
    !------------------------------------------------------------------------------
    ! .true. when the stored basis carries no -i signature states (even
    !  MOCCa block slots all empty): every TR-conserving build, including
    !  parity-broken ones (merged parity blocks, all q(Rz) = +i). Only
    !  such bases have T-image sectors for the retrieval contract.
    !------------------------------------------------------------------------------
    logical :: istr

    istr = (sum(HFBlocks_global((/2, 4, 6, 8/))) .eq. 0)
  end function tr_conserving_basis
!===============================================================================
  subroutine test_retrieval(ifail)
    !------------------------------------------------------------------------------
    ! The 15.8 retrieval checks in the J shape (design 14.3 amendment
    !  of 2026-11-20, 15.9 row 2), exercised in the builds that have
    !  T-image sectors - every TR-conserving build. Two judgements:
    !
    !  a. THE TAU BOTH-WAYS CHECK (mechanism 2): the J_z image-image
    !     blocks served by the accessor's sign flip against the direct
    !     honest computation on materialised T-images (TimeReverse on
    !     spinor and gradient - reversed slice by slice, since T is
    !     spinor-pointwise and commutes with the gradient). This check
    !     ARMS the algebra path for the family.
    !  b. Retrieval agreement of the fold (mechanism 3): the
    !     stored-x-image blocks of J_x/J_y served by the accessor
    !     against the arm-owned twin's diagonal blocks at the matching
    !     indices - the two packagings of the same matrix elements.
    !
    !  In TR-broken builds no sector is an image sector: nothing to
    !  exercise, pass with a report line.
    !------------------------------------------------------------------------------
    integer, intent(out) :: ifail

    integer :: s, sp, i, r, c, n1, n2, k, gi, gj, b, m1, m2
    real(KIND=dp) :: dev, ref, amax
    real(KIND=dp) :: psi_i(mv, 4), psi_j(mv, 4)
    real(KIND=dp) :: derpsi_j(mv, 3, 4)
    real(KIND=dp), allocatable :: A(:, :)
    logical :: checked
    type(sym_matrix) :: M

    ifail = 0
    if (.not. tr_conserving_basis()) then
      print *, '  retrieval: no image sectors in this build - not exercised.'
      return
    endif

    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! a. tau both-ways for J_z: for every stored sector s, the
    !    image-image counterpart block (IEOR(s,1), IEOR(s,1)) is served
    !    by the accessor and recomputed directly on the materialised
    !    images.
    dev = 0.0d0
    amax = 0.0d0
    checked = .false.
    do s = 0, 7
      if (HFBlocks_global(s + 1) .eq. 0) cycle
      sp = ieor(s, 1)
      checked = .true.
      call get_J_block(AX_Z, sp, sp, A)
      if (.not. allocated(A)) then
        ifail = 1
        print *, '  retrieval FAIL: no image-image block served for', sp
        cycle
      endif
      n1 = size(A, 1)
      n2 = size(A, 2)
      do c = 1, n2, max(1, n2/5)
        gj = sum(HFBlocks_global(1:s)) + c
        call Jm_fresh_derive(HFPsi(:, :, gj), gj, derpsi_j)
        psi_j = TimeReverse(HFPsi(:, :, gj))
        do k = 1, 3
          derpsi_j(:, k, :) = TimeReverse(derpsi_j(:, k, :))
        enddo
        do r = 1, n1, max(1, n1/5)
          gi = sum(HFBlocks_global(1:s)) + r
          psi_i = TimeReverse(HFPsi(:, :, gi))
          ref = Jm_oracle_kernel(AX_Z, .false., psi_i, psi_j, derpsi_j)
          dev = max(dev, abs(A(r, c) - ref))
          amax = max(amax, abs(ref))
        enddo
      enddo
      deallocate (A)
    enddo
    if (.not. checked) then
      print *, '  retrieval: tau both-ways: no stored sectors - skip.'
    else if (dev .gt. tol*max(1.0d0, amax)) then
      ifail = 1
      print *, '  retrieval FAIL: tau sign flip contradicts the direct',       &
      &         ' computation, max', dev
    else
      print *, '  retrieval: tau sign flip agrees with the honest kernels.',  &
      &         ' max', dev
    endif

    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! b. Fold agreement for J_x/J_y: for every diagonal twin block (s, s),
    !    the accessor's stored-x-image block (s, IEOR(s,1)) - the image
    !    columns relabelled by their stored partners - must agree with
    !    the twin's own storage.
    dev = 0.0d0
    amax = 0.0d0
    checked = .false.
    do i = AX_X, AX_Y
      call arm_build_twin(i, M)
      do b = 1, M%nblocks
        s = M%blocks(b)%srow
        sp = ieor(s, 1)
        if (HFBlocks_global(s + 1) .eq. 0) cycle
        checked = .true.
        call get_J_block(i, s, sp, A)
        if (.not. allocated(A)) then
          ifail = 1
          print *, '  retrieval FAIL: no stored-x-image block served for',     &
          &         s, sp
          cycle
        endif
        m1 = size(A, 1)
        m2 = size(A, 2)
        if (m1 .ne. size(M%blocks(b)%A, 1) .or.                                 &
        &   m2 .ne. size(M%blocks(b)%A, 2)) then
          ifail = 1
          print *, '  retrieval FAIL: fold block shape mismatch for', s
          cycle
        endif
        do c = 1, m2, max(1, m2/5)
          do r = 1, m1, max(1, m1/5)
            ref = M%blocks(b)%A(r, c)
            dev = max(dev, abs(A(r, c) - ref))
            amax = max(amax, abs(ref))
          enddo
        enddo
        deallocate (A)
      enddo
    enddo
    if (.not. checked) then
      print *, '  retrieval: fold agreement: no twin blocks - skip.'
    else if (dev .gt. tol*max(1.0d0, amax)) then
      ifail = 1
      print *, '  retrieval FAIL: accessor block differs from the twin',      &
      &         ' construction, max', dev
    else
      print *, '  retrieval: accessor agrees with the twin construction.',    &
      &         ' max', dev
    endif
  end subroutine test_retrieval
!===============================================================================

end module angular_momentum_testing
