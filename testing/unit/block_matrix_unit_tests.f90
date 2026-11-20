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
program block_matrix_unit_tests
!==============================================================================
  !
  ! Unit tests for phases 1 and 2 of the symmetry-block matrix module:
  !  construction of the block list from MOCCa's HFBlocks bookkeeping
  !  and population of the stored blocks through callbacks.
  !
  ! Every operator handled by the type is classifiable (design 2a): it is
  !  Hermitian or anti-Hermitian, and it is STY_EVEN or STY_ODD under the
  !  y-time-simplex. init takes both classes as mandatory inputs and
  !  derives the storage kind as kind = herm * sty, so all construction
  !  and population tests loop over the four (herm, sty) combinations.
  !
  ! For a suite of HFBlocks patterns (EV8-like, time-reversal-broken, sparse,
  !  degenerate), every operator symmetry word 0..7 and every physical
  !  class combination, the stored block list is verified against the
  !  design invariants:
  !
  !   *) selection rules: each stored block couples sectors (srow, scol)
  !      with srow = IEOR(scol, symmetry) and nothing else is stored
  !   *) partner skipping: at most one block per sector pair is stored,
  !      the skipped one satisfies scol > srow
  !   *) block shapes, offsets and zero-initialisation
  !   *) sector bounds array st(0:8) and the global size n
  !   *) sector_of mapping, including empty sectors and out-of-range indices
  !
  ! Phase 2 adds the population machinery (design section 5):
  !   *) fill_blocks invokes the callback once per stored block, for the
  !      stored sector pairs only, with the block's sector words and the
  !      right storage shape, and overwrites the entire block
  !   *) the callback hands the unlimited-polymorphic params through
  !      verbatim (an integer offset)
  !   *) an empty block list leaves the callback uncalled
  !
  ! The good-quantum-number mask (design section 2b) adds:
  !   *) the mask is stored and the omitted mask defaults to 15
  !   *) masked pairing: blocks obey srow = IEOR(scol, AND(word, mask))
  !   *) equivalence: a masked construction stores exactly the block list
  !      of the unmasked construction with the effective word
  !   *) the parity-broken pattern with an odd-parity operator: mask 5
  !      merges to diagonal blocks, mask 7 pairs into empty sectors and
  !      stores nothing
  !
  ! The physical classes (design section 2a) add:
  !   *) the characteristic is stored verbatim and the storage kind is
  !      derived as kind = herm * sty
  !   *) no block structure: bounds, block count, sector pairs and
  !      shapes are identical for all four class combinations
  !   *) the mask is stored verbatim too - bit 3 is the consumer's
  !      completeness criterion, not a constructor gate
  !
  ! The partner-reconstructing dense expansion (dense) adds:
  !   *) shape (n, n) of the expansion, including the empty basis
  !   *) every stored block sits at its global position with the exact
  !      callback value; partners are reconstructed with the derived
  !      kind sign at the transposed position; everything else
  !      is exactly zero - verified against a hand-scattered reference
  !      over the full cross product of patterns, symmetry words, class
  !      combinations and good-quantum-number masks
  !
  ! The basis rotations (design section 6a) add:
  !   *) the identity rotation is an exact no-op, built both by init
  !      and extracted from a full identity array by init_from_dense
  !   *) a rotated matrix keeps every metadata field (symmetry word,
  !      kind, characteristic, goodqn, n, st) and its block list
  !      still satisfies the construction invariants
  !   *) the rotated storage agrees with the dense reference
  !      C^T dense(M) C within reassociation tolerance
  !   *) the round trip through the sector-wise transposed rotation
  !      recovers the original numbers
  !   *) the traces of the diagonal blocks are invariant (word 0)
  !   *) for a partner-consistent operand, filled block-wise as a
  !      symmetric or antisymmetric function of the global indices,
  !      the dense expansion of the rotated matrix obeys
  !      D(gj, gi) = kind * D(gi, gj) for both kinds
  ! The structure guards of the rotation constructors are exercised by
  !  the companion program block_matrix_refusal_tests.f90, driven by
  !  the run script - they stop the process and cannot be caught
  !  inside this suite.
  !
  ! The tests are self-contained: no other MOCCa module is used. Compile
  !  together with src_orig/block_matrix_module.f90 and run.
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! The suite consists of eight routines (plus the check helper that
  !  registers every individual condition, and the probed callbacks the
  !  fill and rotation tests use):
  !   test_construction : construction invariants over the full cross
  !                      product of patterns, symmetry words and class
  !                      combinations
  !   test_reinit      : repeated init on the same object
  !   test_sector_of   : global-index-to-sector mapping
  !   test_fill_blocks : block-wise population through the callback
  !   test_goodqn      : the good-quantum-number mask of design 2b
  !   test_sty         : the physical classes of design 2a
  !   test_dense_expansion : partner-reconstructing dense expansion
  !   test_rotation   : basis rotations of design section 6a
  !==============================================================================
  use block_matrix_module

  implicit none

  integer, parameter :: NPAT = 7

  integer :: patterns(NPAT, 8)
  integer :: nfail, ntest

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  ! State shared with the probed callbacks of the fill tests: call counters,
  !  the map of sector pairs the block callback was invoked for, the counts
  !  of the pattern currently under test (the callbacks check the shapes
  !  they receive against it) and the label of the current case.
  integer :: nblockcalls           ! invocations of fill_block_probed
  integer :: seen(0:7, 0:7)        ! sector pairs the block callback served
  integer :: cur_counts(0:7)      ! current pattern, for shape checks
  character(len=160) :: cur_label  ! current case label
  integer :: cur_st(0:8)           ! sector bounds of the current case,
                                  !  for the global-index-value callback
  integer :: cur_knd              ! derived kind of the current case,
                                  !  selects the (anti)symmetric partner
                                  !  callback of the rotation tests

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! HFBlocks patterns: MOCCa ordering, per the table in the module header.
  !   1: EV8-like bookkeeping (only blocks 1, 3, 5, 7 populated)
  !   2: time-reversal broken (all blocks populated)
  !   3: sparse pattern with interior empty sectors
  !   4: empty leading sectors
  !   5: single populated sector
  !   6: completely empty basis
  !   7: parity-broken bookkeeping (blocks 1, 2, 5, 6 populated only)
  patterns(1, :) = [12, 0, 9, 0, 7, 0, 5, 0]
  patterns(2, :) = [12, 11, 9, 8, 7, 6, 5, 4]
  patterns(3, :) = [3, 0, 0, 2, 5, 0, 1, 0]
  patterns(4, :) = [0, 4, 2, 0, 0, 6, 0, 3]
  patterns(5, :) = [0, 0, 0, 0, 0, 0, 9, 0]
  patterns(6, :) = [0, 0, 0, 0, 0, 0, 0, 0]
  patterns(7, :) = [8, 6, 0, 0, 5, 4, 0, 0]

  nfail = 0
  ntest = 0

  call test_construction()
  call test_reinit()
  call test_sector_of()
  call test_fill_blocks()
  call test_goodqn()
  call test_sty()
  call test_dense_expansion()
  call test_rotation()

  print '(''block_matrix unit tests: '', i0, '' checks, '', i0, '' failures'')', &
  &     ntest, nfail
  if (nfail .gt. 0) then
    print *, 'FAILED'
    stop 1
  endif
  print *, 'All block_matrix unit tests passed.'

contains

  !----------------------------------------------------------------------------
  subroutine check(cond, label)
    !---------------------------------------------------------------------------
    ! Register the outcome of a single boolean test condition in the global
    !  counters, printing the label only when the condition fails.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !   cond  : the condition under test
    !   label : identifying text for this check, printed upon failure
    !---------------------------------------------------------------------------
    logical, intent(in) :: cond
    character(len=*), intent(in) :: label

    ntest = ntest + 1
    if (.not. cond) then
      nfail = nfail + 1
      print '(''FAIL: '', a)', label
    endif
  end subroutine check
  !----------------------------------------------------------------------------
  subroutine test_construction()
    !---------------------------------------------------------------------------
    ! Verify every construction invariant of init over the full cross
    !  product of the NPAT HFBlocks patterns, all operator symmetry words
    !  0..7 and all four (herm, sty) class combinations.
    !
    ! For each combination a fresh matrix is built and checked for:
    !   *) stored metadata: symmetry word, derived kind, characteristic,
    !      total size, nblocks <= 8
    !   *) the bounds array st: st(0) = 1, st(8) = n+1, and the sector
    !      sizes reproduce the differences st(s+1) - st(s)
    !   *) the number of stored blocks: one per populated ordered sector
    !      pair, with the s > srow partner always skipped (the derived
    !      kind is never zero)
    !   *) per stored block: srow = IEOR(scol, word), uniqueness of the
    !      (srow, scol) pair, both sectors populated, A allocated with the
    !      sector-pair shape, zero-filled, correct offsets is/js
    !   *) completeness: every ordered sector pair that should be stored
    !      is present in the list
    !   *) hygiene: block slots beyond nblocks carry no allocation
    !---------------------------------------------------------------------------
    type(sym_matrix) :: M
    integer :: p, word, h, s, sc, knd, b, sp, counts(0:7)
    integer :: expected, sumexp
    integer :: herms(2), stys(2)
    logical :: found, occupied(0:7, 0:7)
    character(len=120) :: label

    herms = [HERMITIAN, ANTIHERMITIAN]
    stys = [STY_EVEN, STY_ODD]

    do p = 1, NPAT
      counts(0:7) = patterns(p, 1:8)
      do word = 0, 7
        do h = 1, 2
          do sc = 1, 2
            knd = herms(h)*stys(sc)
            write (label, '(a,i0,a,i0,a,i0,a,i0)') 'pattern ', p,                &
            &  ' word ', word, ' herm ', herms(h), ' sty ', stys(sc)
            call M%init(patterns(p, :), word, herms(h), stys(sc))

            call check(M%symmetry .eq. word, trim(label)//': symmetry word')
            call check(M%kind .eq. knd, trim(label)//': derived kind')
            call check(M%sty .eq. stys(sc), trim(label)//': characteristic')
            call check(M%n .eq. sum(counts), trim(label)//': total size n')
            call check(M%nblocks .ge. 0 .and. M%nblocks .le. 8,                    &
            &          trim(label)//': nblocks in 0..8')

          ! Bounds array: st(s) is the first global index of sector s,
          !  st(0) = 1, st(8) = n + 1 and sector sizes reproduce st.
          call check(M%st(0) .eq. 1, trim(label)//': st(0) = 1')
          call check(M%st(8) .eq. M%n + 1, trim(label)//': st(8) = n+1')
          do s = 0, 7
            call check(M%st(s + 1) - M%st(s) .eq. counts(s),                    &
            &          trim(label)//': sector bounds / sizes')
          enddo

          ! Expected number of stored blocks: one per populated ordered
          !  pair, with partner skipping (the derived kind is always
          !  nonzero, so exactly one block per pair survives).
          sumexp = 0
          do s = 0, 7
            sp = ieor(s, word)
            if (counts(s) .eq. 0) cycle
            if (counts(sp) .eq. 0) cycle
            if (s .gt. sp) cycle
            sumexp = sumexp + 1
          enddo
          call check(M%nblocks .eq. sumexp, trim(label)//': expected nblocks')

          ! Every stored block obeys the XOR selection rule, is unique,
          !  points at the right offset, is zero-filled and carries the
          !  right shape. Partner blocks stay away.
          occupied = .false.
          do b = 1, M%nblocks
            associate (Blk => M%blocks(b))
              write (label, '(a,i0,a,i0,a,i0,a,i0,a,i0)')                       &
              &  'pattern ', p, ' word ', word, ' herm ', herms(h),              &
              &  ' sty ', stys(sc), ' block ', b

              call check(Blk%srow .eq. ieor(Blk%scol, word),                    &
              &          trim(label)//': srow = IEOR(scol, word)')
              call check(Blk%scol .ge. 0 .and. Blk%scol .le. 7,                 &
              &          trim(label)//': scol in 0..7')
              call check(.not. occupied(Blk%srow, Blk%scol),                    &
              &          trim(label)//': unique sector pair')
              occupied(Blk%srow, Blk%scol) = .true.

              ! Partner-skip: a stored block never sits on the skipped
              !  side of its pair.
              call check(Blk%scol .le. Blk%srow,                               &
              &          trim(label)//': partner skipped, scol <= srow')

              call check(counts(Blk%srow) .gt. 0 .and. counts(Blk%scol) .gt. 0,&
              &          trim(label)//': populated sectors')
              call check(allocated(Blk%A), trim(label)//': A allocated')
              if (allocated(Blk%A)) then
                call check(size(Blk%A, 1) .eq. counts(Blk%srow),                &
                &          trim(label)//': rows = size of row sector')
                call check(size(Blk%A, 2) .eq. counts(Blk%scol),                &
                &          trim(label)//': cols = size of col sector')
                call check(all(abs(Blk%A) .eq. 0.0d0),                          &
                &          trim(label)//': zero-filled')
              endif
              call check(Blk%is .eq. M%st(Blk%srow),                             &
              &          trim(label)//': row offset is')
              call check(Blk%js .eq. M%st(Blk%scol),                             &
              &          trim(label)//': col offset js')
            end associate
          enddo

          ! Every populated ordered pair that should be stored is stored.
          do s = 0, 7
            sp = ieor(s, word)
            expected = 0
            if (counts(s) > 0 .and. counts(sp) > 0) then
              if (s .le. sp) expected = 1
            endif
            found = .false.
            do b = 1, M%nblocks
              if (M%blocks(b)%srow .eq. sp .and. M%blocks(b)%scol .eq. s)       &
              &  found = .true.
            enddo
            write (label, '(a,i0,a,i0,a,i0,a,i0,a,i0)') 'pattern ', p,         &
            &  ' word ', word, ' herm ', herms(h), ' sty ', stys(sc),           &
            &  ' pair (row,col) ', sp
            call check(merge(1, 0, found) .eq. expected,                         &
            &          trim(label)//': expected pair present')
          enddo

          ! Unused block slots must be cleaned out.
          do b = M%nblocks + 1, 8
            write (label, '(a,i0,a,i0,a,i0,a,i0,a,i0)') 'pattern ', p,          &
            &  ' word ', word, ' herm ', herms(h), ' sty ', stys(sc),           &
            &  ' unused slot ', b
            call check(.not. allocated(M%blocks(b)%A),                          &
            &          trim(label)//': A not allocated')
          enddo
          enddo
        enddo
      enddo
    enddo
  end subroutine test_construction
  !----------------------------------------------------------------------------
  subroutine test_reinit()
    !---------------------------------------------------------------------------
    ! Repeated init on the same object must discard the previous state
    !  cleanly.
    !
    ! Scenario: start from the TR-broken pattern with a symmetric matrix
    !  of word 0, where the block count equals the number of populated
    !  sectors. Re-init with the EV8 pattern and word 7: that word pairs
    !  every populated sector {0,2,4,6} with an empty partner ({7,5,3,1}),
    !  so the block list must come out empty and no block storage may
    !  survive. Finally rebuild with a word-5 matrix of an anti-Hermitian
    !  operator to verify the list and the symmetry word are restored.
    !---------------------------------------------------------------------------
    type(sym_matrix) :: M
    integer   :: b, nexp
    character(len=120) :: label

    call M%init(patterns(2, :), 0, HERMITIAN, STY_EVEN)
    nexp = count(patterns(2, :) > 0)
    call check(M%nblocks .eq. nexp, 'reinit: symmetric word-0 block count')

    call M%init(patterns(1, :), 7, HERMITIAN, STY_ODD)
    call check(M%n .eq. sum(patterns(1, :)), 'reinit: n updated')
    call check(M%nblocks .ge. 0, 'reinit: nblocks after shrink')
    ! Word 7 couples every populated sector to a different partner; with the
    !  EV8 pattern only sectors {0,2,4,6} are populated, hence pairs
    !  {0,7},{2,5},{4,3},{6,1}: none has both partners populated, so zero
    !  blocks survive.
    do b = 1, 8
      write (label, '(a,i0)') 'reinit block ', b
      call check(.not. allocated(M%blocks(b)%A), trim(label)//': A deallocated')
    enddo
    write (label, '(a,i0,a,i0)') 'reinit: word ', 7, ' EV8 pattern'
    call check(M%nblocks .eq. 0, trim(label)//': no surviving pairs')

    ! And back to the full pattern, checking the list is rebuilt.
    call M%init(patterns(2, :), 5, ANTIHERMITIAN, STY_EVEN)
    call check(M%nblocks .gt. 0, 'reinit: blocks rebuilt')
    call check(M%symmetry .eq. 5, 'reinit: word updated')
    call check(M%kind .eq. KIND_ANTISYMMETRIC, 'reinit: derived kind updated')
  end subroutine test_reinit
  !----------------------------------------------------------------------------
  subroutine test_sector_of()
    !---------------------------------------------------------------------------
    ! Verify sector_of, the global-index-to-sector map, on the sparse
    !  pattern with interior empty sectors.
    !
    ! Coverage: the start, middle and end of every populated sector (the
    !  empty-sector boundary cases included), out-of-range indices (which
    !  must return the sentinel -1) and an exhaustive sweep over all valid
    !  indices, each of which must fall inside the bounds of the sector it
    !  is mapped to.
    !---------------------------------------------------------------------------
    type(sym_matrix) :: M
    integer :: b, s, i
    logical :: inrange

    ! TR-broken pattern with interior empty sectors.
    call M%init(patterns(3, :), 0, HERMITIAN, STY_EVEN)

    do b = 1, M%nblocks
      s = M%blocks(b)%scol
      ! Start, middle and end of the sector must map back to the sector,
      !  also when neighbouring sectors are empty.
      do i = M%st(s), M%st(s) + max(0, M%st(s + 1) - M%st(s)) - 1,              &
      &    max(1, M%st(s + 1) - M%st(s) - 1)
        call check(M%sector_of(i) .eq. s, 'sector_of: interior points')
      enddo
    enddo

    ! Out-of-range indices return the sentinel.
    call check(M%sector_of(-1) .eq. -1, 'sector_of: negative index sentinel')
    if (M%n .gt. 0) then
      call check(M%sector_of(M%n) .ne. -1, 'sector_of: i = n valid')
      call check(M%sector_of(M%n + 1) .eq. -1, 'sector_of: i = n + 1 sentinel')
    endif

    ! Every valid index maps to the sector containing it.
    inrange = .true.
    do i = 1, M%n
      s = M%sector_of(i)
      inrange = inrange .and. s .ge. 0 .and. s .le. 7
      inrange = inrange .and. i .ge. M%st(s) .and. i .lt. M%st(s + 1)
    enddo
    call check(inrange, 'sector_of: all indices in range')
  end subroutine test_sector_of
  !----------------------------------------------------------------------------
  subroutine test_fill_blocks()
    !---------------------------------------------------------------------------
    ! Verify block population through the block-wise callback over the
    !  full cross product of patterns, symmetry words and kinds.
    !
    ! The probed callback overwrites each block with values identifying
    !  the sector pair and the local element, offset by the integer
    !  carried in params. Verified per case:
    !   *) the callback runs exactly once per stored block, and only for
    !      the stored sector pairs, which must match the blocks' srow/scol
    !   *) the storage handed to the callback has the sector-pair shape
    !   *) the stored content equals the callback's values everywhere:
    !      the fill overwrites the pollution written before the call
    !   *) a degenerate empty block list leaves the callback uncalled
    !---------------------------------------------------------------------------
    type(sym_matrix) :: M
    integer :: p, word, h, s, sc, b, r, c, mrow, ncol, sumseen
    integer :: herms(2), stys(2)
    real(8) :: offset
    logical :: match
    character(len=160) :: label, blabel

    herms = [HERMITIAN, ANTIHERMITIAN]
    stys = [STY_EVEN, STY_ODD]
    offset = 7.0d0

    do p = 1, NPAT
      do word = 0, 7
        do h = 1, 2
          do sc = 1, 2
          write (label, '(a,i0,a,i0,a,i0,a,i0)') 'pattern ', p, ' word ', word,  &
          &  ' herm ', herms(h), ' sty ', stys(sc)
          cur_counts(0:7) = patterns(p, 1:8)
          cur_label = trim(label)

          call M%init(patterns(p, :), word, herms(h), stys(sc))

          ! Pollution: every fill must overwrite regardless of prior content.
          do b = 1, M%nblocks
            M%blocks(b)%A = -12345.0d0
          enddo

          nblockcalls = 0
          seen = 0
          call M%fill_blocks(fill_block_probed, int(offset))

          ! Exactly the stored blocks, each once: the counter bounds the
          !  total, the per-block seen entries pin the sector pairs.
          call check(nblockcalls .eq. M%nblocks,                               &
          &          trim(label)//': one call per stored block')
          sumseen = sum(seen)
          call check(sumseen .eq. M%nblocks,                                    &
          &          trim(label)//': no calls outside the list')
          do b = 1, M%nblocks
            write (blabel, '(a,i0)') trim(label)//' block ', b
            call check(seen(M%blocks(b)%srow, M%blocks(b)%scol) .eq. 1,       &
            &          trim(blabel)//': called for its sector pair')
          enddo

          ! Content: every element must carry the callback's value.
          do b = 1, M%nblocks
            mrow = size(M%blocks(b)%A, 1)
            ncol = size(M%blocks(b)%A, 2)
            match = .true.
            do r = 1, mrow
              do c = 1, ncol
                match = match                                                &
                &  .and. (M%blocks(b)%A(r, c) .eq.                            &
                &          block_val(M%blocks(b)%srow, M%blocks(b)%scol,      &
                &                      r, c, offset))
              enddo
            enddo
            write (blabel, '(a,i0)') trim(label)//' block ', b
            call check(match, trim(blabel)//': content overwritten everywhere')
          enddo
          enddo
        enddo
      enddo
    enddo
  end subroutine test_fill_blocks
  !----------------------------------------------------------------------------
  subroutine test_goodqn()
    !---------------------------------------------------------------------------
    ! Verify the good-quantum-number mask of design section 2b.
    !
    ! Over the full cross product of patterns, physical words, kinds and
    !  masks (4 = only Tz good, 5 = Tz and q(Rz) good, 7 = the sector
    !  gradings all good, without any time-simplex bit):
    !   *) the mask is stored on the object and the physical word survives
    !   *) blocks pair as srow = IEOR(scol, AND(word, goodqn)) and the
    !      block count follows the masked pairing as in test_construction
    !   *) a masked construction stores exactly the block list of the
    !      unmasked construction with the effective word
    !
    ! Dedicated scenario - the core of design section 2b: the parity-broken
    !  pattern with an odd-parity operator (word 2). Mask 5 (parity not
    !  good) gives the four diagonal blocks of the merged sectors
    !  0, 1, 4, 5; the all-good mask pairs every populated sector into the
    !  empty odd-parity sectors and stores nothing at all.
    !---------------------------------------------------------------------------
    type(sym_matrix) :: M, Meff
    integer :: p, word, h, s, sc, q, masks(4), b, sp, sumexp
    integer :: herms(2), stys(2)
    character(len=160) :: label

    herms = [HERMITIAN, ANTIHERMITIAN]
    stys = [STY_EVEN, STY_ODD]
    masks = [4, 5, 7, 13]

    do p = 1, NPAT
      do word = 0, 7
        do q = 1, 4
          do h = 1, 2
            do sc = 1, 2
            write (label, '(a,i0,a,i0,a,i0,a,i0,a,i0)') 'pattern ', p,         &
            &  ' word ', word, ' mask ', masks(q), ' herm ', herms(h),          &
            &  ' sty ', stys(sc)
            call M%init(patterns(p, :), word, herms(h), stys(sc),                &
            &          goodqn=masks(q))

            call check(M%goodqn .eq. masks(q), trim(label)//': mask stored')
            call check(M%symmetry .eq. word,                                    &
            &          trim(label)//': physical word stored')

            ! The masked pairing word governs the block list.
            sumexp = 0
            do s = 0, 7
              sp = ieor(s, iand(word, masks(q)))
              if (patterns(p, s + 1) .eq. 0) cycle
              if (patterns(p, sp + 1) .eq. 0) cycle
              if (s .gt. sp) cycle
              sumexp = sumexp + 1
            enddo
            call check(M%nblocks .eq. sumexp,                                   &
            &          trim(label)//': expected nblocks')
            do b = 1, M%nblocks
              call check(M%blocks(b)%srow .eq.                                  &
              &          ieor(M%blocks(b)%scol, iand(word, masks(q))),           &
              &          trim(label)//': srow = IEOR(scol, masked word)')
            enddo

            ! Equivalence: the block list equals that of the unmasked
            !  matrix built with the effective word.
            call Meff%init(patterns(p, :), iand(word, masks(q)), herms(h),      &
            &               stys(sc))
            call check(Meff%nblocks .eq. M%nblocks,                             &
            &          trim(label)//': effective-word nblocks')
            do b = 1, min(M%nblocks, Meff%nblocks)
              call check(Meff%blocks(b)%srow .eq. M%blocks(b)%srow,             &
              &          trim(label)//': effective-word srow')
              call check(Meff%blocks(b)%scol .eq. M%blocks(b)%scol,             &
              &          trim(label)//': effective-word scol')
              call check(size(Meff%blocks(b)%A, 1) .eq. size(M%blocks(b)%A, 1)  &
              &     .and. size(Meff%blocks(b)%A, 2) .eq. size(M%blocks(b)%A, 2), &
              &          trim(label)//': effective-word shapes')
            enddo
            enddo
          enddo
        enddo
      enddo
    enddo

    ! Default: the omitted mask is 15 (all good, incl. the time-simplex
    !  bit of design 2a), and it equals an explicit 15.
    call M%init(patterns(2, :), 5, HERMITIAN, STY_EVEN)
    call check(M%goodqn .eq. 15, 'goodqn: default mask is 15')
    call Meff%init(patterns(2, :), 5, HERMITIAN, STY_EVEN, goodqn=15)
    call check(Meff%nblocks .eq. M%nblocks, 'goodqn: default nblocks equal')

    ! Core scenario of design section 2b: the parity-broken pattern with
    !  an odd-parity operator (word 2). Parity not good (mask 5): diagonal
    !  blocks inside the merged sectors 0, 1, 4, 5.
    call M%init(patterns(7, :), 2, HERMITIAN, STY_EVEN, goodqn=5)
    call check(M%nblocks .eq. 4, 'goodqn: merged diagonal blocks')
    do b = 1, M%nblocks
      call check(M%blocks(b)%srow .eq. M%blocks(b)%scol,                        &
      &          'goodqn: diagonal within merged blocks')
      call check(patterns(7, M%blocks(b)%srow + 1) .gt. 0,                      &
      &          'goodqn: diagonal block is populated')
    enddo

    ! All quantum numbers good (mask 7): the physical word pairs every
    !  populated sector into the empty odd-parity sectors and nothing is
    !  stored - empty, which is where a missing mask would silently go
    !  wrong.
    call M%init(patterns(7, :), 2, HERMITIAN, STY_EVEN, goodqn=7)
    call check(M%nblocks .eq. 0, 'goodqn: unmasked word gives empty list')
  end subroutine test_goodqn
  !----------------------------------------------------------------------------
  subroutine test_sty()
    !---------------------------------------------------------------------------
    ! Verify the physical classes of design section 2a: the mandatory
    !  Hermiticity and y-time-simplex characteristic, from which init
    !  derives the storage kind.
    !
    ! Over the full cross product of patterns, symmetry words and the
    !  four (herm, sty) class combinations (mask 15):
    !   *) the characteristic is stored verbatim
    !   *) the storage kind is derived as kind = herm * sty
    !   *) no block structure: bounds, block count, sector pairs and
    !      shapes of the block list are identical for all four class
    !      combinations - the classes select the stored HALF (the real
    !      part or the imaginary part of the matrix elements), never a
    !      block
    !
    ! Completeness bookkeeping:
    !   *) goodqn bit 3 is stored verbatim - it is the consumer's
    !      completeness criterion (the matrix is the entire operator
    !      iff the bit is set), not a constructor gate
    !---------------------------------------------------------------------------
    type(sym_matrix) :: M, Mref
    integer :: p, word, h, sc, b
    integer :: herms(2), stys(2)
    character(len=160) :: label

    herms = [HERMITIAN, ANTIHERMITIAN]
    stys = [STY_EVEN, STY_ODD]

    do p = 1, NPAT
      do word = 0, 7
        do h = 1, 2
          do sc = 1, 2
            write (label, '(a,i0,a,i0,a,i0,a,i0)') 'pattern ', p,              &
            &  ' word ', word, ' herm ', herms(h), ' sty ', stys(sc)
            call Mref%init(patterns(p, :), word, HERMITIAN, STY_EVEN,           &
            &              goodqn=15)
            call M%init(patterns(p, :), word, herms(h), stys(sc),               &
            &          goodqn=15)

            ! The characteristic is stored verbatim, the kind derived.
            call check(M%sty .eq. stys(sc),                                      &
            &          trim(label)//': characteristic stored')
            call check(M%kind .eq. herms(h)*stys(sc),                            &
            &          trim(label)//': kind = herm * sty')

            ! No block structure: the classes never touch the block
            !  list - they select the stored half, not a block.
            call check(M%n .eq. Mref%n, trim(label)//': n unchanged')
            call check(all(M%st .eq. Mref%st),                                  &
            &          trim(label)//': bounds unchanged')
            call check(M%nblocks .eq. Mref%nblocks,                             &
            &          trim(label)//': nblocks unchanged')
            do b = 1, min(M%nblocks, Mref%nblocks)
              call check(M%blocks(b)%srow .eq. Mref%blocks(b)%srow,             &
              &          trim(label)//': srow unchanged')
              call check(M%blocks(b)%scol .eq. Mref%blocks(b)%scol,             &
              &          trim(label)//': scol unchanged')
              call check(size(M%blocks(b)%A, 1) .eq. size(Mref%blocks(b)%A, 1)  &
              &     .and. size(M%blocks(b)%A, 2) .eq. size(Mref%blocks(b)%A, 2),&
              &          trim(label)//': shapes unchanged')
            enddo
          enddo
        enddo
      enddo
    enddo

    ! The mask is stored verbatim, bit 3 included: the characteristic is
    !  not gated by it (a classifiable operator keeps its classes in
    !  every build; bit 3 only tells consumers whether the stored half
    !  is the complete operator).
    call M%init(patterns(2, :), 5, ANTIHERMITIAN, STY_ODD, goodqn=7)
    call check(M%sty .eq. STY_ODD, 'sty: bit 3 clear keeps characteristic')
    call check(M%kind .eq. KIND_SYMMETRIC,                                      &
    &          'sty: bit 3 clear keeps derived kind')
    call check(M%goodqn .eq. 7, 'sty: bit 3 stored clear (mask 7)')
    call M%init(patterns(2, :), 5, ANTIHERMITIAN, STY_ODD, goodqn=13)
    call check(M%goodqn .eq. 13, 'sty: bit 3 stored set (mask 13)')
  end subroutine test_sty
  !----------------------------------------------------------------------------
  subroutine test_dense_expansion()
    !---------------------------------------------------------------------------
    ! Verify the partner-reconstructing dense expansion over the full
    !  cross product of patterns, symmetry words, class combinations and
    !  masks 4, 5, 7, 13.
    !
    ! Every matrix is populated through the probed block callback with
    !  deterministic values and its dense expansion compared element by
    !  element against a hand-scattered reference:
    !   *) shape (n, n), including the empty basis
    !   *) every stored block sits at its global position carrying the
    !      callback's exact value
    !   *) partners are reconstructed at the transposed position with
    !      the derived kind sign herm * sty (bitwise, since the transpose
    !      reassigns the same values); diagonal blocks are not mirrored
    !   *) every element outside stored blocks and their partners is
    !      exactly zero, so selection-rule holes survive the expansion
    !---------------------------------------------------------------------------
    type(sym_matrix) :: M
    real(8), allocatable :: D(:, :), Rref(:, :)
    integer :: p, word, h, sc, q, masks(4), b, r, c, mrow, ncol
    integer :: i, j, knd
    integer :: herms(2), stys(2)
    real(8) :: offset
    logical :: match
    character(len=160) :: label

    herms = [HERMITIAN, ANTIHERMITIAN]
    stys = [STY_EVEN, STY_ODD]
    masks = [4, 5, 7, 13]
    offset = 11.0d0

    do p = 1, NPAT
      do word = 0, 7
        do h = 1, 2
          do sc = 1, 2
          do q = 1, 4
            knd = herms(h)*stys(sc)
            write (label, '(a,i0,a,i0,a,i0,a,i0,a,i0)') 'pattern ', p,         &
            &  ' word ', word, ' herm ', herms(h), ' sty ', stys(sc),           &
            &  ' mask ', masks(q)
            cur_counts(0:7) = patterns(p, 1:8)
            cur_label = trim(label)

            call M%init(patterns(p, :), word, herms(h), stys(sc),               &
            &          goodqn=masks(q))
            call M%fill_blocks(fill_block_probed, int(offset))

            D = M%dense()

            ! Shape: always (n, n), also for the empty basis.
            call check(size(D, 1) .eq. M%n, trim(label)//': shape (n, n)')
            call check(size(D, 2) .eq. M%n, trim(label)//': shape (n, n)')

            ! Hand-scattered reference: stored blocks at their global
            !  positions, partners with the derived kind sign, zero
            !  everywhere else.
            allocate(Rref(M%n, M%n))
            Rref = 0.0d0
            do b = 1, M%nblocks
              mrow = size(M%blocks(b)%A, 1)
              ncol = size(M%blocks(b)%A, 2)
              do r = 1, mrow
                do c = 1, ncol
                  Rref(M%blocks(b)%is + r - 1, M%blocks(b)%js + c - 1)        &
                  & = block_val(M%blocks(b)%srow, M%blocks(b)%scol, r, c,    &
                  &             offset)
                  if (M%blocks(b)%srow .ne. M%blocks(b)%scol) then
                    Rref(M%blocks(b)%js + c - 1, M%blocks(b)%is + r - 1)      &
                    & = knd*block_val(M%blocks(b)%srow,                        &
                    &                M%blocks(b)%scol, r, c, offset)
                  endif
                enddo
              enddo
            enddo

            ! The expansion equals the reference exactly.
            match = .true.
            do i = 1, M%n
              do j = 1, M%n
                match = match .and. (D(i, j) .eq. Rref(i, j))
              enddo
            enddo
            call check(match, trim(label)//': dense equals reference')
            deallocate(Rref)
          enddo
          enddo
        enddo
      enddo
    enddo
  end subroutine test_dense_expansion
  !----------------------------------------------------------------------------
  subroutine test_rotation()
    !---------------------------------------------------------------------------
    ! Verify the basis rotations of design section 6a over the full cross
    !  product of patterns, symmetry words, class combinations and a
    !  representative set of good-quantum-number masks (15 unmasked, the
    !  parity-masked 5, 7 and 13 of the dense tests - so the unmasked
    !  pairing and the diagonal parity-broken case are both exercised).
    !
    ! Every case fills the matrix through the probed block callback with
    !  deterministic values and builds a random orthogonal rotation (one
    !  hashed Gram-Schmidt block per populated sector, handed through
    !  set_sector). Judgments, per design section 8:
    !   *) the identity rotation is an exact no-op, built both by init
    !      and extracted from a full identity array by init_from_dense
    !   *) a rotated matrix keeps every metadata field and its block
    !      list still satisfies the construction invariants
    !   *) the rotated storage agrees with the dense reference
    !      C^T dense(M) C within reassociation tolerance
    !   *) the round trip through the sector-wise transposed rotation
    !      recovers the original numbers
    !   *) the traces of the diagonal blocks are invariant (word 0)
    !
    ! A closing judgment per case, on a partner-consistent operand
    !  filled through a block callback as a symmetric (kind +1) or antisymmetric
    !  (kind -1) function of the global indices: after rotation the
    !  dense expansion obeys D(gj, gi) = kind * D(gi, gj) - the kind
    !  partner relation survives the rotation, for both kinds.
    !---------------------------------------------------------------------------
    type(sym_matrix) :: M
    type(sym_rotation) :: Rident, Rrand, Rin
    real(8), allocatable :: D0(:, :), Drot(:, :), Dref(:, :), C(:, :)
    real(8), allocatable :: Ublk(:, :), tr_before(:)
    integer :: p, word, h, sc, q, s, b, k, i, knd, sumexp, sp
    integer :: herms(2), stys(2), masks(4)
    real(8), parameter :: tol = 1.0d-8
    real(8) :: offset, dev, tr_after
    character(len=160) :: label

    herms = [HERMITIAN, ANTIHERMITIAN]
    stys = [STY_EVEN, STY_ODD]
    masks = [15, 5, 7, 13]
    offset = 11.0d0

    do p = 1, NPAT
      do word = 0, 7
        do h = 1, 2
          do sc = 1, 2
          do q = 1, 4
            knd = herms(h)*stys(sc)
            write (label, '(a,i0,a,i0,a,i0,a,i0,a,i0)') 'pattern ', p,        &
            &  ' word ', word, ' herm ', herms(h), ' sty ', stys(sc),          &
            &  ' mask ', masks(q)
            cur_counts(0:7) = patterns(p, 1:8)
            cur_label = trim(label)

            call M%init(patterns(p, :), word, herms(h), stys(sc),              &
            &          goodqn=masks(q))
            cur_st = M%st
            cur_knd = knd
            call M%fill_blocks(fill_block_probed, int(offset))
            D0 = M%dense()

            ! Identity rotation through init: an exact no-op.
            call Rident%init(patterns(p, :))
            call M%rotate(Rident)
            Drot = M%dense()
            call check(all(Drot .eq. D0),                                      &
            &          trim(label)//': identity init exact no-op')

            ! Identity rotation extracted from a full identity array.
            if (allocated(C)) deallocate(C)
            allocate(C(M%n, M%n))
            C = 0.0d0
            do i = 1, M%n
              C(i, i) = 1.0d0
            enddo
            call Rident%init_from_dense(patterns(p, :), C)
            deallocate(C)
            call M%rotate(Rident)
            Drot = M%dense()
            call check(all(Drot .eq. D0),                                      &
            &          trim(label)//': identity init_from_dense exact no-op')

            ! A random orthogonal rotation, one block per populated sector.
            call Rrand%init(patterns(p, :))
            do s = 0, 7
              k = Rrand%st(s + 1) - Rrand%st(s)
              if (k .eq. 0) cycle
              Ublk = random_orthogonal(k, 100*s + 7*p)
              call Rrand%set_sector(s, Ublk)
            enddo
            C = dense_of_rotation(Rrand)
            Dref = matmul(transpose(C), matmul(D0, C))

            ! Diagonal-block traces before the rotation (word 0: every
            !  stored block is square and diagonal in sector space).
            if (word .eq. 0) then
              allocate(tr_before(M%nblocks))
              do b = 1, M%nblocks
                tr_before(b) = 0.0d0
                do i = 1, min(size(M%blocks(b)%A, 1),                           &
                &             size(M%blocks(b)%A, 2))
                  tr_before(b) = tr_before(b) + M%blocks(b)%A(i, i)
                enddo
              enddo
            endif

            call M%rotate(Rrand)

            ! Every metadata field is untouched, and the block list still
            !  satisfies the construction invariants of that metadata.
            call check(M%symmetry .eq. word,                                    &
            &          trim(label)//': symmetry word kept')
            call check(M%kind .eq. knd,                                         &
            &          trim(label)//': derived kind kept')
            call check(M%sty .eq. stys(sc),                                     &
            &          trim(label)//': characteristic kept')
            call check(M%goodqn .eq. masks(q),                                  &
            &          trim(label)//': goodqn mask kept')
            call check(M%n .eq. sum(cur_counts),                                &
            &          trim(label)//': total size kept')
            sumexp = 0
            do s = 0, 7
              sp = ieor(s, iand(word, masks(q)))
              if (cur_counts(s) .eq. 0) cycle
              if (cur_counts(sp) .eq. 0) cycle
              if (s .gt. sp) cycle
              sumexp = sumexp + 1
            enddo
            call check(M%nblocks .eq. sumexp,                                   &
            &          trim(label)//': block count kept')
            do b = 1, M%nblocks
              call check(M%blocks(b)%srow                                       &
              &          .eq. ieor(M%blocks(b)%scol, iand(word, masks(q))),     &
              &          trim(label)//': block pairing kept')
              call check(M%blocks(b)%is .eq. M%st(M%blocks(b)%srow),            &
              &          trim(label)//': row offset kept')
              call check(M%blocks(b)%js .eq. M%st(M%blocks(b)%scol),            &
              &          trim(label)//': col offset kept')
              call check(size(M%blocks(b)%A, 1)                                 &
              &          .eq. cur_counts(M%blocks(b)%srow),                     &
              &          trim(label)//': row shape kept')
              call check(size(M%blocks(b)%A, 2)                                 &
              &          .eq. cur_counts(M%blocks(b)%scol),                     &
              &          trim(label)//': col shape kept')
            enddo

            ! The rotated storage equals the dense reference.
            Drot = M%dense()
            dev = 0.0d0
            if (M%n .gt. 0) dev = maxval(abs(Drot - Dref))
            call check(dev .le. tol,                                            &
            &          trim(label)//': rotated equals dense reference')

            ! Diagonal-block traces are invariant.
            if (word .eq. 0) then
              do b = 1, M%nblocks
                tr_after = 0.0d0
                do i = 1, min(size(M%blocks(b)%A, 1),                           &
                &             size(M%blocks(b)%A, 2))
                  tr_after = tr_after + M%blocks(b)%A(i, i)
                enddo
                call check(abs(tr_after - tr_before(b))                        &
                &          .le. tol*max(1.0d0, abs(tr_before(b))),                &
                &          trim(label)//': diagonal-block trace invariant')
              enddo
              deallocate(tr_before)
            endif

            ! Round trip: the sector-wise transposed rotation undoes Rrand.
            call Rin%init(patterns(p, :))
            do s = 0, 7
              if (.not. allocated(Rrand%U(s)%U)) cycle
              call Rin%set_sector(s, transpose(Rrand%U(s)%U))
            enddo
            call M%rotate(Rin)
            Drot = M%dense()
            dev = 0.0d0
            if (M%n .gt. 0) dev = maxval(abs(Drot - D0))
            call check(dev .le. tol,                                            &
            &          trim(label)//': round trip recovers original')

            ! Partner-consistent operand, both kinds: rotating again, the
            !  dense expansion obeys D(gj, gi) = kind * D(gi, gj).
            call M%fill_blocks(fill_block_partner, int(0))
            call M%rotate(Rrand)
            Drot = M%dense()
            dev = 0.0d0
            if (M%n .gt. 0) dev = maxval(abs(Drot - knd*transpose(Drot)))
            call check(dev .le. tol,                                            &
            &          trim(label)//': partner relation survives rotation')
          enddo
          enddo
        enddo
      enddo
    enddo
  end subroutine test_rotation
  !----------------------------------------------------------------------------
  subroutine fill_block_probed(srow, scol, A, params)
    !---------------------------------------------------------------------------
    ! block_filler-shaped test callback: counts its invocations, records
    !  the sector pair it was called for, checks the shape of the block
    !  storage it receives against the current pattern (cur_counts) and
    !  overwrites the block with values from block_val, offset by the
    !  integer carried in params.
    !---------------------------------------------------------------------------
    integer, intent(in)   :: srow, scol
    real(8), intent(out)  :: A(:, :)
    class(*), intent(in)  :: params

    real(8) :: offset
    integer :: r, c, mrow, ncol
    character(len=160) :: label

    offset = 0.0d0
    select type (params)
    type is (integer)
      offset = real(params, 8)
    end select

    nblockcalls = nblockcalls + 1
    seen(srow, scol) = seen(srow, scol) + 1

    label = trim(cur_label)//' callback'
    mrow = size(A, 1)
    ncol = size(A, 2)
    call check(mrow .eq. cur_counts(srow), trim(label)//': row size')
    call check(ncol .eq. cur_counts(scol), trim(label)//': col size')

    do r = 1, mrow
      do c = 1, ncol
        A(r, c) = block_val(srow, scol, r, c, offset)
      enddo
    enddo
  end subroutine fill_block_probed
  !----------------------------------------------------------------------------
  subroutine fill_block_partner(srow, scol, A, params)
    !---------------------------------------------------------------------------
    ! block_filler-shaped callback for the partner-consistency judgment of
    !  the rotation test: the block is filled with a symmetric (kind +1)
    !  or antisymmetric (kind -1, per cur_knd) function of the GLOBAL
    !  index pair, through the cur_st bounds of the current case. Such a
    !  matrix satisfies  A(s', s) = kind * A(s, s')^T  globally, over
    !  stored blocks and unstored partners alike.
    !---------------------------------------------------------------------------
    integer, intent(in)   :: srow, scol
    real(8), intent(out)  :: A(:, :)
    class(*), intent(in)  :: params

    integer :: r, c, gi, gj
    real(8) :: val

    do c = 1, size(A, 2)
      gj = cur_st(scol) + c - 1
      do r = 1, size(A, 1)
        gi = cur_st(srow) + r - 1
        if (cur_knd .gt. 0) then
          val = 1.5d0*gi*gj + 0.5d0*(gi + gj)
        else
          val = 1.5d0*(gi - gj)
        endif
        A(r, c) = val
      enddo
    enddo
  end subroutine fill_block_partner
  !----------------------------------------------------------------------------
  function random_orthogonal(k, seed) result(U)
    !---------------------------------------------------------------------------
    ! Deterministic pseudo-random orthogonal block: hashed columns through
    !  modified Gram-Schmidt, re-orthogonalised once for numerical
    !  quality, with a unit-vector fallback against degenerate hashes.
    !  No library dependency - this suite links nothing but the module.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Input:
    !   k    : the block size
    !   seed : deterministic variation between callers
    ! Result:
    !   U    : a (k, k) orthogonal block
    !---------------------------------------------------------------------------
    integer, intent(in) :: k
    integer, intent(in) :: seed
    real(8) :: U(k, k)

    real(8) :: v(k), nrm, prj
    integer :: c, r, rr, pass, j

    if (k .eq. 0) return

    do c = 1, k
      do r = 1, k
        U(r, c) = cos(0.37d0*(seed + 7*c + 61*r))
      enddo
    enddo

    do pass = 1, 2
      do c = 1, k
        v = U(:, c)
        do rr = 1, c - 1
          prj = dot_product(U(:, rr), v)
          v = v - prj*U(:, rr)
        enddo
        nrm = sqrt(dot_product(v, v))
        if (nrm .lt. 1.0d-8) then
          ! Degenerate hash column: fall back to a unit vector outside the
          !  span built so far - one always exists while c <= k.
          do j = 1, k
            v = 0.0d0
            v(j) = 1.0d0
            do rr = 1, c - 1
              prj = dot_product(U(:, rr), v)
              v = v - prj*U(:, rr)
            enddo
            nrm = sqrt(dot_product(v, v))
            if (nrm .gt. 0.1d0) exit
          enddo
        endif
        U(:, c) = v/nrm
      enddo
    enddo
  end function random_orthogonal
  !----------------------------------------------------------------------------
  function dense_of_rotation(R) result(C)
    !---------------------------------------------------------------------------
    ! The full n x n array of a sym_rotation: the sector-diagonal blocks
    !  at their global positions, zero elsewhere - the array a consumer
    !  would own in the basis_transform convention (columns are the new
    !  basis vectors expressed in the old basis).
    !---------------------------------------------------------------------------
    type(sym_rotation), intent(in) :: R
    real(8), allocatable :: C(:, :)

    integer :: s

    allocate(C(R%n, R%n))
    C = 0.0d0
    do s = 0, 7
      if (.not. allocated(R%U(s)%U)) cycle
      C(R%st(s):R%st(s + 1) - 1, R%st(s):R%st(s + 1) - 1) = R%U(s)%U
    enddo
  end function dense_of_rotation
  !----------------------------------------------------------------------------
  pure function block_val(srow, scol, i, j, offset) result(val)
    !---------------------------------------------------------------------------
    ! Deterministic value identifying a block element: distinct for each
    !  combination of sector pair, local element and offset, so any
    !  misrouting of a fill shows up as a content mismatch. Kept an exact
    !  multiple of 1/2 so equality comparisons are exact.
    !---------------------------------------------------------------------------
    integer, intent(in) :: srow, scol, i, j
    real(8), intent(in) :: offset
    real(8)             :: val

    val = offset + 8.0d0*srow + 64.0d0*scol + i + 512.0d0*j
  end function block_val
  !----------------------------------------------------------------------------

end program block_matrix_unit_tests
