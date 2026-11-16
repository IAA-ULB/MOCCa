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
  ! Unit tests for phase 1 of the symmetry-block matrix module: construction
  !  of the block list from MOCCa's HFBlocks bookkeeping.
  !
  ! For a suite of HFBlocks patterns (EV8-like, time-reversal-broken, sparse,
  !  degenerate), every operator symmetry word 0..7 and every matrix kind,
  !  the stored block list is verified against the design invariants:
  !
  !   *) selection rules: each stored block couples sectors (srow, scol)
  !      with srow = IEOR(scol, symmetry) and nothing else is stored
  !   *) kind-based partner skipping: for (anti-)symmetric matrices at most
  !      one block per sector pair is stored, the skipped one satisfies
  !      scol > srow
  !   *) one block per populated sector pair for general matrices
  !   *) block shapes, offsets and zero-initialisation
  !   *) sector bounds array st(0:8) and the global size n
  !   *) sector_of mapping, including empty sectors and out-of-range indices
  !
  ! The tests are self-contained: no other MOCCa module is used. Compile
  !  together with src_orig/block_matrix_module.f90 and run.
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! The suite consists of three routines (plus the check helper that
  !  registers every individual condition):
  !   test_construction : construction invariants over the full cross
  !                      product of patterns, symmetry words and kinds
  !   test_reinit      : repeated init on the same object
  !   test_sector_of   : global-index-to-sector mapping
  !==============================================================================
  use block_matrix_module

  implicit none

  integer, parameter :: NPAT = 6

  integer :: patterns(NPAT, 8)
  integer :: nfail, ntest

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! HFBlocks patterns: MOCCa ordering, per the table in the module header.
  !   1: EV8-like bookkeeping (only blocks 1, 3, 5, 7 populated)
  !   2: time-reversal broken (all blocks populated)
  !   3: sparse pattern with interior empty sectors
  !   4: empty leading sectors
  !   5: single populated sector
  !   6: completely empty basis
  patterns(1, :) = [12, 0, 9, 0, 7, 0, 5, 0]
  patterns(2, :) = [12, 11, 9, 8, 7, 6, 5, 4]
  patterns(3, :) = [3, 0, 0, 2, 5, 0, 1, 0]
  patterns(4, :) = [0, 4, 2, 0, 0, 6, 0, 3]
  patterns(5, :) = [0, 0, 0, 0, 0, 0, 9, 0]
  patterns(6, :) = [0, 0, 0, 0, 0, 0, 0, 0]

  nfail = 0
  ntest = 0

  call test_construction()
  call test_reinit()
  call test_sector_of()

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
    !  0..7 and all three matrix kinds.
    !
    ! For each combination a fresh matrix is built and checked for:
    !   *) stored metadata: symmetry word, kind, total size, nblocks <= 8
    !   *) the bounds array st: st(0) = 1, st(8) = n+1, and the sector
    !      sizes reproduce the differences st(s+1) - st(s)
    !   *) the number of stored blocks: one per populated ordered sector
    !      pair, with the s > srow partner skipped unless the kind is
    !      KIND_GENERAL
    !   *) per stored block: srow = IEOR(scol, word), uniqueness of the
    !      (srow, scol) pair, both sectors populated, A allocated with the
    !      sector-pair shape, zero-filled, correct offsets is/js and a
    !      cleared transposed flag
    !   *) completeness: every ordered sector pair that should be stored
    !      is present in the list
    !   *) hygiene: block slots beyond nblocks carry no allocation
    !---------------------------------------------------------------------------
    type(sym_matrix) :: M
    integer :: p, word, k, kinds(3), b, s, sp, counts(0:7)
    integer :: expected, sumexp
    logical :: found, occupied(0:7, 0:7)
    character(len=120) :: label

    kinds = [KIND_GENERAL, KIND_SYMMETRIC, KIND_ANTISYMMETRIC]

    do p = 1, NPAT
      counts(0:7) = patterns(p, 1:8)
      do word = 0, 7
        do k = 1, 3
          write (label, '(a,i0,a,i0,a,i0)') 'pattern ', p, ' word ', word, &
          &  ' kind ', kinds(k)
          call M%init(patterns(p, :), word, kinds(k))

          call check(M%symmetry .eq. word, trim(label)//': symmetry word')
          call check(M%kind .eq. kinds(k), trim(label)//': kind')
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

          ! Expected number of stored blocks: one per populated ordered pair,
          !  with kind-based skipping of the s > srow partner for
          !  (anti-)symmetric matrices.
          sumexp = 0
          do s = 0, 7
            sp = ieor(s, word)
            if (counts(s) .eq. 0) cycle
            if (counts(sp) .eq. 0) cycle
            if (kinds(k) .ne. KIND_GENERAL .and. s .gt. sp) cycle
            sumexp = sumexp + 1
          enddo
          call check(M%nblocks .eq. sumexp, trim(label)//': expected nblocks')

          ! Every stored block obeys the XOR selection rule, is unique,
          !  points at the right offset, is zero-filled and carries the right
          !  shape. Partner blocks skipped for non-general kinds stayed away.
          occupied = .false.
          do b = 1, M%nblocks
            associate (Blk => M%blocks(b))
              write (label, '(a,i0,a,i0,a,i0,a,i0)')                            &
              &  'pattern ', p, ' word ', word, ' kind ', kinds(k),             &
              &  ' block ', b

              call check(Blk%srow .eq. ieor(Blk%scol, word),                    &
              &          trim(label)//': srow = IEOR(scol, word)')
              call check(Blk%scol .ge. 0 .and. Blk%scol .le. 7,                 &
              &          trim(label)//': scol in 0..7')
              call check(.not. occupied(Blk%srow, Blk%scol),                    &
              &          trim(label)//': unique sector pair')
              occupied(Blk%srow, Blk%scol) = .true.

              ! Partner-skip: a stored block of an (anti-)symmetric matrix
              !  never sits on the skipped side of its pair.
              if (kinds(k) .ne. KIND_GENERAL) then
                call check(Blk%scol .le. Blk%srow,                             &
                &          trim(label)//': partner skipped, scol <= srow')
              endif

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
              call check(.not. Blk%transposed, trim(label)//': not transposed')
            end associate
          enddo

          ! Every populated ordered pair that should be stored is stored.
          do s = 0, 7
            sp = ieor(s, word)
            expected = 0
            if (counts(s) > 0 .and. counts(sp) > 0) then
              if (kinds(k) .eq. KIND_GENERAL .or. s .le. sp) expected = 1
            endif
            found = .false.
            do b = 1, M%nblocks
              if (M%blocks(b)%srow .eq. sp .and. M%blocks(b)%scol .eq. s)       &
              &  found = .true.
            enddo
            write (label, '(a,i0,a,i0,a,i0,a,i0)') 'pattern ', p,               &
            &  ' word ', word, ' kind ', kinds(k), ' pair (row,col) ', sp
            call check(merge(1, 0, found) .eq. expected,                         &
            &          trim(label)//': expected pair present')
          enddo

          ! Unused block slots must be cleaned out.
          do b = M%nblocks + 1, 8
            write (label, '(a,i0,a,i0,a,i0,a,i0)') 'pattern ', p,               &
            &  ' word ', word, ' kind ', kinds(k), ' unused slot ', b
            call check(.not. allocated(M%blocks(b)%A),                          &
            &          trim(label)//': A not allocated')
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
    !  survive. Finally rebuild with a word-5 antisymmetric matrix to
    !  verify the list and the symmetry word are restored.
    !---------------------------------------------------------------------------
    type(sym_matrix) :: M
    integer   :: b, nexp
    character(len=120) :: label

    call M%init(patterns(2, :), 0, KIND_SYMMETRIC)
    nexp = count(patterns(2, :) > 0)
    call check(M%nblocks .eq. nexp, 'reinit: symmetric word-0 block count')

    call M%init(patterns(1, :), 7, KIND_GENERAL)
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
    call M%init(patterns(2, :), 5, KIND_ANTISYMMETRIC)
    call check(M%nblocks .gt. 0, 'reinit: blocks rebuilt')
    call check(M%symmetry .eq. 5, 'reinit: word updated')
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
    call M%init(patterns(3, :), 0, KIND_GENERAL)

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

end program block_matrix_unit_tests
