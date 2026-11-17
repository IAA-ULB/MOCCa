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
module block_matrix_module
!==============================================================================
  !
  ! This module contains the machinery to efficiently represent matrices
  !  in the configuration space spanned by a finite set of single-particle
  !  wavefunctions on the MOCCa grid.
  !
  ! Implementation of the design summary
  !  'mocca-symmetry-block-matrix-design-summary.md'. The current version
  !  covers phases 1 and 2 of the design: sector encoding, the sym_block
  !  and sym_matrix types, construction of the block list from MOCCa's
  !  HFBlocks bookkeeping, and population of the stored blocks through
  !  dummy-procedure callbacks (the block-wise fill_blocks and the
  !  per-element fallback fill_elements), plus the new_sym_matrix
  !  constructor sugar and the good-quantum-number mask in init (design
  !  section 2b). Later phases add mat-vec/mat-mat products and
  !  the complex-operator wrapper.
  !
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Symmetry model
  ! ---------------
  ! Operators acting on the HF basis carry three discrete quantum numbers,
  !  each of which acts as a Z2 grading on the basis. A symmetry-sector word
  !  is a 3-bit integer s = 0..7 with the encoding
  !
  !   bit 0 : z-signature (Rz) : 0 = q(Rz) = +i , 1 = q(Rz) = -i
  !           (Rz is linear and antihermitian, hence eigenvalues +/- i)
  !   bit 1 : parity (P)       : 0 = even      , 1 = odd
  !   bit 2 : isospin (Tz)     : 0 = neutron   , 1 = proton
  !
  ! An operator carries a symmetry word 'op' with the same encoding: a bit
  !  is set if the operator changes the corresponding quantum number of the
  !  state it acts on. Selection rules are then a single XOR: an operator
  !  with word 'op' couples sector s to sector IEOR(s, op) and to nothing
  !  else. An operator matrix therefore has at most 8 non-vanishing dense
  !  blocks, one per surviving sector pair.
  !
  ! MOCCa's HFBlocks bookkeeping in the wavefunctions module uses exactly
  !  the same ordering: HFBlocks(B) counts the spwfs of sector s = B - 1,
  !
  !    B = 1: neutron, even, Rz = +i   =>  s = 0     B = 5: proton, even, +i  => s = 4
  !    B = 2: neutron, even, Rz = -i   =>  s = 1     B = 6: proton, even, -i  => s = 5
  !    B = 3: neutron, odd , Rz = +i   =>  s = 2     B = 7: proton, odd , +i  => s = 6
  !    B = 4: neutron, odd , Rz = -i   =>  s = 3     B = 8: proton, odd , -i  => s = 7
  !
  !  Typical occupations, all consistent with this single table:
  !    EV8 symmetries ('Rz,T,P,STy'):   all spwfs in blocks 1,3,5,7 (q(Rz)=+i)
  !    time-reversal broken           :  blocks 2,4,6,8 additionally populated
  !    parity broken                  :  spwfs in blocks 1,2,5,6 only, as the
  !                                      parity grading is absent (the P bit
  !                                      of any operator is then vacuous)
  !
  !  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !  Good quantum numbers of the stored basis
  !  -----------------------------------------
  !  The symmetry word of an operator is physical: it tells which quantum
  !  number of the state it acts on is flipped. Not every build labels its
  !  basis states with all three quantum numbers: when a self-consistent
  !  symmetry is broken, the states of both characters are merged into a
  !  single storage block and the corresponding quantum number is not
  !  'good' - its bit is vacuous, and an operator acts diagonally in it.
  !  The classic case is parity: in parity-broken builds an odd-l operator
  !  couples states WITHIN blocks 1, 2, 5, 6 rather than pairing them into
  !  the empty blocks 3, 4, 7, 8 as its physical word suggests.
  !
  !  init therefore takes an optional 3-bit mask 'goodqn' (same bit
  !  assignment as the sector words; default 7 = every quantum number
  !  good) marking the quantum numbers that actually label the basis, and
  !  pairs every sector s with IEOR(s, AND(symmetry, goodqn)). The mask is
  !  a fact of the build scenario, not recoverable from the runtime block
  !  counts: in MOCCa it is the geninfo constant 'goodqn', injected by
  !  Hephaestos ($GOODQN). With all bits good the pairing is the pure XOR
  !  above. Standalone callers (unit tests) omit the mask.
  !
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Conventions
  ! ------------
  !  *) All real work is done in sector space, i.e. all sector-indexed
  !     arrays are declared with lower bound 0 (st(0:8)) such that IEOR
  !     applies directly. No +/- 1 adjustments appear in any logic.
  !  *) MOCCa itself is 1-based. The single conversion from MOCCa's 1-based
  !     HFBlocks ordering to the 0-based sector words happens in init.
  !  *) Global-basis slices follow design option 2 (1-based bounds):
  !     st(s) is the first global index of sector s and st(8) = n + 1,
  !     i.e. sector s occupies the global indices st(s) ... st(s+1) - 1.
  !  *) Each block is a separate contiguous rank-2 array with natural
  !     leading dimension; no flattened storage, no packing.
  !
  ! Storage is single-node by design: rank-dependent logic does not belong
  !  in this module; gathering of whole-basis vectors is the caller's task.
  !
  ! Population contract (design section 5): blocks are populated
  !  exclusively through dummy-procedure callbacks with the abstract
  !  interfaces block_filler (whole blocks) and element_filler (per-element
  !  fallback); procedure pointers are never stored as components:
  !   *) callbacks run for stored blocks only. Blocks skipped by the
  !      kind-based partner rule are not the filler's business: partner
  !      values are resolved at access time (a later phase), never by
  !      writing here.
  !   *) callbacks receive the sector words srow/scol of the block so
  !      they can apply MOCCa's phase conventions for matrix elements.
  !   *) both paths OVERWRITE: fill_blocks hands the block storage out as
  !      an intent(out) dummy, fill_elements assigns element by element.
  !   *) 'params' is unlimited-polymorphic state handed verbatim to the
  !      callback. The signature is frozen: extend params, not the
  !      interface. Callers without state pass any token.
  !   *) fill_elements costs one callback invocation per element and
  !      exists for tests and prototypes only; hot paths fill whole
  !      blocks via fill_blocks.
  !==============================================================================
  implicit none

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Hermiticity class of the matrix. For KIND_SYMMETRIC / KIND_ANTISYMMETRIC
  !  matrices, partner blocks are related by
  !     A(s', s) = sigma * A(s, s')^T      with sigma = +/- 1
  !  and only one block per sector pair is stored; the partner is accessed
  !  with the transposed flag and the sign. A block whose sector is its own
  !  partner is square and internally symmetric/antisymmetric; it is stored
  !  dense and the symmetry serves only as a consistency check.
  integer, parameter :: KIND_GENERAL        =  0
  integer, parameter :: KIND_SYMMETRIC      =  1
  integer, parameter :: KIND_ANTISYMMETRIC  = -1

  !------------------------------------------------------------------------------
  ! One dense block: rows span global indices is : is+m-1 of sector srow,
  !  columns span js : js+n-1 of sector scol, with m = size(A,1), n = size(A,2).
  type :: sym_block
    integer :: srow, scol      ! sector words 0..7 of the row/column sector
    integer :: is, js          ! 1-based global start of the row/column sector
    logical :: transposed      ! .true. if this block is reached as the
                              !  +/- transpose of a stored partner
    real(8), allocatable :: A(:, :)  ! (m, n), contiguous, LDA = m
  end type sym_block

  !------------------------------------------------------------------------------
  ! The full symmetry-adapted matrix: at most 8 blocks in a flat list.
  !  After init, all runtime code loops over blocks(1:nblocks) only;
  !  quantum numbers XOR selection rules never appear at runtime again.
  type :: sym_matrix
    integer :: symmetry       ! operator symmetry word, 0..7 (physical)
    integer :: goodqn         ! good-quantum-number mask of the basis, 0..7
                              ! (design 2b); 7 = all quantum numbers good
    integer :: kind           ! KIND_GENERAL / KIND_SYMMETRIC / KIND_ANTISYMMETRIC
    integer :: n              ! total basis size, sum(HFBlocks)
    integer :: nblocks        ! number of stored blocks, <= 8
    type(sym_block) :: blocks(8)
    integer :: st(0:8)        ! first global (1-based) index of each sector;
                              !  st(8) = n + 1; sector s spans st(s) .. st(s+1)-1
  contains
    procedure :: init => init_sym_matrix
    procedure :: sector_of => sector_of_sym_matrix
    procedure :: fill_blocks => fill_blocks_sym_matrix
    procedure :: fill_elements => fill_elements_sym_matrix
    procedure :: print_dense => print_dense_sym_matrix
  end type sym_matrix

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Population callbacks (design section 5): dummy procedures with an
  !  abstract interface, never stored as components.
  abstract interface
  !----------------------------------------------------------------------------
  subroutine block_filler(srow, scol, A, params)
    !---------------------------------------------------------------------------
    ! Overwrite one whole block. Invoked by fill_blocks once per stored
    !  block, with the sector words srow/scol of the block and a rank-2
    !  view of the block storage itself.
    !
    ! A is an intent(out) dummy: the callback must overwrite the entire
    !  block. 'params' is caller state handed through verbatim.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Input:
    !   srow, scol : sector words 0..7 of the block
    !   params     : caller state, opaque
    ! Output:
    !   A          : the block storage, to be overwritten entirely
    !---------------------------------------------------------------------------
    integer, intent(in)   :: srow, scol
    real(8), intent(out)  :: A(:, :)
    class(*), intent(in)  :: params
  end subroutine block_filler
  !----------------------------------------------------------------------------
  function element_filler(srow, scol, i, j, params) result(val)
    !---------------------------------------------------------------------------
    ! Return one matrix element of the block (srow, scol) at the LOCAL
    !  row/column indices (i, j) of that block, 1-based within the block
    !  (the global indices are st(srow) + i - 1 and st(scol) + j - 1).
    ! Invoked by fill_elements for every element of every stored block.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Input:
    !   srow, scol : sector words 0..7 of the block
    !   params     : caller state, opaque
    ! Result:
    !   val        : the matrix element A(i, j) of the block
    !---------------------------------------------------------------------------
    integer, intent(in)   :: srow, scol, i, j
    class(*), intent(in)  :: params
    real(8)              :: val
  end function element_filler
  end interface

  private
  public :: sym_block, sym_matrix
  public :: KIND_GENERAL, KIND_SYMMETRIC, KIND_ANTISYMMETRIC
  public :: block_filler, element_filler
  public :: new_sym_matrix

contains

  !----------------------------------------------------------------------------
  subroutine init_sym_matrix(this, HFBlocks, symmetry, kind, goodqn)
    !---------------------------------------------------------------------------
    ! Construct the block list of a symmetry-adapted matrix.
    !
    ! For every populated sector s, the operator couples it to sector
    !  sp = IEOR(s, AND(symmetry, goodqn)) and to nothing else - the
    !  physical word is intersected with the good-quantum-number mask of
    !  design section 2b, so bits of quantum numbers that do not label the
    !  basis (broken gradings, e.g. parity in parity-broken builds) act
    !  diagonally. With the default mask 7 this is the pure XOR
    !  sp = IEOR(s, symmetry). A dense block is stored
    !  for each populated sector pair; for KIND_SYMMETRIC / KIND_ANTISYMMETRIC
    !  only one block per pair is stored (partner access through the
    !  transposed flag covers the other), and a sector that is its own
    !  partner gives a single square block.
    !
    ! All blocks are allocated and zero-filled: the matrix is ready to be
    !  populated by a fill callback (later phase) or element-by-element.
    ! init may be called repeatedly on the same object; previous state is
    !  discarded cleanly.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !   HFBlocks(1:8) : number of spwfs per MOCCa symmetry block. This is
    !                   the GLOBAL count array (HFBlocks_global in the
    !                   wavefunctions module); the rank-local twin must not
    !                   be passed, as this type is single-node by design.
    !   symmetry      : physical operator symmetry word, 0..7
    !   kind          : KIND_GENERAL, KIND_SYMMETRIC or KIND_ANTISYMMETRIC
    !   goodqn        : optional 3-bit mask of the quantum numbers that are
    !                   good for this build's basis (design section 2b);
    !                   default 7 = every quantum number good. In MOCCa,
    !                   pass the geninfo constant 'goodqn'.
    !---------------------------------------------------------------------------
    class(sym_matrix), intent(inout) :: this
    integer, intent(in)               :: HFBlocks(8)
    integer, intent(in)               :: symmetry
    integer, intent(in)               :: kind
    integer, intent(in), optional      :: goodqn

    integer :: s, sp, b, qn, word
    integer :: counts(0:7)

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Guards: sector words are 3-bit integers, kinds are the module constants.
    if (symmetry < 0 .or. symmetry > 7) then
      print *, 'init_sym_matrix: symmetry word must be in 0..7, got', symmetry
      stop 1
    endif
    if (kind .ne. KIND_GENERAL .and. kind .ne. KIND_SYMMETRIC            &
    & .and. kind .ne. KIND_ANTISYMMETRIC) then
      print *, 'init_sym_matrix: unknown kind', kind
      stop 1
    endif
    if (any(HFBlocks < 0)) then
      print *, 'init_sym_matrix: HFBlocks entries must be >= 0'
      stop 1
    endif

    ! Good-quantum-number mask (design section 2b): default 7 = every
    !  quantum number good. It is applied exactly once, in the pairing
    !  word below; quantum numbers that do not label the basis act
    !  diagonally.
    qn = 7
    if (present(goodqn)) qn = goodqn
    if (qn < 0 .or. qn > 7) then
      print *, 'init_sym_matrix: goodqn mask must be in 0..7, got', qn
      stop 1
    endif

    ! The only place MOCCa's 1-based block ordering is converted to 0-based
    !  sector words: sector s = MOCCa block s + 1 (see table in module header).
    counts(0:7) = HFBlocks(1:8)

    this%symmetry = symmetry
    this%goodqn = qn
    this%kind = kind
    this%n = sum(counts)

    ! The pairing word: the physical operator word intersected with the
    !  mask - the one place where the good quantum numbers enter the maths.
    word = iand(symmetry, qn)

    ! 1-based bounds: st(s) is the first global index of sector s.
    this%st(0) = 1
    do s = 0, 7
      this%st(s + 1) = this%st(s) + counts(s)
    enddo

    ! Discard any pre-existing state of the object.
    do b = 1, 8
      this%blocks(b)%srow = 0
      this%blocks(b)%scol = 0
      this%blocks(b)%is = 0
      this%blocks(b)%js = 0
      this%blocks(b)%transposed = .false.
      if (allocated(this%blocks(b)%A)) deallocate(this%blocks(b)%A)
    enddo

    ! Append one block per surviving sector pair. The pairing is a single
    !  XOR with the (masked) word; the kind-based skipping keeps at most
    !  one block per pair for (anti-)symmetric matrices.
    b = 0
    do s = 0, 7
      sp = ieor(s, word)
      if (counts(s) .eq. 0) cycle
      if (counts(sp) .eq. 0) cycle
      if (kind .ne. KIND_GENERAL .and. s .gt. sp) cycle

      b = b + 1
      this%blocks(b)%srow = sp
      this%blocks(b)%scol = s
      this%blocks(b)%is = this%st(sp)
      this%blocks(b)%js = this%st(s)
      this%blocks(b)%transposed = .false.
      allocate(this%blocks(b)%A(counts(sp), counts(s)))
      this%blocks(b)%A = 0.0d0
    enddo
    this%nblocks = b
  end subroutine init_sym_matrix
  !----------------------------------------------------------------------------
  pure integer function sector_of_sym_matrix(this, i)
    !---------------------------------------------------------------------------
    ! Map a global spwf index to its symmetry sector. Returns -1 if the index
    !  is out of range. Exploits that the bounds array is non-decreasing.
    !---------------------------------------------------------------------------
    class(sym_matrix), intent(in) :: this
    integer, intent(in)            :: i

    integer :: s

    sector_of_sym_matrix = -1
    if (i < 1 .or. i > this%n) return
    do s = 7, 0, -1
      if (i .ge. this%st(s)) then
        sector_of_sym_matrix = s
        exit
      endif
    enddo
  end function sector_of_sym_matrix
  !----------------------------------------------------------------------------
  subroutine fill_blocks_sym_matrix(this, fill, params)
    !---------------------------------------------------------------------------
    ! Populate all stored blocks through a block-wise callback.
    !
    ! For every stored block the callback 'fill' is invoked once, receiving
    !  the block's sector words and a direct rank-2 view of the block
    !  storage as an intent(out) dummy: the callback OVERWRITES the entire
    !  block, typically via whole-array assignments or dense kernels.
    !
    ! The callback runs for stored blocks only; partners skipped by the
    !  kind-based rule are untouched, as their values are resolved at
    !  access time (a later phase) and never stored here. Distinct blocks
    !  occupy disjoint storage, so the loop is safe to parallelise over
    !  blocks; the OpenMP treatment is deferred to the design's step 7.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Input/output:
    !   this   : the matrix whose stored blocks are populated
    !   fill   : callback overwriting one block (interface block_filler)
    !   params : state handed to the callback, verbatim; callers without
    !            state pass any token
    !---------------------------------------------------------------------------
    class(sym_matrix), intent(inout) :: this
    procedure(block_filler)           :: fill
    class(*), intent(in)              :: params

    integer :: b

    do b = 1, this%nblocks
      call fill(this%blocks(b)%srow, this%blocks(b)%scol,                        &
      &         this%blocks(b)%A, params)
    enddo
  end subroutine fill_blocks_sym_matrix
  !----------------------------------------------------------------------------
  subroutine fill_elements_sym_matrix(this, fill, params)
    !---------------------------------------------------------------------------
    ! Populate all stored blocks through a per-element callback.
    !
    ! Element-wise fallback of the design (section 5): for every stored
    !  block and every local element (i, j) of it, the callback is queried
    !  for the value A_ij, which is assigned into the block storage,
    !  overwriting whatever init or an earlier fill left there.
    !
    ! This costs one callback invocation per element and exists for tests
    !  and prototypes only; hot paths fill whole blocks via fill_blocks.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Input/output:
    !   this   : the matrix whose stored blocks are populated
    !   fill   : callback returning one element (interface element_filler)
    !   params : state handed to the callback, verbatim
    !---------------------------------------------------------------------------
    class(sym_matrix), intent(inout) :: this
    procedure(element_filler)         :: fill
    class(*), intent(in)              :: params

    integer :: b, i, j, m, ncol

    do b = 1, this%nblocks
      m = size(this%blocks(b)%A, 1)
      ncol = size(this%blocks(b)%A, 2)
      do i = 1, m
        do j = 1, ncol
          this%blocks(b)%A(i, j) = fill(this%blocks(b)%srow,                     &
          &                          this%blocks(b)%scol, i, j, params)
        enddo
      enddo
    enddo
  end subroutine fill_elements_sym_matrix
  !----------------------------------------------------------------------------
  function new_sym_matrix(HFBlocks, symmetry, kind, fill, params, goodqn)     &
  & result(M)
    !---------------------------------------------------------------------------
    ! Constructor sugar (design section 5): build and fill in one statement.
    !
    ! Equivalent to init followed by fill_blocks with the same arguments.
    !  A named function rather than a generic 'sym_matrix': in Fortran a
    !  name cannot serve as both a derived-type name and a function name.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Input:
    !   HFBlocks(1:8) : GLOBAL MOCCa count array, see init
    !   symmetry      : physical operator symmetry word, 0..7
    !   kind          : KIND_GENERAL, KIND_SYMMETRIC or KIND_ANTISYMMETRIC
    !   fill          : callback overwriting each stored block
    !   params        : state handed to the callback, verbatim
    !   goodqn        : optional good-quantum-number mask (design 2b),
    !                   handed to init
    ! Output:
    !   M             : the freshly built and filled matrix
    !---------------------------------------------------------------------------
    integer, intent(in)     :: HFBlocks(8)
    integer, intent(in)     :: symmetry
    integer, intent(in)     :: kind
    procedure(block_filler) :: fill
    class(*), intent(in)    :: params
    integer, intent(in), optional :: goodqn
    type(sym_matrix)        :: M

    call M%init(HFBlocks, symmetry, kind, goodqn)
    call M%fill_blocks(fill, params)
  end function new_sym_matrix
  !----------------------------------------------------------------------------
  subroutine print_dense_sym_matrix(this, lun, title)
    !---------------------------------------------------------------------------
    ! Human-friendly dump of a sym_matrix as one dense n x n matrix.
    !
    ! Debug and test helper for SMALL matrices ONLY: the full dense matrix
    !  is materialised, so memory and output scale as n**2 and this routine
    !  does not transfer to realistic basis sizes.
    !
    ! (Anti-)symmetric matrices are printed complete: partner blocks that
    !  init does not store are reconstructed from their stored twin as
    !      A(s', s) = kind * A(s, s')^T ,
    !  the same expansion the dense reference of the design's validation
    !  tests (summary section 9) is built from. The transposed flag is
    !  ignored here: init never sets it and only later phases assign it.
    !
    ! Layout: a header decoding the symmetry word, the kind, the stored
    !  block list and the populated sectors, then the matrix itself, with
    !  '|' and '-' lines drawn at sector boundaries and '.' for exact
    !  zeros so the block structure can be read at a glance. Entries that
    !  do not fit the fixed F10.4 field print as asterisks.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !   this  : the matrix to dump, as built by init
    !   lun   : optional output unit, defaults to stdout (unit 6)
    !   title : optional line printed verbatim above the dump
    !---------------------------------------------------------------------------
    class(sym_matrix), intent(in) :: this
    integer, intent(in), optional :: lun
    character(*), intent(in), optional :: title

    character(2) :: sig
    character(4) :: par
    character(7) :: iso
    character(6) :: rowlab
    character(10) :: cell
    character(len=:), allocatable :: line
    logical, allocatable :: boundary(:)
    real(8), allocatable :: dense(:, :)
    integer :: u, i, j, b, m, ncol, s

    if (present(lun)) then
      u = lun
    else
      u = 6
    endif

    if (present(title)) write(u, '(A)') trim(title)
    if (this%n <= 0) then
      write(u, '(A)') ' sym_matrix is empty or not initialised, nothing to dump'
      return
    endif

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Header: symmetry word decoded into the bits it flips, the kind, the
    !  stored block list and the populated sectors with index ranges.
    write(u, '(A,I0)') ' sym_matrix: n = ', this%n

    if (this%symmetry == 0) then
      line = 'sector-diagonal operator, flips nothing'
    else
      line = 'flips:'
      if (btest(this%symmetry, 0)) line = line//' Rz'
      if (btest(this%symmetry, 1)) line = line//' P'
      if (btest(this%symmetry, 2)) line = line//' Tz'
    endif
    write(u, '(A,I0,A,A)') '  symmetry word ', this%symmetry, ': ', trim(line)
    if (this%goodqn .ne. 7) then
      write(u, '(A,I0,A)') '  good quantum numbers mask ', this%goodqn,     &
      &                    ' (design 2b: bits set = quantum numbers good)'
    endif

    select case (this%kind)
    case (KIND_GENERAL)
      line = 'general'
    case (KIND_SYMMETRIC)
      line = 'symmetric'
    case (KIND_ANTISYMMETRIC)
      line = 'antisymmetric'
    case default
      line = 'UNKNOWN'
    end select
    write(u, '(A,A,A,I0)') '  kind: ', trim(line), ', stored blocks: ', this%nblocks

    do b = 1, this%nblocks
      m = size(this%blocks(b)%A, 1)
      ncol = size(this%blocks(b)%A, 2)
      write(u, '(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0)')              &
      & '  block ', b, ': sector ', this%blocks(b)%srow, ' rows ',                    &
      & this%blocks(b)%is, ':', this%blocks(b)%is + m - 1, ' x sector ',              &
      & this%blocks(b)%scol, ' cols ', this%blocks(b)%js, ':',                         &
      & this%blocks(b)%js + ncol - 1
    enddo

    do s = 0, 7
      if (this%st(s + 1) <= this%st(s)) cycle
      if (btest(s, 0)) then
        sig = '-i'
      else
        sig = '+i'
      endif
      if (btest(s, 1)) then
        par = 'odd'
      else
        par = 'even'
      endif
      if (btest(s, 2)) then
        iso = 'proton'
      else
        iso = 'neutron'
      endif
      write(u, '(A,I0,A,I0,A,I0,A,A,A,A,A,A,A,A)')                                     &
      & '  sector ', s, ': indices ', this%st(s), ':', this%st(s + 1) - 1,             &
      & ' [Rz ', trim(sig), ', P ', trim(par), ', ', trim(iso), ']'
    enddo

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! The dense picture: expand all blocks into one n x n matrix,
    !  reconstructing unstored partners of (anti-)symmetric matrices.
    !  Sector boundaries: column/row i starts a sector boundary if it is
    !  the first index of sectors 1..7 (empty sectors leave no gap, so two
    !  successive boundaries can coincide; boundary 1 is never marked).
    allocate(boundary(this%n))
    boundary = .false.
    do s = 1, 7
      if (this%st(s) >= 2 .and. this%st(s) <= this%n) boundary(this%st(s)) = .true.
    enddo

    allocate(dense(this%n, this%n))
    dense = 0.0d0
    do b = 1, this%nblocks
      if (.not. allocated(this%blocks(b)%A)) cycle
      m = size(this%blocks(b)%A, 1)
      ncol = size(this%blocks(b)%A, 2)
      dense(this%blocks(b)%is:this%blocks(b)%is + m - 1,                             &
      &     this%blocks(b)%js:this%blocks(b)%js + ncol - 1) = this%blocks(b)%A
      if (this%kind .ne. KIND_GENERAL                                                &
      & .and. this%blocks(b)%srow .ne. this%blocks(b)%scol) then
        dense(this%blocks(b)%js:this%blocks(b)%js + ncol - 1,                          &
        &     this%blocks(b)%is:this%blocks(b)%is + m - 1)                             &
        &   = this%kind * transpose(this%blocks(b)%A)
      endif
    enddo

    ! Column ruler with a '(s<sector>)' tag at each sector's first column.
    write(u, '(A)') ''
    line = repeat(' ', 6)
    do j = 1, this%n
      if (boundary(j)) line = line//'|'
      cell = ' '
      if (this%st(this%sector_of(j)) == j) then
        write(cell, '("(s",I0,")")') this%sector_of(j)
        cell = adjustl(cell)
      endif
      line = line//cell
    enddo
    line = line//'|'
    write(u, '(A)') line

    ! One line per global row, '-' lines between row sectors.
    do i = 1, this%n
      if (boundary(i)) then
        line = repeat(' ', 6)
        do j = 1, this%n
          if (boundary(j)) line = line//'|'
          line = line//repeat('-', 10)
        enddo
        line = line//'|'
        write(u, '(A)') line
      endif
      write(rowlab, '(I6)') i
      line = rowlab
      do j = 1, this%n
        if (boundary(j)) line = line//'|'
        cell = ' '
        if (dense(i, j) == 0.0d0) then
          cell(6:6) = '.'
        else
          write(cell, '(F10.4)') dense(i, j)
        endif
        line = line//cell
      enddo
      line = line//'|'
      write(u, '(A)') line
    enddo

    ! Closing line under the matrix.
    line = repeat(' ', 6)
    do j = 1, this%n
      if (boundary(j)) line = line//'|'
      line = line//repeat('-', 10)
    enddo
    line = line//'|'
    write(u, '(A)') line
    write(u, '(A)') ''
  end subroutine print_dense_sym_matrix
  !----------------------------------------------------------------------------

end module block_matrix_module
