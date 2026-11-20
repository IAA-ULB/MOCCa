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
  ! This module contains the machinery to efficiently handle the matrix
  !  representation of (generalized) single-particle operators in the
  !  configuration space spanned by a finite set of spwfs on a MOCCa grid.
  !
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Single-particle operators acting on spwfs can be classified with the
  !  following properties
  !
  !    - z-signature (Rz): even or odd
  !    - parity (P)      : even or odd
  !    - isospin (Tz)    : even or odd
  !    - hermitianness   : hermitian or antihermitian
  !    - y-time-simplex  : even or odd
  !    - time-reversal   : even or odd
  !
  ! Single-particle quantum numbers
  ! --------------------------------
  !
  ! Only three (Rz,P,Tz) of these correspond to quantum numbers of the
  !  single-particle states; hence any set of spwfs can be divided into
  !  sectors, labelled by a single 'symmetry word' or a 3bit integer
  !  s = 0..7 with the encoding
  !
  !   bit 0 : z-signature (Rz) : 0 = q(Rz) = +i, 1 = q(Rz) = -i
  !   bit 1 : parity (P)       : 0 = even      , 1 = odd
  !   bit 2 : isospin (Tz)     : 0 = neutron   , 1 = proton
  !
  ! The matrix-representation of a single-particle operator carries
  !  exactly this encoding: it sets these bits if the operator changes
  !  the corresponding quantum number. This means that selection rules
  !  can be implemented as a single XOR operator:
  !
  !      an operator couples sector s <-> sector IEOR(s,op)
  !
  !  Because of this property,
  !   (a) there is essentially no need for any if-statements that depend
  !       on Hephaestos symmetry-decisions
  !   (b) the matrix representation of an operator has AT MOST 8
  !       dense subblocks.
  !
  ! The resulting encoding is exactly the one assumed by MOCCa everywhere
  !  else, as accounted for in HFBlocks - except for a shift by 1:
  !
  !  B = 1, s = 0;  n, P=+1, Rz = +i   B = 5: p, P=+1, +i  => s = 4
  !  B = 2, s = 1;  n, P=+1, Rz = -i   B = 6: p, P=+1, -i  => s = 5
  !  B = 3, s = 2;  n, P=+1, Rz = +i   B = 7: p, P=-1, +i  => s = 6
  !  B = 4, s = 3;  n, P=+1, Rz = -i   B = 8: p, P=-1, -i  => s = 7
  !
  !  Many-body quantum numbers, broken symmetries, hermitianness
  !  ------------------------------------------------------------
  !
  !  The many-body state can carry two additional symmetries:
  !
  !  - y-time-simplex, S^T_y
  !  - time-reversal, T
  !
  !  These do not result in quantum numbers and have no associated
  !  selection rule and block structure. Operators however CAN anticommute
  !  with these, hence why the matrix representation of a one-body
  !  operator carries the logical "sty" and "tau".
  !
  !  Broken symmetries in a MOCCa calculation do NOT change the properties
  !  of any operator - the mathematical object does not change. What does
  !  change is the division of spwfs into blocks: the calculation loses
  !  a quantum number and the associated bit in the operator should have
  !  no more effect. This is why an operator carries also the set
  !   goodqns = "good quantum numbers."
  !
  ! Finally, operators can be either hermitian or antihermitian
  !
  !      O^\dagger = \pm O
  !
  ! The (anti)Hermitian and S^T_y classification can be combined
  !  to determine that the matrix can be built with purely
  !  real or purely imaginary elements and with definite symmetry.
  !
  !   hermitian, STY-even     -> real and symmetric     -> kind = +1
  !   hermitian, STY-odd      -> imag and antisymmetric -> kind = -1
  !   antihermitian, STY-odd  -> real and symmetric     -> kind = +1
  !   antihermitian, STY-even -> imag and antisymmetric -> kind = -1
  !
  ! This classification is computed and stored in the integer "KIND"
  ! of each block_matrix object.
  !
  ! Caveats - things that took me a while to internalize
  ! -----------------------------------------------------
  ! - the product of two (anti)hermitian matrices need not be
  !    either hermitian or antihermitian, it will generally be
  !    of mixed classification. This is one reason why this
  !    module does not include a routine to do matrix-matrix
  !    multiplication.
  !
  ! Additional remarks
  ! -------------------------------------------------------
  ! - this module is stand-alone: it can be compiled alone
  !   and unit-tested easily and possibly reused in other codes.
  ! - there is no MPI implementation here (for now), all
  !   storage is rank-local
  !==============================================================================
  implicit none

  !------------------------------------------------------------------------------
  ! The classification possibilities for an operator
  !  We store these as parameters: writing e.g. HERMITIAN is less error-prone
  !   and more easy to debug than a bunch of '+1' and '-1's.
  ! Anti-hermitianness
  integer, parameter :: HERMITIAN       =  1
  integer, parameter :: ANTIHERMITIAN   = -1
  ! Characteristic under S^T_y
  integer, parameter :: STY_EVEN =  1
  integer, parameter :: STY_ODD  = -1
  ! Storage kind of the stored real numbers
  !     A(s', s) = kind * A(s, s')^T
  integer, parameter :: KIND_SYMMETRIC      =  1
  integer, parameter :: KIND_ANTISYMMETRIC  = -1
  ! Time-reversal character
  integer, parameter :: TAU_EVEN =  1
  integer, parameter :: TAU_ODD  = -1

  !------------------------------------------------------------------------------
  ! Type for a single symmetry block in a block_matrix
  !  - rows    span global indices ( is : is+m-1 ) of sector srow,
  !  - columns span global indices ( js : js+n-1 ) of sector scol,
  ! with m = size(A,1), n = size(A,2), the sizes of the sectors.
  type :: sym_block
    integer :: srow, scol      ! sectors to which this submatrix corresponds
                               ! i.e. quantum numbers of the spwfs associated
                               ! with rows and columns.
    integer :: is, js          ! 1-based global start of the row/column sector
    real(8), allocatable :: A(:, :)  ! actual storage of real numbers
  end type sym_block

  !------------------------------------------------------------------------------
  ! The block_matrix: a flat list of at most 8 sym_blocks.
  type :: block_matrix
    integer :: symmetry          ! operator symmetry word, 0..7
    integer :: goodqn            ! good-quantum-number mask of the basis,
                                 ! 0..31; bits 0-2 as the sector words,
                                 ! bit 3 = STy conserved, bit 4 = T conserved
    integer :: kind              ! real-symmetric or imaginary-antisymmetric
    integer :: sty               ! characteristic under the y-time-simplex
    integer :: tau               ! time parity of the operator
    logical :: folded            ! .true. = the TR-folded composite <a|O T|b>
    integer :: n                 ! total basis size, sum(HFBlocks)
    integer :: nblocks           ! number of stored blocks, <= 8
    type(sym_block) :: blocks(8) ! actual storage of sym_blocks
    integer :: st(0:8)           ! first global index of each sector;
                                 ! st(8) = n + 1
  contains
    procedure :: init           => init_block_matrix
    procedure :: sector_of      => sector_of_block_matrix
    procedure :: fill_blocks    => fill_blocks_block_matrix
    procedure :: dense          => dense_block_matrix
    procedure :: print_dense    => print_dense_block_matrix
    procedure :: full_trace     => full_trace_block_matrix
    procedure :: rotate         => rotate_block_matrix
  end type block_matrix

  !----------------------------------------------------------------------------
  ! A single block of a matrix representing a unitary transformation
  type :: sector_rotation
    real(8), allocatable :: U(:, :)   ! (k, k), orthogonal: U^T U = I
  end type sector_rotation

  !----------------------------------------------------------------------------
  ! A basis rotation represented as a unitary transformation.
  ! I.e. a change of spwfs specified as
  !
  !                  psi' = C^T psi
  !
  ! where C is a unitary matrix that respects all symmetries; i.e. it
  ! only has real entries and is sector-diagonal.
  !
  ! Note: this transpose convention is exactly that of basis_transform.f90
  type :: sym_rotation
    integer :: n = 0                  ! total basis size, sum(HFBlocks)
    integer :: st(0:8) = 0            ! sector bounds; st(0) = 1 marks an
                                      !  initialised object, default 0 not
    type(sector_rotation) :: U(0:7)   ! populated sectors only
  contains
    procedure :: init            => init_sym_rotation
    procedure :: init_from_dense => init_from_dense_sym_rotation
    procedure :: set_sector      => set_sector_sym_rotation
  end type sym_rotation

  abstract interface
      subroutine block_filler(srow, scol, A, params)
      !---------------------------------------------------------------------------
      ! Overwrite one whole block of a block matrix.
      !
      ! This routine is the abstract interface for operator-specific routines
      !  that implement the matrix element calculation. Such routines should
      !  depend on
      !   (a) the sector words (= quantum numbers) srow/scol of the block.
      !   (b) any other parameters that can be passed in through a type "params"
      !
      ! Note that A is an intent(out) dummy: the routine must overwrite the
      !  entire block and cannot write individual variables.
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! Input:
      !   srow, scol : sector words 0..7 of the block
      !   params     : additional parameters
      ! Output:
      !   A          : the block storage, to be overwritten entirely
      !---------------------------------------------------------------------------
      integer, intent(in)   :: srow, scol
      real(8), intent(out)  :: A(:, :)
      class(*), intent(in)  :: params
    end subroutine block_filler
  end interface

  private

  public :: sym_block, block_matrix
  public :: sector_rotation, sym_rotation
  public :: HERMITIAN, ANTIHERMITIAN
  public :: KIND_SYMMETRIC, KIND_ANTISYMMETRIC
  public :: STY_EVEN, STY_ODD
  public :: TAU_EVEN, TAU_ODD
  public :: block_filler

contains

  subroutine init_block_matrix(this, HFBlocks, symmetry, herm, sty, tau,          &
  &                           folded, goodqn)
    !---------------------------------------------------------------------------
    ! Construct a block_matrix based on the classification of the operator.
    !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !   HFBlocks(1:8) : GLOBAL spwf counts per MOCCa block (header table);
    !                   the rank-local twin must not be passed.
    !   symmetry      : operator symmetry word, 0..7
    !   herm          : HERMITIAN or ANTIHERMITIAN, i.e. +1/-1
    !   sty           : STY_EVEN or STY_ODD       , i.e. +1/-1
    !   tau           : TAU_EVEN or TAU_ODD       , i.e. +1/-1
    !   folded        : optional, default .false. - .true. marks the object
    !                   as the TR-folded composite <a|O T|b> of design 2d;
    !                   identity metadata, no fill logic reads it
    !   goodqn        : optional good-quantum-number mask of the basis
    !---------------------------------------------------------------------------
    class(block_matrix), intent(inout) :: this
    integer, intent(in)               :: HFBlocks(8)
    integer, intent(in)               :: symmetry
    integer, intent(in)               :: herm
    integer, intent(in)               :: sty
    integer, intent(in)               :: tau
    logical, intent(in), optional     :: folded
    integer, intent(in), optional      :: goodqn

    integer :: s, sp, b, qn, word, kind
    integer :: counts(0:7)

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Guards: sector words are 3-bit integers, kinds are the module constants.
    if (symmetry < 0 .or. symmetry > 7) then
      print *, 'init_block_matrix: symmetry word must be in 0..7, got', symmetry
      stop 1
    endif
    if (herm .ne. HERMITIAN .and. herm .ne. ANTIHERMITIAN) then
      print *, 'init_block_matrix: herm must be HERMITIAN or ANTIHERMITIAN,',     &
      &         ' got', herm
      stop 1
    endif
    if (sty .ne. STY_EVEN .and. sty .ne. STY_ODD) then
      print *, 'init_block_matrix: sty must be STY_EVEN or STY_ODD, got', sty
      stop 1
    endif
    if (tau .ne. TAU_EVEN .and. tau .ne. TAU_ODD) then
      print *, 'init_block_matrix: tau must be TAU_EVEN or TAU_ODD, got', tau
      stop 1
    endif
    if (any(HFBlocks < 0)) then
      print *, 'init_block_matrix: HFBlocks entries must be >= 0'
      stop 1
    endif

    ! Bits 3 and 4 of the mask are stored verbatim: the consumer's
    !  completeness criterion (STy) and doubling scenario (T), not
    !  constructor gates.
    qn = 15
    if (present(goodqn)) qn = goodqn
    if (qn < 0 .or. qn > 31) then
      print *, 'init_block_matrix: goodqn mask must be in 0..31, got', qn
      stop 1
    endif

    ! The only place MOCCa's 1-based block ordering is converted to 0-based
    !  sector words: sector s = MOCCa block s + 1 (see table in module header).
    counts(0:7) = HFBlocks(1:8)

    this%symmetry = symmetry
    this%goodqn = qn
    this%n = sum(counts)

    ! kind fixes the partner sign at access time; the block list built
    !  below depends only on its being nonzero.
    this%sty = sty
    kind = herm*sty
    this%kind = kind

    ! Descriptive TR metadata (D-7, 2026-11-20), stored verbatim: no
    !  constructor logic derives from either.
    this%tau = tau
    this%folded = .false.
    if (present(folded)) this%folded = folded

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
      if (s .gt. sp) cycle

      b = b + 1
      this%blocks(b)%srow = sp
      this%blocks(b)%scol = s
      this%blocks(b)%is = this%st(sp)
      this%blocks(b)%js = this%st(s)
      allocate(this%blocks(b)%A(counts(sp), counts(s)))
      this%blocks(b)%A = 0.0d0
    enddo
    this%nblocks = b
  end subroutine init_block_matrix
  !----------------------------------------------------------------------------
  pure integer function sector_of_block_matrix(this, i)
    !---------------------------------------------------------------------------
    ! Map a global spwf index to its symmetry sector. Returns -1 if the index
    !  is out of range. Exploits that the bounds array is non-decreasing.
    !---------------------------------------------------------------------------
    class(block_matrix), intent(in) :: this
    integer, intent(in)            :: i

    integer :: s

    sector_of_block_matrix = -1
    if (i < 1 .or. i > this%n) return
    do s = 7, 0, -1
      if (i .ge. this%st(s)) then
        sector_of_block_matrix = s
        exit
      endif
    enddo
  end function sector_of_block_matrix
  !----------------------------------------------------------------------------
  subroutine fill_blocks_block_matrix(this, fill, params)
    !---------------------------------------------------------------------------
    ! Populate all stored blocks through the callback: one call per stored
    !  block, overwriting the storage handed out as an intent(out) dummy.
    !  Unstored partners are never touched (dense reconstructs them);
    !  distinct blocks occupy disjoint storage, so the loop is safe to
    !  parallelise over blocks.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Input/output:
    !   this   : the matrix whose stored blocks are populated
    !   fill   : callback overwriting one block (interface block_filler)
    !   params : state handed to the callback, verbatim; callers without
    !            state pass any token
    !---------------------------------------------------------------------------
    class(block_matrix), intent(inout) :: this
    procedure(block_filler)           :: fill
    class(*), intent(in)              :: params

    integer :: b

    do b = 1, this%nblocks
      call fill(this%blocks(b)%srow, this%blocks(b)%scol,                        &
      &         this%blocks(b)%A, params)
    enddo
  end subroutine fill_blocks_block_matrix
  !----------------------------------------------------------------------------
  function dense_block_matrix(this) result(D)
    !---------------------------------------------------------------------------
    ! Partner-reconstructing dense expansion: the full n x n matrix with
    !  all stored blocks at their global positions, unstored partners as
    !      A(s', s) = kind * A(s, s')^T ,
    !  zero elsewhere. The return shape of dense wrappers such as
    !  Qlm_spme_fm in moments.f90 and of the rotation unit tests.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Input:
    !   this : the matrix whose stored blocks are expanded
    ! Output:
    !   D    : the dense n x n expansion, freshly allocated
    !---------------------------------------------------------------------------
    class(block_matrix), intent(in) :: this
    real(8), allocatable :: D(:, :)

    integer :: b, m, ncol

    allocate(D(this%n, this%n))
    D = 0.0d0
    do b = 1, this%nblocks
      if (.not. allocated(this%blocks(b)%A)) cycle
      m = size(this%blocks(b)%A, 1)
      ncol = size(this%blocks(b)%A, 2)
      D(this%blocks(b)%is:this%blocks(b)%is + m - 1,                              &
      &   this%blocks(b)%js:this%blocks(b)%js + ncol - 1)                           &
      &   = this%blocks(b)%A
      if (this%blocks(b)%srow .ne. this%blocks(b)%scol) then
        D(this%blocks(b)%js:this%blocks(b)%js + ncol - 1,                            &
        &   this%blocks(b)%is:this%blocks(b)%is + m - 1)                             &
        &   = this%kind * transpose(this%blocks(b)%A)
      endif
    enddo
  end function dense_block_matrix
  !----------------------------------------------------------------------------
  subroutine print_dense_block_matrix(this, lun, title)
    !---------------------------------------------------------------------------
    ! Human-friendly dump of a block_matrix as one dense n x n matrix,
    !  partners reconstructed by the dense expansion. Debug/test helper
    !  for SMALL matrices ONLY - memory and output scale as n**2.
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
    class(block_matrix), intent(in) :: this
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
      write(u, '(A)') ' block_matrix is empty or not initialised, nothing to dump'
      return
    endif

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Header: symmetry word decoded into the bits it flips, the kind, the
    !  stored block list and the populated sectors with index ranges.
    write(u, '(A,I0)') ' block_matrix: n = ', this%n

    if (this%symmetry == 0) then
      line = 'sector-diagonal operator, flips nothing'
    else
      line = 'flips:'
      if (btest(this%symmetry, 0)) line = line//' Rz'
      if (btest(this%symmetry, 1)) line = line//' P'
      if (btest(this%symmetry, 2)) line = line//' Tz'
    endif
    write(u, '(A,I0,A,A)') '  symmetry word ', this%symmetry, ': ', trim(line)
    if (this%goodqn .ne. 15) then
      write(u, '(A,I0,A)') '  good quantum numbers mask ', this%goodqn,     &
      &                    ' (design 2b: bits set = quantum numbers good)'
    endif
    if (this%sty .eq. STY_EVEN) then
      write(u, '(A)') '  y-time-simplex even: stored numbers are the real part'
    else
      write(u, '(A)')                                                        &
      & '  y-time-simplex odd: stored numbers are the imaginary part'
    endif
    if (.not. btest(this%goodqn, 3)) then
      write(u, '(A)')                                                        &
      & '  STy not conserved in this build: one half of a complex operator'
    endif
    if (this%folded) then
      write(u, '(A)')                                                        &
      & '  TR-folded composite <a|O T|b>: full_trace refuses this object'
    endif
    if (this%tau .eq. TAU_EVEN) then
      line = 'time-even'
    else
      line = 'time-odd'
    endif
    if (btest(this%goodqn, 4)) then
      write(u, '(A,A,A)') '  operator ', trim(line),                         &
      &   ': full trace over the T-doubled basis (goodqn bit 4 set)'
    else
      write(u, '(A,A,A)') '  operator ', trim(line),                          &
      &   ': full trace is the plain numerical trace (goodqn bit 4 clear)'
    endif

    select case (this%kind)
    case (KIND_SYMMETRIC)
      line = 'symmetric'
    case (KIND_ANTISYMMETRIC)
      line = 'antisymmetric'
    case default
      line = 'UNKNOWN'
    end select
    write(u, '(A,A,A,A,I0)') '  kind: ', trim(line), ' (derived, design 2a),',   &
    &                      ' stored blocks: ', this%nblocks

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
    ! The dense picture: the partner-reconstructing dense expansion of the
    !  dense helper - (anti-)symmetric matrices are printed complete since
    !  init does not store their partner blocks.
    !  Sector boundaries: column/row i starts a sector boundary if it is
    !  the first index of sectors 1..7 (empty sectors leave no gap, so two
    !  successive boundaries can coincide; boundary 1 is never marked).
    allocate(boundary(this%n))
    boundary = .false.
    do s = 1, 7
      if (this%st(s) >= 2 .and. this%st(s) <= this%n) boundary(this%st(s)) = .true.
    enddo

    dense = this%dense()

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
  end subroutine print_dense_block_matrix
  !----------------------------------------------------------------------------
  function full_trace_block_matrix(this) result(t)
    !---------------------------------------------------------------------------
    ! The honest operator's physical trace as one complex number (D-7,
    !  2026-11-20; module header for the conventions).
    !  With goodqn bit 4 set - the stored basis is the q(Rz) = +i half of
    !  a T-symmetric set - the trace over the T-doubled basis (stored
    !  states plus their time-reversal images) is
    !      Tr_full = Tr_stored + tau * conj(Tr_stored)
    !  by mechanism 1b of design 15.8, <Tu|O|Tv> = tau * conj(<u|O|v>):
    !  2 Re, 2 i Im, or 0, by (tau, sty). With bit 4 clear the stored
    !  basis is the whole physical basis and the plain numerical trace
    !  is the physical one. Signature-flipping words contribute nothing
    !  either way: no diagonal blocks exist, the diagonal vanishes by
    !  the XOR selection rules.
    !  Refuses folded composites: their diagonal entries are fold
    !  diagonals <a|O T|a>, a different physical statement that the
    !  design keeps D-5-gated (DECISIONS.md).
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !   this : the honest (folded = .false.) matrix whose trace is asked
    ! Output:
    !   t    : the physical trace, complex
    !---------------------------------------------------------------------------
    class(block_matrix), intent(in) :: this
    complex(8) :: t, ts

    integer :: b, i, m
    real(8) :: sdiag

    if (this%folded) then
      print *, 'full_trace: refusal - the object is the TR-folded composite',  &
      &         ' <a|O T|b>: its diagonals are fold diagonals, not the trace'
      stop 1
    endif

    ! The complex trace of the stored half: the stored numbers are the
    !  real (STY_EVEN) or imaginary (STY_ODD) part of each element, so the
    !  diagonal sum acquires the factor 1 or i once, not per element.
    sdiag = 0.0d0
    do b = 1, this%nblocks
      if (this%blocks(b)%srow .ne. this%blocks(b)%scol) cycle
      if (.not. allocated(this%blocks(b)%A)) cycle
      m = size(this%blocks(b)%A, 1)
      do i = 1, m
        sdiag = sdiag + this%blocks(b)%A(i, i)
      enddo
    enddo
    ts = cmplx(sdiag, 0.0d0, 8)
    if (this%sty .eq. STY_ODD) ts = (0.0d0, 1.0d0)*ts

    if (btest(this%goodqn, 4)) then
      t = ts + this%tau*conjg(ts)
    else
      t = ts
    endif
  end function full_trace_block_matrix

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! BASIS ROTATIONS (design section 6a)
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  subroutine rotate_block_matrix(this, R)
    !---------------------------------------------------------------------------
    ! Rotate a block_matrix into the basis defined by the sector-diagonal
    !  orthonormal blocks of R: every stored block transforms independently
    !  as  A_b <- U_{srow}^T A_b U_{scol} , the same operator in the rotated
    !  basis. Unstored partners are never touched - the kind relation
    !  survives, so dense() reconstructs them as before - and every
    !  metadata field is kept (goodqn plays no role: rotate acts on storage
    !  sectors, never on quantum numbers). A similarity transform, not an
    !  operator product: products of classifiable operators are
    !  unclassifiable and stay out of the type.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input/output:
    !   this : the matrix whose stored blocks are rotated, in place
    !   R    : the rotation; must share this matrix's global HFBlocks
    !---------------------------------------------------------------------------
    class(block_matrix), intent(inout) :: this
    class(sym_rotation), intent(in) :: R

    integer :: b, srow, scol
    real(8), allocatable :: tmp(:, :)

    if (R%n .ne. this%n .or. any(R%st .ne. this%st)) then
      print *, 'rotate: rotation and matrix must share the same global',       &
      &         ' HFBlocks (n or st differ)'
      stop 1
    endif

    do b = 1, this%nblocks
      if (.not. allocated(this%blocks(b)%A)) cycle
      srow = this%blocks(b)%srow
      scol = this%blocks(b)%scol
      if (.not. allocated(R%U(srow)%U) .or.                                     &
      &   .not. allocated(R%U(scol)%U)) then
        print *, 'rotate: rotation lacks a block of a populated sector'
        stop 1
      endif
      tmp = matmul(this%blocks(b)%A, R%U(scol)%U)
      this%blocks(b)%A = matmul(transpose(R%U(srow)%U), tmp)
    enddo
  end subroutine rotate_block_matrix
  !----------------------------------------------------------------------------
  subroutine init_sym_rotation(this, HFBlocks)
    !---------------------------------------------------------------------------
    ! Construct the IDENTITY rotation of the basis described by HFBlocks:
    !  one k x k identity block per populated sector, nothing for empty
    !  sectors. Rotating with the result of init is an exact no-op, and
    !  set_sector/init_from_dense can overwrite any of its blocks later.
    !  init may be called repeatedly; previous state is discarded cleanly.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !   HFBlocks(1:8) : number of spwfs per MOCCa symmetry block, the
    !                   GLOBAL count array (HFBlocks_global); see init
    !                   of block_matrix for the block-to-sector table
    !---------------------------------------------------------------------------
    class(sym_rotation), intent(inout) :: this
    integer, intent(in) :: HFBlocks(8)

    integer :: s, k, i
    integer :: counts(0:7)

    if (any(HFBlocks < 0)) then
      print *, 'init_sym_rotation: HFBlocks entries must be >= 0'
      stop 1
    endif

    counts(0:7) = HFBlocks(1:8)
    this%n = sum(counts)

    this%st(0) = 1
    do s = 0, 7
      this%st(s + 1) = this%st(s) + counts(s)
    enddo

    do s = 0, 7
      if (allocated(this%U(s)%U)) deallocate(this%U(s)%U)
      k = counts(s)
      if (k .eq. 0) cycle
      allocate(this%U(s)%U(k, k))
      this%U(s)%U = 0.0d0
      do i = 1, k
        this%U(s)%U(i, i) = 1.0d0
      enddo
    enddo
  end subroutine init_sym_rotation
  !----------------------------------------------------------------------------
  subroutine init_from_dense_sym_rotation(this, HFBlocks, C, tol)
    !---------------------------------------------------------------------------
    ! Construct a rotation from a FULL n x n rotation array C, in the
    !  basis_transform convention: each COLUMN of C is a new basis vector
    !  expressed in the old basis (psi' = C^T psi). Only the
    !  sector-diagonal parts of C are kept, one square block per
    !  populated sector.
    !
    ! Loud structure checks, both against the tolerance 'tol'
    !  (default 1e-12 absolute): C must be sector-diagonal - every
    !  off-sector submatrix vanishes within tol - and every extracted
    !  block must be orthonormal within tol. A C that mixes sectors is
    !  REFUSED: the rotated basis would have no definite good quantum
    !  numbers and no block_matrix could represent operators over it.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !   HFBlocks(1:8) : global MOCCa count array, see init_sym_rotation
    !   C             : the full n x n rotation array
    !   tol           : optional structure tolerance, default 1e-12
    !---------------------------------------------------------------------------
    class(sym_rotation), intent(inout) :: this
    integer, intent(in) :: HFBlocks(8)
    real(8), intent(in) :: C(:, :)
    real(8), intent(in), optional :: tol

    integer :: s, sp
    real(8) :: tol_, dev

    call this%init(HFBlocks)

    if (size(C, 1) .ne. this%n .or. size(C, 2) .ne. this%n) then
      print *, 'init_from_dense: C must be the full n x n rotation, got',     &
      &         size(C, 1), 'x', size(C, 2), ' for n =', this%n
      stop 1
    endif

    tol_ = 1.0d-12
    if (present(tol)) tol_ = tol
    if (tol_ <= 0.0d0) then
      print *, 'init_from_dense: tolerance must be positive, got', tol_
      stop 1
    endif

    do s = 0, 7
      if (this%st(s + 1) <= this%st(s)) cycle

      do sp = 0, 7
        if (sp .eq. s) cycle
        if (this%st(sp + 1) <= this%st(sp)) cycle
        dev = maxval(abs(C(this%st(s):this%st(s + 1) - 1,                       &
        &                    this%st(sp):this%st(sp + 1) - 1)))
        if (dev .gt. tol_) then
          print *, 'init_from_dense: C is not sector-diagonal: the',           &
          &         ' off-sector block of sectors ', s, ' and ', sp,             &
          &         ' has max abs', dev
          stop 1
        endif
      enddo

      this%U(s)%U = C(this%st(s):this%st(s + 1) - 1,                             &
      &                this%st(s):this%st(s + 1) - 1)
      dev = deviation_orthonormal(this%U(s)%U)
      if (dev .gt. tol_) then
        print *, 'init_from_dense: the block of sector', s, 'is not',           &
        &         ' orthonormal: max |U^T U - I| =', dev
        stop 1
      endif
    enddo
  end subroutine init_from_dense_sym_rotation
  !----------------------------------------------------------------------------
  subroutine set_sector_sym_rotation(this, s, U, tol)
    !---------------------------------------------------------------------------
    ! Hand the orthonormal rotation block of a single sector directly -
    !  the natural entry for unit tests and for callers that own
    !  per-sector rotations (e.g. diagonalising one sector at a time).
    !
    ! Guards, all loud: the object must have been initialised (init or
    !  init_from_dense), s must be a populated sector of the basis, U
    !  must be square of the sector's size, and U must be orthonormal
    !  within the tolerance 'tol' (default 1e-12 absolute).
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !   s   : the sector the block belongs to, 0..7
    !   U   : the square orthonormal block of that sector
    !   tol : optional orthonormality tolerance, default 1e-12
    !---------------------------------------------------------------------------
    class(sym_rotation), intent(inout) :: this
    integer, intent(in) :: s
    real(8), intent(in) :: U(:, :)
    real(8), intent(in), optional :: tol

    integer :: k
    real(8) :: tol_, dev

    if (this%st(0) .ne. 1 .or. this%n .lt. 0) then
      print *, 'set_sector: call init or init_from_dense first'
      stop 1
    endif
    if (s < 0 .or. s > 7) then
      print *, 'set_sector: sector word must be in 0..7, got', s
      stop 1
    endif
    k = this%st(s + 1) - this%st(s)
    if (k .eq. 0) then
      print *, 'set_sector: sector', s, 'is empty in this basis'
      stop 1
    endif
    if (size(U, 1) .ne. k .or. size(U, 2) .ne. k) then
      print *, 'set_sector: the block of sector', s, 'must be square of',    &
      &         'size', k, ', got', size(U, 1), 'x', size(U, 2)
      stop 1
    endif

    tol_ = 1.0d-12
    if (present(tol)) tol_ = tol
    dev = deviation_orthonormal(U)
    if (dev .gt. tol_) then
      print *, 'set_sector: the block of sector', s, 'is not orthonormal:',   &
      &         ' max |U^T U - I| =', dev
      stop 1
    endif

    if (allocated(this%U(s)%U)) deallocate(this%U(s)%U)
    allocate(this%U(s)%U(k, k))
    this%U(s)%U = U
  end subroutine set_sector_sym_rotation
  !----------------------------------------------------------------------------
  pure function deviation_orthonormal(U) result(dev)
    !---------------------------------------------------------------------------
    ! Maximal deviation of a square block from orthonormality,
    !      dev = max |U^T U - I| ,
    ! zero for an empty block. Callers guarantee squareness; this helper
    ! is the single place the orthonormality checks of the rotation
    ! constructors are computed.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !   U   : the square block to measure
    ! Result:
    !   dev : the maximal deviation from the identity of U^T U
    !---------------------------------------------------------------------------
    real(8), intent(in) :: U(:, :)
    real(8) :: dev

    real(8), allocatable :: G(:, :)
    integer :: k, i

    dev = 0.0d0
    k = size(U, 1)
    if (k .eq. 0) return

    G = matmul(transpose(U), U)
    do i = 1, k
      G(i, i) = G(i, i) - 1.0d0
    enddo
    dev = maxval(abs(G))
  end function deviation_orthonormal
  !----------------------------------------------------------------------------

end module block_matrix_module
