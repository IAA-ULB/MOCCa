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
  !  wavefunctions on the MOCCa grid. It implements the design summary
  !  'mocca-symmetry-block-matrix-design-summary.md', whose sections are
  !  referenced below: sector encoding, the sym_block / block_matrix types,
  !  the block-list construction from MOCCa's HFBlocks bookkeeping, block
  !  population through the fill_blocks callback, the partner-reconstructing
  !  dense expansion, and the sector-diagonal basis rotations (section 6a)
  !  - real orthogonal only, one square orthonormal block per populated
  !  sector, the shape MOCCa's own producers make (HFB's Canonical
  !  diagonalises rho per MOCCa block, so cantransfo is exactly
  !  sector-diagonal). The mat-vec / mat-mat products once planned for
  !  later phases are dropped: no consumer exists (section 6).
  !
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Symmetry model
  ! ---------------
  ! Operators acting on the HF basis carry three discrete quantum numbers,
  !  each of which labels the basis states. A symmetry-sector word
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
  !    parity broken                  :  spwfs in blocks 1,2,5,6 only, as
  !                                      parity is not a good quantum number
  !                                      (the P bit of any operator is then
  !                                      vacuous)
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
  !  init therefore takes an optional 5-bit mask 'goodqn' (bits 0-2 as
  !  the sector words, bit 3 = the y-time-simplex conserved, see the next
  !  section, bit 4 = time reversal conserved, the scenario bit of the
  !  doubling section below; default 15 = every sector grading and the
  !  simplex good, doubling off) marking the
  !  quantum numbers that actually label the basis, and pairs every
  !  sector s with IEOR(s, AND(symmetry, goodqn)) - bits 0-2 enter the
  !  pairing; bit 3 is stored verbatim for consumers of the characteristic.
  !  The mask is a fact of the build scenario, not recoverable from the
  !  runtime block counts: in MOCCa it comes from the geninfo constant
  !  'goodqn', injected by Hephaestos ($GOODQN). With all bits good the
  !  pairing is the pure XOR above. Standalone callers (unit tests) omit
  !  the mask.
  !
  !  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !  Characteristic under the y-time-simplex
  !  ----------------------------------------
  !  The y-time-simplex  STy = P . T . Ry  is a hermitian antilinear
  !  involution (STy**2 = +1), and the spwfs MOCCa stores are its
  !  invariants,
  !       STy psi = psi^*(x,-y,z,sigma) = psi .
  !  An operator O of definite characteristic,
  !       STy O STy^{-1} = sigma * O        (sigma = +/- 1),
  !  therefore has matrix elements over the stored basis that obey
  !       < i | O | j > = sigma * < i | O | j >*  :
  !  STY_EVEN operators have purely REAL, STY_ODD operators purely
  !  IMAGINARY matrix elements.
  !
  !  STy carries NO block structure: every stored state is an invariant,
  !  so there is no good quantum number to pair and the characteristic
  !  never touches the block list. It fixes how matrix elements are
  !  calculated (the real or the imaginary part of the spinor product)
  !  and how the stored numbers must be interpreted by consumers.
  !
  !  Every operator handled by this type is classifiable (design summary
  !  section 2a): it is Hermitian or anti-Hermitian, and it is even or
  !  odd under the y-time-simplex. init therefore takes both physical
  !  classes as MANDATORY inputs, 'herm' (HERMITIAN / ANTIHERMITIAN)
  !  and 'sty' (STY_EVEN / STY_ODD), and DERIVES the storage kind from
  !  them:
  !
  !      kind = herm * sty
  !
  !  i.e. the stored numbers - the real part (STY_EVEN) or the imaginary
  !  part (STY_ODD) of the matrix elements - form a symmetric matrix for
  !  a Hermitian STY_EVEN or an anti-Hermitian STY_ODD operator, and an
  !  antisymmetric one for a Hermitian STY_ODD or an anti-Hermitian
  !  STY_EVEN operator.
  !
  !  In a build that conserves the time simplex (goodqn bit 3, the case
  !  for every SYMSTRING in the tree) the characteristic means more: the
  !  omitted half of the complex matrix elements is then identically
  !  zero and the single stored matrix IS the complete operator. With
  !  bit 3 clear the same inputs still define a consistent single matrix
  !  - the selected half of the (re, im) pair model of design section
  !  12, the other half a separate matrix - but no completeness claim
  !  can be made; consumers select on the stored goodqn mask themselves.
  !  init consumes bits 0-2 of the mask for the pairing and stores the
  !  mask verbatim.
  !
  !  Producers hand 'sty' to their fill callbacks through params (the
  !  callback interface is frozen): the callback computes the real or
  !  the imaginary part of the spinor product accordingly. This presumes
  !  the +1 invariance convention above; a basis with STy psi = -psi
  !  states would need per-state signs, which MOCCa's conventions do
  !  not produce.
  !
  !  Typical characteristics: r^l Re Y_lm is
  !  STY_EVEN and r^l Im Y_lm is STY_ODD for every m, and J_y is
  !  STY_ODD. Character tables per operator family are producer-side
  !  knowledge; for operators acting on spin (e.g. spin densities) they
  !  depend on MOCCa's spin convention for STy and must be taken from
  !  the symmetry machinery, not re-derived here.
  !
  !  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !  Time parity and the doubling scenario (D-7, 2026-11-20)
  !  -----------------------------------------------------
  !  Two descriptive metadata components carry what the generator
  !  entry-point names used to be the only trace of:
  !   *) tau - the operator's time parity T O T^-1 = tau O, a mandatory
  !      init class like herm and sty (TAU_EVEN / TAU_ODD); the folded
  !      composite <a|O T|b> inherits its parent's. No fill logic reads
  !      it - full_trace is the one in-module consumer.
  !   *) folded - .true. for the TR-folded composite (design 2d): the
  !      honest-vs-composite identity, stored verbatim.
  !  goodqn bit 4 ("T good") is the scenario bit: set = the stored basis
  !  is the q(Rz) = +i half of a T-symmetric set and the doubled-basis
  !  contracts are defined; clear = the basis is the whole physical
  !  basis and they degenerate. Marker-free code reads the scenario off
  !  the object instead of the $TR/$NTR marker.
  !
  !  full_trace returns the honest operator's physical trace as one
  !  complex number: with bit 4 set, over the T-doubled basis,
  !      Tr_full = Tr_stored + tau * conj(Tr_stored)
  !  (design 15.8 mechanism 1b: <Tu|O|Tv> = tau * conj(<u|O|v>)) - i.e.
  !  2 Re, 2 i Im, or 0 by (tau, sty); with bit 4 clear the plain
  !  numerical trace. Signature-flipping words contribute nothing: no
  !  diagonal blocks exist, the diagonal vanishes by selection rules.
  !  FOLDED composites are refused: their diagonals are fold diagonals
  !  <a|O T|a>, a different statement, kept D-5-gated by design.
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
  !  exclusively through the dummy-procedure callback of the abstract
  !  interface block_filler (whole blocks); procedure pointers are
  !  never stored as components:
  !   *) callbacks run for stored blocks only. Blocks skipped by the
  !      kind-based partner rule are not the filler's business: partner
  !      values are resolved at access time (a later phase), never by
  !      writing here.
  !   *) callbacks receive the sector words srow/scol of the block so
  !      they can apply MOCCa's phase conventions for matrix elements.
  !   *) the callback OVERWRITES: fill_blocks hands the block storage
  !      out as an intent(out) dummy.
  !   *) 'params' is unlimited-polymorphic state handed verbatim to the
  !      callback. The signature is frozen: extend params, not the
  !      interface. Callers without state pass any token.
  !==============================================================================
  implicit none

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Physical classes of an operator (design summary section 2a). Every
  !  operator handled by this type is classifiable: it is Hermitian or
  !  anti-Hermitian under dagger conjugation, and it is even or odd
  !  under the y-time-simplex. Both classes are MANDATORY init inputs;
  !  the storage kind below is derived from them as kind = herm * sty.
  integer, parameter :: HERMITIAN       =  1
  integer, parameter :: ANTIHERMITIAN   = -1

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Characteristic of an operator under the y-time-simplex STy = P.T.Ry
  !  (design summary section 2a, and the module header above for the full
  !  treatment): STY_EVEN operators have purely real, STY_ODD operators
  !  purely imaginary matrix elements over the stored basis (MOCCa's
  !  spwfs are STy invariants). Mandatory init input selecting which
  !  part of the complex matrix elements the stored numbers are.
  integer, parameter :: STY_EVEN =  1
  integer, parameter :: STY_ODD  = -1

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Storage kind of the stored real numbers: DERIVED in init from the
  !  physical classes,  kind = herm * sty . Partner blocks are related by
  !     A(s', s) = kind * A(s, s')^T
  !  and only one block per sector pair is stored; consumers reconstruct
  !  partners through this same relation (dense). A block whose sector is
  !  its own partner is square and internally symmetric/antisymmetric; it
  !  is stored dense and the symmetry serves only as a consistency check.
  integer, parameter :: KIND_SYMMETRIC      =  1
  integer, parameter :: KIND_ANTISYMMETRIC  = -1

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Time parity of an operator (D-7, 2026-11-20; header section above for
  !  the doubling semantics): T O T^-1 = tau * O, a physical class like
  !  herm and sty; the folded composite OT inherits its parent's.
  !  Mandatory init input, descriptive metadata: init stores it verbatim
  !  and full_trace is the one in-module reader.
  integer, parameter :: TAU_EVEN =  1
  integer, parameter :: TAU_ODD  = -1

  !------------------------------------------------------------------------------
  ! One dense block: rows span global indices is : is+m-1 of sector srow,
  !  columns span js : js+n-1 of sector scol, with m = size(A,1), n = size(A,2).
  type :: sym_block
    integer :: srow, scol      ! sector words 0..7 of the row/column sector
    integer :: is, js          ! 1-based global start of the row/column sector
    real(8), allocatable :: A(:, :)  ! (m, n), contiguous, LDA = m
  end type sym_block

  !------------------------------------------------------------------------------
  ! The full symmetry-adapted matrix: at most 8 blocks in a flat list.
  !  After init, all runtime code loops over blocks(1:nblocks) only;
  !  quantum numbers XOR selection rules never appear at runtime again.
  type :: block_matrix
    integer :: symmetry       ! operator symmetry word, 0..7 (physical)
    integer :: goodqn         ! good-quantum-number mask of the basis,
                              ! 0..31 (design 2b; bits 0-2 as the sector
                              ! words, bit 3 = STy conserved, design 2a,
                              ! bit 4 = T conserved, the scenario bit of
                              ! the header's doubling section);
                              ! 15 = every grading/simplex good,
                              ! doubling off (the standalone default)
    integer :: kind           ! storage kind of the stored numbers,
                              ! DERIVED in init from the physical classes:
                              ! KIND_SYMMETRIC / KIND_ANTISYMMETRIC,
                              ! kind = herm * sty
    integer :: sty            ! characteristic under the y-time-simplex
                              ! (design 2a): STY_EVEN / STY_ODD as handed
                              ! to init. Carries no block structure: tells
                              ! producers and consumers which part of the
                              ! complex matrix elements the stored real
                              ! numbers are - the real part (even) or the
                              ! imaginary part (odd). This matrix is the
                              ! complete operator iff goodqn bit 3 is set.
    integer :: tau            ! time parity of the operator (D-7, 2026-11-
                              ! 20): TAU_EVEN / TAU_ODD as handed to init;
                              ! composites inherit the parent's. Descriptive
                              ! metadata: no fill logic reads it; full_trace
                              ! is the one in-module consumer
    logical :: folded         ! .true. = the TR-folded composite <a|O T|b>
                              ! (design 2d), not the honest operator: the
                              ! identity metadata of D-7, stored verbatim;
                              ! the doubling contracts refuse such objects
    integer :: n              ! total basis size, sum(HFBlocks)
    integer :: nblocks        ! number of stored blocks, <= 8
    type(sym_block) :: blocks(8)
    integer :: st(0:8)        ! first global (1-based) index of each sector;
                              !  st(8) = n + 1; sector s spans st(s) .. st(s+1)-1
  contains
    procedure :: init => init_block_matrix
    procedure :: sector_of => sector_of_block_matrix
    procedure :: fill_blocks => fill_blocks_block_matrix
    procedure :: dense => dense_block_matrix
    procedure :: print_dense => print_dense_block_matrix
    procedure :: full_trace => full_trace_block_matrix
    procedure :: rotate => rotate_block_matrix
  end type block_matrix

  !----------------------------------------------------------------------------
  ! One square orthonormal rotation block belonging to a single sector
  !  (design section 6a). Left unallocated for empty sectors and for
  !  objects that have not been through init yet.
  type :: sector_rotation
    real(8), allocatable :: U(:, :)   ! (k, k), orthogonal: U^T U = I
  end type sector_rotation

  !----------------------------------------------------------------------------
  ! A basis rotation: one real orthonormal block per populated sector,
  !  i.e. the sector-diagonal change of basis  psi' = C^T psi  whose
  !  columns are the new basis vectors expressed in the old basis (the
  !  basis_transform convention). Applied to a block_matrix through rotate,
  !  it expresses the SAME operator in the rotated basis and preserves
  !  every metadata field (design section 6a).
  type :: sym_rotation
    integer :: n = 0                  ! total basis size, sum(HFBlocks)
    integer :: st(0:8) = 0            ! 1-based sector bounds, scheme 2:
                                      !  sector s spans st(s) .. st(s+1)-1;
                                      !  st(0) = 1 marks an initialised
                                      !  object, the default 0 does not
    type(sector_rotation) :: U(0:7)   ! populated sectors only
  contains
    procedure :: init => init_sym_rotation
    procedure :: init_from_dense => init_from_dense_sym_rotation
    procedure :: set_sector => set_sector_sym_rotation
  end type sym_rotation

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

  !----------------------------------------------------------------------------
  subroutine init_block_matrix(this, HFBlocks, symmetry, herm, sty, tau,          &
  &                           folded, goodqn)
    !---------------------------------------------------------------------------
    ! Construct the block list of a symmetry-adapted matrix: every populated
    !  sector s is paired with sp = IEOR(s, AND(symmetry, goodqn)) - quantum
    !  numbers that do not label the basis act diagonally - and each
    !  surviving pair gets one zero-filled dense block, with only one side
    !  per pair stored (dense reconstructs the partner); a sector that is
    !  its own partner gives a single square block. The storage kind is
    !  derived as kind = herm * sty (header 'Characteristic' section).
    !  init may be called repeatedly; previous state is discarded.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !   HFBlocks(1:8) : GLOBAL spwf counts per MOCCa block (header table);
    !                   the rank-local twin must not be passed.
    !   symmetry      : operator symmetry word, 0..7
    !   herm          : HERMITIAN or ANTIHERMITIAN
    !   sty           : STY_EVEN or STY_ODD - selects the stored half
    !   tau           : TAU_EVEN or TAU_ODD - the operator's time parity,
    !                   descriptive metadata (header 'Time parity' section)
    !   folded        : optional, default .false. - .true. marks the object
    !                   as the TR-folded composite <a|O T|b> of design 2d;
    !                   identity metadata, no fill logic reads it
    !   goodqn        : optional good-quantum-number mask of the basis,
    !                   default 15 = every grading/simplex good, doubling
    !                   off; in MOCCa pass the geninfo constant 'goodqn'
    !                   ($GOODQN).
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
    !  2026-11-20; header 'Time parity' section for the conventions).
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
