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
  !  covers phase 1 of the design: sector encoding, the sym_block and
  !  sym_matrix types and construction of the block list from MOCCa's
  !  HFBlocks bookkeeping. Later phases add block filling, mat-vec/mat-mat
  !  products and the complex-operator wrapper.
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
    integer :: symmetry       ! operator symmetry word, 0..7
    integer :: kind           ! KIND_GENERAL / KIND_SYMMETRIC / KIND_ANTISYMMETRIC
    integer :: n              ! total basis size, sum(HFBlocks)
    integer :: nblocks        ! number of stored blocks, <= 8
    type(sym_block) :: blocks(8)
    integer :: st(0:8)        ! first global (1-based) index of each sector;
                              !  st(8) = n + 1; sector s spans st(s) .. st(s+1)-1
  contains
    procedure :: init => init_sym_matrix
    procedure :: sector_of => sector_of_sym_matrix
  end type sym_matrix

  private
  public :: sym_block, sym_matrix
  public :: KIND_GENERAL, KIND_SYMMETRIC, KIND_ANTISYMMETRIC

contains

  !----------------------------------------------------------------------------
  subroutine init_sym_matrix(this, HFBlocks, symmetry, kind)
    !---------------------------------------------------------------------------
    ! Construct the block list of a symmetry-adapted matrix.
    !
    ! For every populated sector s, the operator couples it to sector
    !  sp = IEOR(s, symmetry) and to nothing else. A dense block is stored
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
    !   symmetry      : operator symmetry word, 0..7
    !   kind          : KIND_GENERAL, KIND_SYMMETRIC or KIND_ANTISYMMETRIC
    !---------------------------------------------------------------------------
    class(sym_matrix), intent(inout) :: this
    integer, intent(in)               :: HFBlocks(8)
    integer, intent(in)               :: symmetry
    integer, intent(in)               :: kind

    integer :: s, sp, b
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

    ! The only place MOCCa's 1-based block ordering is converted to 0-based
    !  sector words: sector s = MOCCa block s + 1 (see table in module header).
    counts(0:7) = HFBlocks(1:8)

    this%symmetry = symmetry
    this%kind = kind
    this%n = sum(counts)

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
    !  XOR; the kind-based skipping keeps at most one block per pair for
    !  (anti-)symmetric matrices.
    b = 0
    do s = 0, 7
      sp = ieor(s, symmetry)
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

end module block_matrix_module
