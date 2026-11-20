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
module angular_momentum
!==============================================================================
  !
  ! Single-particle matrix elements of the angular-momentum operators on
  !  the symmetry-block matrix type (design summary section 14).
  !
  !  One interface, three objects (D-7 of 2026-11-20, design 14.3
  !  amendment / 15.9 row 2): the public generators build the honest trio
  !
  !   J_x, J_y, J_z  - the one-body operators over the stored basis.
  !  J_x/J_y flip the z-signature (word 1), J_z conserves it (word 0);
  !  all three are Hermitian, parity-even and isoscalar; J_y is the
  !  STY_ODD member of the triad. In TR-conserving builds the word-1
  !  partners of the stored sectors are empty, so honest J_x/J_y come
  !  out as structural zeros there - the odd-m multipole precedent, no
  !  special case.
  !
  !  The TR-folded content <a|J_mu T|b> is served by the get_J_block
  !  accessor: the 15.8 retrieval contract of the gradient family in the
  !  J shape - any block of the complete doubled-basis matrix on demand,
  !  the fold materialised by the honest kernels on the TimeReverse
  !  images of image-side states, the J_z image-image blocks by the tau
  !  sign flip. The eager twin entry point of the first landing has
  !  retired; the t-kernels of wavefunctions.f90 remain the owner of the
  !  fold conventions, exercised by the unit arm as the reference the
  !  accessor packaging is asserted against. The composites are
  !  HERMITIAN as operators,
  !      (J.T)^dagger = T^-1 J = +(J.T)   for MOCCa's FERMIONIC T**2=-1,
  !  not a corollary of the linear-operator algebra of design 2a; their
  !  class tags are pinned to the kernel-transpose oracle of the
  !  angular_momentum_testing module (design 14.1): the oracle's mirror
  !  measurement is the authority the tags have to agree with,
  !  re-verified by the unit arm at every run - never re-derive the
  !  kind signs.
  !
  ! The generators wrap the angmom kernels of wavefunctions.f90 VERBATIM -
  !  one kernel invocation per stored pair, exactly the arithmetic of the
  !  legacy loops of calcJ2andBelyaev_HFB, so those loops remain the A/B
  !  oracle of the consumer migration (the step-3 precedent).
  !
  ! Scope decisions (design 14 and 14.3):
  !  - this module is free of $-markers; consumers gate by symmetry
  !    scenario downstream - only $TR-marked call sites may request
  !    image-involving blocks from the accessor.
  !  - single-node, serial semantics like the multipole generator: the
  !    basis is indexed globally through HFpsi, matching the matrix's own
  !    HFBlocks_global bookkeeping.
  !
  ! Conventions
  ! ------------
  !  *) The composite order is J.T acting on the ket: the kernels compute
  !     < wf2 | J_mu T | wf1 > with the T-fold baked into their per-
  !     component sign tables.
  !  *) dwf1 is the complete gradient of wf1 only, exactly as the kernels
  !     receive it in the legacy path; stored derivatives (HFdpsi) are used
  !     when available, otherwise the Derive_X/Y/Z_spwf primitives with the
  !     sx/sy/sz per-component reflection signs - the serial branch of
  !     Transfer_derpsi_complete.
  !==============================================================================
  use compilation,       only: dp
  use geninfo
  use wavefunctions,     only: HFPsi, HFdpsi, sx, sy, sz, nwt
  use wavefunctions,     only: Derive_X_spwf, Derive_Y_spwf, Derive_Z_spwf
  use wavefunctions,     only: angmom_x_real, angmom_x_imag
  use wavefunctions,     only: angmom_y_real, angmom_y_imag
  use wavefunctions,     only: angmom_z_real, angmom_z_imag
  use wavefunctions,     only: HFBlocks_global
  use wavefunctions,     only: TimeReverse
  use block_matrix_module, only: sym_matrix, block_filler
  use block_matrix_module, only: HERMITIAN
  use block_matrix_module, only: STY_EVEN, STY_ODD

  implicit none

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  ! Component labels of the generators: 1 = x, 2 = y, 3 = z (the honest trio).
  integer, parameter :: AX_X = 1
  integer, parameter :: AX_Y = 2
  integer, parameter :: AX_Z = 3

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
  ! Time parity of the family (design 14.1, D-7 of 2026-11-20):
  !  T J_mu T^-1 = -J_mu; the composite inherits it
  !  (T(J.T)T^-1 = tau(J) (J.T)). Feeds the objects' tau metadata and the
  !  accessor's image-image sign flip - the both-ways check of the unit
  !  arm is the arbiter, never intuition.
  integer, parameter :: TAU_ANGFAM = -1

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
  ! State handed through the block_filler callback (design section 5): the
  !  sector bounds of the matrix being filled and the object selector.
  type :: Jm_params
    integer :: st(0:8)     ! sector bounds of the matrix under construction
    integer :: mu          ! component: AX_X / AX_Y / AX_Z
  end type Jm_params

  private
  public :: AX_X, AX_Y, AX_Z
  public :: TAU_ANGFAM
  public :: Jm_params
  public :: Jmu_spme_matrix, get_J_block

contains

  !----------------------------------------------------------------------------
  subroutine Jmu_spme_matrix(mu, me)
    !---------------------------------------------------------------------------
    ! Stored half of the matrix elements of the honest J_mu, as one
    !  symmetry-block matrix: Re for x/z (STY_EVEN, kind symmetric), Im
    !  for y (STY_ODD, kind antisymmetric), all Hermitian. Words: 1 for
    !  x/y (z-signature flipping), 0 for z; parity- and isospin-even.
    !  Symmetry-forbidden objects are not an error: they store zero
    !  matrices through empty block lists - guarding user input is the
    !  caller's business.
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !   mu : component selector, AX_X / AX_Y / AX_Z
    ! Output:
    !   me : the stored half of the matrix elements (dimensionless)
    !---------------------------------------------------------------------------
    integer, intent(in)           :: mu
    type(sym_matrix), intent(out) :: me

    call Jm_build(mu, me)
  end subroutine Jmu_spme_matrix
  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  subroutine Jm_build(mu, me)
    !---------------------------------------------------------------------------
    ! Construction core of the honest generators: resolve the physical
    !  classes of the 14.1 table, init from the global HFBlocks and
    !  fill through the block_filler callback.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !   mu     : component selector, AX_X / AX_Y / AX_Z
    ! Output:
    !   me     : the built and filled matrix
    !---------------------------------------------------------------------------
    integer, intent(in)           :: mu
    type(sym_matrix), intent(out) :: me

    type(Jm_params) :: params
    integer :: word, sty_class

    if (mu .lt. AX_X .or. mu .gt. AX_Z) then
      print *, 'angular_momentum: component selector must be AX_X, AX_Y or',  &
      &         ' AX_Z, got', mu
      call stp('Invalid component in the angular-momentum generators.')
    endif

    ! The honest trio's own physical classes (the 14.1 table).
    word      = 1
    sty_class = STY_EVEN
    if (mu .eq. AX_Y) sty_class = STY_ODD
    if (mu .eq. AX_Z) word = 0
    call me%init(HFBlocks_global, word, HERMITIAN, sty_class,                  &
    &            TAU_ANGFAM, goodqn=goodqn)

    ! Symmetry-forbidden objects (or an empty basis) store zero matrices
    !  through their empty block lists: nothing left to do.
    if (me%nblocks .eq. 0) return

    params%st     = me%st
    params%mu     = mu
    call me%fill_blocks(Jm_fill_block, params)
  end subroutine Jm_build
  !----------------------------------------------------------------------------
  subroutine Jm_fill_block(srow, scol, A, params)
    !---------------------------------------------------------------------------
    ! block_filler callback of the generators: overwrite one stored block
    !  with the matrix elements of the selected object, one angmom kernel
    !  invocation per stored pair - the verbatim arithmetic of the legacy
    !  loops of calcJ2andBelyaev_HFB.
    !
    ! The block spans the global row range st(srow) : st(srow+1)-1 and the
    !  global column range st(scol) : st(scol+1)-1 of the HF basis. Every
    !  element of the dense rectangle is computed from its own pair of
    !  spinors: A(r, c) = kernel(psi_i, psi_j, derpsi_j) with i the global
    !  row index and j the global column index.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !   srow, scol : sector words 0..7 of the block as handed by fill_blocks
    !   params     : type(Jm_params): sector bounds + object selector
    ! Output:
    !   A          : the block storage, overwritten entirely
    !---------------------------------------------------------------------------
    integer, intent(in)   :: srow, scol
    real(KIND=dp), intent(out) :: A(:, :)
    class(*), intent(in)  :: params

    integer :: si, sj, m, n, r, c, gi, gj
    real(KIND=dp) :: psi_i(mv, 4)
    real(KIND=dp) :: psi_j(mv, 4)
    real(KIND=dp) :: derpsi_j(mv, 3, 4)

    select type (params)
    type is (Jm_params)

      ! Global index ranges of the block's sectors, cross-checked against
      !  the shape of the storage handed over.
      si = params%st(srow)
      sj = params%st(scol)
      m = params%st(srow + 1) - si
      n = params%st(scol + 1) - sj
      if (m .ne. size(A, 1) .or. n .ne. size(A, 2)) then
        print *, 'Jm_fill_block: params%st inconsistent with the block.'
        call stp('Sector bounds inconsistent in Jm_fill_block.')
      endif

      do c = 1, n
        gj = sj + c - 1
        psi_j = HFPsi(:, :, gj)
        call Jm_derive(psi_j, gj, derpsi_j)
        do r = 1, m
          gi = si + r - 1
          if (gi .eq. gj) then
            psi_i = psi_j
          else
            psi_i = HFPsi(:, :, gi)
          endif
          A(r, c) = Jm_element(params%mu, psi_i, psi_j, derpsi_j)
        enddo
      enddo

    end select
  end subroutine Jm_fill_block
  !----------------------------------------------------------------------------
  function Jm_element(mu, psi_i, psi_j, derpsi_j) result(val)
    !---------------------------------------------------------------------------
    ! The matrix element < i | J_mu | j > of the honest operator: one
    !  dispatch to the angmom kernel of the right component and half -
    !  the verbatim kernel arithmetic, no reimplementation.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !   mu       : component selector
    !   psi_i    : the row spinor  (mv, 4)
    !   psi_j    : the column spinor (mv, 4)
    !   derpsi_j : the complete gradient of the column spinor (mv, 3, 4)
    ! Result:
    !   val      : the matrix element
    !---------------------------------------------------------------------------
    integer, intent(in)  :: mu
    real(KIND=dp), intent(in) :: psi_i(:, :), psi_j(:, :), derpsi_j(:, :, :)
    real(KIND=dp)        :: val

    select case (mu)
    case (AX_X)
      val = angmom_x_real(psi_i, psi_j, derpsi_j)
    case (AX_Y)
      val = angmom_y_imag(psi_i, psi_j, derpsi_j)
    case (AX_Z)
      val = angmom_z_real(psi_i, psi_j, derpsi_j)
    case default
      print *, 'Jm_element: component must be AX_X, AX_Y or AX_Z.'
      call stp('Invalid component for the honest object.')
      val = 0.0d0
    end select
  end function Jm_element
  !----------------------------------------------------------------------------
  subroutine Jm_derive(psi, glob, derpsi)
    !---------------------------------------------------------------------------
    ! Complete gradient of the stored spinor with the global index glob,
    !  exactly the serial branch of Transfer_derpsi_complete of the
    !  wavefunctions module: the stored derivatives when the run keeps
    !  them, otherwise the Derive_X/Y/Z_spwf primitives with the per-
    !  component reflection signs of the spwf.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !   psi    : the spinor, (mv, 4)
    !   glob   : its global index
    ! Output:
    !   derpsi : the complete gradient, (mv, 3, 4)
    !---------------------------------------------------------------------------
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
  end subroutine Jm_derive
  !----------------------------------------------------------------------------
  subroutine get_J_block(mu, srow, scol, A)
    !---------------------------------------------------------------------------
    ! The 15.8 retrieval contract in the J shape (design 14.3 amendment,
    !  15.9 row 2): serve any block (srow, scol) of the COMPLETE
    !  doubled-basis matrix of J_mu on demand - stored blocks from the
    !  honest generator, image-involving blocks materialised, without
    !  any object ever storing them.
    !
    ! Sector classification follows the global block counts: a populated
    !  sector holds stored states; an EMPTY sector whose T-partner
    !  IEOR(s,1) is populated holds the time-reversal images of those
    !  partner states, one image in the same order as its stored original
    !  (the stored-k <-> image-k correspondence). The block is served by
    !  the cheapest valid mechanism (design 15.8):
    !    1. both sectors stored: the honest generator's own block - a
    !       zero rectangle when the word does not pair the two sectors.
    !       The whole matrix is (re)built per call in this first version:
    !       straightforwardness first, cache at the call site when a
    !       consumer needs many blocks of one object.
    !    2. both sectors image: J_z is signature-conserving, its
    !       image-image block is A = TAU_ANGFAM * the stored counterpart
    !       block, by <Tu|O|Tv> = tau * conj(<u|O|v>) and the real
    !       single-half storage - the both-ways check of the unit arm
    !       verifies this relation against a direct honest computation
    !       before this path may serve. The signature-flipping J_x/J_y
    !       have structurally zero image-image blocks - both images
    !       carry the same z-signature, which the operator flips.
    !    3. one sector stored, one image: the honest per-pair kernels
    !       with the materialised T-image on the image side - TimeReverse
    !       is spinor-pointwise and commutes with the gradient, so the
    !       image's derivative is the reversal of the stored original's.
    !       The stored-x-image blocks ARE the fold content the retired
    !  Jt entry used to serve, packaged by the stored-k <->
    !  image-k correspondence.
    !  A block pairing two sectors that are neither stored nor image is
    !  structurally empty: A is returned unallocated.
    !  The word gate uses the EFFECTIVE masked word, init's iand rule
    !  (design 2c): parity grading bits the build has merged act
    !  diagonally also here; Rz is good in every build, so the doubled
    !  basis stays signature-graded everywhere.
    !
    !  Only $TR-marked call sites may request image-involving blocks (the
    !  caller-gates-scenarios rule): image sectors carry the T-images of
    !  the stored states exactly while the build conserves time reversal.
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !   mu        : component selector, AX_X / AX_Y / AX_Z
    !   srow, scol : sector words 0..7, rows / columns of the block
    ! Output:
    !   A         : the block, (n_row, n_col) allocated and filled; zero
    !               when the word does not pair the sectors; unallocated
    !               when the pairing is structurally empty
    !---------------------------------------------------------------------------
    integer, intent(in)          :: mu
    integer, intent(in)          :: srow, scol
    real(KIND=dp), allocatable, intent(out) :: A(:, :)

    integer :: word, word_eff, nrow, ncol, r, c, gi, gj, k
    integer :: rrow, rcol
    logical :: row_img, col_img
    real(KIND=dp) :: psi_i(mv, 4)
    real(KIND=dp) :: psi_j(mv, 4)
    real(KIND=dp) :: derpsi_j(mv, 3, 4)
    type(sym_matrix) :: M
    real(KIND=dp), allocatable :: S(:, :)

    ! Sizes: an image sector counts its stored partner's states; a sector
    !  that is neither stored nor an image is structurally empty.
    nrow = Jm_sector_count(srow)
    ncol = Jm_sector_count(scol)
    if (nrow .eq. 0 .or. ncol .eq. 0) then
      if (allocated(A)) deallocate (A)
      return
    endif

    row_img = (.not. Jm_sector_stored(srow))
    col_img = (.not. Jm_sector_stored(scol))

    word = 1
    if (mu .eq. AX_Z) word = 0
    word_eff = iand(word, iand(goodqn, 7))
    if (ieor(srow, scol) .ne. word_eff) then
      allocate (A(nrow, ncol))
      A = 0.0d0
      return
    endif

    if ((.not. row_img) .and. (.not. col_img)) then
      ! 1. Stored-stored: the honest generator's own block.
      call Jmu_spme_matrix(mu, M)
      call Jm_copy_block(M, srow, scol, A)
      return
    endif

    if (row_img .and. col_img) then
      ! 2. Image-image: zero for the signature-flipping components; the
      !    tau sign flip of the stored counterpart for J_z.
      if (mu .eq. AX_Z) then
        rrow = ieor(srow, 1)
        rcol = ieor(scol, 1)
        call Jmu_spme_matrix(mu, M)
        call Jm_copy_block(M, rrow, rcol, S)
        if (.not. allocated(S)) then
          allocate (A(nrow, ncol))
          A = 0.0d0
          return
        endif
        allocate (A(nrow, ncol))
        A = TAU_ANGFAM * S
        deallocate (S)
      else
        allocate (A(nrow, ncol))
        A = 0.0d0
      endif
      return
    endif

    ! 3. One stored sector, one image sector: honest per-pair kernels with
    !    the T-image materialised on the image side.
    allocate (A(nrow, ncol))
    do c = 1, ncol
      gj = Jm_image_index(scol, c)
      psi_j = HFPsi(:, :, gj)
      call Jm_derive(psi_j, gj, derpsi_j)
      if (col_img) then
        psi_j = TimeReverse(psi_j)
        do k = 1, 3
          derpsi_j(:, k, :) = TimeReverse(derpsi_j(:, k, :))
        enddo
      endif
      do r = 1, nrow
        gi = Jm_image_index(srow, r)
        psi_i = HFPsi(:, :, gi)
        if (row_img) psi_i = TimeReverse(psi_i)
        A(r, c) = Jm_element(mu, psi_i, psi_j, derpsi_j)
      enddo
    enddo
  end subroutine get_J_block
  !----------------------------------------------------------------------------
  function Jm_sector_stored(sector) result(stored)
    !---------------------------------------------------------------------------
    ! .true. when the sector carries stored spwfs (count > 0) in the
    !  global block count bookkeeping, .false. for empty sectors - the
    !  T-image slots in TR-conserving builds.
    !---------------------------------------------------------------------------
    integer, intent(in) :: sector
    logical :: stored

    stored = (HFBlocks_global(sector + 1) .gt. 0)
  end function Jm_sector_stored
  !----------------------------------------------------------------------------
  function Jm_sector_count(sector) result(n)
    !---------------------------------------------------------------------------
    ! Number of states of the sector in the doubled-basis view: the stored
    !  count for stored sectors, the count of the T-partner sector for
    !  image sectors, zero when the sector and its partner are both empty.
    !---------------------------------------------------------------------------
    integer, intent(in) :: sector
    integer :: n

    if (Jm_sector_stored(sector)) then
      n = HFBlocks_global(sector + 1)
    else
      n = HFBlocks_global(ieor(sector, 1) + 1)
    endif
  end function Jm_sector_count
  !----------------------------------------------------------------------------
  function Jm_image_index(sector, k) result(glob)
    !---------------------------------------------------------------------------
    ! The global spwf index of the k-th state of the sector in the
    !  doubled-basis view: its own index for stored sectors; the index of
    !  the k-th state of the T-partner sector for image sectors (the
    !  stored-k <-> image-k correspondence).
    !---------------------------------------------------------------------------
    integer, intent(in) :: sector, k
    integer :: glob
    integer :: partner

    if (Jm_sector_stored(sector)) then
      glob = sum(HFBlocks_global(1:sector)) + k
    else
      partner = ieor(sector, 1)
      glob = sum(HFBlocks_global(1:partner)) + k
    endif
  end function Jm_image_index
  !----------------------------------------------------------------------------
  subroutine Jm_copy_block(M, srow, scol, A)
    !---------------------------------------------------------------------------
    ! Copy the block (srow, scol) - or its kind-weighted partner - out of
    !  the sym_matrix M into a freshly allocated rectangle; unallocated A
    !  when no part of M covers the requested pair.
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !   M         : the matrix whose block list is searched
    !   srow, scol: sector words 0..7 of the requested block
    ! Output:
    !   A         : the copied block, allocated when found
    !---------------------------------------------------------------------------
    type(sym_matrix), intent(in) :: M
    integer, intent(in) :: srow, scol
    real(KIND=dp), allocatable, intent(out) :: A(:, :)

    integer :: b

    do b = 1, M%nblocks
      if (M%blocks(b)%srow .eq. srow .and. M%blocks(b)%scol .eq. scol) then
        allocate (A(size(M%blocks(b)%A, 1), size(M%blocks(b)%A, 2)))
        A = M%blocks(b)%A
        return
      endif
      if (M%blocks(b)%srow .eq. scol .and. M%blocks(b)%scol .eq. srow) then
        allocate (A(size(M%blocks(b)%A, 2), size(M%blocks(b)%A, 1)))
        A = M%kind * transpose(M%blocks(b)%A)
        return
      endif
    enddo
    if (allocated(A)) deallocate (A)
  end subroutine Jm_copy_block

end module angular_momentum
