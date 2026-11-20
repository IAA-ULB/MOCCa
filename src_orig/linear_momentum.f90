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
module linear_momentum
!==============================================================================
  !
  ! Single-particle matrix elements of the gradient operators on the
  !  symmetry-block matrix type (design summary section 15).
  !
  ! Five objects, two families:
  !
  !   (a) the honest trio \nabla_x, \nabla_y, \nabla_z - the one-body
  !       operators over the stored basis, as many of them as the build
  !       can represent. All three are ANTIHERMITIAN, parity-odd and
  !       isoscalar; \nabla_x/\nabla_y flip the z-signature (word 3),
  !       \nabla_z conserves it (word 2). The stored halves follow the
  !       15.2 class table: Re for \nabla_x/\nabla_z (STY_EVEN), Im for
  !       \nabla_y (STY_ODD), kinds -1/-1/+1. In TR-conserving builds
  !       the partner sectors of honest \nabla_x/\nabla_y are empty, so
  !       they come out as structural zeros there - the odd-m multipole
  !       precedent, no special case.
  !
  !   (b) the TR-folded pair \nabla_x.T, \nabla_y.T - the matrix elements
  !       with the time-reversal images of the stored states as kets, the
  !       objects the hidden_TR branches of the legacy CompNablaMelements
  !       compute (design sections 2d and 15): T conserves parity and
  !       isospin and flips the z-signature, so the composites have word
  !       IEOR(word, 1) = 2, within reach of the stored sectors. The fold
  !       bakes in exactly as the legacy does: TimeReverse applied to
  !       the transferred derivative of the ket. Their class tags are
  !       pinned to the mirror oracle of the linear_momentum_testing
  !       module, which measures both composites SYMMETRIC:
  !            \nabla_x.T : HERMITIAN x STY_EVEN    -> kind +1
  !            \nabla_y.T : ANTIHERMITIAN x STY_ODD -> kind +1
  !       the per-component shape the J family landed with for K_x/K_y.
  !       The legacy's in-loop mirror fill contradicts the measured
  !       composite for the folded x - the angular-momentum 14.2b wart
  !       in reverse, consumer-invisible because every consumer
  !       contracts squares; the record lives in design 15.2.
  !
  ! D-7 / design 15.9 row 2 (2026-11-20): the five-object public
  !  interface collapses. The public face is the honest generator, the
  !  Nabla_half_fm packaging (its folded branch is the private twin
  !  engine, the CompCOMCorrection production lane) and the accessor;
  !  the eager twin sibling entry is private since this row.
  !
  ! The generators wrap the legacy per-pair arithmetic VERBATIM - one
  !  dv*sum per stored pair, the exact component contractions of
  !  CompNablaMelements - so that routine remains the A/B oracle of the
  !  consumer migration (the step-3 precedent). Densified block products
  !  (the design's 13.1 pattern) are deliberately NOT used: per-pair
  !  fidelity first, performance later.
  !
  ! Retrieval contract (design section 15.8): the additional matrix
  !  elements - those involving the T-image states - are carried
  !  IMPLICITLY. get_nabla_block serves any block of the complete
  !  doubled-basis matrix on demand: stored blocks from the honest
  !  generators, image-image blocks by the tau sign flip of their
  !  stored counterparts (the algebra route - armed per family by the
  !  both-ways check of the unit arm), image-row/image-column blocks by
  !  the honest kernels on materialised T-images, and everything not
  !  paired by the operator word as structural zero. block_matrix_module
  !  itself learns nothing about time reversal; the tau value lives here
  !  as a compile-time constant because the accessor is per-family.
  !
  ! Scope decisions (design section 15):
  !  - this module is free of $-markers; the folded generators exist in
  !    every build, their CALLERS gate by symmetry scenario downstream
  !    (only $TR-marked call sites ever need them).
  !  - single-node, serial semantics like the angular-momentum and
  !    multipole generators: the basis is indexed globally through
  !    the spwf arrays, matching the matrix's own HFBlocks_global
  !    bookkeeping. CompNablaMelements keeps its MPI machinery until
  !    the deletion gate.
  !  - the basis argument mirrors the legacy routine: 'HF', 'CAN' or
  !    'DEN' select HFPsi/canpsi/denpsi; the DERivative companions are
  !    taken from storage when the run keeps them, otherwise derived
  !    on the fly - the serial branch of Transfer_derpsi.
  !
  ! Conventions
  ! ------------
  !  *) MOCCa's time reversal is FERMIONIC: T**2 = -1 (design 14).
  !  *) tau, the gradient family's time parity (design 15.8): T is
  !     antilinear, so the i of p = -i(\nabla) is what makes the
  !     momentum time-odd - \nabla itself is time-EVEN, tau = +1. The
  !     composite \nabla_mu.T inherits the parent's tau:
  !     T(\nabla_mu T)T^-1 = tau (\nabla_mu T).
  !==============================================================================
  use compilation,       only: dp
  use geninfo
  use wavefunctions,     only: HFPsi, canpsi, denpsi, HFdpsi, candpsi, dendpsi
  use wavefunctions,     only: sx, sy, sz, nwt
  use wavefunctions,     only: Derive_X_spwf, Derive_Y_spwf, Derive_Z_spwf
  use wavefunctions,     only: TimeReverse
  use wavefunctions,     only: HFBlocks_global
  use block_matrix_module, only: sym_matrix, block_filler
  use block_matrix_module, only: HERMITIAN, ANTIHERMITIAN
  use block_matrix_module, only: STY_EVEN, STY_ODD

  implicit none

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Time parity of the gradient family (design 15.8):
  !  T (\nabla_mu) T^-1 = +(\nabla_mu). Compile-time constant of the
  !  accessor. Feeds the objects' tau metadata (D-7, 2026-11-20); the
  !  folded composite inherits it (T(O.T)T^-1 = tau(O) (O.T)).
  integer, parameter :: TAU_NABLAFAM = +1

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Component labels of the generators: 1 = x, 2 = y, 3 = z (the honest trio).
  integer, parameter :: AX_X = 1
  integer, parameter :: AX_Y = 2
  integer, parameter :: AX_Z = 3

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! State handed through the block_filler callback (design section 5): the
  !  sector bounds of the matrix being filled and the object selector.
  type :: Nabla_params
    integer :: st(0:8)     ! sector bounds of the matrix under construction
    integer :: mu           ! component: AX_X / AX_Y / AX_Z
    logical :: folded       ! .true.: the \nabla_mu.T composite
    integer :: basis        ! 1 = 'HF', 2 = 'CAN', 3 = 'DEN'
  end type Nabla_params

  private
  public :: AX_X, AX_Y, AX_Z
  public :: TAU_NABLAFAM
  public :: Nabla_params
  public :: Nabla_spme_matrix
  public :: Nabla_half_fm
  public :: get_nabla_block

contains

  !----------------------------------------------------------------------------
  subroutine Nabla_spme_matrix(mu, spwf_basis, me)
    !---------------------------------------------------------------------------
    ! Stored half of the matrix elements of the honest gradient operator
    !  \nabla_mu, as one sym_matrix: Re for x/z (STY_EVEN), Im for y
    !  (STY_ODD); all three are ANTIHERMITIAN, so init derives the kinds
    !  -1/-1/+1 - exactly the partner-fill signs of the legacy
    !  CompNablaMelements. Words: 3 for x/y (z-signature flipping,
    !  parity-odd), 2 for z (signature conserving, parity-odd); isoscalar.
    !  The goodqn-masked init reproduces every $PBROKEN branch structure
    !  of the legacy routine without any marker here. Symmetry-forbidden
    !  objects are not an error: they store zero matrices through empty
    !  block lists - guarding user input is the caller's business.
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !   mu        : component selector, AX_X / AX_Y / AX_Z
    !   spwf_basis: 'HF' / 'CAN' / 'DEN' - the basis of the stored MEs,
    !               case-insensitive; the legacy default 'DEN' is not
    !               assumed here, every caller states its basis.
    ! Output:
    !   me        : the stored half of the matrix elements (fm^-1)
    !---------------------------------------------------------------------------
    integer, intent(in)           :: mu
    character(len=*), intent(in) :: spwf_basis
    type(sym_matrix), intent(out) :: me

    call Nabla_build(mu, .false., spwf_basis, me)
  end subroutine Nabla_spme_matrix
  !----------------------------------------------------------------------------
  subroutine NablaT_spme_matrix(mu, spwf_basis, me)
    !---------------------------------------------------------------------------
    ! Stored half of the matrix elements of the TR-folded composite
    !  \nabla_mu.T - the honest matrix elements with the time-reversal
    !  images of the stored states as kets (design 2d): Re for x, Im for
    !  y. Word IEOR(word, 1) = 2; there is no folded \nabla_z - the
    !  honest operator serves both TR scenarios. Class tags follow the
    !  mirror oracle (module header). PRIVATE since D-7 / design 15.9
    !  row 2 (2026-11-20): demoted from the public interface behind
    !  get_nabla_block, surviving as the engine of Nabla_half_fm's
    !  folded branch (the CompCOMCorrection production lane). Consumers
    !  request fold content through the accessor; the unit arm owns the
    !  eager-reference arithmetic.
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !   mu        : component selector, AX_X / AX_Y
    !   spwf_basis: 'HF' / 'CAN' / 'DEN', see Nabla_spme_matrix
    ! Output:
    !   me        : the stored half of the composite's matrix elements
    !---------------------------------------------------------------------------
    integer, intent(in)           :: mu
    character(len=*), intent(in)  :: spwf_basis
    type(sym_matrix), intent(out) :: me

    if (mu .eq. AX_Z) then
      print *, 'NablaT_spme_matrix: no folded twin exists for \nabla_z',        &
      &         ' (word 2),'
      call stp('use Nabla_spme_matrix for \nabla_z.')
    endif

    call Nabla_build(mu, .true., spwf_basis, me)
  end subroutine NablaT_spme_matrix
  !----------------------------------------------------------------------------
  subroutine Nabla_build(mu, folded, spwf_basis, me)
    !---------------------------------------------------------------------------
    ! Shared construction core of the two generators: resolve the word
    !  and classes, init from the global HFBlocks and fill through the
    !  block_filler callback, which loops the legacy per-pair arithmetic.
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !   mu        : component selector, AX_X / AX_Y (/ AX_Z if not folded)
    !   folded    : .true. for the \nabla_mu.T composite
    !   spwf_basis: 'HF' / 'CAN' / 'DEN'
    ! Output:
    !   me        : the built and filled matrix
    !---------------------------------------------------------------------------
    integer, intent(in)           :: mu
    logical, intent(in)           :: folded
    character(len=*), intent(in)  :: spwf_basis
    type(sym_matrix), intent(out) :: me

    type(Nabla_params) :: params
    integer :: word, sty_class, ibasis

    if (mu .lt. AX_X .or. mu .gt. AX_Z) then
      print *, 'linear_momentum: component selector must be AX_X, AX_Y or',    &
      &         ' AX_Z, got', mu
      call stp('Invalid component in the gradient generators.')
    endif

    ! Class selection per the 15.2 table, oracle-arbitered. The honest
    !  kinds (-1/-1/+1) match the legacy partner fills (CompNablaMelements
    !  slots (1,1), (2,2), (3,1)). The folded tags follow the mirror
    !  oracle of linear_momentum_testing - the design's arbiter of
    !  composite classes:
    !       \nabla_x.T : HERMITIAN x STY_EVEN    -> kind +1 symmetric
    !       \nabla_y.T : ANTIHERMITIAN x STY_ODD -> kind +1 symmetric
    !  the same per-component (herm, sty) shape the angular-momentum
    !  family landed with for K_x and K_y; antilinear composites do not
    !  inherit the linear-operator argument of design 2a (the JxT/JyT
    !  14.2b note). The oracle contradiction of the first-drafted tags
    !  is the record of design 15.2.
    word = 3
    if (mu .eq. AX_Z) word = 2
    if (folded) word = ieor(word, 1)
    sty_class = STY_EVEN
    if (mu .eq. AX_Y) sty_class = STY_ODD

    if (folded .and. mu .eq. AX_X) then
      call me%init(HFBlocks_global, word, HERMITIAN, sty_class,                &
      &            TAU_NABLAFAM, folded=folded, goodqn=goodqn)
    else
      call me%init(HFBlocks_global, word, ANTIHERMITIAN, sty_class,            &
      &            TAU_NABLAFAM, folded=folded, goodqn=goodqn)
    endif

    ! Symmetry-forbidden objects (or an empty basis) store zero matrices
    !  through their empty block lists: nothing left to do.
    if (me%nblocks .eq. 0) return

    ibasis = Nabla_basis_code(spwf_basis)
    params%st     = me%st
    params%mu     = mu
    params%folded = folded
    params%basis  = ibasis
    call me%fill_blocks(Nabla_fill_block, params)
  end subroutine Nabla_build
  !----------------------------------------------------------------------------
  integer function Nabla_basis_code(spwf_basis) result(code)
    !---------------------------------------------------------------------------
    ! Map the legacy basis strings of CompNablaMelements onto an integer
    !  selector for the fill machinery - the same bases the wavefunctions
    !  module stores (HFPsi / canpsi / denpsi and their derivatives).
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !   spwf_basis : 'HF' / 'CAN' / 'DEN', case-insensitive
    ! Output:
    !   code       : 1 = 'HF', 2 = 'CAN', 3 = 'DEN'
    !---------------------------------------------------------------------------
    character(len=*), intent(in) :: spwf_basis

    select case (to_upper(adjustl(spwf_basis)))
    case ('HF')
      code = 1
    case ('CAN')
      code = 2
    case ('DEN')
      code = 3
    case default
      print *, 'linear_momentum: spwf_basis must be HF, CAN or DEN, got:',      &
      &         spwf_basis
      call stp('Invalid basis in the gradient generators.')
      code = 0
    end select
  end function Nabla_basis_code
  !----------------------------------------------------------------------------
  subroutine Nabla_fill_block(srow, scol, A, params)
    !---------------------------------------------------------------------------
    ! block_filler callback of the generators: overwrite one stored block
    !  with the matrix elements of the selected object, one dv*sum per
    !  stored pair - the verbatim arithmetic of the legacy
    !  CompNablaMelements loops, so the per-element oracle comparison of
    !  linear_momentum_testing can be bit-identical.
    !
    ! The block spans the global row range st(srow) : st(srow+1)-1 and the
    !  global column range st(scol) : st(scol+1)-1 of the selected basis.
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !   srow, scol : sector words 0..7 of the block as handed by fill_blocks
    !   params     : type(Nabla_params): sector bounds + object selector
    ! Output:
    !   A          : the block storage, overwritten entirely
    !---------------------------------------------------------------------------
    integer, intent(in)   :: srow, scol
    real(KIND=dp), intent(out) :: A(:, :)
    class(*), intent(in)  :: params

    integer :: si, sj, m, n, r, c, gi, gj
    real(KIND=dp) :: psi_i(mv, 4)
    real(KIND=dp) :: psi_j(mv, 4)
    real(KIND=dp) :: der_j(mv, 4)

    select type (params)
    type is (Nabla_params)

      ! Global index ranges of the block's sectors, cross-checked against
      !  the shape of the storage handed over.
      si = params%st(srow)
      sj = params%st(scol)
      m = params%st(srow + 1) - si
      n = params%st(scol + 1) - sj
      if (m .ne. size(A, 1) .or. n .ne. size(A, 2)) then
        print *, 'Nabla_fill_block: params%st inconsistent with the block.'
        call stp('Sector bounds inconsistent in Nabla_fill_block.')
      endif

      do c = 1, n
        gj = sj + c - 1
        call Nabla_basis_spinor(gj, params%basis, psi_j)
        call Nabla_basis_derive(gj, params%mu, params%basis, der_j)
        if (params%folded) der_j = TimeReverse(der_j)
        do r = 1, m
          gi = si + r - 1
          if (gi .eq. gj) then
            psi_i = psi_j
          else
            call Nabla_basis_spinor(gi, params%basis, psi_i)
          endif
          A(r, c) = Nabla_element(params%mu, psi_i, der_j)
        enddo
      enddo

    class default
      call stp('Nabla_fill_block: params is not of type Nabla_params.')
    end select
  end subroutine Nabla_fill_block
  !----------------------------------------------------------------------------
  function Nabla_element(mu, psi_i, der_j) result(val)
    !---------------------------------------------------------------------------
    ! One matrix element < i | \nabla_mu | j > for a ket whose \nabla_mu
    !  derivative der_j is already in hand (and possibly already
    !  time-reversed: the fold is applied by the caller, exactly where the
    !  legacy applies its hidden_TR): the verbatim legacy contraction,
    !       mu = AX_X / AX_Z (Re half):
    !           val = dv * sum_c der_j(:,c) psi_i(:,c)
    !       mu = AX_Y (Im half):
    !           val = dv * ( der_j(:,2) psi_i(:,1) - der_j(:,1) psi_i(:,2)
    !                       - der_j(:,3) psi_i(:,4) + der_j(:,4) psi_i(:,3) )
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !   mu     : component selector
    !   psi_i  : the row spinor (mv, 4)
    !   der_j  : the (possibly T-folded) mu-derivative of the ket (mv, 4)
    ! Result:
    !   val    : the matrix element, the stored half per the STy class
    !---------------------------------------------------------------------------
    integer, intent(in)  :: mu
    real(KIND=dp), intent(in) :: psi_i(:, :), der_j(:, :)
    real(KIND=dp)       :: val
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
  end function Nabla_element
  !----------------------------------------------------------------------------
  subroutine Nabla_basis_spinor(glob, basis, psi)
    !---------------------------------------------------------------------------
    ! The spinor with the global index glob of the selected basis.
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !   glob  : global spwf index, 1..nwt
    !   basis : 1 = 'HF', 2 = 'CAN', 3 = 'DEN'
    ! Output:
    !   psi   : the spinor, (mv, 4)
    !---------------------------------------------------------------------------
    integer, intent(in)        :: glob, basis
    real(KIND=dp), intent(out) :: psi(:, :)

    select case (basis)
    case (1)
      psi = HFPsi(:, :, glob)
    case (2)
      psi = canpsi(:, :, glob)
    case (3)
      psi = denpsi(:, :, glob)
    case default
      call stp('Nabla_basis_spinor: invalid basis selector.')
    end select
  end subroutine Nabla_basis_spinor
  !----------------------------------------------------------------------------
  subroutine Nabla_basis_derive(glob, mu, basis, derpsi)
    !---------------------------------------------------------------------------
    ! The mu-derivative of the spinor with the global index glob of the
    !  selected basis: the stored derivatives when the run keeps them,
    !  otherwise the Derive_X/Y/Z_spwf primitives with the per-component
    !  reflection signs of the spwf - the serial branch of the legacy
    !  Transfer_derpsi, whose pointer association of the basis arrays
    !  is made here.
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !   glob  : global spwf index, 1..nwt
    !   mu    : direction of the derivative, AX_X / AX_Y / AX_Z
    !   basis : 1 = 'HF', 2 = 'CAN', 3 = 'DEN'
    ! Output:
    !   derpsi : the mu-derivative of the spinor, (mv, 4)
    !---------------------------------------------------------------------------
    integer, intent(in)        :: glob, mu, basis
    real(KIND=dp), intent(out) :: derpsi(:, :)

    real(KIND=dp), pointer :: dpsis(:, :, :, :)
    real(KIND=dp), pointer :: psis(:, :, :)

    nullify (dpsis, psis)
    select case (basis)
    case (1)
      psis  => HFPsi
      dpsis => HFdpsi
    case (2)
      psis  => canpsi
      dpsis => candpsi
    case (3)
      psis  => denpsi
      dpsis => dendpsi
    case default
      call stp('Nabla_basis_derive: invalid basis selector.')
    end select

    if (store_derivatives) then
      derpsi = dpsis(:, mu, :, glob)
    else
      select case (mu)
      case (AX_X)
        call Derive_X_spwf(psis(:, :, glob), sx(:, glob), derpsi)
      case (AX_Y)
        call Derive_Y_spwf(psis(:, :, glob), sy(:, glob), derpsi)
      case (AX_Z)
        call Derive_Z_spwf(psis(:, :, glob), sz(:, glob), derpsi)
      case default
        call stp('Nabla_basis_derive: invalid component selector.')
      end select
    endif
  end subroutine Nabla_basis_derive
  !----------------------------------------------------------------------------
  function Nabla_half_fm(mu, folded, spwf_basis) result(me)
    !---------------------------------------------------------------------------
    ! Dense legacy shape of the stored half of the selected object, the
    !  Qlm_spme_fm precedent: the partner-reconstructing dense expansion of
    !  the sym_matrix built by the generators, real(nwt, nwt).
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !   mu        : component selector, AX_X / AX_Y (/ AX_Z if not folded)
    !   folded    : .true. for the \nabla_mu.T composite
    !   spwf_basis: 'HF' / 'CAN' / 'DEN'
    ! Output:
    !   me        : the dense stored half, real(nwt, nwt)
    !---------------------------------------------------------------------------
    integer, intent(in)          :: mu
    logical, intent(in)          :: folded
    character(len=*), intent(in) :: spwf_basis
    real(KIND=dp), allocatable   :: me(:, :)

    type(sym_matrix) :: M

    if (folded) then
      call NablaT_spme_matrix(mu, spwf_basis, M)
    else
      call Nabla_spme_matrix(mu, spwf_basis, M)
    endif
    me = M%dense()
  end function Nabla_half_fm
  !----------------------------------------------------------------------------
  subroutine get_nabla_block(mu, srow, scol, spwf_basis, A)
    !---------------------------------------------------------------------------
    ! The 15.8 retrieval contract: serve any block (srow, scol) of the
    !  COMPLETE doubled-basis matrix of \nabla_mu on demand - stored blocks
    !  from the honest generators, image-involving blocks materialised,
    !  without any object ever storing them.
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
    !    2. both sectors image, signature-conserving operator (\nabla_z,
    !       word 2): A = TAU_NABLAFAM * the stored counterpart block, by
    !       <Tu|O|Tv> = tau <u|O|v>* and the real single-half storage.
    !       The sign flip is armed per family: the unit arm's both-ways
    !       check verifies this relation against a direct honest
    !       computation before this path may serve.
    !       For the signature-flipping operators (\nabla_x/\nabla_y,
    !       word 3) the image-image blocks are STRUCTURALLY ZERO - both
    !       images carry the same z-signature, which the operator flips -
    !       so no counterpart exists and none is needed.
    !    3. one sector stored, one image: the honest per-pair kernels
    !       with the materialised T-image on the image side - bra-side
    !       time reversal for image rows, the legacy hidden_TR fold on
    !       the ket derivative for image columns, exactly the arithmetic
    !       the folded generators store columnwise.
    !  A block pairing two sectors that are neither stored nor image is
    !  structurally empty: A is returned unallocated.
    !
    !  Only $TR-marked call sites may request image-involving blocks (the
    !  caller-gates-scenarios rule): image sectors carry the T-images of
    !  the stored states exactly while the build conserves time reversal.
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !   mu        : component selector, AX_X / AX_Y / AX_Z
    !   srow, scol : sector words 0..7, rows / columns of the block
    !   spwf_basis: 'HF' / 'CAN' / 'DEN'
    ! Output:
    !   A         : the block, (n_row, n_col) allocated and filled; zero
    !               when the word does not pair the sectors; unallocated
    !               when the pairing is structurally empty
    !---------------------------------------------------------------------------
    integer, intent(in)          :: mu
    integer, intent(in)          :: srow, scol
    character(len=*), intent(in) :: spwf_basis
    real(KIND=dp), allocatable, intent(out) :: A(:, :)

    integer :: word, nrow, ncol, r, c, gi, gj, basis
    integer :: rrow, rcol
    logical :: row_img, col_img
    real(KIND=dp) :: der(mv, 4)
    real(KIND=dp) :: bra(mv, 4)
    type(sym_matrix) :: M
    real(KIND=dp), allocatable :: S(:, :)

    basis = Nabla_basis_code(spwf_basis)

    ! Sizes: an image sector counts its stored partner's states; a sector
    !  that is neither stored nor an image is structurally empty.
    nrow = Nabla_sector_count(srow)
    ncol = Nabla_sector_count(scol)
    if (nrow .eq. 0 .or. ncol .eq. 0) then
      if (allocated(A)) deallocate (A)
      return
    endif

    row_img = (.not. Nabla_sector_stored(srow))
    col_img = (.not. Nabla_sector_stored(scol))

    ! The effective pairing word of the honest operator - init's iand
    !  rule (design 2c): parity grading bits the build has merged act
    !  diagonally also here, and the served blocks match the generator's
    !  own masked init. Rz is good in every build, so the doubled basis
    !  carries definite signature everywhere (stored states and images
    !  alike). In `$NTR` builds no sector is an image (T is not
    !  conserved) and only mechanism 1 is reachable.
    word = 3
    if (mu .eq. AX_Z) word = 2
    word = iand(word, iand(goodqn, 7))
    if (ieor(srow, scol) .ne. word) then
      allocate (A(nrow, ncol))
      A = 0.0d0
      return
    endif

    if ((.not. row_img) .and. (.not. col_img)) then
      ! 1. Stored-stored: the honest generator's own block.
      call Nabla_spme_matrix(mu, spwf_basis, M)
      call Nabla_copy_block(M, srow, scol, A)
      return
    endif

    if (row_img .and. col_img) then
      ! 2. Image-image: zero for the signature-flipping operators; the
      !    tau sign flip of the stored counterpart for \nabla_z.
      if (mu .eq. AX_Z) then
        rrow = ieor(srow, 1)
        rcol = ieor(scol, 1)
        call Nabla_spme_matrix(mu, spwf_basis, M)
        call Nabla_copy_block(M, rrow, rcol, S)
        if (.not. allocated(S)) then
          allocate (A(nrow, ncol))
          A = 0.0d0
          return
        endif
        allocate (A(nrow, ncol))
        A = TAU_NABLAFAM * S
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
      gj = Nabla_image_index(scol, c)
      call Nabla_basis_derive(gj, mu, basis, der)
      if (col_img) then
        ! The ket is the T-image: the legacy hidden_TR convention folds
        !  the derivative of the stored original.
        der = TimeReverse(der)
      endif
      do r = 1, nrow
        gi = Nabla_image_index(srow, r)
        call Nabla_basis_spinor(gi, basis, bra)
        if (row_img) then
          ! The bra is the T-image of the stored original.
          bra = TimeReverse(bra)
        endif
        A(r, c) = Nabla_element(mu, bra, der)
      enddo
    enddo
  end subroutine get_nabla_block
  !----------------------------------------------------------------------------
  function Nabla_sector_stored(sector) result(stored)
    !---------------------------------------------------------------------------
    ! .true. when the sector carries stored spwfs (count > 0) in the
    !  global block count bookkeeping, .false. for empty sectors - the
    !  T-image slots in TR-conserving builds.
    !---------------------------------------------------------------------------
    integer, intent(in) :: sector
    logical :: stored

    stored = (HFBlocks_global(sector + 1) .gt. 0)
  end function Nabla_sector_stored
  !----------------------------------------------------------------------------
  function Nabla_sector_count(sector) result(n)
    !---------------------------------------------------------------------------
    ! Number of states of the sector in the doubled-basis view: the stored
    !  count for stored sectors, the count of the T-partner sector for
    !  image sectors, zero when the sector and its partner are both empty.
    !---------------------------------------------------------------------------
    integer, intent(in) :: sector
    integer :: n

    if (Nabla_sector_stored(sector)) then
      n = HFBlocks_global(sector + 1)
    else
      n = HFBlocks_global(ieor(sector, 1) + 1)
    endif
  end function Nabla_sector_count
  !----------------------------------------------------------------------------
  function Nabla_image_index(sector, k) result(glob)
    !---------------------------------------------------------------------------
    ! The global spwf index of the k-th state of the sector in the
    !  doubled-basis view: its own index for stored sectors; the index of
    !  the k-th state of the T-partner sector for image sectors (the
    !  stored-k <-> image-k correspondence).
    !---------------------------------------------------------------------------
    integer, intent(in) :: sector, k
    integer :: glob
    integer :: partner

    if (Nabla_sector_stored(sector)) then
      glob = sum(HFBlocks_global(1:sector)) + k
    else
      partner = ieor(sector, 1)
      glob = sum(HFBlocks_global(1:partner)) + k
    endif
  end function Nabla_image_index
  !----------------------------------------------------------------------------
  subroutine Nabla_copy_block(M, srow, scol, A)
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
  end subroutine Nabla_copy_block
  !----------------------------------------------------------------------------

end module linear_momentum
