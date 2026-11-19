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
  ! Five objects, two families:
  !
  !   (a) the honest trio  J_x, J_y, J_z  - the one-body operators over the
  !       stored basis. J_x/J_y flip the z-signature (word 1), J_z conserves
  !       it (word 0); all three are Hermitian, parity-even and isoscalar;
  !       J_y is the STY_ODD member of the triad (design section 2a). In
  !       TR-conserving builds the word-1 partners of the stored sectors are
  !       empty, so honest J_x/J_y come out as structural zeros there - the
  !       odd-m multipole precedent, no special case.
  !
  !   (b) the TR-folded pair  J_x.T, J_y.T  - the antilinear composites the
  !       t-kernels of wavefunctions.f90 implement, with MOCCa's FERMIONIC
  !       time reversal T**2 = -1. Both preserve the z-signature (T flips
  !       it, J_mu flips it back): word 0, within the stored sectors. The
  !       honest word-1 operators being structurally empty in TR-conserving
  !       builds is exactly why those builds' consumers need the folded
  !       objects. The composites are Hermitian as operators,
  !            (J.T)^dagger = T^-1 J = -(T**2) (J.T) = +(J.T)  for T**2 = -1 ,
  !       an assertion recorded with the user 2026-11-19 - NOT a corollary
  !       of the linear-operator algebra of design section 2a. Their class
  !       tags below are therefore pinned to the empirical arbiter of
  !       design section 14.3: the kernel-transpose oracle of the
  !       angular_momentum_testing module measures the mirror character of
  !       the stored numbers directly, and its printout is the authority
  !       the tags have to agree with.
  !
  ! The generators wrap the angmom kernels of wavefunctions.f90 VERBATIM -
  !  one kernel invocation per stored pair, exactly the arithmetic of the
  !  legacy loops of calcJ2andBelyaev_HFB, so those loops remain the A/B
  !  oracle of the consumer migration (the step-3 precedent).
  !
  ! Scope decisions (design 14 and 14.3):
  !  - this module is free of $-markers; the folded generators exist in
  !    every build, their CALLERS gate by symmetry scenario downstream.
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
  use wavefunctions,     only: angmom_xt_real
  use wavefunctions,     only: angmom_yt_imag
  use wavefunctions,     only: HFBlocks_global
  use block_matrix_module, only: sym_matrix, block_filler
  use block_matrix_module, only: HERMITIAN, ANTIHERMITIAN
  use block_matrix_module, only: STY_EVEN, STY_ODD

  implicit none

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Component labels of the generators: 1 = x, 2 = y, 3 = z (the honest trio).
  integer, parameter :: AX_X = 1
  integer, parameter :: AX_Y = 2
  integer, parameter :: AX_Z = 3

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! State handed through the block_filler callback (design section 5): the
  !  sector bounds of the matrix being filled and the object selector.
  type :: Jm_params
    integer :: st(0:8)     ! sector bounds of the matrix under construction
    integer :: mu          ! component: AX_X / AX_Y / AX_Z
    logical :: folded      ! .true.: the J_mu.T composite
  end type Jm_params

  private
  public :: AX_X, AX_Y, AX_Z
  public :: Jm_params
  public :: Jmu_spme_matrix, Jt_spme_matrix

contains

  !----------------------------------------------------------------------------
  subroutine Jmu_spme_matrix(mu, me)
    !---------------------------------------------------------------------------
    ! Build the stored half of the single-particle matrix elements of the
    !  honest angular-momentum operator J_mu in the HF basis, as one
    !  symmetry-block matrix:
    !       mu = AX_X :  me = Re < i | J_x | j >   - STY_EVEN, kind symmetric
    !       mu = AX_Y :  me = Im < i | J_y | j >   - STY_ODD, kind antisymm.
    !       mu = AX_Z :  me = Re < i | J_z | j >   - STY_EVEN, kind symmetric
    !  All three operators are Hermitian (design section 14.1); init derives
    !  the storage kinds from the classes. In a build that conserves the time
    !  simplex the stored half is the complete operator.
    !
    !  The symmetry words are physical: J_x/J_y flip the z-signature (word 1),
    !  J_z conserves it (word 0); parity and isospin bits are clear. The
    !  good-quantum-number mask of geninfo merges the words with the basis as
    !  everywhere in the design.
    !
    !  Symmetry-forbidden combinations are not an error here: they store zero
    !  matrices through empty block lists (e.g. honest J_x/J_y in a
    !  TR-conserving build) - guarding user input is the caller's business.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !   mu : component selector, AX_X / AX_Y / AX_Z
    ! Output:
    !   me : the stored half of the matrix elements (fm**0 - dimensionless)
    !---------------------------------------------------------------------------
    integer, intent(in)           :: mu
    type(sym_matrix), intent(out) :: me

    call Jm_build(mu, .false., me)
  end subroutine Jmu_spme_matrix
  !----------------------------------------------------------------------------
  subroutine Jt_spme_matrix(mu, me)
    !---------------------------------------------------------------------------
    ! Build the stored half of the single-particle matrix elements of the
    !  TR-folded composite J_mu.T in the HF basis, as one symmetry-block
    !  matrix - the object the t-kernels of wavefunctions.f90 implement:
    !       mu = AX_X :  me = Re < i | J_x T | j >
    !       mu = AX_Y :  me = Im < i | J_y T | j >
    !  (there is no folded J_z: J_z preserves the z-signature, so the honest
    !  operator serves both TR-conserving and TR-broken builds.)
    !
    !  The composites are Hermitian as operators (design section 14.1, the
    !  fermionic T**2 = -1 chain); their storage kinds are pinned to the
    !  kernel-transpose oracle of the angular_momentum_testing module - the
    !  class tags here are the oracle's measured facts, re-verified by the
    !  unit arm at every run.
    !
    !  Word 0, within the stored sectors: only $TR-marked consumers have a
    !  use for these objects (design section 14.3's caller-gates-scenarios
    !  rule) - the generators themselves exist in every build.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !   mu : component selector, AX_X / AX_Y
    ! Output:
    !   me : the stored half of the composite's matrix elements
    !---------------------------------------------------------------------------
    integer, intent(in)           :: mu
    type(sym_matrix), intent(out) :: me

    if (mu .eq. AX_Z) then
      print *, 'Jt_spme_matrix: no folded twin exists for J_z (word 0),'
      call stp('use Jmu_spme_matrix for J_z.')
    endif

    call Jm_build(mu, .true., me)
  end subroutine Jt_spme_matrix
  !----------------------------------------------------------------------------
  subroutine Jm_build(mu, folded, me)
    !---------------------------------------------------------------------------
    ! Shared construction core of the two generators: resolve the physical
    !  classes of the selected object, init the block list from the global
    !  HFBlocks bookkeeping and fill every stored block through the
    !  block_filler callback, which loops the angmom kernels over the
    !  stored pairs.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !   mu     : component selector, AX_X / AX_Y (/ AX_Z if not folded)
    !   folded : .true. for the J_mu.T composite
    ! Output:
    !   me     : the built and filled matrix
    !---------------------------------------------------------------------------
    integer, intent(in)           :: mu
    logical, intent(in)           :: folded
    type(sym_matrix), intent(out) :: me

    type(Jm_params) :: params
    integer :: word, sty_class

    if (mu .lt. AX_X .or. mu .gt. AX_Z) then
      print *, 'angular_momentum: component selector must be AX_X, AX_Y or',  &
      &         ' AX_Z, got', mu
      call stp('Invalid component in the angular-momentum generators.')
    endif

    ! Physical classes per the 14.1 table. The honest operators are
    !  Hermitian, their STy classes matching the stored halves the
    !  production kernels compute (J_x/J_z real, J_y imaginary - the
    !  legacy storage conventions of calcJ2andBelyaev_HFB). For the
    !  folded composites the (herm, sty) pair does not carry its
    !  linear-operator meaning (design sections 2a and 14): sty is the
    !  STORAGE selector - which half the t-kernel computes, real for
    !  J_x.T, imaginary for J_y.T - and herm is the kind-completing
    !  sign, herm = kind * sty, with kind the mirror character the
    !  kernel-transpose oracle MEASURES: symmetric, +1, for BOTH folded
    !  objects (angular_momentum_testing, 2026-11-19: max |a_ij - a_ji|
    !  ~ 1e-14 for both components at the O16/SLy4 unit state). This
    !  matches the production fills of calcJ2andBelyaev_HFB exactly.
    if (folded) then
      word = 0
      if (mu .eq. AX_Y) then
        sty_class = STY_ODD
        call me%init(HFBlocks_global, word, ANTIHERMITIAN, sty_class,          &
        &            goodqn=goodqn)
      else
        sty_class = STY_EVEN
        call me%init(HFBlocks_global, word, HERMITIAN, sty_class,              &
        &            goodqn=goodqn)
      endif
    else
      word      = 1
      sty_class = STY_EVEN
      if (mu .eq. AX_Y) sty_class = STY_ODD
      if (mu .eq. AX_Z) word = 0
      call me%init(HFBlocks_global, word, HERMITIAN, sty_class,                &
      &            goodqn=goodqn)
    endif

    ! Symmetry-forbidden objects (or an empty basis) store zero matrices
    !  through their empty block lists: nothing left to do.
    if (me%nblocks .eq. 0) return

    params%st     = me%st
    params%mu     = mu
    params%folded = folded
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
          A(r, c) = Jm_element(params%mu, params%folded,                     &
          &                   psi_i, psi_j, derpsi_j)
        enddo
      enddo

    end select
  end subroutine Jm_fill_block
  !----------------------------------------------------------------------------
  function Jm_element(mu, folded, psi_i, psi_j, derpsi_j) result(val)
    !---------------------------------------------------------------------------
    ! The matrix element < i | J_mu(.T) | j > of the selected object: one
    !  dispatch to the angmom kernel of the right component and half -
    !  the verbatim kernel arithmetic, no reimplementation.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !   mu       : component selector
    !   folded   : .true. for the J_mu.T composite
    !   psi_i    : the row spinor  (mv, 4)
    !   psi_j    : the column spinor (mv, 4)
    !   derpsi_j : the complete gradient of the column spinor (mv, 3, 4)
    ! Result:
    !   val      : the matrix element
    !---------------------------------------------------------------------------
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
        print *, 'Jm_element: the folded family covers AX_X and AX_Y.'
        call stp('Invalid component for the folded object.')
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
        print *, 'Jm_element: component must be AX_X, AX_Y or AX_Z.'
        call stp('Invalid component for the honest object.')
        val = 0.0d0
      end select
    endif
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

end module angular_momentum
