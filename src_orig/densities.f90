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
module densities
!===============================================================================
! Module that defines, calculates and generally deals with all densities. 
! 
!=============================================================================== 
! Hephaestos keywords
!
! DECLARATION                      : [WAY too long to include here]
! INITIALIZATION                   : [WAY too long to include here]
! ZEROING                          : [WAY too long to include here]
! 
! EXPRESSION                       : [WAY too long to include here]
! EXPRESSION_OFFDIAG_SYMMETRIC     : [WAY too long to include here]
! EXPRESSION_OFFDIAG_ANTISYMMETRIC : [WAY too long to include here]
! 
! BCSEXPRESSION                    : [WAY too long to include here]
! HFBEXPRESSION                    : [WAY too long to include here]
! DERIVATION                       : [WAY too long to include here]
!
! ISOSPINCOUPL                     : [WAY too long to include here]
! ISOSPINCOUPL_SYMMETRIC           : [WAY too long to include here]
! ISOSPINCOUPL_ANTISYMMETRIC       : [WAY too long to include here]
! MPIDEN                           : [WAY too long to include here]
!
! TR              : $TR
! NTR             : $NTR 
! 
! SX_RHO          : $SX_RHO 
! SY_RHO          : $SY_RHO
! SZ_RHO          : $SZ_RHO
! 
! SX_RHO_ANTISYM  : $SX_RHO_ANTISYM 
! SY_RHO_ANTISYM  : $SY_RHO_ANTISYM
! SZ_RHO_ANTISYM  : $SZ_RHO_ANTISYM
!
! SX_SX/SY/SZ     : $SX_SX, $SX_SY, $SX_SZ
! SY_SX/SY/SZ     : $SY_SX, $SY_SY, $SY_SZ
! SY_SX/SY/SZ     : $SZ_SX, $SZ_SY, $SZ_SZ
!
! PBROKEN         : $PBROKEN
!
! WRITEDENSITIES_HDF5 : [WAY too long to include here]
!
! DISABLED KEYWORD, still present in Hephaestos
! CLEANING        : [WAY too long to include here]
!===============================================================================
!
! A density D_L_R is stored as
!
!      D_L_R (mv, [cartesian indices], [isospin indices])
!             |     |                            |
!             > spatial indices                  |
!                   |                            |
!                   > all cartesian indices      |
!                                                > isospin indices
!                                                  4 for normal densities
!                                                     (n, p, 0, 1)
!                                                  2 for pairing densities
!                                                     (n,p)
!
!=============================================================================== 
! Some technical notes:
!
! a) The densities are represented as vectors on the mesh, instead of 3D arrays. 
!    Advantages: 
!      *) Two less indices, which is needed to naively store densities
!         and their derivatives for N2/3LO. Fortran 90 was limited to 
!         arrays of rank 7, while Fortran 2008 lifted this limit to 15
!         many compilers (i.e. gfortran) do not yet support it.
!      *) Less involved coding in the functional.f90 file, since there
!         are less ':' to be put in the sums, and all of the spatial
!         indices are combined into one.
!      *) It MIGHT result in faster summing of densities, but this has never 
!         been a problem on the mean-field level.
!
!    Disadvantages: 
!      *) This needed a wrapper routine for the derivative functions Derive_grad
!         and derive_lap. Currently everything is accomplished by defining 
!         3D pointers, which is hopefully more efficient than a call to RESHAPE.
!         If this ever takes up a significant fraction of computation time, one
!         can think of explicitly writing the 1D derivative code.
!
! b) Note that the storage scheme for derivatives is not yet implemented on the 
!    level of densities, only on the level of derivatives of densities.
!    Thus
!           D_N_N is fully stored with indices (nx*ny*nz,3,3,2)
!    But 
!           Der_Der_D_I_I is stored as (nx*ny*nz,7,2)
!
!===============================================================================
! NOTE: the use statements below carry only:-lists. gfortran (13) refuses to
! attach a submodule to a parent whose .mod re-exports procedure
! characteristics from bare `use` statements - see the interface bodies below
! for which submodules these are.
use compilation, only: dp
use geninfo, only: mv, nx, ny, nz, mpi_rank
use wavefunctions,                            &
               only: nwt, nwn, nwt_local, hfpsi,               &
                     hftransfo, canpsi, candpsi, canddpsi,     &
                     candddpsi, hfdpsi, hfddpsi, hfdddpsi,     &
                     DenPsi, DendPsi, DenDDPsi, DenDPsi,       &
                     DendddPsi, diagsphamil, Pauli,            &
                     ImagMultiplySpinor, hfblocks,             &
                     hfblocks_global, rank_map, spwf_inverse,  &
                     spwf_map, momentum_updates, derivecan,    &
                     efficienthfb, canenergies, sphamil,     &
                     sx, sy, sz,                            &
                     transform_spwfs_inplace, DeriveHF,     &
                     Transfer_psi, Transfer_derpsi
use pairing, only: pairingtype, rho_pairing, hfbgaps, cantransfo, cancuttransfo, bogoliubov, Pcutoffs, kappa_can
use HFB, only: Canonical
use derivatives
use preconditioning
use basis_transform, only: transform_mat, transform_bogo, transform_diag, transform_vec, transform_spwfs
use timing

use vectors, only: DensityVector, PotentialVector, memory, construct_charge_density

implicit none

    interface operator (+)
      !Overloading "+" to be used to add density vectors.
      module procedure Add_densityvector
    end interface
    
    interface operator (*)
      !Overloading "*" to be used to multiply density vectors with scalars
      module procedure multiply_densityvector
    end interface
    !---------------------------------------------------------------------------
    ! Different density-vectors that can be used throughout the code
    type(DensityVector), target :: Density
    type(DensityVector), target :: Density_out
    type(DensityVector), target, allocatable :: DensityHistory(:)
    !---------------------------------------------------------------------------
    ! As several other modules deal with the density D_I_I and its derivatives
    ! in various forms,  Hephaestos fills in here the appropriate symmetries.
    integer, parameter :: sx_rho = $SX_RHO, sx_rho_antisym = $SX_RHO_ANTISYM
    integer, parameter :: sy_rho = $SY_RHO, sy_rho_antisym = $SY_RHO_ANTISYM
    integer, parameter :: sz_rho = $SZ_RHO, sz_rho_antisym = $SZ_RHO_ANTISYM
    ! and similar for the vector spin density s, which is needed in the 
    ! preconditioning of the functionals
    integer, parameter :: sx_s(3) = (/$SX_SX,$SX_SY,$SX_SZ/)
    integer, parameter :: sy_s(3) = (/$SY_SX,$SY_SY,$SY_SZ/)
    integer, parameter :: sz_s(3) = (/$SZ_SX,$SZ_SY,$SZ_SZ/)
    
    interface mixup_rhokappa
      module procedure mixup_rhokappa_real
      module procedure mixup_rhokappa_complex
    end interface
    
  interface
    !-------------------------------------------------------------------------
    ! The following procedures are implemented in submodules (files
    ! src_orig/densities_{diag,offdiag,perturbed,sphamil}.f90).
    ! Only their interfaces live here.
    !-------------------------------------------------------------------------
    module function densit(rho, kappa) result(R)
      real(KIND=dp), intent(in) :: rho(:), kappa(:,:)
      type(DensityVector)       :: R
    end function densit

    module subroutine densit_offdiag(rho, kappa_plus, kappa_minus, Rs, Ra, R_pp_plus, R_pp_minus)
      complex(KIND=dp), intent(in)     :: rho(:,:), kappa_plus(:,:), kappa_minus(:,:)
      type(DensityVector), intent(out) :: Rs, Ra, R_pp_plus, R_pp_minus
    end subroutine densit_offdiag

    module function densit_offdiag_restricted(rho, kappa) result(R)
      real(KIND=dp), intent(in) :: rho(:,:), kappa(:,:)
      type(DensityVector)       :: R
    end function densit_offdiag_restricted

    module function densit_offdiag_pp(kappa) result(R_total)
      COMPLEX(KIND=dp), intent(in)       :: kappa(:,:)
      type(DensityVector)                :: R_total
    end function densit_offdiag_pp

    module function densit_offdiag_ph_symmetric(rho) result(R_total)
      complex(KIND=dp), intent(in) :: rho(:,:)
      type(DensityVector)          :: R_total
    end function densit_offdiag_ph_symmetric

    module function densit_offdiag_ph_antisymmetric(rho) result(R_total)
      complex(KIND=dp), intent(in) :: rho(:,:)
      type(DensityVector)          :: R_total
    end function densit_offdiag_ph_antisymmetric

    module function calc_sphamil_me(denpsi, dendpsi, denddpsi, Fs, Fa, onthefly) result(sphamil_me)
      real(KIND=dp), intent(in)         :: denpsi(:,:,:), dendpsi(:,:,:,:), denddpsi(:,:,:,:)
      logical, intent(in)               :: onthefly
      type(PotentialVector), intent(in) :: Fs, Fa
      complex(KIND=dp), allocatable     :: sphamil_me(:,:)
    end function calc_sphamil_me

    module function calc_sphamil_me_sym(denpsi, dendpsi, denddpsi, F,  onthefly) result(sphamil_me)
      real(KIND=dp), intent(in)         :: denpsi(:,:,:), dendpsi(:,:,:,:), denddpsi(:,:,:,:)
      logical, intent(in)               :: onthefly
      type(PotentialVector), intent(in) :: F
      complex(KIND=dp), allocatable     :: sphamil_me(:,:)
    end function calc_sphamil_me_sym

    module function calc_sphamil_me_antisym(denpsi, dendpsi, denddpsi, F,  onthefly) result(sphamil_me)
      real(KIND=dp), intent(in)         :: denpsi(:,:,:), dendpsi(:,:,:,:), denddpsi(:,:,:,:)
      logical, intent(in)               :: onthefly
      type(PotentialVector), intent(in) :: F
      complex(KIND=dp), allocatable     :: sphamil_me(:,:)
    end function calc_sphamil_me_antisym

    module function calc_delta_me(denpsi, dendpsi, denddpsi, F,  onthefly) result(delta_me)
      real(KIND=dp), intent(in)         :: denpsi(:,:,:), dendpsi(:,:,:,:), denddpsi(:,:,:,:)
      logical, intent(in)               :: onthefly
      type(PotentialVector), intent(in) :: F
      complex(KIND=dp), allocatable     :: delta_me(:,:)
    end function calc_delta_me

  end interface

contains
 
 subroutine save_density_history(R)
  !-----------------------------------------------------------------------------
  ! Update the history of the mean-field densities with R.
  !-----------------------------------------------------------------------------
  integer :: i
  type(DensityVector), intent(in) :: R
  
  if(.not.allocated(DensityHistory)) allocate(DensityHistory(memory))
    
  do i=1,memory-1
      DensityHistory(memory-i+1) = DensityHistory(memory-i)
  enddo
  DensityHistory(1) = R

end subroutine save_density_history

function Add_densityvector(R1, R2) result(R)
  !-----------------------------------------------------------------------------
  ! Add two density vectors together
  !-----------------------------------------------------------------------------
  type(DensityVector), intent(in) :: R1, R2
  type(DensityVector)             :: R

$INITIALIZATION
$ADD

  allocate(R%chargedensity(nx,ny,nz))
  R%chargedensity = R1%chargedensity + R2%chargedensity

end function Add_densityvector

function multiply_densityvector(a, R1) result(R)
  !-----------------------------------------------------------------------------
  ! Add multiply a density vector with a real scalar.
  !-----------------------------------------------------------------------------
  type(DensityVector), intent(in) :: R1
  real(KIND=dp), intent(in)       :: a
  type(DensityVector)             :: R
!   real(KIND=dp) :: stor ! Commented for now: required for memory estimation
!                         ! through Hephaestos

$INITIALIZATION
$MULTIPLY

  allocate(R%chargedensity(nx,ny,nz))
  R%chargedensity = a * R1%chargedensity 

end function multiply_densityvector
 
subroutine construct_canonical_basis(rho, kappa, rho_c, kappa_c)
  !-----------------------------------------------------------------------------
  ! High-level routine that constructs the canonical basis of spwfs.
  ! Does nothing for a HF or BCS calculation.
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Input :
  !  rho   : matrix elements of the density in the basis spanned by the spwfs 
  !          currently in memory
  !  kappa : matrix elements of the anomalous density in spanned by the spwfs 
  !          currently in memory
  ! Output:
  !  rho_can   : diagonal matrix elements of the density in the canonical basis
  !  kappa_can : "diagonal" matrix elements of kappa in the canonical basis
  !
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Side effects:
  !  *) The spwfs of canonical basis get constructed. Depending on whether 
  !     efficientHFB is active, these get stored in CanPsi or overwrite the 
  !     existing HFPsi array.
  !  *) Cantransfo is saved; this is the transformation matrix between 
  !     the HF-basis and the canonical basis. It is a trivial transformation if
  !     efficientHFB is active.
  !-----------------------------------------------------------------------------
  real(KIND=dp), intent(inout) :: rho(:,:), kappa(:,:)
  real(KIND=dp), intent(out)   :: rho_c(:), kappa_c(:)
  integer :: ifail, wave

  if(pairingtype.ne.2) return

  call start_timer(T_den_can)

    if(.not.allocated(Canenergies)) allocate(Canenergies(nwt)) 

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! a. construct the transformation CanTransfo that brings us into the 
    !    canonical basis by diagonalizing rho
    print *, 'LOL', size(rho,1), size(kappa,1), size(rho_c), size(kappa_c), size(cantransfo), size(cancuttransfo)
    call Canonical(rho, kappa, rho_c, kappa_c, cantransfo,cancuttransfo,ifail)
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! b. use this transformation to construct the physical wavefunctions
    if(efficientHFB) then
       ! Saving memory
       call transform_spwfs_inplace(hfpsi, cantransfo)
       if(allocated(momentum_updates)) then
         call transform_spwfs_inplace(momentum_updates, cantransfo)
       endif
    else
       ! Full-on transformation, we have memory to burn
       call transform_spwfs(hfpsi, canpsi, cantransfo)
    endif
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! c. Transform all relevant matrices into the new basis
    if(efficientHFB) then
      rho_pairing   = transform_mat(rho, cantransfo)
      ! Kappa transforms differently from rho, but in case of real 
      ! matrices this is largely irrelevant
      kappa = transform_mat(kappa, cantransfo)

      sphamil = transform_mat(sphamil, cantransfo)
      do wave=1,nwt
        canenergies(wave) = sphamil(wave,wave)
      enddo 
      HFBgaps    = transform_mat(HFBgaps, cantransfo)

      HFtransfo   = transform_vec(HFtransfo  , cantransfo)
      Bogoliubov  = transform_bogo(Bogoliubov, cantransfo)

      ! The new basis IS the canonical basis, hence the canonical transformation
      ! is trivial.
      cantransfo = 0.0d0
      do wave=1,nwt
        cantransfo(wave,wave) = 1.0d0
      enddo 
    else
      canenergies = transform_diag(sphamil, cantransfo) 
    endif
    call stop_timer(T_den_can)

    ! If we are not doing HFB efficiently, we should rederive the spwfs
    if((.not. efficientHFB) .and. store_derivatives) call derivecan()

 end subroutine construct_canonical_basis

 function gen_unitary_transform() result(transfo)
  !-----------------------------------------------------------------------------
  ! Generate a random unitary transformation in the space of the single-particle
  ! states that respects the symmetries of the calculation.
  !
  ! Input
  ! Output:
  !  - transfo : the unitary transformation
  !-----------------------------------------------------------------------------
  integer                    :: si, B, i, j, N
  real(KIND=dp), allocatable :: transfo(:,:)
!  real(KIND=dp), allocatable :: check(:,:)
  real(KIND=dp)              :: fac

  allocate(transfo(nwt,nwt))
  ! Generate a (symmetry-respecting) random unitary transformation
  transfo = 0.0d0
  si = 0
  do B=1,8
    N = HFBlocks(B)
    ! Generate random numbers
    call random_number(transfo(si+1:si+N, si+1:si+N))
    ! Orthonormalize the columns
    do i=1,N
      ! Normalize
      fac =  sqrt(sum(transfo(si+1:si+N,si+i)**2))
      transfo(si+1:si+N,si+i) = transfo(si+1:si+N,si+i)/ fac
      do j=i+1,N
        ! Orthogonalize the rest
        fac =  sum(transfo(si+1:si+N,si+j)*transfo(si+1:si+N,si+i))
        transfo(si+1:si+N,si+j) = transfo(si+1:si+N,si+j) &
        &                                       - fac *  transfo(si+1:si+N,si+i)
      enddo
    enddo
!     print *, 'B=', B
!     do i=1,N
!       print ('(99f10.3)'), transfo(i,1:N)
!     enddo
    si = si + N
  enddo


 end function gen_unitary_transform

subroutine mixup_rhokappa_real(rho, kappa, transfo)
  !-----------------------------------------------------------------------------
  ! Construct a random symmetry-respecting unitary transformation and apply
  ! it to (i) the single-particle wavefunctions in the Hartree-Fock basis
  ! and (ii) the matrices rho and kappa.
  !-----------------------------------------------------------------------------

  real(KIND=dp), allocatable, intent(in) :: transfo(:,:)
  real(KIND=dp), intent(inout) :: rho(:,:), kappa(:,:)

  ! Perform the transformation of the spwfs
  call transform_spwfs_inplace(hfpsi, transfo)
  ! reperform derivatives((
  call deriveHF()
  ! Transform the density matrix and canonical kappa
  rho   = transform_mat(rho, transfo)
  ! Kappa transforms differently from rho, but in case of real
  ! matrices this is largely irrelevant
  ! TODO: enable!
  !kappa = transform_mat(kappa, transfo)

end subroutine mixup_rhokappa_real

subroutine mixup_rhokappa_complex(rho, kappa, transfo)
  !-----------------------------------------------------------------------------
  ! Construct a random symmetry-respecting unitary transformation and apply
  ! it to (i) the single-particle wavefunctions in the Hartree-Fock basis
  ! and (ii) the matrices rho and kappa.
  !-----------------------------------------------------------------------------

  real(KIND=dp), allocatable, intent(in) :: transfo(:,:)
  complex(KIND=dp), intent(inout) :: rho(:,:), kappa(:,:)

  ! Perform the transformation of the spwfs
  call transform_spwfs_inplace(hfpsi, transfo)
  ! reperform derivatives((
  call deriveHF()
  ! Transform the density matrix and canonical kappa
  rho   = transform_mat(rho, transfo)
  ! Kappa transforms differently from rho, but in case of real
  ! matrices this is largely irrelevant
  ! TODO: enable!
  !kappa = transform_mat(kappa, transfo)

end subroutine mixup_rhokappa_complex




subroutine print_maxval(name, den_sym, den_asym)
  !
  !
  !
  !
  complex(KIND=dp), intent(in) :: den_sym(:,:), den_asym(:,:)
  character(len=*), intent(in) :: name
  print *, '--------------------------------------------'
  print *, name
  print *, '--------------------------------------------'
  print *, ' sym  Re', maxval(abs(DBLE(den_sym)))
  print *, ' sym  Im', maxval(abs(IMAG(den_sym)))
  print *, ' asym Re', maxval(abs(DBLE(den_asym)))
  print *, ' asym Im', maxval(abs(IMAG(den_asym)))
  print *
end subroutine print_maxval









! subroutine MassageDensity()
!     !---------------------------------------------------------------------------
!     ! Operate on the density before feeding it into the rest of the program.
!     !---------------------------------------------------------------------------
!     real(KIND=dp), target  :: resid(nx*ny*nz,4)
! $NTR real(KIND=dp), target :: sresid(nx*ny*nz,3,4)
!     if(all(D_I_I_hist.eq.0.0)) return
!     !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     ! Compute the residual
!     resid = D_I_I - D_I_I_hist(:,:,1)
! $NTR    sresid = D_I_S - D_I_S_hist(:,:,:,1)
!     !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     ! Perform mixing
!     select case(densitymixing)
!     case(0)
!         !-----------------------------------------------------------------------
!         ! Precondition the potentials instead of the densities.
!         ! So do nothing to the densities.
!     case(1)
!         !-----------------------------------------------------------------------
!         ! Simple linear mixing at the moment.
!         D_I_I = D_I_I_hist(:,:,1) + (1-denmix) * resid
! $NTR    D_I_S = D_I_S_hist(:,:,:,1) + (1-denmix) * sresid
!     end select
!     !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     ! Safeguard
!     where(D_I_I.lt.1d-10) D_I_I = 0
! end subroutine MassageDensity

function couple_iso(density, iso) result(coupled)
    !---------------------------------------------------------------------------
    ! All densities are calculated and stored in proton-neutron format, but 
    ! terms in the EDF are generally calculated in isospin formalism. 
    !
    ! This function takes as argument a density (with ONLY an proton/neutron)
    ! index and produces the isoscalar or isovector combination for on-the-fly
    ! resummation and easy Hephaestos code generation.
    !
    !---------------------------------------------------------------------------
    integer, intent(in)       ::  iso
    real(KIND=dp), intent(in) ::  density(mv,2)
    real(KIND=dp)             ::  coupled(mv)
    
    select case(iso)
    case(0)
      ! Isoscalar = neutron + proton
      coupled = sum(density(:,2))
    case(1)
      ! Isovector = neutron - proton
      coupled = density(:,1) - density(:,2)
    case DEFAULT
      call stp('Invalid iso argument to couple_iso.')
    end select

end function couple_iso


#if($FAM == 0) 
! We don't define this routine for FAM calculations because the pointer remapping
! does not play nice with the complex-valued densities.
subroutine print_boxsize_check(R)
  !-----------------------------------------------------------------------------
  ! Print the maximum values of the densities D_I_I and DP_I_I at the edges of 
  ! the box, to check that we are not dealing with a too small box.
  !-----------------------------------------------------------------------------
  type(DensityVector), intent(in), target :: R
  real(KIND=dp), pointer                  :: rho3D(:,:,:,:)
  real(KIND=dp), pointer                  :: rhoP_3D(:,:,:,:)

  1 format ('------------------------- Box Size Check ---------------------------')
 11 format (' Normal  density        rho  =  D_I_I')
 12 format (' Pairing density \tilde{rho} = DP_I_I')
  2 format (' Xmax = (nx+0.5)dx = ', f10.3, ' fm,  max(rho(X=Xmax)) = ', es12.3 )
! 21 format (' Xmin =-(nx+0.5)dx = ', f6.3, 'fm,  max(rho(X=Xmin)) = ', e12.3 )

  3 format (' Ymax = (ny+0.5)dx = ', f10.3, ' fm,  max(rho(Y=Ymax)) = ', es12.3 )
! 31 format (' Ymin =-(ny+0.5)dx = ', f5.3, 'fm,  max(rho(Y=Ymax)) = ', e12.3 )

  4 format (' Zmax = (nz+0.5)dx = ', f10.3, ' fm,  max(rho(Z=Zmax)) = ', es12.3 )
$PBROKEN 41 format (' Zmin =-(nz+0.5)dx = ', f10.3, ' fm,  max(rho(Z=Zmin)) = ', es12.3 )

  rho3D(1:nx,1:ny,1:nz,1:2)   => R%D_I_I
  
  print 1
  print 11
  print 2, meshX(nx) , maxval(sum(rho3D(nx,:,:,:),3))
  print 3 , meshY(ny), maxval(sum(rho3D(:,ny,:,:),3))
  print 4 , meshZ(nz), maxval(sum(rho3D(:,:,nz,:),3))
$PBROKEN  print 41, meshZ(1) , maxval(sum(rho3D(:,:,1,:),3))

  if(pairingtype .ne. 0) then
    rhoP_3D(1:nx,1:ny,1:nz,1:2) => R%DP_I_I

    print 12
    print 2, meshX(nx) , maxval(abs(sum(rhoP_3D(nx,:,:,:),3)))
    print 3 , meshY(ny), maxval(abs(sum(rhoP_3D(:,ny,:,:),3)))
    print 4 , meshZ(nz), maxval(abs(sum(rhoP_3D(:,:,nz,:),3)))
$PBROKEN  print 41, meshZ(1) , maxval(abs(sum(rhoP_3D(:,:,1,:),3)))
  endif
end subroutine print_boxsize_check
#endif 

#if(USE_HDF5>0 && $FAM == 0) 
! There is no need to write densities to file in a FAM code
subroutine write_hdf5_densities(file_id, R)
  !---------------------------------------------------------------------------
  ! Write the densities to a HDF5 file 
  !
  ! Input:
  !   file_id : HDF5 file identifier
  !   R       : DensityVector type containing the densities to be written
  !---------------------------------------------------------------------------
  use HDF5
  use HDF5_auxiliary 
  integer(HID_T), INTENT(IN) :: file_id 
  type(DensityVector), intent(in) :: R

    ! Hephaestos fills in a call to hdf5_write_dataset_1d for every density
$WRITEDENSITIES_HDF5

end subroutine write_hdf5_densities
#endif
function divJ_spwf(der_index)
    !---------------------------------------------------------------------------
    ! Calculate the
    !            nabla cdot J
    ! where J is the vector component of the spin-current density J_munu.
    ! The contribution from a single spwf is
    !
    !  sum_{mu nu kappa} eps_munukappa
    !            rho_ii  Im [ \nabla_mu Psi_i^* \nabla_kappa \sigma_nu \Psi_ii]
    !
    ! where eps_munukappa is a Levi-Civita symbol.
    !
    !
    ! Input:
    !    der_index:  index of the derivative of the wavefunction in the array
    ! Output:
    !    divJ_spwf: unweighted contribution of this spwf to the divJ
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! The reason this thing is calculated separately from all other densities
    ! is because the sum of spwfs is represented much more accurately on the
    ! mesh than the sum of erivatives of J on the mesh. We use this summation
    ! for the calculation of multipole moments of J, but it should in W.R.'s
    ! opinion NOT be used in the calculation of any energy for consistency reasons.
    !
    ! Originally, this routine included the sums over spwfs; which was
    ! easier to understand. Unfortunately, compatibility with
    !       store_derivatives = .true.
    ! flag requires that this function call be within a loop.
    !
    !---------------------------------------------------------------------------
    integer, intent(in) :: der_index

    real(KIND=dp) :: divJ_spwf(nx*ny*nz)
    ! MB 24/12/14 comment use of Pauli back in now that the memory leak is fixed
    real(KIND=dp) :: temp(nx*ny*nz,4)

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! I code this with calls to the Pauli and ImagMultiplySpinor functions
    ! to make no mistakes
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! x y z
    temp      = Pauli(DenDPsi(:,2,:,der_index), 3)
    divJ_spwf =             ImagMultiplySpinor(DendPsi(:,1,:,der_index), temp)
    ! y x z
    temp = Pauli(DenDPsi(:,1,:,der_index), 3)
    divJ_spwf = divJ_spwf - ImagMultiplySpinor(DendPsi(:,2,:,der_index), temp)
    ! x z y
    temp = Pauli(DenDPsi(:,3,:,der_index), 2)
    divJ_spwf = divJ_spwf - ImagMultiplySpinor(DendPsi(:,1,:,der_index), temp)
    ! z x y
    temp = Pauli(DenDPsi(:,1,:,der_index), 2)
    divJ_spwf = divJ_spwf + ImagMultiplySpinor(DendPsi(:,3,:,der_index), temp)
    ! y z x
    temp = Pauli(DenDPsi(:,3,:,der_index), 1)
    divJ_spwf = divJ_spwf + ImagMultiplySpinor(DendPsi(:,2,:,der_index), temp)
    ! z y x
    temp = Pauli(DenDPsi(:,2,:,der_index), 1)
    divJ_spwf = divJ_spwf - ImagMultiplySpinor(DendPsi(:,3,:,der_index), temp)

end function divJ_spwf

function CompNablaMelements(spwf_basis) result(NablaMelements)
    !---------------------------------------------------------------------------
    ! Computes the matrix elements of Nabla
    !
    !   NablaMElements(mu,1,i,j) = Re < Psi_i | \nabla_{\mu} | \Psi_j >
    !   NablaMElements(mu,2,i,j) = Im < Psi_i | \nabla_{\mu} | \Psi_j >
    !    (mu = x/1 ,y/2, z/3)
    ! 
    ! for a given single-particle spwf_basis = HF/CAN/DEN. The corresponding 
    ! argument is optional, and defaults to 'DEN'.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! The current MPI parallelization strategy to calculate a given matrix 
    ! element is simple: transfer all relevant spwfs to an indicated rank.
    !    1) find out the MPI ranks that store psi_i and psi_j, rank_i and rank_j
    !    2) send both to the designated rank
    !    3) which integrates psi_i \nabla_mu psi_j 
    !    4) at this point we use the symmetry of the operator, i.e.
    !          <psi_i |\nabla_mu |psi_j> = - <psi_j| \nabla_mu| psi_i>^*
    !       Do NOT forget that \nabla is an antihermitian operator!
    ! 
    ! This does not load balance very well, but at least we take the effort of
    ! load balancing the protons and neutrons. 
    !    Sequential : everything done by rank 0 
    !    Parallel   : everything done by rank 0 (neutrons)
    !                 and the last rank. These should in all realistic cases 
    !                 only store neutron/proton spwfs respectively and should
    !                 hence never communicate.
    !    5) MPI_BCAST creates the complete array for all ranks
    !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! TODO: this routine can be refactored by
    !       - adding a bunch of small functions for repeated instructions!
    !       - removing the MPI, as this is functionality that is not 
    !         currently used.
    ! .... but what should really be done is abstracting the concept of a 
    !      single-particle operator with different symmetry properties.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Some notes:
    ! 1) These matrix elements are, in general, complex numbers!
    ! 2) Symmetries restrict them in weird ways:
    !    Same p             => <Nx> = <Ny> = <Nz> = 0
    !    Same signature     => <Nx> = <Ny> = 0
    !    Different signature=> <Nz> = 0
    !    Different isospin  => <Nx> = <Ny> = <Nz> = 0
    ! 3) We must not forget the time-reversed states when time reversal is
    !    conserved!
    ! 4) This correction is a nightmare to efficiently parallelize using MPI.
    !    The reason is that all matrix elements that are not zero involve 
    !    spwfs with different parity and often opposite signature. These spwfs
    !    are unlikely to be found on the same MPI rank => massive communication 
    !    costs. I have at the moment not made any effort to be clever. The 
    !    implementation below is naive but works and can serve as benchmark
    !    to future implementations that are more clever.
    ! 5) This routine is NOT ready for generalization beyond a CR4 (PT-broken)
    !    context.
    !
    ! Luckily, these matrix elements are not necessary in the most demanding 
    ! case with respect to scaling with the number of spwfs: nuclear pasta 
    ! calculations describe essentially infinite systems which are well-located
    ! in space.
    !---------------------------------------------------------------------------
    character(len=*), intent(in), optional :: spwf_basis
    character(len=3) :: basis

    integer       :: i, j, B, N, si, N2, N3, N4, wave, wave2
    integer       :: ranki, rankj, wave_global, wave2_global
    integer       :: designated_rank(2), calc_rank, T
    real(KIND=dp) :: NablaMElements(3,2,nwt,nwt), psi(mv,4)
    real(KIND=dp) :: derx(mv,4), dery(mv,4), derz(mv,4)
    logical       :: hidden_TR = .false.
#if(USE_MPI>0)
    integer       :: mpi_err
#endif

    call start_timer(T_NablaMElements)
    if(present(spwf_basis)) then
      if(      (spwf_basis.ne.'DEN') &
      &  .and. (spwf_basis.ne.'HF' ) &
      &  .and. (spwf_basis.ne.'CAN')) then ! always sanity check
        call stp('NablaMElements needs a single-particle basis to work with.')
      endif
      basis = spwf_basis
   else 
      basis = 'DEN'
    endif

    NablaMElements= 0.0_dp

    designated_rank(1) = 0
    designated_rank(2) = NPROCS - 1

    if(NPROCS .gt. 1) then
      ! Verify that both designated ranks should have no communications; if they
      ! do this calculation will take very long.

      if((MPI_RANK.eq.designated_rank(1)) &
      &          .and.                    &
      &        (sum(HFBlocks(5:8)).gt. 0)) then
        call stp('This load balancing cannot be used in CompNablaMelements.')
      else if((MPI_RANK.eq.designated_rank(2)) &
      &          .and.                         &
      &            (sum(HFBlocks(1:4)).gt. 0)) then
        call stp('This load balancing cannot be used in CompNablaMelements.') 
      endif
    endif

    si = 0 
    do B=1,8,4 ! This loop is essentially over protons vs neutrons
      N = HFblocks_global(B)  ;  N2= HFBlocks_global(B+1)
      N3= HFBlocks_global(B+2);  N4= HFBlocks_global(B+3)

      calc_rank = designated_rank(mod(B,5) + 1)

      !-------------------------------------------------------------------------
      ! Re < p_z > : only non-zero matrix elements when equal signature
      !              and opposite parity (if it is conserved)
      !-------------------------------------------------------------------------
      ! Block 1 with block 3
      do j=1,N3
        wave2_global = si+N+N2+j
        rankj        = rank_map(wave2_global)
        wave2        = spwf_inverse(wave2_global)

        ! sending, receing, copying etc to the array psi for calc_rank only
#if(USE_MPI>0)
        call transfer_derpsi(derz, wave2, 3,basis,.false., rankj, calc_rank)
#else
        call transfer_derpsi(derz, wave2, 3,basis,.false.)
#endif
        do i=1,N
          wave_global = si+i  ! this is the global index of the spwf
          ! .... and we have to find the MPI RANK and the index it is on
          ranki = rank_map(wave_global)
          wave  = spwf_inverse(wave_global)
          ! sending, receing, copying etc to the array psi for calc_rank only
#if(USE_MPI>0)
          call transfer_psi(psi, wave, basis,ranki, calc_rank )
#else
          call transfer_psi(psi, wave, basis)
#endif
          if(MPI_RANK .eq. calc_rank) then
            ! Perform the calculation with the rank designated to calculate
            NablaMElements(3,1,wave_global,wave2_global) = dv*                 &
            & sum(                  derz(:,1) * psi(:,1) + derz(:,2) * psi(:,2)&
            &                    +  derz(:,3) * psi(:,3) + derz(:,4) * psi(:,4))

            NablaMElements(3,1,wave2_global,wave_global) =&
            &                     - NablaMElements(3,1,wave_global,wave2_global)
          endif
        enddo
      enddo     

      ! Block 1 with block 1 (parity broken)
$PBROKEN   do j=1,N
$PBROKEN     wave2_global = si+j
$PBROKEN     rankj        = rank_map(wave2_global)
$PBROKEN     wave2        = spwf_inverse(wave2_global)
$PBROKEN
$PBROKEN     ! sending, receing, copying etc to the array psi for calc_rank only
#if(USE_MPI>0)
$PBROKEN     call transfer_derpsi(derz,wave2, 3,basis,.false.,rankj,calc_rank)
#else
$PBROKEN     call transfer_derpsi(derz,wave2, 3,basis,.false.)
#endif
$PBROKEN     do i=1,N
$PBROKEN       wave_global = si+i  ! this is the global index of the spwf
$PBROKEN       ! .... and we have to find the MPI RANK and the index it is on
$PBROKEN       ranki = rank_map(wave_global)
$PBROKEN       wave  = spwf_inverse(wave_global)
#if(USE_MPI>0)
$PBROKEN       call transfer_psi(psi, wave, basis, ranki, calc_rank)
#else
$PBROKEN       call transfer_psi(psi, wave, basis)
#endif
$PBROKEN
$PBROKEN       if(MPI_RANK .eq. calc_rank) then
$PBROKEN         NablaMElements(3,1,wave_global,wave2_global) = dv*            &
$PBROKEN         & sum(             derz(:,1) * psi(:,1) + derz(:,2) * psi(:,2)&
$PBROKEN         &               +  derz(:,3) * psi(:,3) + derz(:,4) * psi(:,4))
$PBROKEN
$PBROKEN         NablaMElements(3,1,wave2_global,wave_global) = &
$PBROKEN         &               -  NablaMElements(3,1,wave_global,wave2_global)
$PBROKEN        endif
$PBROKEN     enddo
$PBROKEN   enddo

      ! Block 2 with block 4
      do j=1,N4
        wave2_global = si+N+N2+N3+j
        rankj        = rank_map(wave2_global)
        wave2        = spwf_inverse(wave2_global)

        ! sending, receing, copying etc to the array psi for calc_rank only
#if(USE_MPI>0)
        call transfer_derpsi(derz, wave2, 3, basis,.false., rankj, calc_rank)
#else
        call transfer_derpsi(derz, wave2, 3, basis,.false.)
#endif
        do i=1,N2
          wave_global  = si+N+i
          ranki        = rank_map(wave_global)
          wave         = spwf_inverse(wave_global)
          ! sending, receing, copying etc to the array psi for calc_rank only
#if(USE_MPI>0)
          call transfer_psi(psi, wave, basis, ranki, calc_rank)
#else
          call transfer_psi(psi, wave, basis)
#endif
          if(MPI_RANK.eq.calc_rank) then
          ! Re < p_z > 
            NablaMElements(3,1,wave_global,wave2_global) = dv*                 &
            & sum(                  derz(:,1) * psi(:,1) + derz(:,2) * psi(:,2)&
            &                    +  derz(:,3) * psi(:,3) + derz(:,4) * psi(:,4))

            NablaMElements(3,1,wave2_global,wave_global) = &
            &                    -  NablaMElements(3,1,wave_global,wave2_global)
          endif
        enddo
      enddo
      
      ! Block 2 with block 2 (parity broken)
$PBROKEN   do j=1,N2
$PBROKEN     wave2_global = si+N+j
$PBROKEN     rankj        = rank_map(wave2_global)
$PBROKEN     wave2        = spwf_inverse(wave2_global)
$PBROKEN
#if(USE_MPI>0)
$PBROKEN     call transfer_derpsi(derz,wave2,3,basis,.false.,rankj,calc_rank)
#else
$PBROKEN     call transfer_derpsi(derz,wave2,3,basis,.false.)
#endif
$PBROKEN     do i=1,N2
$PBROKEN       wave_global = si+N+i  ! this is the global index of the spwf
$PBROKEN       ! .... and we have to find the MPI RANK and the index it is on
$PBROKEN       ranki = rank_map(wave_global)
$PBROKEN       wave  = spwf_inverse(wave_global)
#if(USE_MPI>0)
$PBROKEN       call transfer_psi(psi, wave, basis, ranki, calc_rank)
#else
$PBROKEN       call transfer_psi(psi, wave, basis)
#endif
$PBROKEN 
$PBROKEN       if(MPI_RANK.eq.calc_rank) then
$PBROKEN         ! Re < p_z > 
$PBROKEN         NablaMElements(3,1,wave_global,wave2_global) = dv*            &
$PBROKEN         & sum(             derz(:,1) * psi(:,1) + derz(:,2) * psi(:,2)&
$PBROKEN         &               +  derz(:,3) * psi(:,3) + derz(:,4) * psi(:,4))
$PBROKEN
$PBROKEN         NablaMElements(3,1,wave2_global,wave_global) = &
$PBROKEN         &               -  NablaMElements(3,1,wave_global,wave2_global)
$PBROKEN       endif
$PBROKEN     enddo
$PBROKEN   enddo

      !----------------------------------------------------------------------------
      ! Re < p_x >, Im <p_y> 
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! hidden_TR = .true. => Include a time-reversal operation in the calculation
      !                       of matrix elements of these operators
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
$NTR  hidden_TR = .false.
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! hidden_TR = .true. => Do not include a time-reversal operation in the
      !                       calculation of matrix elements
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
$TR   hidden_TR = .true.
      !
      ! Note: I define these logicals here in order to cut down on the preprocessor
      !       directives needed just below.
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      ! Block 1 with block 4  (T-broken)
$NTR  do j=1, N4
$NTR     wave2_global = si+N+N2+N3+j
      ! Block 1 with block 3  (T-conserved)
$TR   do j=1, N3
$TR     wave2_global = si+N+j
        rankj        = rank_map(wave2_global)
        wave2        = spwf_inverse(wave2_global)
#if(USE_MPI>0)
        call transfer_derpsi(derx, wave2, 1,basis,hidden_TR,rankj,calc_rank)
        call transfer_derpsi(dery, wave2, 2,basis,hidden_TR,rankj,calc_rank)
#else
        call transfer_derpsi(derx, wave2, 1,basis,hidden_TR)
        call transfer_derpsi(dery, wave2, 2,basis,hidden_TR)
#endif
        do i=1,N
          wave_global  = si+i
          ranki        = rank_map(wave_global)
          wave         = spwf_inverse(wave_global)
          ! sending, receing, copying etc to the array psi for calc_rank only
#if(USE_MPI>0)
          call transfer_psi(psi, wave, basis, ranki, calc_rank)
#else
          call transfer_psi(psi, wave, basis)
#endif
          if(MPI_RANK.eq.calc_rank) then

            NablaMElements(1,1,wave_global,wave2_global) = dv*                 &
            & sum(                 derx(:,1) * psi(:,1) + derx(:,2) * psi(:,2) &
            &                   +  derx(:,3) * psi(:,3) + derx(:,4) * psi(:,4))
  
            NablaMElements(2,2,wave_global,wave2_global) = dv*                 &
            & sum(                  dery(:,2) * psi(:,1) - dery(:,1) * psi(:,2)&
            &                    -  dery(:,3) * psi(:,4) + dery(:,4) * psi(:,3))

            ! Attention to the signs here
            NablaMElements(1,1,wave2_global,wave_global) = &
            &                     - NablaMElements(1,1,wave_global,wave2_global)
            NablaMElements(2,2,wave2_global,wave_global) = &
            &                       NablaMElements(2,2,wave_global,wave2_global)
         endif
        enddo
      enddo
      ! Block 1 with block 2  (T-broken, P-broken)
$PBROKEN $NTR do j=1, N2
$PBROKEN $NTR   wave2_global = si+N+j
$PBROKEN $NTR   rankj        = rank_map(wave2_global)
$PBROKEN $NTR   wave2        = spwf_inverse(wave2_global)
#if(USE_MPI>0)
$PBROKEN $NTR   call transfer_derpsi(derx,wave2,1,basis,hidden_TR,rankj,calc_rank)
$PBROKEN $NTR   call transfer_derpsi(dery,wave2,2,basis,hidden_TR,rankj,calc_rank)
#else
$PBROKEN $NTR   call transfer_derpsi(derx,wave2,1,basis,hidden_TR)
$PBROKEN $NTR   call transfer_derpsi(dery,wave2,2,basis,hidden_TR)
#endif

              ! Block 1 with block 1  (T-conserved, P-broken)
$PBROKEN $TR  do j=1, N
$PBROKEN $TR    wave2_global = si+j
$PBROKEN $TR    rankj        = rank_map(wave2_global)
$PBROKEN $TR    wave2        = spwf_inverse(wave2_global)
#if(USE_MPI>0)
$PBROKEN $TR    call transfer_derpsi(derx,wave2,1,basis,hidden_TR,rankj,calc_rank)
$PBROKEN $TR    call transfer_derpsi(dery,wave2,2,basis,hidden_TR,rankj,calc_rank)
#else
$PBROKEN $TR    call transfer_derpsi(derx,wave2,1,basis,hidden_TR)
$PBROKEN $TR    call transfer_derpsi(dery,wave2,2,basis,hidden_TR)
#endif
$PBROKEN        do i=1,N
$PBROKEN          wave_global  = si+i
$PBROKEN          ranki        = rank_map(wave_global)
$PBROKEN          wave         = spwf_inverse(wave_global)
$PBROKEN          ! sending, receing, copying etc to the array psi for calc_rank only
#if(USE_MPI>0)
$PBROKEN          call transfer_psi(psi, wave, basis, ranki, calc_rank)
#else
$PBROKEN          call transfer_psi(psi, wave, basis)
#endif

$PBROKEN        if(MPI_RANK .eq. calc_rank) then
$PBROKEN          NablaMElements(1,1,wave_global,wave2_global) = dv*           &
$PBROKEN          & sum(         derx(:,1) * psi(:,1) + derx(:,2) * psi(:,2)   &
$PBROKEN          &           +  derx(:,3) * psi(:,3) + derx(:,4) * psi(:,4))

$PBROKEN          NablaMElements(2,2,wave_global,wave2_global) = dv*           &
$PBROKEN          & sum(         dery(:,2) * psi(:,1) - dery(:,1) * psi(:,2)   &
$PBROKEN          &           -  dery(:,3) * psi(:,4) + dery(:,4) * psi(:,3))

$PBROKEN          NablaMElements(1,1,wave2_global,wave_global) =               &
$PBROKEN          &               - NablaMElements(1,1,wave_global,wave2_global)
$PBROKEN          NablaMElements(2,2,wave2_global,wave_global) =               &
$PBROKEN          &                 NablaMElements(2,2,wave_global,wave2_global)
$PBROKEN        endif
$PBROKEN     enddo
$PBROKEN   enddo

     ! Block 2 with block 3  (T-broken)
$NTR do j=1, N3
$NTR    wave2_global = si+N+N2+j
     ! Block 2 with block 4  (T-conserved)
$TR  do j=1, N4
$TR     wave2_global = si+N+N2+N3+j
        rankj        = rank_map(wave2_global)
        wave2        = spwf_inverse(wave2_global)
#if(USE_MPI>0)
        call transfer_derpsi(derx, wave2, 1,basis,hidden_TR,rankj,calc_rank)
        call transfer_derpsi(dery, wave2, 2,basis,hidden_TR,rankj,calc_rank)
#else
        call transfer_derpsi(derx, wave2, 1,basis,hidden_TR)
        call transfer_derpsi(dery, wave2, 2,basis,hidden_TR)
#endif
        do i=1,N2
          wave_global  = si+N+i
          ranki        = rank_map(wave_global)
          wave         = spwf_inverse(wave_global)
          ! sending, receing, copying etc to the array psi for calc_rank only
#if(USE_MPI>0)
          call transfer_psi(psi, wave, basis, ranki, calc_rank)
#else
          call transfer_psi(psi, wave, basis)
#endif
          if(MPI_RANK.eq.calc_rank) then
            NablaMElements(1,1,wave_global,wave2_global) = dv*                 &
            & sum(                 derx(:,1) * psi(:,1) + derx(:,2) * psi(:,2) &
            &                   +  derx(:,3) * psi(:,3) + derx(:,4) * psi(:,4))!

            NablaMElements(2,2,wave_global,wave2_global) = dv*                 &
            & sum(                  dery(:,2) * psi(:,1) - dery(:,1) * psi(:,2)&
            &                    -  dery(:,3) * psi(:,4) + dery(:,4) * psi(:,3))
            ! Attention to the signs below
            NablaMElements(1,1,wave2_global,wave_global) = &
            &                     - NablaMElements(1,1,wave_global,wave2_global)
            NablaMElements(2,2,wave2_global,wave_global) = &
            &                       NablaMElements(2,2,wave_global,wave2_global)
          endif
        enddo
      enddo

              ! Block 2 with block 2 (P-broken, T-conserved)
              ! Note: there is no "time-reversal broken" version of the following
              ! code block since "block 1 with 2" is already taken into account
              ! above..
$PBROKEN $TR  do j=1, N2
$PBROKEN $TR   wave2_global = si+N+j
$PBROKEN $TR   rankj        = rank_map(wave2_global)
$PBROKEN $TR   wave2        = spwf_inverse(wave2_global)
#if(USE_MPI>0)
$PBROKEN $TR   call transfer_derpsi(derx,wave2,1,basis,hidden_TR,rankj,calc_rank)
$PBROKEN $TR   call transfer_derpsi(dery,wave2,2,basis,hidden_TR,rankj,calc_rank)
#else
$PBROKEN $TR   call transfer_derpsi(derx,wave2,1,basis,hidden_TR)
$PBROKEN $TR   call transfer_derpsi(dery,wave2,2,basis,hidden_TR)
#endif

$PBROKEN $TR   do i=1,N2
$PBROKEN $TR      wave_global  = si+N+i
$PBROKEN $TR      ranki        = rank_map(wave_global)
$PBROKEN $TR      wave         = spwf_inverse(wave_global)
#if(USE_MPI>0)
$PBROKEN $TR      call transfer_psi(psi, wave, basis, ranki, calc_rank)
#else
$PBROKEN $TR      call transfer_psi(psi, wave, basis)
#endif
$PBROKEN $TR      if(MPI_RANK.eq.calc_rank) then
$PBROKEN $TR        NablaMElements(1,1,wave_global,wave2_global) = dv*           &
$PBROKEN $TR        & sum(         derx(:,1) * psi(:,1) + derx(:,2) * psi(:,2)   &
$PBROKEN $TR        &           +  derx(:,3) * psi(:,3) + derx(:,4) * psi(:,4))

$PBROKEN $TR        NablaMElements(2,2,wave_global,wave2_global) = dv*           &
$PBROKEN $TR        & sum(          dery(:,2) * psi(:,1) - dery(:,1) * psi(:,2)  &
$PBROKEN $TR        &            -  dery(:,3) * psi(:,4) + dery(:,4) * psi(:,3))

$PBROKEN $TR        NablaMElements(1,1,wave2_global,wave_global) = &
$PBROKEN $TR        &               - NablaMElements(1,1,wave_global,wave2_global)
$PBROKEN $TR        NablaMElements(2,2,wave2_global,wave_global) = &
$PBROKEN $TR        &                 NablaMElements(2,2,wave_global,wave2_global)
$PBROKEN $TR      endif
$PBROKEN $TR   enddo
$PBROKEN $TR  enddo
      T = N + N2 + N3 + N4
      si = si + T 
    enddo

#if(USE_MPI>0)
    ! Broadcast all calculated matrix elements to all ranks. 
    ! Note that this is an explicit, separate, loop such that the parallel 
    ! structure of the loop above does not get disturbed
    si = 0 
    do B=1,8,4 ! This loop is essentially over protons vs neutrons
      N = HFblocks_global(B)  ;  N2= HFBlocks_global(B+1)
      N3= HFBlocks_global(B+2);  N4= HFBlocks_global(B+3)
      T = N + N2 + N3 + N4
      calc_rank = designated_rank(mod(B,5) + 1)

      ! Broadcast all matrix elements from the calculating rank to all the rest
      call MPI_BCAST(NablaMElements(:,:,si+1:si+T, si+1:si+T), 6*T**2, &
      &              MPI_REAL8, calc_rank, MPI_COMM_WORLD, mpi_err)
      si = si + T 
    enddo
#endif
    !---------------------------------------------------------------------------
    call stop_timer(T_NablaMElements)

end function CompNablaMelements
end module densities
