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
module vectors
 !==============================================================================
 !  Module containing type declarations for both density- and potential-vectors.
 !
 !  All of these declations are separated from the densities.f90 and 
 !  functional.f90 modules such that they can be employed by lower-level 
 !  modules without much fuss.
 !
 !==============================================================================
 ! Hephaestos keywords:
 !
 ! DECLARATION         : [WAY TOO LONG]
 ! DECLARATION_FIELDS  : [WAY TOO LONG]
 ! COULOMB_REAL        : $COULOMB_REAL 
 ! COULOMB_COMPLEX     : $COULOMB_COMPLEX
 !==============================================================================

 use geninfo

 implicit none

 type DensityVector
    !-------------------------------------------------------------------------
    ! Custom type for containing all Skyrme densities
    !
    !   \mathcal{R} = (D_I_I, D_N_N, ....)
    !
    ! but also their derivatives, as well as the charge density!
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Note: chargedensity is defined on a (nx,ny,nz) mesh while all the others
    !       are defined as a (nx*ny*nz) vector for storage efficiency
    !-------------------------------------------------------------------------
$DECLARATION
$COULOMB_REAL    real(KIND=dp), allocatable    :: chargedensity(:,:,:)
$COULOMB_COMPLEX complex(KIND=dp), allocatable :: chargedensity(:,:,:)
    !---------------------------------------------------------------------------
    ! Separate, manual, declaration of Div.J(r) as calculated from the spwfs
    ! for the more accurate calculation of its multipole moments
    real(KIND=dp), allocatable :: divJ(:,:)
  end type DensityVector

  type PotentialVector
    !-------------------------------------------------------------------------
    ! Custom type for containing all potentials
    !
    !   \mathcal{F} = (F_I_I, F_N_N, ....)
    !-------------------------------------------------------------------------
$DECLARATION_POTENTIALS
    ! Coulomb potentials; both unfolded and folded.
$COULOMB_REAL    real(KIND=dp), allocatable :: CoulombPotential(:,:,:)
$COULOMB_REAL    real(KIND=dp), allocatable :: ExchangePotential(:,:,:)
$COULOMB_REAL    real(KIND=dp), allocatable :: FoldedCoul(:,:,:,:), FoldedExchange(:,:,:,:)
$COULOMB_COMPLEX complex(KIND=dp), allocatable :: CoulombPotential(:,:,:)
$COULOMB_COMPLEX complex(KIND=dp), allocatable :: ExchangePotential(:,:,:)
$COULOMB_COMPLEX complex(KIND=dp), allocatable :: FoldedCoul(:,:,:,:), FoldedExchange(:,:,:,:)
 end type PotentialVector

 !---------------------------------------------------------------------------
 ! The amount of iterations to keep in memory for the density and/or potential
 ! mixing and estimation of the convergence rate
#if(PASTA == 0)
 integer            :: memory = 3
#else
 ! We squeeze out every drop of memory we can
 integer            :: memory = 0
#endif

contains

subroutine construct_charge_density(R, sx, sy, sz)
    !---------------------------------------------------------------------------
    ! Construct the charge density from the proton and neutron densities,
    ! using nucleonic
    !
    ! Input:
    !     R : the DensityVector for which we need to construct the charge density
    !     sx, sy, sz : symmetries of the density, explicitly passed in because
    !                  R might not be a density obtained in a static mean-field
    !                  calculation
    ! Output:
    !     R%chargedensity : the charge density, defined on a (nx,ny,nz) mesh
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! TODO: document what this routine does precisely
    !       -> this includes a transfer from 1D to 3D
    !---------------------------------------------------------------------------
    use timing, only: start_timer, stop_timer, T_chargedensity
    use Folding, only: fold_form_factor

#if(PASTA==1)
    4 format ('--------------------------------------------------------------------')
    1 format (' Warning: the charge in your system is not equal to the desired one.')
    2 format (' Number of protons - \int charge density = ', es10.3 )
    3 format (' The electron density is compensating.')
#endif

    type(DensityVector),intent(inout) :: R
    integer, intent(in)        :: sx, sy, sz
    integer                    :: i,j,k, it
#if(PASTA==1)
    real(KIND=dp)              :: rho_el, volume
#endif

$COULOMB_COMPLEX    complex(KIND=dp)              :: temp(nx,ny,nz,2)
$COULOMB_COMPLEX    complex(KIND=dp), allocatable :: folded(:,:,:)
$COULOMB_REAL       real(KIND=dp)                 :: temp(nx,ny,nz,2)
$COULOMB_REAL       real(KIND=dp), allocatable    :: folded(:,:,:)

    call start_timer(T_chargedensity)

    ! Deallocation such that rho_charge does not have the wrong dimensions
    if(allocated(R%chargedensity))      deallocate(R%chargedensity)
    if(.not.allocated(R%chargedensity)) allocate(R%chargedensity(nx,ny,nz))
    !---------------------------------------------------------------------------
    ! Transfer the density from a (nx*ny*nz,2) to a (nx,ny,nz,2) array
    do it=1,2
      do k=1,nz
        do j=1,ny
          do i=1,nx
              temp(i,j,k,it) = R%D_I_I(meshindex(i,j,k),it)
          enddo
        enddo
      enddo
    enddo
    !---------------------------------------------------------------------------
    ! Perform the actual folding with form factors
    call fold_form_factor(temp, folded, nx,ny,nz, sx, sy, sz)
    if( .not. allocated(folded)) then
      ! if no folding is applied, just use the proton density 
      R%chargedensity = temp(:,:,:,2)
    else
      ! if folding is applied, use the folded density
      R%chargedensity = folded
    endif
    !---------------------------------------------------------------------------
    ! When performing simulations for nuclear pasta, one assumes the entire 
    ! volume is charge neutral: a constant background of electrons floods the 
    ! entire simulation volume. We subtract this backound here.
#if(PASTA==1)
    volume=nx*ny*nz*dv   ! simplification by WR: the physical volume simulated
                         ! can just be gotten by the volume element...
    rho_el=sum(R%chargedensity)*dv/volume
    ! Note: it is CRUCIAL to put here the integral of the charge density as
    !       opposed to just the number of protons. If, for whatever reason,
    !       the code fails to build a proton + neutron charge density that
    !       does not integrate perfectly to tthe number of protons, then
    !       putting the number of protons here will lead to a small amount
    !       of charge; this will blow up the Coulomb solver if periodic
    !       boundary conditions are applied.

    if(abs(sum(R%chargedensity)*dv - protons) > 1e-7) then
      print 4
      print 1
      print 2, protons - sum(R%chargedensity)*dv
      print 4
    endif
    R%chargedensity = R%chargedensity -rho_el
#endif
    call stop_timer(T_chargedensity)
end subroutine construct_charge_density

function memory_for_densities() result (stor)
  !-----------------------------------------------------------------------------
  ! Return an estimation for the total memory required to store all densities.
  !-----------------------------------------------------------------------------
  use iso_fortran_env, only: int64
  integer(KIND=int64) :: stor
  stor = 0
$MEMORY_DENSITIES

end function memory_for_densities

 end module vectors


