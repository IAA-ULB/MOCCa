module vectors
 !==============================================================================
 !_________ _______  _       _________ _______  _                 _______
 !\__   __/(  ___  )( (    /|\__   __/(  ___  )( \      |\     /|(  ____ \
 !   ) (   | (   ) ||  \  ( |   ) (   | (   ) || (      | )   ( || (    \/
 !   | |   | (___) ||   \ | |   | |   | (___) || |      | |   | || (_____
 !   | |   |  ___  || (\ \) |   | |   |  ___  || |      | |   | |(_____  )
 !   | |   | (   ) || | \   |   | |   | (   ) || |      | |   | |      ) |
 !   | |   | )   ( || )  \  |   | |   | )   ( || (____/\| (___) |/\____) |
 !   )_(   |/     \||/    )_)   )_(   |/     \|(_______/(_______)\_______)
 !
 !  Copyright W. Ryssens & M. Bender
 !
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

subroutine ConstructChargeDensity(R, sx, sy, sz, proton_size, neutron_size)
    !---------------------------------------------------------------------------
    ! Construct the charge density from the proton and neutron densities,
    ! using various effective forms
    !
    ! Input:
    !     TODO:
    ! Output:
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! TODO: document what this routine does precisely
    !       -> this includes a transfer from 1D to 3D
    !---------------------------------------------------------------------------
    use timing, only: start_timer, stop_timer, T_chargedensity
    use Folding, only: ConstructFoldingMatrices, FoldGaussian, Gaussx, Gaussy, Gaussz

#if(PASTA==1)
    4 format ('--------------------------------------------------------------------')
    1 format (' Warning: the charge in your system is not equal to the desired one.')
    2 format (' Number of protons - \int charge density = ', es10.3 )
    3 format (' The electron density is compensating.')
#endif

    type(DensityVector),intent(inout) :: R
    real(KIND=dp), intent(in)  :: proton_size(2), neutron_size(2)
    integer, intent(in)        :: sx, sy, sz
    integer                    :: i,j,k
#if(PASTA==1)
    real(KIND=dp)              :: rho_el, volume
#endif

$COULOMB_COMPLEX    complex(KIND=dp)        :: temp(nx,ny,nz)
$COULOMB_REAL       real(KIND=dp)           :: temp(nx,ny,nz)

    call start_timer(T_chargedensity)

    ! Deallocation such that rho_charge does not have the wrong dimensions
    if(allocated(R%chargedensity))      deallocate(R%chargedensity)
    if(.not.allocated(R%chargedensity)) allocate(R%chargedensity(nx,ny,nz))
    !---------------------------------------------------------------------------
    ! If we account for the finite extent of the charge of the nucleus, then
    ! we need to fold densities and potentials with gaussians. This sets up the
    ! required matrices.
    !
    ! Note: this little piece of code is duplicated, since in different
    !       runmodes of the code different Coulomb routines get called in
    !       different order; this makes sure we get no segfaults.
    !---------------------------------------------------------------------------
    if(any(proton_size .ne. 0.0_dp) .or. any(neutron_size.ne.0.0_dp)) then
      if(.not.allocated(Gaussx)) then
          allocate(Gaussx(nx,nx,2,2), Gaussy(ny,ny,2,2), Gaussz(nz,nz,2,2))
          Gaussx = 0.0 ;  Gaussy = 0.0 ; Gaussz = 0.0
      endif
      call ConstructFoldingMatrices(Gaussx,Gaussy,Gaussz,sx, sy, sz)
    endif
    !---------------------------------------------------------------------------
    ! Proton contributions to the charge density.
    ! We start from the proton point density
    do k=1,nz
      do j=1,ny
        do i=1,nx
            temp(i,j,k) = R%D_I_I(meshindex(i,j,k),2)
            print *, 'temp(i,j,k) = ', temp(i,j,k)
        enddo
      enddo
    enddo
    stop
    if(proton_size(1).gt.0.0) then
        ! Fold the source with a Gaussian
        R%chargedensity = &
        & FoldGaussian(temp, GaussX(:,:,1,2), GaussY(:,:,1,2), GaussZ(:,:,1,2),&
        &                                                            nx, ny, nz)
    endif
    if(proton_size(2).gt.0.0) then
        ! Fold the source with another Gaussian, this time with minus sign.
        R%chargedensity = R%chargedensity + &
        & FoldGaussian(temp, GaussX(:,:,2,2), GaussY(:,:,2,2), GaussZ(:,:,2,2),&
        &                                                            nx, ny, nz)
    endif

    if(all(proton_size.eq.0.0)) then
        R%chargedensity = temp
    endif
    !---------------------------------------------------------------------------
    ! Neutron contributions to the charge density.
    if(any(neutron_size .gt. 0.0d0)) then 
      do k=1,nz
        do j=1,ny
          do i=1,nx
             temp(i,j,k) = R%D_I_I(meshindex(i,j,k),1)
          enddo
        enddo
      enddo

      if(neutron_size(1).gt.0.0) then
          ! Fold the source with a Gaussian
          R%chargedensity = R%chargedensity + &
          & FoldGaussian(temp, GaussX(:,:,1,1), GaussY(:,:,1,1), GaussZ(:,:,1,1),&
          &                                                            nx, ny, nz)
      endif
      if(neutron_size(2).gt.0.0) then
          ! Fold the source with a Gaussian, minus sign this time
          R%chargedensity = R%chargedensity - &
          & FoldGaussian(temp, GaussX(:,:,2,1), GaussY(:,:,2,1), GaussZ(:,:,2,1),&
          &                                                            nx, ny, nz)
      endif
    endif
    !---------------------------------------------------------------------------
    ! When performing simulations for nuclear pasta, one assumes the entire 
    ! volume is charge neutral: a constant background of electrons floods the 
    ! entire simulation volume. We subtract thi s backrgound here.
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
end subroutine ConstructChargeDensity

function memory_for_densities() result (stor)
  !-----------------------------------------------------------------------------
  ! Return an estimation for the total memory required to store all densities.
  !-----------------------------------------------------------------------------
  integer(KIND=LargeInt) :: stor
  stor = 0
$MEMORY_DENSITIES

end function memory_for_densities

 end module vectors


