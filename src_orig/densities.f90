module densities
!===============================================================================
!  #######   ##   #    # #####   ##   #      #    #  ####
!     #     #  #  ##   #   #    #  #  #      #    # #
!     #    #    # # #  #   #   #    # #      #    #  ####
!     #    ###### #  # #   #   ###### #      #    #      #
!     #    #    # #   ##   #   #    # #      #    # #    #
!     #    #    # #    #   #   #    # ######  ####   ####
!
!  Copyright W. Ryssens & M. Bender
!
!===============================================================================
!
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
!===============================================================================
use compilation
use geninfo
use wavefunctions, only: HFPsi, HFdPsi, HFddPsi, HFdddpsi
use wavefunctions, only: occupations, HFBlocks, blocks
use pairing
use derivatives 
use preconditioning 

implicit none

    !---------------------------------------------------------------------------
    ! Type declaration of the various densities
$DECLARATION   

    !---------------------------------------------------------------------------
    ! Density-mixing parameter default value.
    ! This can be set in the scfiteration namelist in the scfiteration model.
    real(KIND=dp) :: denmix = 0.75_dp
    !---------------------------------------------------------------------------
    ! Type of density mixing to perform. (Default = None)
    ! This can be set in the scfiteration namelist in the scfiteration model.
    integer       :: densitymixing = 0
    !---------------------------------------------------------------------------
    ! Previous value(s) of the density rho.
    real(KIND=dp), allocatable :: D_I_I_hist(:,:,:)
    !---------------------------------------------------------------------------
    ! The amount of iterations to keep in memory for the density mixing
    integer           :: memory = 1
    
    !---------------------------------------------------------------------------
    ! Pointer to which basis is supposed to be used to calculate the densities
    ! Based on pairingtype
    !  (0) HF  => use the HF basis
    !  (1) BCS => use the HF basis
    !  (2) HFB => Use the canonical basis
    real(KIND=dp), pointer ::      DenPsi(:,:,:)
    real(KIND=dp), pointer ::   DendPsi(:,:,:,:)
    real(KIND=dp), pointer ::  DenddPsi(:,:,:,:)
    real(KIND=dp), pointer :: DendddPsi(:,:,:,:)
    !---------------------------------------------------------------------------
    
contains

subroutine densit(SaveRho)
    !---------------------------------------------------------------------------
    ! Calculate all of the densities. 
    ! If SaveRho=.false., do not save the previous values to history!
    !---------------------------------------------------------------------------
    integer      :: i, it, wave
    real(KIND=dp):: weight
    logical      :: SaveRho
    
    select case(PairingType)
    case(0,1)
      ! HF or BCS Calculation
      DenPsi => HFPsi ; DenDPsi => HFDPsi ; DenddPsi => HFddPsi 
      DendddPsi => HFdddpsi
    case(2)
      DenPsi => CanPsi ; DenDPsi => CanDPsi ; DenddPsi => CanddPsi 
      DendddPsi => Candddpsi
    end select
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Allocation and initialization
$INITIALIZATION

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Save old density for next iteration and mixing.
    ! Note that this is only necessary at the moment for the ordinary rho
    ! density, it is the one that can make calculations unstable.
    if(.not. allocated(D_I_I_hist)) then
        allocate(D_I_I_hist(nx*ny*nz,2,memory)) ; D_I_I_hist = 0.0_dp
    endif   
    if(SaveRho) then
      do i=1,memory-1
          D_I_I_hist(:,:,memory-i+1) = D_I_I_hist(:,:,memory-i)
      enddo
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
      ! Saving the input density for mixing    
      D_I_I_hist(:,:,1) = D_I_I
    endif
    
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Zero the current density
$ZEROING

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! PARTICLE-HOLE DENSITIES
    do wave=1,nwt
        ! Isospin is neutron in the first half of blocks, proton in the rest
        it = 2
        if(wave.le.sum(HFBlocks(1:Blocks/2))) it = 1
        
        ! For ordinary densities
        weight  = occupations(wave) 

        do i=1,mv
$EXPRESSION
        enddo
    enddo
    
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! PAIRING DENSITIES
    if(PairingType.eq. 0) then
        ! Make sure the cutoffs are calculated
        do wave=1,nwt
            ! Isospin is neutron in the first half of blocks, proton in the rest
            it = 2
            if(wave.le.sum(HFBlocks(1:Blocks/2))) it = 1
            
            ! For ordinary densities
            ! Currently only suitable for BCS pairing with T conserved
            weight  = kappa_pairing(wave,wave) * Pcutoffs(wave)**2
           
            do i=1,mv
    $PAIREXPRESSION
            enddo
        enddo
    endif

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! The asked for mixing+preconditioning scheme.
    call MassageDensity()
    
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Calculation of the 'derived' densities, densities obtainable by 
    ! deriving other ones. 
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    do it=1,2
$DERIVATION  
    enddo  
end subroutine densit

subroutine MassageDensity()
    !---------------------------------------------------------------------------
    ! Operate on the density before feeding it into the rest of the program.
    !---------------------------------------------------------------------------
    real(KIND=dp), target :: resid(nx*ny*nz,2)
    
    if(all(D_I_I_hist.eq.0.0)) return
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Compute the residual
    resid = D_I_I - D_I_I_hist(:,:,1)

    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Perform mixing
    select case(densitymixing) 
    case(0)
        !-----------------------------------------------------------------------
        ! Precondition the potentials instead of the densities. 
        ! So do nothing to the densities.
    case(1)
        !-----------------------------------------------------------------------
        ! Simple linear mixing at the moment.
        D_I_I = D_I_I_hist(:,:,1) + (1-denmix) * resid
    end select
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Safeguard
    where(D_I_I.lt.1d-10) D_I_I = 0 
end subroutine MassageDensity
    
end module densities
