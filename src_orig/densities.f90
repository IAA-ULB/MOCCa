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
! a) Unlike most stuff, the densities are represented as vectors on the mesh, 
!    instead of 3D arrays. 
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
!===============================================================================
use compilation
use geninfo
use wavefunctions, only: HFPsi, HFdPsi, HFddPsi, HFdddpsi
use wavefunctions, only: occupations, HFBlocks, blocks
use wavefunctions, only: upairing, vpairing
use derivatives  

implicit none

    !---------------------------------------------------------------------------
    ! Type declaration of the various densities
$DECLARATION   

    !---------------------------------------------------------------------------
    ! Density-mixing parameter
    real(KIND=dp) :: denmix = 0.65_dp
    
contains

subroutine readdensit


    namelist /densit/ denmix

    read(unit=*, nml = densit)

end subroutine readdensit

subroutine densit(iteration)
    !---------------------------------------------------------------------------
    ! Calculate all of the densities
    !---------------------------------------------------------------------------
    integer      :: i, it, wave
    integer, intent(in) :: iteration
    real(KIND=dp):: weight
    
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Allocation and initialization
$INITIALIZATION
    
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Calculation by summing the densities
    do wave=1,nwt
        ! Isospin is neutron in the first half of blocks, proton in the rest
        it = 2
        if(wave.le.sum(HFBlocks(1:Blocks/2))) it = 1
        weight = occupations(wave)
        if( iteration.ne.0)  then
             weight = weight * (1-denmix)
        endif
        do i=1,mv
$EXPRESSION
        enddo
    enddo

    ! Calculation by summing of the pairing densities
    do wave=1,nwt
        ! Isospin is neutron in the first half of blocks, proton in the rest
        it = 2
        if(wave.le.sum(HFBlocks(1:Blocks/2))) it = 1
        weight = upairing(wave) * vpairing(wave)
        if( iteration.ne.0)  then
             weight = weight * (1-denmix)
        endif
        do i=1,mv
!   PAIRINGEXPRESSION
        enddo
    enddo       

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Calculation of the 'derived' densities, densities obtainable by 
    ! deriving other ones. 
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    do it=1,2
$DERIVATION  
    enddo  
end subroutine densit
    
end module densities
