module densities
!=======================================================================
!  #######   ##   #    # #####   ##   #      #    #  ####
!     #     #  #  ##   #   #    #  #  #      #    # #
!     #    #    # # #  #   #   #    # #      #    #  ####
!     #    ###### #  # #   #   ###### #      #    #      #
!     #    #    # #   ##   #   #    # #      #    # #    #
!     #    #    # #    #   #   #    # ######  ####   ####
!
!  Copyright W. Ryssens & M. Bender
!
!=======================================================================

use compilation
use wavefunctions

implicit none

    !-------------------------------------------------------------------
    ! Type declaration of the various densities
$DECLARATION   

contains

subroutine densit
    !---------------------------------------------------------------------------
    ! Calculate all of the densities
    !---------------------------------------------------------------------------
    integer :: i, j, k, it, wave
    
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Allocation and initialization
$INITIALIZATION
    
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Calculation by summing the densities
    do wave=1,nwt
        ! Isospin is neutron in the first half of blocks, proton in the rest
        it = 2
        if(wave.le.sum(HFBlocks(1:Blocks/2))) it = 1
        do k=1,nz
            do j=1,ny
                do i=1,nx
$EXPRESSION
                enddo
            enddo
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
