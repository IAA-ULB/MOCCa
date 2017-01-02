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
    !---------------------------------------------------------------
    ! Calculate all of the densities
    !---------------------------------------------------------------
    integer :: i, j, k, it
    
$INITIALIZATION
    
    do j=1,nwt
        ! Isospin is neutron in the first half of blocks, proton in the rest
        it = 2
        if(j.le.sum(HFBlocks(1:Blocks/2))) it = 1
        
        do k=1,nz
            do j=1,ny
                do i=1,nx
$EXPRESSION
                enddo
            enddo
        enddo
    enddo

    print *, sum(rho)*dv*8
end subroutine densit
    
end module densities
