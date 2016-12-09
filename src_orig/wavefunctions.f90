module wavefunctions
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
 !
 ! Module containing the single-particle wave-functions (spwfs for short)
 ! for the Tantalus program.
 !
 !
 !
 !
 !=======================================================================
 
 use compilation
 use derivatives
 
 implicit none
 
 !-----------------------------------------------------------------------
 ! Array containing the values of the spwfs in the Hartree-Fock basis
 ! Dimensions (nx,ny,nz,4,nwt)
 real(KIND=dp), allocatable :: HFBasis(:,:,:,:,:) 
 
 !-----------------------------------------------------------------------
 ! Number of the blocks with the same quantum numbers that divide up the 
 ! HFBasis. Any possibility has a maximum of two spatial operators that 
 ! introduce a quantum number, while proton-neutron symmetry adds another one. 
 ! The number of blocks thus needs to be decided on compile time by a
 ! replacement script. 
 ! HFBlocks contains the sizes of the various blocks.  
 !-----------------------------------------------------------------------
 ! Examples:
 !    * EV8-like calculation:     8 blocks (P,Rz,T3)
 !    * EV8-like calculation:     4 blocks (  Rz,T3)
 !    * Rx broken, Sx conserved : 4 blocks (  Sx,T3)
 integer, parameter   :: Blocks          =$BLOCKS
 integer              :: HFBlocks(Blocks)=0

contains 

 subroutine iniwavefunctions   
    !--------------------------------------------------------------------
    ! Placeholder subroutine: just allocates the HFBasis and fills it with
    ! random numbers.
    !---------------------------------------------------------------------
    
    allocate(HFBasis(nx,ny,nz,4,nwt))
    call RANDOM_NUMBER(HFBasis)
 end subroutine iniwavefunctions
 
end module wavefunctions
