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
 ! Hephaestos:
 ! BLOCKS $BLOCKS
 !=======================================================================
 use compilation
 use derivatives
 use nil8
 
 implicit none
 
 !-----------------------------------------------------------------------
 ! Array containing the values of the spwfs in the Hartree-Fock basis
 ! and their derivatives
 ! Dimensions (nx,ny,nz,4,nwt)
 real(KIND=dp), allocatable :: HFBasis(:,:,:,:,:)
 real(KIND=dp), allocatable :: HFDeriv(:,:,:,:,:,:) ! Derivatives
 real(KIND=dp), allocatable :: HFLapla(:,:,:,:,:)   ! Laplacian
 !-----------------------------------------------------------------------
 ! Density matrix rho and anomalous density matrix kappa
 ! Dimensions (nwt, nwt) (although many are zero when symmetries are conserved)
 !real(KIND=dp), allocatable :: rho(:,:), kappa(:,:)
 !-----------------------------------------------------------------------
 ! Occupations of the single-particle wave-functions, i.e. the eigenvalues
 ! of rho. 
 real(KIND=dp), allocatable :: occupations(:)
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
 !    * EV4-like calculation:     4 blocks (  Rz,T3)
 !    * Rx broken, Sx conserved : 4 blocks (  Sx,T3)
 integer, parameter   :: Blocks          =$BLOCKS
 integer              :: HFBlocks(Blocks)=0

contains 

  subroutine iniwavefunctions()   
    !--------------------------------------------------------------------
    !
    !
    !--------------------------------------------------------------------   
    
    real(KIND=dp)        :: homegax, homegay,homegaz, alpha,qqq
    integer              :: meven = 5, modd = 4, i
    integer, allocatable :: kparz(:)
    !--------------------------------------------------------------------
    ! Build harmonic oscillator eigenfunctions by constructing them in  
    ! an EV8-like box and then expanding them to the entire box. 
    !--------------------------------------------------------------------
        
    alpha = 0.2    
    qqq   = 1.0    
    homegaz  = alpha*qqq**(-2.0/3.0)
    homegax  = alpha*qqq**(-2*cos(-2*pi/3)/3)
    homegay  = alpha*qqq**(-2*cos(+2*pi/3)/3)
    
    call nilsson (HFBasis,kparz,meven,modd,20,10,10, 10, 10,nx,ny,nz,0.8d0,homegax,homegay,homegaz)
    
    do i=1,10
        if(kparz(i) .gt. 0) HFBlocks(1) = HFBlocks(1) +1
        if(kparz(i) .lt. 0) HFBlocks(3) = HFBlocks(3) +1
    enddo
    do i=11,20
        if(kparz(i) .gt. 0) HFBlocks(5) = HFBlocks(5) +1
        if(kparz(i) .lt. 0) HFBlocks(7) = HFBlocks(7) +1
    enddo
  end subroutine iniwavefunctions
 
end module wavefunctions
