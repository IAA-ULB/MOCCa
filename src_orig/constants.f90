module constants
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
! Module governing the constant that go into a mean-field program. 
! Static (values of SLy4) for the moment, but will get changeable at 
! runtime rather quickly. 
!=======================================================================
    use compilation

    implicit none

    real(KIND=dp) :: t0=-2488.913 
    real(KIND=dp) :: x0=0.834
    real(KIND=dp) :: t1=486.818
    real(KIND=dp) :: x1=-0.344
    real(KIND=dp) :: t2=-546.395
    real(KIND=dp) :: x2=-1.0
    real(KIND=dp) :: t3a=13777.0
    real(KIND=dp) :: x3a=1.354
    real(KIND=dp) :: yt3a=0.166666666666666666667 
    real(KIND=dp) :: t3b=0.0 
    real(KIND=dp) :: x3b=0.0
    real(KIND=dp) :: yt3b=0
    real(KIND=dp) :: te=0.0
    real(KIND=dp) :: to=0.0
    real(KIND=dp) :: wso=123.0
    real(KIND=dp) :: wsoq=123.0

    !Functional parameters in the BFH representation
    real(KIND=dp) :: B1,B2,B3,B4,B5,B6,B7,B8,Byt3,B9,B9q
    real(KIND=dp) :: B10,B11,B12,B13,B14
    real(KIND=dp) :: B15,B16, B17, B18,B19,B20,B21
    !------------------------------------------------------------------------------
    ! Physical constants.
    real(KIND=dp):: e2             =  1.43996446_dp
    real(KIND=dp):: hbar           =  6.58211928_dp
    real(KIND=dp):: clum           = 29.9792458_dp
    real(KIND=dp):: nucleonmass(2) = (/939.565379_dp , 938.272046_dp /)
    real(KIND=dp):: hbm(2)         = 20.73551910_dp

contains 

    subroutine calcedfcoefs()
        real(KIND=dp), parameter :: rhosat=0.16_dp
        !---------------------------------------------------------------------------
        !Calculating the B-s
        !---------------------------------------------------------------------------
        B1   = t0*(1.0_dp   + 1/2.0_dp*x0)/2.0_dp
        B2   =-t0*(1/2.0_dp + x0         )/2.0_dp

        B3   = (t1*( 1      + 1/2.0_dp*x1) + t2*(1      + 1/2.0_dp*x2))/4.0_dp
        B4   =-(t1*( 1/2.0_dp + x1       ) - t2*(1/2.0_dp +        x2))/4.0_dp

        B5   =-(3.0_dp*t1*(1  + x1/2.0_dp)-t2*(1  + x2/2.0_dp))/16.0_dp
        B6   = (3.0_dp*t1*(x1 +  1/2.0_dp)+t2*(x2 +  1/2.0_dp))/16.0_dp

        B7  = t3a*( 1   + x3a/2.0_dp)/12.0_dp
        B8  =-t3a*( x3a +   1/2.0_dp)/12.0_dp
        !---------------------------------------------------------------------------
        ! B9 and B9q are the constants related to the spin-orbit interaction.
        ! Starting from a Skyrme Force, they should in principle be equal for a
        ! true force, but apparently there are some physics cases where a difference
        ! between the isoscalar and isovector coupling is wanted.
        !---------------------------------------------------------------------------
        B9   = -wso /2.0_dp
        B9q  = -wsoq/2.0_dp

        Byt3= yt3a
        
        return
    end subroutine calcedfcoefs
    
    

end module 
