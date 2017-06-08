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

!    real(KIND=dp) :: t0=-2488.913 
!    real(KIND=dp) :: x0=0.834
!    real(KIND=dp) :: t1=486.818
!    real(KIND=dp) :: x1=-0.344
!    real(KIND=dp) :: t2=-546.395
!    real(KIND=dp) :: x2=-1.0
!    real(KIND=dp) :: t3a=13777.0
!    real(KIND=dp) :: x3a=1.354
!    real(KIND=dp) :: yt3a=0.166666666666666666667 
!    real(KIND=dp) :: t3b=0.0 
!    real(KIND=dp) :: x3b=0.0
!    real(KIND=dp) :: yt3b=0
!    real(KIND=dp) :: te=0.0
!    real(KIND=dp) :: to=0.0
!    real(KIND=dp) :: wso=123.0
!    real(KIND=dp) :: wsoq=123.0
!    
!    real(KIND=dp) :: t1n2=24.3409
!    real(KIND=dp) :: t2n2=-27.31975
!    real(KIND=dp) :: x1n2=-0.344
!    real(KIND=dp) :: x2n2=-1.0   
    

    !---------------------------------------------------------------------------
    ! Physical constants.
    real(KIND=dp):: e2             =  1.43996446_dp
    real(KIND=dp):: clum           = 29.9792458_dp
    real(KIND=dp):: nucleonmass(2) = (/939.565379_dp , 938.272046_dp /)
    
    
contains 

end module 
