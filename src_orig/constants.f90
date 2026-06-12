!===============================================================================
!     __  __  ___   ____ ____
!    |  \/  |/ _ \ / ___/ ___|__ _
!    | |\/| | | | | |  | |   / _` |
!    | |  | | |_| | |__| |__| (_| |
!    |_|  |_|\___/ \____\____\__,_|
!
! Written mainly by W. Ryssens & M. Bender
!
! Opensource software distributed under the GNU AGPLv3 licence, see the
!  LICENCE file in the root of this project.
!===============================================================================
module constants
    use compilation

    implicit none
    !---------------------------------------------------------------------------
    ! Physical constants.
    real(KIND=dp):: e2             =  1.43996446_dp
    real(KIND=dp):: clum           = 29.9792458_dp
    real(KIND=dp):: nucleonmass(2) = (/939.565379_dp , 938.272046_dp /)
    
contains 

end module 
