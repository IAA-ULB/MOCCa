module compilation 

  implicit none
  
  integer, parameter :: dp = selected_real_kind(15,307) 
  integer, parameter :: sp = selected_real_kind(6,37)  
  !---------------------------------------------------------------------------
  ! Pi is always practical (delicious) to have.
  real(KIND=dp), parameter  :: pi=4.0_dp*atan2(1.0_dp,1.0_dp)

end module compilation 
