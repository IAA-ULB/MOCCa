module compilation
  !-----------------------------------------------------------------------------
  ! Module to propagate 
  ! a) the definition of a double and single precision real
  ! b) MPI utilities
  ! throughout the entire code.
  !-----------------------------------------------------------------------------
  
  implicit none

#if (USE_MPI > 0)
  include 'mpif.h'
#endif

  integer, parameter :: dp = selected_real_kind(15,307) 
  integer, parameter :: sp = selected_real_kind(6,37)  

end module compilation

