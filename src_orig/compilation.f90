module compilation
  !-----------------------------------------------------------------------------
  ! Module to propagate 
  !
  ! a) the definition of a double and single precision real
  ! b) MPI utilities (i.e. the USE MPI statement)
  !
  ! throughout the entire code.
  !-----------------------------------------------------------------------------

#if(USE_MPI > 0)
  use MPI
#endif
  
  implicit none

  integer, parameter :: dp = selected_real_kind(15,307) 
  integer, parameter :: sp = selected_real_kind(6,37)  
  integer, parameter :: LargeInt = selected_int_kind (12)
  ! Large integers are required for the calculation of memory requirements

end module compilation

