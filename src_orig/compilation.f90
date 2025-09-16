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

  implicit none (external)

  public

  integer, parameter :: dp = selected_real_kind(15,307)
  integer, parameter :: sp = selected_real_kind(6,37)

end module compilation

