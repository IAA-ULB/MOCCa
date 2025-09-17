module compilation
  !-----------------------------------------------------------------------------
  ! Module to propagat the definition of a double and single precision real
  ! MPI utilities (i.e. the USE MPI statement) throughout the entire code.
  !-----------------------------------------------------------------------------


  implicit none (external)

  public

  integer, parameter :: dp = selected_real_kind(15,307)
  integer, parameter :: sp = selected_real_kind(6,37)

end module compilation

