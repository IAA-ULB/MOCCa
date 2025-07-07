! file helloworld.f90
subroutine sayhello(comm)
  use mpi
  use iso_c_binding
  implicit none
  integer(c_int), intent(in) :: comm
  ! integer ::comm
  integer :: rank, size, ierr
  call MPI_Comm_size(comm, size, ierr)
  call MPI_Comm_rank(comm, rank, ierr)
  print *, 'Hello, World! I am process ',rank,' of ',size,'.'
end subroutine sayhello
