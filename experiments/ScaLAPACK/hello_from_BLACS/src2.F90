program helloworld
  use mpi
  implicit none
  integer :: rank, comsize, ierr,ICONTEXT,NPROW,NPCOL, ROWID, COLID
  
  print *,'src2.F90'
  call MPI_Init(ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, comsize, ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  print *,'MPI_COMM_WORLD=',MPI_COMM_WORLD
  print *,'Hello World, from task ', rank, 'of', comsize
  
  CALL BLACS_GET( -1, 0, ICONTEXT )
  print *, 'context=',ICONTEXT
  ! CALL BLACS_PINFO(rank,comsize)
  print *,'Hello World2 from task ', rank, 'of', comsize
  CALL BLACS_GRIDINIT(ICONTEXT, 'R', 3, 2)
  CALL BLACS_GRIDINFO(ICONTEXT, NPROW, NPCOL, ROWID, COLID)

  call MPI_Finalize(ierr)
end program helloworld

! integer :: NPROC,IERR,MPIME,comm

! comm = MPI_COMM_WORLD
! !
! NPROW = 3 ! cartesian direction 0
! NPCOL = 2 ! cartesian direction 1
! ! Get a default BLACS context
! !
! ICONTEXT = comm
! ! Initialize a default BLACS context
! end program src