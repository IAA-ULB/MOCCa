! #define __MPI__ 1

program hello_from_BLACS
  use mpi
  implicit none
      
! local variables
  integer info, nproc, nprow, npcol, myid, myrow, mycol, ctxt, ctxt_sys, ctxt_all,comm, size, rank

  print *,'src.F90'
  
#ifdef __MPI__
  write(*,*) "__MPI__"
  call mpi_init(info)
  call MPI_Comm_size(MPI_COMM_WORLD, size, info)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, info)
  print *, 'Hello, World! I am process ',rank,' of ',size,'MPI_COMM_WORLD=',MPI_COMM_WORLD
  comm = MPI_COMM_WORLD
!   ctxt_sys = comm
! ??? Can we call BLACS_PINFO after mpi_init ???
  write(*,*) "BLACS_PINFO after mpi_init"
! determine rank of calling process and the size of processor set
  call BLACS_PINFO(rank,size)
#else
  write(*,*) "BLACS_PINFO only"
! determine rank of calling process and the size of processor set
  call BLACS_PINFO(rank,size)
  print *, 'Hello, World2 I am process ',rank,' of ',size,'.'
#endif

! get the internal default context  
  call BLACS_GET( -1, 0, ctxt_sys )
  print *, 'rank',rank,' of ',size,'BLACS context=',ctxt_sys,'MPI_COMM_WORLD=',MPI_COMM_WORLD

! Set up a process grid for the process set
  ctxt_all = ctxt_sys
  call BLACS_GRIDINIT( ctxt_all, 'C', size, 1)
  call BLACS_BARRIER(ctxt_all,'A')
! Set up a process grid of size 3*2
  ctxt = ctxt_sys
  call BLACS_GRIDINIT( ctxt, 'C', 3, 2)
! All processes not belonging to ctxt jump to the end of the program
  if (ctxt.ge.0) then
    ! Get the process coordinates in the grid
    call BLACS_GRIDINFO( ctxt, nprow, npcol, myrow, mycol )
    if (myid.eq.0) write(6,*) 'hello from process        rank       myrow       mycol       nprow       npcol'
    write(6,*) 'hello from process', rank, myrow, mycol, nprow, npcol
    ! ...
    !  now ScaLAPACK or PBLAS procedures can be used
    ! ...
  endif
  call BLACS_EXIT(0)
end
