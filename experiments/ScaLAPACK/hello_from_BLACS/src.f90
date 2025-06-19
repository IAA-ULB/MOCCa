      program hello_from_BLACS
      implicit none
!  local variables
      integer info, nproc, nprow, npcol, myid, myrow, mycol, ctxt, ctxt_sys, ctxt_all

! determine rank of calling process and the size of processor set
      call BLACS_PINFO(myid,nproc)
! get the internal default context
      call BLACS_GET( 0, 0, ctxt_sys )
! Set up a process grid for the process set
      ctxt_all = ctxt_sys
      call BLACS_GRIDINIT( ctxt_all, 'C', nproc, 1)
      call BLACS_BARRIER(ctxt_all,'A')
! Set up a process grid of size 3*2
      ctxt = ctxt_sys
      call BLACS_GRIDINIT( ctxt, 'C', 3, 2)
! All processes not belonging to ctxt jump to the end of the program
      if (ctxt.lt.0) go to 1000
! Get the process coordinates in the grid
      call BLACS_GRIDINFO( ctxt, nprow, npcol, myrow, mycol )
      if (myid.eq.0) write(6,*) 'hello from process        rank       myrow       mycol       nprow       npcol'
      write(6,*) 'hello from process', myid, myrow, mycol, nprow, npcol
! ...
!  now ScaLAPACK or PBLAS procedures can be used
! ...
 1000  continue
! return all BLACS contexts
      call BLACS_EXIT(0)
      stop
      end
