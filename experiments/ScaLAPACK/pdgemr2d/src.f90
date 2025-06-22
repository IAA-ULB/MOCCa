      program prace_tut
!     =================

      implicit none 
      integer i,j

      integer,parameter :: m=7, n=5
      real*8  A0(m,n) ! de hele matrix 
      integer A0_dsc(9) ! zijn descriptor

      integer Al_m,Al_n
      real*8, allocatable :: Al(:,:) ! de lokale delen van A0
      integer Al_dsc(9)
      
      integer myid, nproc
      integer nprow, npcol, myrow, mycol

      integer ctxt, ctxt_sys, ctxt_all, info, ctxt_0
      integer NUMROC, INDXG2L, INDXG2P

      integer llda,jl,il,ipr,ipc
      real*8  matval,Aij
      

      call BLACS_PINFO(myid,nproc)
      call BLACS_GET( 0, 0, ctxt_sys ) ! construct blacs context
      ctxt_all = ctxt_sys           
      call BLACS_GRIDINIT( ctxt_all, 'R', 2, 2) ! initialize blacs context
      call BLACS_GRIDINFO( ctxt_all, nprow,npcol,myrow,mycol)
      print '("[",i0,"/",i0,"] myrow="i0" mycol="i0" nprow="i0" npcol="i0"")', myid, nproc, myrow, mycol, nprow, npcol

    ! Distribute matrix A0 (M x N) from root node to all processes in context ictxt.
    !
    ! call SL_INIT(ICTXT, NPROW, NPCOL)
    ! Create a 1 node context for the full A matrix
    ! for loading matrices
      ctxt_0 = ctxt_sys
      call BLACS_GRIDINIT( ctxt_0,  'R', 1, 1)
    
    
    ! LOAD MATRIX ON ROOT NODE AND CREATE DESC FOR IT
      if (myrow == 0 .and. mycol == 0) then
        ! initialiseer de hele matrix op rank 0
        do i = 1, m
            do j = 1, n
                A0(i,j) = matval(i,j)
            end do
        end do
        llda = NUMROC( m, m, myrow, 0, nprow )
        print '("[",i0,"/",i0,"] myrow="i0" mycol="i0" llda="i0"")', myid, nproc, myrow, mycol, llda
        call DESCINIT( A0_dsc, m, n, m, n, 0, 0, ctxt_0, max(1, llda), info )
      else
        ! zet de context van A0 op -1 voor alle andere ranks dan rank 0.
        A0_dsc(1:9) = 0
        A0_dsc(2) = -1
      end if
    
    ! CREATE DESC FOR DISTRIBUTED MATRIX
    llda = NUMROC( m, 2, myrow, 0, nprow )
    print '("[",i0,"/",i0,"] myrow="i0" mycol="i0" llda="i0"")', myid, nproc, myrow, mycol, llda
    CALL DESCINIT( Al_dsc, m, n, 2, 2, 0, 0, ctxt_all, max(1, llda), info )
    
    ! ! DISTRIBUTE DATA
    write(*,*) "node r=", myrow, "c=", mycol, "m=", m, "n=", n
    call PDGEMR2D( m, n, A0, 1, 1, A0_dsc, Al, 1, 1, Al_dsc, A0_dsc( 2 ) )

1000 continue
      call BLACS_EXIT(0)

      stop
      end

    function matval(i,j)
      real*8 matval, aij
      integer i,j 
      aij = 10.0*i + j
      ! write (*,*) i,j,aij
      matval = aij
    end function