      program prace_tut
!     =================

      implicit none 
      integer i,j

      integer,parameter :: m=8, n=8
      real*8  A_ful(m,n) ! de hele matrix 
      integer A_ful_dsc(9) ! zijn descriptor

      real*8, allocatable :: A_sub(:,:) ! de lokale delen van A0
      integer A_sub_dsc(9)
      integer lwork
      real*8, allocatable :: work(:) 
      real*8, allocatable :: eigenvalues(:)
      real*8, allocatable :: eigenvectors(:,:)
      
      integer myid, nproc
      integer nprow, npcol, myrow, mycol

      integer ctxt, ctxt_sys, ctxt_all, info, ctxt_0
      integer NUMROC, INDXG2L, INDXG2P

      ! integer llda,jl,il,ipr,ipc
      

      call BLACS_PINFO(myid,nproc)
      call BLACS_GET( 0, 0, ctxt_sys ) ! construct blacs context
      ctxt_all = ctxt_sys           
      call BLACS_GRIDINIT( ctxt_all, 'R', 2, 2) ! initialize blacs context
      call BLACS_GRIDINFO( ctxt_all, nprow,npcol,myrow,mycol)
      print '("ctxt_all [",i0,"/",i0,"] -> ("i0","i0")")', myid, nproc, myrow, mycol

      ctxt_0 = ctxt_sys
      call BLACS_GRIDINIT( ctxt_0, 'R', 1, 1)
      call BLACS_GRIDINFO( ctxt_all, nprow,npcol,myrow,mycol)
      print '("ctxt_0   [",i0,"/",i0,"] -> ("i0","i0")")', myid, nproc, myrow, mycol

      ! initialize full matrix on proces (0,0)
      if (myrow == 0 .and. mycol == 0) then
        write (*,*) 'hello', myid
      
        do i = 1,m
          do j = 1,n
            if (i.le.j) then
              A_ful(i,j) = 10.0*i + j
            else
              A_ful(i,j) = 10.0*j + i
            endif
          enddo
        enddo
        do i = 1,m
          do j = 1,n
            write (*,'(f6.0)', advance='no') A_ful(i,j)
          enddo
          write (*,*)
        enddo
        call DESCINIT( A_ful_dsc, m, n, m, n, 0, 0, ctxt_0, n, info )
        write (*,*) 'ok'
      else
        write (*,*) 'hello', myid

        ! zet de context van A0 op -1 voor alle andere ranks dan rank 0.
        A_ful_dsc(1:9) = 0
        A_ful_dsc(2) = -1
      end if
    
      allocate(A_sub(4,4))
      ! CREATE DESC FOR DISTRIBUTED MATRIX
      ! llda = NUMROC( m, 2, myrow, 0, nprow )
      ! print '("[",i0,"/",i0,"] myrow="i0" mycol="i0" llda="i0"")', myid, nproc, myrow, mycol, llda
      CALL DESCINIT( A_sub_dsc, m, n, 2, 2, 0, 0, ctxt_all, 4, info )
      
      ! ! DISTRIBUTE DATA
      ! write(*,*) "node r=", myrow, "c=", mycol, "m=", m, "n=", n
      call PDGEMR2D( m, n, A_ful, 1, 1, A_ful_dsc, A_sub, 1, 1, A_sub_dsc, A_sub_dsc( 2 ) )

      ! print '("[",i0,"/",i0,"] -> ("i0","i0")")', myid, nproc, myrow, mycol
      ! write (*,'(f6.0)', advance='no') A_sub
      do i = 1,4
        write (*, '("[",i0,"/",i0,"] -> ("i0","i0")")',advance='no'), myid, nproc, myrow, mycol
        do j = 1,4
          write (*,'(f6.0)', advance='no') A_sub(i,j)
        enddo
        write (*,*)
      enddo

      allocate(eigenvalues(n))
      allocate(eigenvectors(n,n))
      allocate(work(1))
      CALL PDSYEV ('V','L',n,A_sub,1,1,A_sub_dsc,eigenvalues,eigenvectors,1,1,A_sub_dsc, work,-1,info)

      lwork=int(work(1))
      deallocate(work)
      allocate(work(lwork))
      ! .... and now do the actual work
      CALL PDSYEV ('V','L',n,A_sub,1,1,A_sub_dsc,eigenvalues,eigenvectors,1,1,A_sub_dsc, work,lwork,info)

      write (*,*) 'A_sub', myid, eigenvalues
      
      if (myrow == 0 .and. mycol == 0) then
        CALL PDSYEV ('V','L',n,A_ful,1,1,A_ful_dsc,eigenvalues,eigenvectors,1,1,A_ful_dsc, work,-1,info)

        lwork=int(work(1))
        deallocate(work)
        allocate(work(lwork))
        ! .... and now do the actual work
        CALL PDSYEV ('V','L',n,A_ful,1,1,A_ful_dsc,eigenvalues,eigenvectors,1,1,A_ful_dsc, work,lwork,info)
        write (*,*) 'A_ful', myid, eigenvalues

      endif

1000 continue
      call BLACS_EXIT(0)

      stop
      end
