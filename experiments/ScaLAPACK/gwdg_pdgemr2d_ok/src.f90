      program prace_tut
!     =================
!     

      implicit none 
      integer i,j
      integer,parameter :: m=5,n=7        ! aantal rijen/kolommen in de volledige matrix A_full
      integer,parameter :: mb=2,nb=2      ! blocking factor voor rijen/kolommen
      real*8, allocatable :: A_full(:,:)  ! de volledige matrix
      real*8, allocatable :: A_sub(:,:)   ! submatrix, lokaal
      integer ms,ns                       ! dimensies van de submatrix
      integer A_full_dsc(9), A_sub_dsc(9) ! descriptors

      integer myid,nproc                  ! rank/aantal processen
      integer nprow, npcol                ! dimensies van het processor grid
      integer myrow, mycol                ! coordinaten van deze rank in het processor grid

      integer llda,is,js,ipr,ipc

      integer ctxt, ctxt_sys, ctxt_all, info, ctxt_0
      integer NUMROC, INDXG2L, INDXG2P
      real*8  matval,Aij
      

      call BLACS_PINFO(myid,nproc)
      call BLACS_GET( 0, 0, ctxt_sys )
      ctxt_all = ctxt_sys
      call BLACS_GRIDINIT( ctxt_all, 'C', 2, 2)
      call BLACS_GRIDINFO( ctxt_all, nprow,npcol,myrow,mycol)
      print '("[",i0,"/",i0,"] myrow="i0" mycol="i0" nprow="i0" npcol="i0"")', myid, nproc, myrow, mycol, nprow, npcol
      
      call SL_INIT(ctxt_0, 1, 1) ! create 1 node context for loading matrices
      ! descriptor for A_full§
      if (myid.eq.0) then
            llda = NUMROC( m, m, myrow, 0, nprow ) ! blocking factor equals dimension, because it is not distributed. llda = 7
            ! print '("llda=",i0)', llda
            call DESCINIT( A_full_dsc, m, n, mb, mb, 0, 0, ctxt_0, max(1, llda), info )
            ! initialize the full matrix
            allocate(A_full(m,n))
            do i = 1,m 
               do j = 1,n 
                  A_full(i,j) = matval(i,j)
               enddo
            enddo
            do i = 1,m
               print '(f0.0," ",f0.0," ",f0.0," ",f0.0," ",f0.0," ",f0.0," ",f0.0)', A_full(i,1), A_full(i,2),A_full(i,3),A_full(i,4),A_full(i,5),A_full(i,6),A_full(i,7)
            enddo
                  
      else
            A_full_dsc(1:9) = 0
            A_full_dsc(2)   = -1
      endif

      ! allocate local matrix A_sub
      ms = NUMROC(m,2,myrow,0,nprow)
      ns = NUMROC(n,2,mycol,0,npcol)
      allocate(A_sub(ms,ns))
      print '("rank[",i0,"("i0","i0")/",i0,"] A_sub: "i0"x"i0" (allocated)")', myid, myrow, mycol, nproc, ms, ns
      ! initialize the descriptor for the local matrix:
      call DESCINIT( A_sub_dsc, m, n, 2, 2, 0, 0, ctxt, ms, info )
      
      ! copy a_0 on process 0 block-cyclically to the local arrays al in the process grid
      call PDGEMR2D(m, n, A_full, 1, 1, A_full_dsc, A_sub, 1, 1, A_sub_dsc, A_sub_dsc(2))
      ! end if      
      
      print '("[",i0,"/",i0,"]=("i0","i0") A_sub:")', myid, nproc, myrow,mycol      
      do is = 1,ms
            do js = 1,ns
                  write (*,'(f6.0)', advance="no") A_sub(is,js)
            end do
            write (*,*)
      end do

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