      program prace_tut
!     =================
!     

      implicit none 
      integer i,j
      real*8 Ag(5,5)
      integer myid,nproc
      integer ctxt, ctxt_sys, ctxt_all
      integer loc_m,loc_n
      integer NUMROC, INDXG2L, INDXG2P
      integer nprow, npcol, myrow, mycol
      real*8, allocatable :: Al(:,:)
      integer llda,jl,il,ipr,ipc
      real*8  matval,Aij

      ! Initialize A
      do i=1,5
            do j=1,5
                  Ag(i,j) = matval(i,j)
            end do          
      end do          
      
      call BLACS_PINFO(myid,nproc)
      
      if (myid.eq.0) then
            print '("[",i0,"/",i0,"] Ag:")', myid, nproc
            do i=1,5
                  print '(f4.0,f4.0,f4.0,f4.0,f4.0)', Ag(i,1),Ag(i,2),Ag(i,3),Ag(i,4),Ag(i,5)
            end do          
      end if 
      call BLACS_GET( 0, 0, ctxt_sys )
      ctxt_all = ctxt_sys
      call BLACS_GRIDINIT( ctxt_all, 'C', 2, 2)
      call BLACS_GRIDINFO( ctxt_all, nprow,npcol,myrow,mycol)
      print '("[",i0,"/",i0,"] myrow="i0" mycol="i0"")', myid, nproc, myrow, mycol

      ! allocate local matrix Al
      loc_m = NUMROC(5,2,myrow,0,nprow)
      loc_n = NUMROC(5,2,mycol,0,npcol)
      allocate(Al(loc_m,loc_n))
      print '("[",i0,"/",i0,"] loc_m="i0" loc_n="i0"")', myid, nproc, loc_m, loc_n

      ! copy the values of the global matrix Ag
      ! to the local parts of the distributed matrix
      llda =  max(1,NUMROC( 5, 2, myrow, 0, nprow ))
      print '("[",i0,"/",i0,"] llda="i0"")', myid, nproc, llda

      do i = 1,5
         il  = INDXG2L(i,2,0,0,nprow)
         ipr = INDXG2P(i,2,0,0,nprow)
         do j = 1,5
            jl  = INDXG2L(j,2,0,0,npcol)
            ipc = INDXG2P(j,2,0,0,npcol)
            ! print '("[",i0,"/",i0,"] j="i0" jl="i0" ipc="i0"")', myid, nproc, j, jl, ipc
            ! print '("[",i0,"/",i0,"] i="i0" il="i0" ipr="i0"")', myid, nproc, i, il, ipr
            Aij = matval(i,j) 
            ! write (*,*) Aij
            if (ipr.eq.myrow.and.ipc.eq.mycol) then
                  print '("[",i0,"/",i0,"] i,j=("i0","i0") il,jl=("i0","i0") Aij="f6.0"")', myid, nproc, i,j,il,jl,Aij
                  Al(il,jl) = Aij
            else
                  print '("[",i0,"/",i0,"] i,j=("i0","i0") il,jl=("i0","i0") Aij=------")', myid, nproc, i,j,il,jl                  
            end if
         end do
      end do


      print '("[",i0,"/",i0,"] Al")', myid, nproc      
      do il = 1,loc_m
            do jl = 1,loc_n
                  ! print '("[",i0,"/",i0,"] Al("i0","i0")="f6.0"")', myid, nproc, il,jl, Al(il,jl)

                  write (*,'(f6.0)', advance="no") Al(il,jl)
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