      program use_PDSYEV
!     ==================
!  +-------------------------------------------------------------------+
!  |  Example program for using the ScaLAPACK driver routine PDSYEV    |
!  |  for calculating eigenvalues and eigenvectors                     |
!  |  of a real symmetric matrix.                                      |
!  |  The matrix is given as a submatrix sub(A)                        |
!  |  of size nxn of a global distributed matrix A,                    |
!  |  the eigenvalues are returned in a global vector,                 |
!  |  the eigenvectors are stored as columns in a submatrix sub(Z)     |
!  |  of size nxn of a global distributed matrix Z.                    |
!  |                                                                   |
!  |  sub(A)(1:n,1:n) = A(i_A:i_A+n-1,j_A:j_A+n-1)                     |
!  |  sub(Z)(1:n,1:n) = Z(i_Z:i_Z+n-1,j_Z:j_Z+n-1)                     |
!  |  The global matrix A is chosen as                                 |
!  |  A(i,j) = 1./(mati + matj*abs(i-j))                               |
!  +-------------------------------------------------------------------+

      implicit none 

!  parameter constants
      integer           max_al, max_n, max_w
      parameter        (max_al=300000, max_n=2000,  max_w=300000)
!  variables for definition of matrix A 
      real*8 mati,matj
      common/matdat/mati,matj

!  local variables
      real*8    al(max_al), al_s(max_al), zl(max_al), rl(max_al), w(max_n), work(max_w), matval, bval, err
      integer   ctxt, ctxt_sys, ctxt_all, &
     &          nproc, nprow, npcol, myid, myrow, mycol, &
     &          m_A, n_A, m_Z, n_Z, rsrc_A, csrc_A, rsrc_Z, csrc_Z, &
     &          mb_A, nb_A, mb_Z, nb_Z, nb, lel1, lel2, &
     &          m_, n_, rsrc_, csrc_, ra, ca, &
     &          m_al, n_al, m_zl, n_zl,  &
     &          llda, lldz, desc_A(9), desc_Z(9), &
     &          n, ncheck, i_A, j_A, i_Z, j_Z, i_s, &
     &          ipr, ipc, il, jl, i, j, k, info, lwork, mll, nll, &
     &          NUMROC, INDXL2G, INDXG2L, INDXG2P

!  calling BLACS_PINFO(myid,nproc) initializes MPI and returns size and rank
      call BLACS_PINFO(myid,nproc)
      ! print '("[",i0,"/",i0,"]")', myid, nproc
      call BLACS_GET( 0, 0, ctxt_sys )
      ctxt_all = ctxt_sys
      call BLACS_GRIDINIT( ctxt_all, 'C', nproc, 1)
      ! Assigns available processes into BLACS process grid.
      ! Create a nproc x 1 process grid
      ! Use Column major ordering (as use by default by Fortran)

! input parameters
! Parameters for distributed global matrices A, Z
! Restrictions are:
!   mb_A = nb_A = mb_Z = nb_Z
!   m_A = m_Z, n_A = n_Z
      m_ =10; n_ = 10; nb = 4; rsrc_ = 0; csrc_=0
      m_A=m_; n_A=n_; mb_A=nb; nb_A=nb
      m_Z=m_; n_Z=n_; mb_Z=nb; nb_Z=nb
      rsrc_A=rsrc_; csrc_A=csrc_; rsrc_Z=rsrc_; csrc_Z=csrc_
      
! Parameters for submatrices sub(A) and sub(Z)
! Restrictions are:
!   i_A+n-1 <= m_A, j_A+n-1<=n_A, i_Z+n-1<=m_Z, j_Z+n-1<=n_Z
!   mod(j_A-1,nb) =  mod(i_A-1,nb) = mod(i_Z-1,nb) = 0
!   the process row containing the row i_A of A
!   must also contain the row i_Z of Z
!   The documenting comments in the code for PDSYEV do not mention
!   the following restriction, which, as experience shows, has to be met:
!    i_A = j_A = i_Z = j_Z
      n = 4; i_s = 1 
      i_A=1+i_s*nprow*nb; j_A=1+i_s*nb*nprow
      i_Z=1+i_s*nprow*nb; j_Z=1+i_s*nb*nprow
! Parameters for processor grid
      nprow=1; npcol=1
! Set up a process grid of size nprow*npcol
      if (nprow*npcol.gt.nproc) then
         write(6,*) 'nproc = ',nproc,' less then nprow*npcol = '
      end if
      ctxt = ctxt_sys; call BLACS_GRIDINIT( ctxt, 'C', nprow, npcol)
! Processes not belonging to the grid jump to the end of program
      if (ctxt.lt.0) go to 1000
! Get the process coordinates in the grid
      call BLACS_GRIDINFO( ctxt, nprow, npcol, myrow, mycol )
! number of rows and columns of local parts al, zl for A, Z
      m_al = NUMROC( m_A, mb_A, myrow, rsrc_A, nprow )
      n_al = NUMROC( n_A, nb_A, mycol, csrc_A, npcol )
      m_zl = NUMROC( m_Z, mb_Z, myrow, rsrc_Z, nprow )
      n_zl = NUMROC( n_Z, nb_Z, mycol, csrc_Z, npcol )
! size of workarray
      lel1 = nb**2*(n/(nb*nprow)+n/(nb*npcol)+3)
      lel2 = nb*((n-1)/(nb*nprow*npcol)+1)
      lwork = n*5 + max(2*n,lel1) + n*lel2 + 1

! Test for sufficient memory
      if (m_al*n_al.gt.max_al .or. &
     &    m_zl*n_zl.gt.max_al .or. &
     &            n.gt.max_n  .or. &
     &        lwork.gt.max_w) then
       write(6,*)'not enough memory:  al ',m_al*n_al,max_al, &
     & ' zl', m_zl*n_zl,max_al,' n',n,max_n,' lwork',lwork,max_w 
       end if
       
! set data for A and B; A(i,j) = 1./(mati+matj*(abs(i-j))
      mati = 1.; matj = 5.

! initializing descriptors for the distributed matrices A, B:
      llda = max(1,m_al); lldz = max(1,m_zl)
      call DESCINIT( desc_A, m_A, n_A, mb_A, nb_A, rsrc_A, csrc_A, ctxt, llda, info )
      call DESCINIT( desc_Z, m_Z, n_Z, mb_Z, nb_Z, rsrc_Z, csrc_Z, ctxt, lldz, info )

! initialize in parallel the local parts of A in al and in al_s
      do jl = 1 , n_al
         j = INDXL2G(jl,nb_A,mycol,csrc_A,npcol)
         do il = 1 , m_al
            i = INDXL2G(il,mb_A,myrow,rsrc_A,nprow)
            al(il+m_al*(jl-1))= matval(i,j)
            al_s(il+m_al*(jl-1))= matval(i,j)
         end do
      end do

! calculate eigenvalues and eigenvectors 
! setting lw = -1, no calculation is performed, only the needed
! workspace is returned in work(1)
      call PDSYEV( 'V', 'U', n, al, i_A, j_A, desc_A, w, zl, i_Z, j_Z, desc_Z, work, lwork, info )
!      call dgamx2d( ctxt, 'A', ' ', 1, 1, work(1), 1, ra, ca, -1, 0, 0 )
!       if (myid.eq.0) 
!     &    write(6,*)int(work(1)), lwork

! Inspect the elements of the distributed matrix z with eigenvectors
! on each process
!      call inspect(zl,'zl   ',desc_Z)

! check the result
! matrix vector multiplication of coefficient matrix in al_s
! times eigenvector k, result is stored in global array R with 
! descriptor desc_Z
      err = 0.d0
      do k = 1 , n
         call PDGEMV( 'N',n,n,  1.d0,al_s,i_A,j_A,desc_A,zl,i_Z,j_Z+k-1,desc_Z,1,0.d0,rl,i_Z,j_Z+k-1,desc_Z,1, )
         do i = 1 , n
            jl = INDXG2L(j_Z+k-1,nb,0,0,npcol)
            ipc = INDXG2P(j_Z+k-1,nb,0,csrc_,npcol)
            il = INDXG2L(i_Z+i-1,nb,0,0,nprow)
            ipr = INDXG2P(i_Z+i-1,nb,0,rsrc_,nprow)
            if (ipr.eq.myrow.and.ipc.eq.mycol) then
              if(i.le.ncheck) write(6,*) k,w(k),i, rl(il+lldz*(jl-1)),zl(il+lldz*(jl-1))
              err = err + abs(rl(il+lldz*(jl-1)) - w(k)*zl(il+lldz*(jl-1)))
            end if
         end do
      end do
      call dgsum2d( ctxt, 'A', ' ', 1, 1, err, 1, 0, 0 )
      if (myid.eq.0) then
            write(6,*)'error in eigensolution is',err/real(n**2)
      end if
 1000 continue
      call BLACS_EXIT(0)

      stop
      end

      real*8 function matval(i,j)
      implicit none
      real*8 mati,matj
      common/matdat/mati,matj
      integer i,j
!      matval = i+100*j
      matval = 1./(mati + matj*abs(i-j))
      return
      end
