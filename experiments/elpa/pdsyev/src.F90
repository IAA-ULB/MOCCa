program PDGEMR2D_PDSYEV

  use mpi
  use elpa
  
  implicit none 

  class(elpa_t), pointer :: lpa

  integer :: i,j,ms,ns,is,js


  integer,parameter :: n=8, nblk=2
  integer :: A_ful_dsc(9) ! zijn descriptor

  integer :: A_sub_dsc(9)
  integer :: lwork
  real*8 :: A_ful(n,n) ! de hele matrix 
  real*8, allocatable :: A_sub(:,:) ! de lokale delen van A0
  real*8, allocatable :: work(:) 
  real*8, allocatable :: eigenvalues(:)
  real*8, allocatable :: eigenvectors(:,:)
  
  integer rank, nranks 
  integer nprow, npcol, myrow, mycol

  integer ctxt, ctxt_sys, ctxt_all, info, ctxt_0, success
  integer NUMROC, INDXG2L, INDXG2P
    
  call mpi_init(info)
  call mpi_comm_rank(mpi_comm_world,rank,info)
  call mpi_comm_size(mpi_comm_world,nranks,info)  
  print *, "mpi_comm_world=", mpi_comm_world

! Transfer the MPI communicator to blacs
  ctxt_all = mpi_comm_world
  call BLACS_GET( -1, 0, ctxt_all ) ! construct blacs context

  ! ctxt_all = ctxt_sys           
  call BLACS_GRIDINIT( ctxt_all, 'R', 2, 2) ! initialize blacs context
  call BLACS_GRIDINFO( ctxt_all, nprow,npcol,myrow,mycol)
  print '("ctxt_all [",i0,"/",i0,"] -> ("i0","i0")/("i0","i0")")', rank, nranks, myrow, mycol, nprow,npcol

! Transfer the MPI communicator to blacs
  ctxt_0 = mpi_comm_world
  call BLACS_GET( -1, 0, ctxt_0 ) ! construct blacs context
  call BLACS_GRIDINIT( ctxt_0, 'R', 1, 1)
  
  ! initialize full matrix on proces (0,0)
  if (myrow == 0 .and. mycol == 0) then      
    do i = 1,n
      do j = 1,n
        if (i.eq.j) then
            A_ful(i,j) = 1.5
        else 
            A_ful(i,j) = 1.0/sqrt(float(abs(i-j)))
        endif
      enddo
    enddo
    do i = 1,n
      do j = 1,n
        write (*,'(f6.2)', advance='no') A_ful(i,j)
      enddo
      write (*,*)
    enddo
    call DESCINIT( A_ful_dsc, n, n, n, n, 0, 0, ctxt_0, n, info )
    write (*,*) 'ok'
  else
    write (*,*) 'hello', rank

    ! zet de context van A0 op -1 voor alle andere ranks dan rank 0.
    A_ful_dsc(1:9) = 0
    A_ful_dsc(2) = -1
  end if
  
  ms = NUMROC(n,nblk,myrow,0,nprow)
  ns = NUMROC(n,nblk,mycol,0,npcol)
  print *,"ms = NUMROC(", n,nblk,myrow,0,nranks,") = ", ms
  print *,"ns = NUMROC(", n,nblk,mycol,0,nranks,") = ", ns
  print *, n, "ms=",ms," ns",ns
  call DESCINIT( A_sub_dsc, n, n, nblk, nblk, 0, 0, ctxt_all, ms, info )
  
  ! allocate local distributed matrix ad
  allocate(A_sub(ms,ns))
  print '("rank[",i0,"("i0","i0")/",i0,"] ad: "i0"x"i0" (allocated)")', rank, myrow, mycol, nranks, ms, ns

  
  ! ! DISTRIBUTE DATA
  ! write(*,*) "node r=", myrow, "c=", mycol, "n=", n, "n=", n
  call PDGEMR2D( n, n, A_ful, 1, 1, A_ful_dsc, A_sub, 1, 1, A_sub_dsc, A_sub_dsc( 2 ) )

  ! print '("[",i0,"/",i0,"] -> ("i0","i0")")', rank, nranks, myrow, mycol
  ! write (*,'(f6.0)', advance='no') A_sub
  do i = 1,4
    write (*, '("[",i0,"/",i0,"] -> ("i0","i0") row ",i0)',advance='no'), rank, nranks, myrow, mycol, i
    do j = 1,4
      write (*,'(f6.2)', advance='no') A_sub(i,j)
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

  write (*,*) 'pdsyev: eigenvalues', rank, eigenvalues
  
! UP to this point this where all the prerequisites which have to be done in the
! application if you have a distributed eigenvalue problem to be solved, independent of
! whether you want to use ELPA, Scalapack, EigenEXA or alike

! ! DISTRIBUTE DATA
  ! write(*,*) "node r=", myrow, "c=", mycol, "n=", n, "n=", n
  call PDGEMR2D( n, n, A_ful, 1, 1, A_ful_dsc, A_sub, 1, 1, A_sub_dsc, A_sub_dsc( 2 ) )
  do i = 1,4
    write (*, '("[",i0,"/",i0,"] -> ("i0","i0") row ",i0)',advance='no'), rank, nranks, myrow, mycol, i
    do j = 1,4
      write (*,'(f6.2)', advance='no') A_sub(i,j)
    enddo
    write (*,*)
  enddo

! Now you can start using ELPA

!   if (elpa_init(20250131) /= ELPA_OK) then        ! put here the API version that you are using
!     print *, "ELPA API version not supported"
!     stop 1
!   endif

!   lpa => elpa_allocate(success)
!   if (success /= ELPA_OK) then
!     ! react on the error
!     ! we urge every user to always check the error codes
!     ! of all ELPA functions
!   endif

! ! set parameters decribing the matrix and it's MPI distribution
!   call lpa%set("n", n, success)                          ! size of the n x n matrix
!   call lpa%set("nev", n, success)                        ! number of eigenvectors that should be computed ( 1<= nev <= n)
!   call lpa%set("local_nrows", ms, success)            ! number of local rows of the distributed matrix on this MPI task 
!   call lpa%set("local_ncols", ns, success)            ! number of local columns of the distributed matrix on this MPI task
!   call lpa%set("nblk", nblk, success)                      ! size of the BLACS block cyclic distribution
!   call lpa%set("mpi_comm_parent", MPI_COMM_WORLD, success) ! the global MPI communicator
!   call lpa%set("process_row", myrow, success)            ! row coordinate of MPI process
!   call lpa%set("process_col", mycol, success)            ! column coordinate of MPI process

!   success = lpa%setup()
!   if (success.ne.0) then
!     print *,"lpa%setup() error"
!     stop
!   endif

! ! if desired, set any number of tunable run-time options
! ! look at the list of possible options as detailed later in
! ! USERS_GUIDE.md
!   call lpa%set("solver", ELPA_SOLVER_2STAGE, success)
!   if (success.ne.0) then
!     print *,"lpa%set('solver', ELPA_SOLVER_2STAGE, success) error"
!     stop
!   endif

! ! set the AVX BLOCK2 kernel, otherwise ELPA_2STAGE_REAL_DEFAULT will
! ! be used
!   ! call lpa%set("real_kernel", ELPA_2STAGE_REAL_AVX_BLOCK2, success)

! ! use method solve to solve the eigenvalue problem to obtain eigenvalues
! ! and eigenvectors
! ! other possible methods are desribed in USERS_GUIDE.md
!   call lpa%eigenvalues(A_sub, eigenvalues, success)
!   ! call lpa%eigenvectors(A_sub, eigenvalues, eigenvectors, success)
!   print*,success,elpa_strerr(3)
! ! cleanup
!   call elpa_deallocate(lpa)

!   call elpa_uninit()

  call BLACS_EXIT(0)

end program PDGEMR2D_PDSYEV
