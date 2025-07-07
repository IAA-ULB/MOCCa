program src2
use mpi
use elpa

implicit none

class(elpa_t), pointer :: elpa_
integer :: mpierr, rank, nprocs
integer :: npcol, nprow, mycol, myrow
integer :: my_blacs_ctxt, info
integer :: na = 8 ! global dimension of the matrix to be solved
integer :: nblk = 2 ! the block size of the scalapack block cyclic distribution
real*8, allocatable :: a(:,:), ev(:), ad(:,:), z(:,:)
integer :: i,j,ms,ns,is,js
integer :: success
real*8  :: NUMROC
integer :: a_desc(9), ad_desc(9)

!-------------------------------------------------------------------------------
!  MPI Initialization

call mpi_init(mpierr)
call mpi_comm_rank(mpi_comm_world,rank,mpierr)
call mpi_comm_size(mpi_comm_world,nprocs,mpierr)  

!-------------------------------------------------------------------------------
! Selection of number of processor rows/columns
! the application has to decide how the matrix should be distributed
npcol = 2
nprow = 2

!-------------------------------------------------------------------------------
! Set up BLACS context and MPI communicators
!
! The BLACS context is only necessary for using Scalapack.
!
! For ELPA, the MPI communicators along rows/cols are sufficient,
! and the grid setup may be done in an arbitrary way as long as it is
! consistent (i.e. 0<=myrow<nprow, 0<=my_pcol<npcol and every
! process has a unique (myrow,my_pcol) pair).
! For details look at the documentation of  BLACS_Gridinit and
! BLACS_Gridinfo of your BLACS installation

my_blacs_ctxt = mpi_comm_world
call BLACS_Gridinit( my_blacs_ctxt, 'C', nprow, npcol )
call BLACS_Gridinfo( my_blacs_ctxt, nprow, npcol, myrow, mycol )

! compute for your distributed matrix the number of local rows and columns 
! per MPI task, e.g. with
! the Scalapack tools routine NUMROC 

! Set up a scalapack descriptor for the checks below.
! For ELPA the following restrictions hold:
! - block sizes in both directions must be identical (args 4+5)
! - first row and column of the distributed matrix must be on row/col 0/0 (args 6+7)

call descinit( a_desc, na, na, nblk, nblk, 0, 0, my_blacs_ctxt, na, info )
if (info .ne. 0) then
  print *,"Failed: descinit( a_desc, na, na, nblk, nblk, 0, 0, my_blacs_ctxt, na, info )"
  stop 1
endif
if (rank == 0) then
    ! Allocate matrices 
    allocate(a (na,na))
    allocate(ev(na))

    ! fill the matrix with resonable values
    do i=1,na
        do j=1,na
            if (i.eq.j) then
                a(i,j) = 1.5
            else 
                a(i,j) = 1.0/sqrt(float(abs(i-j)))
            endif
            print *,advance-'no' a(i,j)
        enddo
        print *
    enddo
endif

print *,"na   =",na,"nblk =",nblk,"myrow=",myrow,"nprow=",nprow,"mycol=",mycol,"npcol=",npcol,"nprocs=",nprocs
ms = NUMROC(na,nblk,myrow,0,nprocs)
ns = NUMROC(na,nblk,mycol,0,nprocs)
print *,"ms=",ms," ns",ns
call DESCINIT( ad_desc, na, na, 2, 2, 0, 0, my_blacs_ctxt, ms, info )
! call DESCINIT( zd_desc, na, na, 2, 2, 0, 0, my_blacs_ctxt, ms, info )
! allocate local distributed matrix ad
allocate(ad(ms,ns))
allocate(z (ms,ns))
print '("rank[",i0,"("i0","i0")/",i0,"] ad: "i0"x"i0" (allocated)")', rank, myrow, mycol, nprocs, ms, ns
! initialize the descriptor for the loc§al matrix:

! copy a_0 on process 0 block-cyclically to the local arrays al in the process grid
call PDGEMR2D(na, na, a, 1, 1, a_desc, ad, 1, 1, ad_desc, ad_desc(2))
! end if      

print '("[",i0,"/",i0,"]=("i0","i0") ad:")', rank, nprocs, myrow,mycol      
do is = 1,ms
    do js = 1,ns
            write (*,'(f6.0)', advance="no") ad(is,js)
    end do
    write (*,*)
end do


! UP to this point this where all the prerequisites which have to be done in the
! application if you have a distributed eigenvalue problem to be solved, independent of
! whether you want to use ELPA, Scalapack, EigenEXA or alike

! Now you can start using ELPA

if (elpa_init(20250131) /= ELPA_OK) then        ! put here the API version that you are using
   print *, "ELPA API version not supported"
   stop 1
endif

 elpa_ => elpa_allocate(success)
 if (success /= ELPA_OK) then
   ! react on the error
   ! we urge every user to always check the error codes
   ! of all ELPA functions
 endif

 ! set parameters decribing the matrix and it's MPI distribution
 call elpa_%set("na", na, success)                          ! size of the na x na matrix
 call elpa_%set("nev", na, success)                        ! number of eigenvectors that should be computed ( 1<= nev <= na)
 call elpa_%set("local_nrows", ns, success)            ! number of local rows of the distributed matrix on this MPI task 
 call elpa_%set("local_ncols", ns, success)            ! number of local columns of the distributed matrix on this MPI task
 call elpa_%set("nblk", nblk, success)                      ! size of the BLACS block cyclic distribution
 call elpa_%set("mpi_comm_parent", MPI_COMM_WORLD, success) ! the global MPI communicator
 call elpa_%set("process_row", myrow, success)            ! row coordinate of MPI process
 call elpa_%set("process_col", mycol, success)            ! column coordinate of MPI process

 success = elpa_%setup()

 ! if desired, set any number of tunable run-time options
 ! look at the list of possible options as detailed later in
 ! USERS_GUIDE.md
 call elpa_%set("solver", ELPA_SOLVER_2STAGE, success)

 ! set the AVX BLOCK2 kernel, otherwise ELPA_2STAGE_REAL_DEFAULT will
 ! be used
 call elpa_%set("real_kernel", ELPA_2STAGE_REAL_AVX_BLOCK2, success)

 ! use method solve to solve the eigenvalue problem to obtain eigenvalues
 ! and eigenvectors
 ! other possible methods are desribed in USERS_GUIDE.md
 call elpa_%eigenvectors(a, ev, z, success)

 ! cleanup
 call elpa_deallocate(elpa_)

 call elpa_uninit()

end program src2
