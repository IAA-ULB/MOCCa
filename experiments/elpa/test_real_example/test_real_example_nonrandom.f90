!   This file is part of ELPA.
!
!   The ELPA library was originally created by the ELPA consortium,
!   consisting of the following organizations:
!
!   - Max Planck Computing and Data Facility (MPCDF), formerly known as
!    Rechenzentrum Garching der Max-Planck-Gesellschaft (RZG),
!   - Bergische Universität Wuppertal, Lehrstuhl für angewandte
!    Informatik,
!   - Technische Universität München, Lehrstuhl für Informatik mit
!    Schwerpunkt Wissenschaftliches Rechnen ,
!   - Fritz-Haber-Institut, Berlin, Abt. Theorie,
!   - Max-Plack-Institut für Mathematik in den Naturwissenschaften,
!    Leipzig, Abt. Komplexe Strukutren in Biologie und Kognition,
!    and
!   - IBM Deutschland GmbH
!
!
!   More information can be found here:
!   http://elpa.mpcdf.mpg.de/
!
!   ELPA is free software: you can redistribute it and/or modify
!   it under the terms of the version 3 of the license of the
!   GNU Lesser General Public License as published by the Free
!   Software Foundation.
!
!   ELPA is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU Lesser General Public License for more details.
!
!   You should have received a copy of the GNU Lesser General Public License
!   along with ELPA.  If not, see <http://www.gnu.org/licenses/>
!
!   ELPA reflects a substantial effort on the part of the original
!   ELPA consortium, and we ask you to respect the spirit of the
!   license that we chose: i.e., please contribute any changes you
!   may have back to the original ELPA library distribution, and keep
!   any derivatives of ELPA under the same license that we chose for
!   the original distribution, the GNU Lesser General Public License.
!
!
!>
!> Fortran test programm to demonstrates the use of
!> ELPA 1 real case library.
!> If "HAVE_REDIRECT" was defined at build time
!> the stdout and stderr output of each MPI task
!> can be redirected to files if the environment
!> variable "REDIRECT_ELPA_TEST_OUTPUT" is set
!> to "true".
!>
!> By calling executable [arg1] [arg2] [arg3] [arg4]
!> one can define the size (arg1), the number of
!> Eigenvectors to compute (arg2), and the blocking (arg3).
!> If these values are not set default values (4000, 1500, 16)
!> are choosen.
!> If these values are set the 4th argument can be
!> "output", which specifies that the EV's are written to
!> an ascii file.
!>
program test_real_example

!-------------------------------------------------------------------------------
! Standard eigenvalue problem - REAL version
!
! This program demonstrates the use of the ELPA module
! together with standard scalapack routines
!
! Copyright of the original code rests with the authors inside the ELPA
! consortium. The copyright of any additional modifications shall rest
! with their original authors, but shall adhere to the licensing terms
! distributed along with the original code in the file "COPYING".
!
!-------------------------------------------------------------------------------

  use iso_c_binding
  use elpa
  use mpi
  implicit none

  !-------------------------------------------------------------------------------
  ! Please set system size parameters below!
  ! na:  System size
  ! nev:  Number of eigenvectors to be calculated
  ! nblk: Blocking factor in block cyclic distribution
  !-------------------------------------------------------------------------------

  integer                  :: nblk
  integer                  :: na, nev

  integer                  :: np_rows, np_cols, na_rows, na_cols

  integer                  :: myrank, nranks, my_prow, my_pcol, mpi_comm_rows, mpi_comm_cols
  integer                  :: i, mpierr, my_blacs_ctxt, sc_desc(9), info, nprow, npcol

  integer, external        :: numroc

  real(kind=c_double), allocatable :: a(:,:), z(:,:), ev(:)
  real(kind=c_double)      :: aij

  integer                  :: iseed(4096) ! Random seed, size should be sufficient for every generator

  integer                  :: STATUS
  integer                  :: nargs, command_argument_count
  integer                  :: success
  character(len=8)         :: task_suffix
  integer                  :: j

  integer, parameter       :: error_units = 0
  character(100)           :: fmt1
  integer                  :: il,jl,indxl2g
  
  class(elpa_t), pointer   :: e
  
!-------------------------------------------------------------------------------
! default parameters (may be overwritten by command line arguments)
  na   = 1000
  nblk =   16

!-------------------------------------------------------------------------------
! command line arguments
! 1 -> matrix size
! 2 -> blocking factor
  nargs = command_argument_count()
  if (nargs.ge.1) then
    call get_command_argument(1, fmt1)
    read(fmt1,*,iostat=STATUS) na
  endif
  if (nargs.ge.2) then
    call get_command_argument(2, fmt1)
    read(fmt1,*,iostat=STATUS) nblk
  endif

!-------------------------------------------------------------------------------
! Fix nev (we might want another command line argument for that)
  nev = na
  
  call mpi_init(mpierr)
  call mpi_comm_rank(mpi_comm_world,myrank,mpierr)
  call mpi_comm_size(mpi_comm_world,nranks,mpierr)

  if (myrank.eq.0) then
    print '(a)', 'Source file : experiments/elpa/test_real_example/test_real_example.f90'
    write (*,"('matrix size     :', i8)") na
    write (*,"('blocking factor :', i8)") nblk
    write (*,"('#eigenvalues.   :', i8)") nev
  endif

  do np_cols = NINT(SQRT(REAL(nranks))),2,-1
    if(mod(nranks,np_cols) == 0 ) exit
  enddo
  ! at the end of the above loop, nranks is always divisible by np_cols
  np_rows = nranks/np_cols


  ! initialise BLACS
  my_blacs_ctxt = mpi_comm_world
  call BLACS_Gridinit(my_blacs_ctxt, 'C', np_rows, np_cols)
  call BLACS_Gridinfo(my_blacs_ctxt, nprow, npcol, my_prow, my_pcol)

  fmt1 = "('[', i0, '/', i0, '] -> ', i0, 'x', i0, ' (', i0, ',', i0, ')')"
  write (*,fmt1) myrank, nranks, np_rows, np_cols, my_prow, my_pcol

  if (myrank==0) then
    print '(a)','| Past BLACS_Gridinfo.'
  end if
  ! determine the neccessary size of the distributed matrices,
  ! we use the scalapack tools routine NUMROC

  na_rows = numroc(na, nblk, my_prow, 0, np_rows)
  na_cols = numroc(na, nblk, my_pcol, 0, np_cols)

  ! set up the scalapack descriptor for the checks below
  ! For ELPA the following restrictions hold:
  ! - block sizes in both directions must be identical (args 4 a. 5)
  ! - first row and column of the distributed matrix must be on
  !  row/col 0/0 (arg 6 and 7)

  call descinit(sc_desc, na, na, nblk, nblk, 0, 0, my_blacs_ctxt, na_rows, info)

  if (info .ne. 0) then
    write(error_units,*) 'Error in BLACS descinit! info=',info
    write(error_units,*) 'Most likely this happend since you want to use'
    write(error_units,*) 'more MPI tasks than are possible for your'
    write(error_units,*) 'problem size (matrix size and blocksize)!'
    write(error_units,*) 'The blacsgrid can not be set up properly'
    write(error_units,*) 'Try reducing the number of MPI tasks...'
    call MPI_ABORT(mpi_comm_world, 1, mpierr)
  endif

  if (myrank==0) then
    print '(a)','| Past scalapack descriptor setup.'
  end if

  allocate(a (na_rows,na_cols)) ! local (sub-)matrix
  allocate(z (na_rows,na_cols)) ! local (sub-)matrix

  allocate(ev(na)) ! global array for all eigenvalues, exists on all MPI ranks.

  ! we want different random numbers on every process
  ! (otherwise A might get rank deficient):

  do il=1,na_rows
    i = indxl2g(il, nblk, my_prow, 0, np_rows)
    if (.not.((i.ge.1).and.(i.le.na))) then
      print "('ERROR : indxl2g maps il=', i0, ' to i=', i0)", il, i
    endif
    do jl=1,na_cols
      j = indxl2g(jl, nblk, my_pcol, 0, np_cols)
      if (.not.((j.ge.1).and.(j.le.na))) then
        print "('ERROR : indxl2g maps jl=', i0, ' to j=', i0)", il, i
      endif
      a(il,jl) = aij(i,j)
    enddo
    write (*,'(i1)', advance='no') myrank
    do jl=1,na_cols 
      write (*,'(f6.2)', advance='no') a(il,jl)
    enddo
    write (*,*)
  enddo

  ! print '(a)','| Symmetric matrix block has been set up.'
  write (*,"('| rank ', i0, ' : Symmetric matrix block has been set up.')"), myrank
  
  !-------------------------------------------------------------------------------

  if (elpa_init(20240501) /= elpa_ok) then 
    print *, "ERROR : ELPA API version not supported"
    stop 1
  endif

  e => elpa_allocate(success)
  if (success /= ELPA_OK) then
    print *,"ERROR : e => elpa_allocate(success) failed." 
    stop 1
  endif

  ! set parameters decribing the matrix and it's MPI distribution
  call e%set("na" , na , success)
  call e%set("nev", nev, success)
  call e%set("local_nrows", na_rows, success)
  call e%set("local_ncols", na_cols, success)
  call e%set("nblk", nblk, success)
  call e%set("mpi_comm_parent", mpi_comm_world, success)
  call e%set("process_row", my_prow, success)
  call e%set("process_col", my_pcol, success)

  success = e%setup()
  if (success /= ELPA_OK) then
    print *,"ERROR : success = e%setup() failed." 
    stop 1
  endif

  call e%set("solver", elpa_solver_1stage, success)

  ! Calculate eigenvalues/eigenvectors

  if (myrank==0) then
    print '(a)','| Entering one-step ELPA solver '
    print '(a)','| ... '
  end if

  call mpi_barrier(mpi_comm_world, mpierr) ! for correct timings only
  call e%eigenvectors(a, ev, z, success)

  if (myrank.eq.0) then
    write (*,*) "Eigenvalues:"
    do i=1,nev
      write (*,*) i,ev(i)
    enddo
  endif

  if (myrank==0) then
    print '(a)','| One-step ELPA solver complete.'
    print '(a)','| cleaning up.'
  end if

  call elpa_deallocate(e)
  call elpa_uninit()
  call blacs_gridexit(my_blacs_ctxt)
  call mpi_finalize(mpierr)

end

function aij(i,j) 
  use iso_c_binding
  integer, intent(in) :: i,j ! input
  real(kind=c_double) :: aij ! output
  integer, parameter  :: symmetric = 1

  if (i.lt.1) then 
    write(*,*) "ERROR aij : i.lt.1"
    stop 1
  endif
  if (j.lt.1) then 
    write(*,*) "ERROR aij : j.lt.1"
    stop 1
  endif

  if (symmetric.eq.1) then
    if (i.eq.j) then
      aij = 1.5
    else
      aij = 1.0 / sqrt(real(abs(i-j)))
    endif
  else
    aij = 10*i + j 
  endif
end function
