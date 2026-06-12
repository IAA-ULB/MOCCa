!===============================================================================
!     __  __  ___   ____ ____
!    |  \/  |/ _ \ / ___/ ___|__ _
!    | |\/| | | | | |  | |   / _` |
!    | |  | | |_| | |__| |__| (_| |
!    |_|  |_|\___/ \____\____\__,_|
!
! Written mainly by W. Ryssens & M. Bender
!
! Opensource software distributed under the GNU AGPLv3 licence, see the
!  LICENCE file in the root of this project.
!===============================================================================
module timing
!-------------------------------------------------------------------------------
! Timing module for Tantalus, allowing for the definition of times in multiple
! (possibly nested) contexts.
!
! W.R. 2019, but heavily based on earlier work by Christopher Gilbreth on the
! Yale SMMC code in 2010/2014.
!-------------------------------------------------------------------------------
  use compilation,     only : dp
  use geninfo,         only : MPI_RANK, NPROCS, MPI_BLOCK_ASSIGNMENTS
  use iso_fortran_env, only : int64, real64
#if(USE_MPI > 0)
  use MPI            , only : MPI_COMM_WORLD, MPI_BARRIER
#endif

  implicit none (external)

  public

  !-----------------------------------------------------------------------------
  ! Timer IDs. These are set to values by add_timer().
  integer :: T_wfini, T_wfoutput, T_wfinput
  integer :: T_derivatives, T_derivatives_can, T_evolution, T_ortho, T_tantalus
  integer :: T_subspace_rotation, T_subrot_transfo, T_subrot_diag
  integer :: T_calc_sph, T_calc_sph_me, T_update_sph
  integer :: T_densities, T_potentials, T_energy, T_pairing, T_den_ph, T_den_pp
  integer :: T_den_der, T_sphamil, T_coulomb, T_den_can, T_MOI
  integer :: T_COM, T_COM1, T_COM2, T_gaps, T_moments, T_feasible, T_spwfangmom
  integer :: T_chargedensity, T_collective_moi, T_microscopic_pairing
  integer :: T_HFdiag, T_Hortho, T_moment_cutoff, T_norm_ortho, T_diag_ortho
  integer :: T_NablaMElements, T_basistransfo, T_COM2_summation, T_pot_precon
  integer :: T_fam, T_ZPE
  integer :: T_den_perturbed, T_den_perturbed_sym, T_den_perturbed_asym
  integer :: T_den_perturbed_pp, T_spme_perturbed_pp
  integer :: T_spme_perturbed, T_spme_perturbed_sym, T_spme_perturbed_asym
#if(USE_MPI > 0)
  integer :: T_transfer_psi_1to2,T_transfer_psi_2to1, T_allreduce
#endif
  !-----------------------------------------------------------------------------
  ! There are two ways to record the time:
  !  1. cpu_time measures CPU time (excludes time spent in other programs)
  !  2. system_clock measures walltime
  integer, parameter :: CPU_TIME_ = 1, SYSTEM_CLOCK_ = 2
  integer, parameter :: TIME_TYPE = SYSTEM_CLOCK_
  integer(int64)     :: count_rate  ! Conversion for system_clock()

  !-----------------------------------------------------------------------------
  ! Timing info for a timer in a particular context.
  ! When a timer is started, there are usually other timers going already.
  ! The set of which timers are already going is called a "context", and
  ! is represented by a set of binary on/off flags (the key).
  ! We record times in each context separately.
  type context
     ! Binary description of the context
     integer(int64) :: key
     ! # of times this timer has been started in this context
     integer(int64) :: ncalls
     ! for cpu_time
     real(real64)    :: tstart, tstop
     ! for system_clock
     integer(int64) :: itstart, itstop
     ! Total time so far in this context
     real(real64) :: tsum
  end type context

  !-----------------------------------------------------------------------------
  ! Structure to identify a timer and the contexts in which it's been called
  type timer
     integer :: id
     character(len=64) :: name
     ! Array of context-specific timer structures
     type(context), allocatable :: contexts(:)
     ! When this timer is running, this is the index of the relevant context in
     ! contexts(:). When not running, this is the index of the last context in
     ! which it was executed.
     integer :: context_index
  end type timer

  integer                          :: ntimers = 0
  type(timer), allocatable, target :: timers(:)
  integer(int64)                       :: current_context = 0


contains

  elemental subroutine init_context(c)
    type(context), intent(inout) :: c
    c%key = -1
    c%ncalls  = 0
    c%tstart  = 0
    c%tstop   = 0
    c%itstart = 0
    c%itstop  = 0
    c%tsum    = 0
  end subroutine init_context

  subroutine add_timer(name,id)
    !---------------------------------------------------------------------------
    ! Start a new timer with a name and an ID. Contexts are automatically
    ! assigned.
    !---------------------------------------------------------------------------
    character(len=*), intent(in) :: name
    integer, intent(out) :: id
    type(timer), pointer :: t

    integer(int64) :: count

    call allocate_timers(1)
    if (count_rate == 0) then
       call system_clock(count,count_rate)
    end if

    ntimers = ntimers + 1
    id = ntimers
!    if (id > bit_size(current_context)-1) then
!       ! Can't store a context for > 63 timers in an 8-byte integer.
!       stop "Error: too many timers!"
!    end if

    t => timers(id)
    t%id = id
    t%name = name
    allocate(t%contexts(4))
    call init_context(t%contexts)
    t%context_index = 1
  end subroutine add_timer

  subroutine clear_timers()
    ntimers = 0
  end subroutine clear_timers

  subroutine allocate_timers(n)
    !---------------------------------------------------------------------------
    ! Allocate new timers (if necessary).
    !---------------------------------------------------------------------------
    integer, intent(in)      :: n
    integer                  :: new_size
    type(timer), allocatable :: timers_tmp(:)

    if (.not. allocated(timers)) then
       allocate(timers(n+15))
       return
    end if

    if (ntimers + n > size(timers)) then
       new_size = size(timers)
       do
          new_size = new_size*2
          if (new_size >= size(timers) + n) exit
       end do
       allocate(timers_tmp(size(timers)))
       timers_tmp = timers
       deallocate(timers)
       allocate(timers(new_size))
       timers(1:ntimers) = timers_tmp
       deallocate(timers_tmp)
    end if
  end subroutine allocate_timers

  subroutine start_timer(id)
    !---------------------------------------------------------------------------
    ! Start the timer with given id.
    !---------------------------------------------------------------------------
    integer, intent(in) :: id

    type(timer),   pointer :: t
    type(context), pointer :: contexts(:), c
    integer :: i, nc

    if (id < 1 .or. id > ntimers) then
      print *, "Error: invalid timer id"
      print *, id
      stop
    endif
    t => timers(id)
    if (btest(current_context,id-1)) then
       write (0,*) "Error: timer ", trim(adjustl(t%name)), " Already running"
       stop
    end if

    ! Find the context we're running in, or add a new one if necessary
    contexts => t%contexts
    nc = size(contexts)
    do i=1,nc
       c => contexts(i)
       if (c%key < 0 .or. c%key == current_context) exit
    end do
    if (i > nc) then
       call allocate_more_contexts(t)
       c => t%contexts(i)
    end if

    if (c%key < 0) then
       c%key = current_context
    end if
    current_context = ibset(current_context, id-1)

    if (TIME_TYPE == CPU_TIME_) then
       call cpu_time(c%tstart)
    else if (TIME_TYPE == SYSTEM_CLOCK_) then
       call system_clock(c%itstart)
    else
       stop "Error in timers.f90: invalid time_type"
    end if

    t%context_index = i
    c%ncalls = c%ncalls + 1
  end subroutine start_timer

  subroutine stop_timer(id)
    !---------------------------------------------------------------------------
    ! Stop the timer with given id.
    !---------------------------------------------------------------------------
    integer, intent(in)    :: id
    type(timer), pointer   :: t
    type(context), pointer :: contexts(:), c

    if (id < 1 .or. id > ntimers) stop "Error: invalid timer id"
    t => timers(id)
    if (.not. btest(current_context,id-1)) stop "Error: timer already stopped"
    contexts => t%contexts

    current_context = ibclr(current_context, id-1)
    c => contexts(t%context_index)
    if (c%key /= current_context) then
       write (0,*) "Error: tried to stop timer ", trim(t%name), " in wrong context!"
       write (0,*) "Timer's context: ", c%key
       write (0,*) "Expected context: ", current_context
       stop
    end if

    if (TIME_TYPE == CPU_TIME_) then
       call cpu_time(c%tstop)
       c%tsum = c%tsum + c%tstop-c%tstart
    else if (TIME_TYPE == SYSTEM_CLOCK_) then
       call system_clock(c%itstop)
       c%tsum = c%tsum + (c%itstop-c%itstart)/real(count_rate,kind(1.d0))
    else
       stop "Error in timers.f90: invalid time_type"
    end if
  end subroutine stop_timer

  real(KIND=dp) function get_time(id) result(time)
    !---------------------------------------------------------------------------
    ! Get time elapsed for the timer indicated by 'id', allowing for the
    ! possibility that the timer is currently running.
    !---------------------------------------------------------------------------
    integer, intent(in) :: id

    type(timer),   pointer :: t
    type(context), pointer :: c
    real(KIND=dp)          :: tstop
    integer(int64)             :: itstop
    integer                :: i
    logical                :: running

    if (id < 1 .or. id > ntimers) stop "Error: get_time(): invalid timer id"
    t => timers(id)

    time = 0.d0
    do i=1,size(t%contexts)
       c => t%contexts(i)
       if (c%key >= 0) then
          time = time + c%tsum
       end if
    end do

    running = btest(current_context,id-1)
    if (running) then
       c => t%contexts(t%context_index)
       if (TIME_TYPE == CPU_TIME_) then
          call cpu_time(tstop)
          time = time + tstop-c%tstart
       else if (TIME_TYPE == SYSTEM_CLOCK_) then
          call system_clock(itstop)
          time = time + (itstop-c%itstart)/real(count_rate,kind(1.d0))
       else
          stop "Error in timers.f90: invalid time_type"
       end if
    end if
  end function get_time

  subroutine calc_total_time(total)
    !---------------------------------------------------------------------------
    ! Sum time in all outermost contexts. We'll call this the total.
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(out) :: total
    integer                    :: id, ic
    type(timer), pointer       :: t
    type(context), pointer     :: c

    total = 0.d0
    do id=1,ntimers
       t => timers(id)
       do ic=1,size(t%contexts)
          c => t%contexts(ic)
          if (c%key == 0) total = total + c%tsum
       end do
    end do
  end subroutine calc_total_time

  subroutine print_all_timers_flat()
    !---------------------------------------------------------------------------
    !
    !---------------------------------------------------------------------------
    integer :: id, ic, r
    type(timer), pointer   :: t
    type(context), pointer :: c
    real(KIND=DP)          :: time, total
    integer(int64)             :: ncalls
#if(USE_MPI>0)
    integer :: mpi_err
#endif
    do r=1, NPROCS
      if(MPI_RANK == r) then
        if (current_context /= 0) then
           write (*,*) "WARNING: There are timers still running. They should be &
                &stopped"
           write (*,*) "before printing out the timers. Some numbers may be&
                & incorrect."
        end if

        call calc_total_time(total)

        print *, 'Timers of MPI-rank ', r
        write (*,'(a2,tr1,a46,tr2,a10,tr2,a12,tr2,a7)') &
             "id", "Timer name                                    ", &
             "# of calls", "time (s)", "% total"
        write (*,'(2("*"),tr1,46("*"),tr2,10("*"),tr2,12("*"),tr2,7("*"))')
        write (*,'(2x,tr1,a,tr41,tr2,a10,tr2,f12.4,tr2,f6.2,"%")') &
                "TOTAL", "-", total, 100.d0

        do id=1,ntimers
           t => timers(id)
           time = 0.d0
           ncalls = 0
           do ic=1,size(t%contexts)
              c => t%contexts(ic)
              if (c%key >= 0) then
                 time = time + c%tsum
                 ncalls = ncalls + c%ncalls
              end if
           end do
           write (*,'(i2,tr1,a46,tr2,i10,tr2,f12.4,tr2,f6.2,"%")') &
                t%id, t%name, ncalls, time, time/total * 100.d0
        end do
      endif
#if(USE_MPI>0)
      call MPI_BARRIER(MPI_COMM_WORLD,mpi_err)
#endif
    enddo
  end subroutine print_all_timers_flat

  recursive subroutine print_all_timers_aux(icontext,depth,nsub,tsub,&
       total,prnt)
    !---------------------------------------------------------------------------
    ! Recursively print all timers in a given context, and their subtimers.
    !---------------------------------------------------------------------------
    integer,    intent(in)        :: depth
    integer(int64), intent(in)        :: icontext
    integer,    intent(out)       :: nsub
    real(KIND=DP),    intent(out) :: tsub
    real(KIND=DP),    intent(in)  :: total
    logical,    intent(in)        :: prnt

    integer(int64)                    :: isubcontext
    integer                       :: id, ic, nsub1
    real(KIND=DP)                 :: tsub1, tinternal, tpercent
    type(timer),   pointer        :: t
    type(context), pointer        :: c

    character(len=60) :: spaces, str

    nsub = 0; tsub = 0
    spaces = ' '
    do id=1,ntimers
       t => timers(id)
       do ic=1,size(t%contexts)
          c => t%contexts(ic)
          if (c%key == icontext) then
             nsub = nsub + 1
             tsub = tsub + c%tsum
             isubcontext = ibset(icontext,id-1)
             if (prnt) then
                write (str,'(a,a)') spaces(1:depth*2), trim(t%name)
                write (*,'(i2,tr1,a40,tr2,i10,tr2,f12.4,tr2,f6.2,"%")') &
                     id, str, c%ncalls, c%tsum, c%tsum/total*100.d0
             end if
             ! Add up the amount of time spent in subtimers
             call print_all_timers_aux(isubcontext,depth+1,nsub1,tsub1,&
                  total,.false.)
             if (nsub1 > 0 .and. prnt) then
                ! Print amount of time spent in this timer, and not in subtimers
                tinternal = c%tsum - tsub1
                tpercent = tinternal/total*100.d0
                if (tpercent >= .01d0) then
                   write (str,'(a,a)') spaces(1:(depth+1)*2), '(internal)'
                   write (*,'(a2,tr1,a40,tr2,a10,tr2,f12.4,tr2,f6.2,"%")') &
                        '', str, '-', tinternal, tpercent
                end if
             end if
             ! Print all the subtimers
             call print_all_timers_aux(isubcontext,depth+1,nsub1,tsub1,&
                  total,prnt)
          end if
       enddo
    enddo

  end subroutine print_all_timers_aux

  subroutine print_all_timers()
    !--------------------------------------------------------------------------
    ! Print the timing information on a run.
    !
    ! For serial execution, this prints a single table reflecting all calculations.
    ! For parallel execution, this prints one table for every symmetry block
    ! with timing information taken from the lowest rank assigned to that
    ! symmetry block.
    !--------------------------------------------------------------------------
    integer :: nsub, r
    real(real64) :: tsub, total

#if(USE_MPI>0)
    integer :: mpi_err, B
#endif

    call calc_total_time(total)

#if(USE_MPI>0)
    do B=1,8
      call MPI_BARRIER(MPI_COMM_WORLD,mpi_err)
      do r = 0, NPROCS-1
         ! Searching for the first rank assigned to this symmetry block
         if(MPI_BLOCK_ASSIGNMENTS(r+1) .eq. B ) then
            if(MPI_RANK .eq. r) then
               print *,'------------------------------------------------------------------------------'
               print *, 'Timers of rank ', r, ' which deals with block B = ', B
               print *,'------------------------------------------------------------------------------'
               write (*,'(a2,tr1,a40,tr2,a10,tr2,a12,tr2,a7)') &
                     "id", "Timer name                                    ", &
                     "# of calls", "time (s)", "% total"
               write (*,'(2("_"),tr1,40("_"),tr2,10("_"),tr2,12("_"),tr2,7("_"))')
               write (*,'(2x,tr1,a,tr35,tr2,a10,tr2,f12.4,tr2,f6.2,"%")') &
                        "TOTAL", "-", total, 100.d0
               call print_all_timers_aux(0_int64,0,nsub,tsub,total,.true.)
               write (*,'(2("_"),tr1,40("_"),tr2,10("_"),tr2,12("_"),tr2,7("_"))')
            endif
            exit
         endif
      enddo
      call MPI_BARRIER(MPI_COMM_WORLD,mpi_err)
   enddo
#else
   r = 0
   print *,'------------------------------------------------------------------------------'
   print *, 'Timers of rank ', r
   print *,'------------------------------------------------------------------------------'
   write (*,'(a2,tr1,a40,tr2,a10,tr2,a12,tr2,a7)') &
         "id", "Timer name                                    ", &
         "# of calls", "time (s)", "% total"
   write (*,'(2("_"),tr1,40("_"),tr2,10("_"),tr2,12("_"),tr2,7("_"))')
   write (*,'(2x,tr1,a,tr35,tr2,a10,tr2,f12.4,tr2,f6.2,"%")') &
            "TOTAL", "-", total, 100.d0
   call print_all_timers_aux(0_8,0,nsub,tsub,total,.true.)
   write (*,'(2("_"),tr1,40("_"),tr2,10("_"),tr2,12("_"),tr2,7("_"))')
#endif

  end subroutine print_all_timers

  subroutine allocate_more_contexts(t)
    type(timer), intent(inout) :: t

    type(context), allocatable :: tmp(:)
    integer :: nc

    nc = size(t%contexts)
    allocate(tmp(nc))
    tmp = t%contexts
    deallocate(t%contexts)
    allocate(t%contexts(nc*2))
    t%contexts(1:nc) = tmp(1:nc)
    call init_context(t%contexts(nc+1:))
  end subroutine allocate_more_contexts

end module timing
