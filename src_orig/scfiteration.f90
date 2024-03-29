module SCFiteration
!===============================================================================
!  #######   ##   #    # #####   ##   #      #    #  ####
!     #     #  #  ##   #   #    #  #  #      #    # #
!     #    #    # # #  #   #   #    # #      #    #  ####
!     #    ###### #  # #   #   ###### #      #    #      #
!     #    #    # #   ##   #   #    # #      #    # #    #
!     #    #    # #    #   #   #    # ######  ####   ####
!
!  Copyright W. Ryssens & M. Bender
!
!===============================================================================
! 
! This module is the central control center for the treatment of the evolution 
! of the mean-field potentials and densities from one iteration to the next. 
!
!===============================================================================

  use functional
  use densities

  implicit none
  
  !-----------------------------------------------------------------------------
  ! Determine the SCF-evolution scheme
  !  (0) => Preconditioning of necessary potentials  (here F_I_I)
  !  (1) => Linear mixing of the necessary densities (here D_I_I)
  !
  ! Currently not changeable; only (0) is working.
  integer, parameter :: scfscheme = 0
  !-----------------------------------------------------------------------------
  ! Determine what to do with mixing of the potentials
  integer       :: mixingscheme = 0
  real(KIND=dp) :: mixstepsize  = 1.0d0

contains

  subroutine readscfiteration(file_number)
    !---------------------------------------------------------------------------
    ! Read the namelist determining the SCF-update.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !   file_number : optional integer. If present, read from (open) channel
    !                 with this number. If absent, read from STDIN.
    !---------------------------------------------------------------------------
    integer(dp), intent(in), optional   :: file_number   
#if(USE_MPI>0)
    integer                             :: mpi_err
#endif

    namelist /scfiteration/ preconfactor, mixingscheme,mixstepsize, memory

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Only the very first MPI rank reads the input
    if(MPI_RANK.eq.0) then
      if(present(file_number)) then
        read (unit=file_number, nml=scfiteration)
      else
        read (unit=*, nml=scfiteration)
      endif

      ! Sanity checks
      if((scfscheme .ne. 0) .and. (scfscheme.ne.1)) then
        call stp('Invalid scfscheme value.')
      endif
    endif
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Broadcasting the MPI information
#if(USE_MPI > 0)
    !call MPI_BCAST(scfscheme   , 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_BCAST(preconfactor, 1, MPI_REAL8  , 0, MPI_COMM_WORLD, mpi_err)
    call MPI_BCAST(mixingscheme, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_BCAST(mixstepsize , 1, MPI_REAL8  , 0, MPI_COMM_WORLD, mpi_err)
    call MPI_BCAST(memory      , 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)
#endif
  
  end subroutine readscfiteration

  subroutine printscfiteration
    !---------------------------------------------------------------------------
    ! Print some info on the SCF-update.
    !---------------------------------------------------------------------------
    
    1 format(80('-'))
    2 format(' SCF iteration strategy: ',/, 2x, a30 )
    !3 format('   denmix= '            , f7.4)        
    4 format('   Preconfactor= '      , f7.4)
    6 format(' Potential mixing active!', /,     &  
    &        '                  memory:' 2x, i4, &
    &        '                stepsize:',2x, f7.4)    
        
    print 1
    select case(scfscheme)
    case(0)
      print 2, 'Potential preconditioning'
      print 4, preconfactor
    !case(1)
    !  print 2, 'Linear mixing of densities'
    !  print 3, denmix
    end select
    if(mixingscheme.eq.1) then
      print 6, memory, mixstepsize
    endif
  end subroutine printscfiteration
  
  function AndersonMixPotentials( iterates, updates, stepsize, Nsaved) result(F)
    !---------------------------------------------------------------------------
    ! Perform an Anderson Mixing step on the potential vectors to accelerate
    ! convergence. This particular implementation is based on the description
    ! in 
    ! 
    ! M. F. Herbst and A. Levitt, 
    ! A robust and efficient line search for self-consistent field iterations
    ! arXiv:2109.14018 (unpublished at the time of typing)
    ! 
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! For every self-consistent field iteration we write the starting = input
    ! potentialvector as F_n. The output of the iteration is then 
    !
    !                F(R(F_n))
    !
    ! where R(F_n) is the density-vector generated from an iteration, and F is  
    ! the potentialvector calculated from the density-vector R(F_n).
    ! Denoting further 
    !
    !         delta F_n = P^{-1} [ F(R(F_n)) - F_n ] 
    !
    ! where P^{-1} is the preconditioner acting on (the difference of) the 
    ! functional vectors.
    !
    ! This routine produces a new potential vector, built out of a linear 
    ! combinations of the past n previous updates.
    !
    !  F    = M_n + a^{-1} \sum_{i=1}^{n-1} beta_i M_i
    !  M_i  = F_i + alpha delta F_i 
    !
    ! where alpha is a step-size parameter, that MOCCa generally puts to 1.
    ! The beta_i are determined by an optimisation: 
    !
    !    || delta F_n + \sum_{i=1}^{n-1} beta_i [ delta F_(n-i) - delta F_n] ||
    !
    ! which yields the system of linear equations
    !
    !      A beta = b
    !
    ! with matrix 
    !   A_ij = <  delta F_(n-i) - delta F_n |  delta F_(n-j) - delta F_n >
    ! and rhs 
    !    b_i = - < delta F_n | delta F_(n-i) - delta F_n > 
    !
    ! Key to making this thing work is making sure the conditioning number of 
    ! the least-square problems. This means (i) not too much iterates and (ii)
    ! actively checking the conditioning number.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Input:
    !   iterates: the potentialvectors that were the inputs at every iteration
    !             F_i, including the current one (F_n).
    !   updates : the updates produced by previous iterations,
    !              delta F_i, including the current one (delta F_n).
    !              NOTE: these are the proposed updates without any kind of
    !                    Anderson mixing, but WITH preconditioning applied
    !                    for ease of use. 
    !   stepsize: parameter for the size of the update, alpha in the equations
    !             above. Typically set to one.
    !   Nsaved  : number of updates already in memory
    ! Output: 
    !   F       : a new potentialvector, suited for starting the next iteration.
    !---------------------------------------------------------------------------
    
    type(Potentialvector), intent(in)  :: iterates(:), updates(:)
    real(KIND=dp), intent(in)          :: stepsize
    type(Potentialvector)              :: F
    integer, intent(in)                :: Nsaved

    type(PotentialVector), allocatable :: X(:)
    integer                    :: N, info, i, j, lwork
    integer, allocatable       :: ipiv(:)
    real(KIND=dp), allocatable :: beta(:), A(:,:), Acopy(:,:), eigval(:), rhs(:)
    real(KIND=dp), allocatable :: work(:)
    real(KIND=dp)              :: cond
    
    N = min(size(iterates), Nsaved)
    
    if(N .eq. 1) then
      ! If the history is just a single iteration, just take that update
      F = iterates(1) + stepsize * updates(1)
      return
    endif
     
    allocate(beta(N-1)) ; beta = 0.0d0
    allocate(A(N-1,N-1)); A    = 0.0d0
    allocate(X(N-1))
    
    !---------------------------------------------------------------------------
    ! We introduce the shorthand 
    ! X(i) = delta F(n-i) - delta F(n)
    do i=1,N-1
      X(i) = updates(i+1) + (-1.0d0) * updates(1)
    enddo    

    !---------------------------------------------------------------------------
    ! We build the complete matrix of the linear system
    do i=1,N-1
      do j=i,N-1
        A(i,j) = PVectorInproduct(X(i), X(j))
        A(j,i) = A(i,j)
      enddo
    enddo  
    Acopy = A

    !---------------------------------------------------------------------------
    ! Then, we investigate the condition number of the matrix A.
    ! Since this matrix is not big (N ~ 10 at most) and symmetric, we don't mess
    ! around with sophisticated methods 
    allocate(eigval(N-1))    
    lwork = 100 !3*(N-1) - 1
    allocate(Work(lwork))
    ! We simply diagonalize the matrix...
    call DSYEV( 'V', 'U', N-1, Acopy, N-1, eigval, work, lwork,info)

    if(info .ne. 0) then
      print *, 'Problem with call to DSYEV in AndersonMixPotentials'
      print *, 'INFO = ', info      
      stop
    endif
    deallocate(work)
    !... and calculate the condition number directly
    cond = sqrt(eigval(N-1))/sqrt(eigval(1))
    print *, 'condition number of A', cond, eigval(N-1), eigval(1)
    print *, '2x2 submatrix', A(1,1:2)
    print *, '2x2 submatrix', A(2,1:2)
    !
    ! MAYBE DO SOMETHING WITH THIS INFORMATION?
    !
    !---------------------------------------------------------------------------
    ! And then we solve the linear system: 
    allocate(rhs(N-1)) ; rhs = 0.0d0  
    do i=1,N-1
      rhs(i) = - PVectorInproduct(X(i), updates(1))
    enddo
    
    allocate(Ipiv(N-1))
    allocate(work(1))
    call dsysv ('U', N-1, 1, A, N-1, ipiv, rhs, N-1, work, -1,info)
    lwork = int(work(1)) ; deallocate(work) ; allocate(work(lwork))
    call dsysv ('U', N-1, 1, A, N-1, ipiv, rhs, N-1, work, lwork,info)
    
    if(info .ne. 0) then
      print *, 'Problem with call to DSYSV in AndersonMixPotentials'
      print *, 'INFO = ', info      
      stop
    endif
    
    print *, 'MIXCOEFFS', rhs
    
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! rhs now contains the solution to the linear system, i.e. the mixing 
    ! values beta_i.
    F = iterates(1) +  stepsize * updates(1)
    do i=1,N-1
      F = F +               rhs(i) * (iterates(i+1) + (-1.0d0) * iterates(1)) &
       &    +   stepsize * (rhs(i) * (updates(i+1)  + (-1.0d0) * updates(1)))
    enddo

    return 
  end function AndersonMixPotentials

end module
