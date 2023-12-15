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
  integer :: scfscheme = 0
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

    namelist /scfiteration/ scfscheme, denmix, preconfactor

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
    call MPI_BCAST(scfscheme   , 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_BCAST(denmix      , 1, MPI_REAL8  , 0, MPI_COMM_WORLD, mpi_err)
    call MPI_BCAST(preconfactor, 1, MPI_REAL8  , 0, MPI_COMM_WORLD, mpi_err)
#endif
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Bookkeeping for all MPI ranks
    ! Interpreting the scfscheme choice in terms of densities and potentials
    select case(scfscheme)
    case(0)
      ! Potential preconditioning
      densitymixing            = 0
      potentialpreconditioning = 1
    case(1)
      ! Linear mixing of the densities
      densitymixing            = 1
      potentialpreconditioning = 0
    end select

  end subroutine readscfiteration

  subroutine printscfiteration
    !---------------------------------------------------------------------------
    ! Print some info on the SCF-update.
    !---------------------------------------------------------------------------
    
    1 format(80('-'))
    2 format(' SCF iteration strategy: ',/, 2x, a30 )
    3 format('   denmix= '            , f7.4)        
    4 format('   Preconfactor= '      , f7.4)
    
    print 1
    select case(scfscheme)
    case(0)
      print 2, 'Potential preconditioning'
      print 4, preconfactor
    case(1)
      print 2, 'Linear mixing of densities'
      print 3, denmix
    end select

  end subroutine printscfiteration

end module
