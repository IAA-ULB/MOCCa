module version

 !==============================================================================
 !_________ _______  _       _________ _______  _                 _______ 
 !\__   __/(  ___  )( (    /|\__   __/(  ___  )( \      |\     /|(  ____ \
 !   ) (   | (   ) ||  \  ( |   ) (   | (   ) || (      | )   ( || (    \/
 !   | |   | (___) ||   \ | |   | |   | (___) || |      | |   | || (_____ 
 !   | |   |  ___  || (\ \) |   | |   |  ___  || |      | |   | |(_____  )
 !   | |   | (   ) || | \   |   | |   | (   ) || |      | |   | |      ) |
 !   | |   | )   ( || )  \  |   | |   | )   ( || (____/\| (___) |/\____) |
 !   )_(   |/     \||/    )_)   )_(   |/     \|(_______/(_______)\_______)
 !                                                                       
 !  Copyright W. Ryssens & M. Bender
 !==============================================================================
 ! Module that contains the routine to print the header of the code, which 
 ! contains information about the code version, symmetry, compilation options,
 ! and environment as obtained in the compilation process from the Makefile.
 !==============================================================================

 implicit none 

contains 

subroutine print_header(fam)
    !------------------------------------------------------------------------------
    ! Print the header to STDOUT.
    !
    ! Input
    !  fam : logical, indicating whether a regular header or FAM header is printed
    !------------------------------------------------------------------------------
    use GenInfo, only: SYMSTRING, MPI_RANK, NPROCS, reduX, reduY, reduZ
    use IO,      only : SYM_CODE, TRANS_CODE

    logical, intent(in) :: fam

    ! Information gleaned from git and the Makefile
    character(len=44), parameter        :: versiontag=VTAG
    character(len=58), parameter        :: version1  =VERSION1
    character(len=58), parameter        :: version2  =VERSION2
    character(len=58), parameter        :: version3  =VERSION3
    !  character(len=58), parameter        :: version4  =VERSION4
    character(len=58), parameter        :: compiler  =COMPCOMP
    character(len=58), parameter        :: cflags    =CFLAGS
    character(len=58), parameter        :: optflags  =OPTFLAGS

    ! Formatting statements for printing a nice header
    character(len=1000), parameter :: header = "(                                 &
        &     8x,' ____________________________________________________________', &
        &   /,8x,'|                                                           |'  &
        &   /,8x,'|                                                           |', &
        &   /,8x,'|  #######   ##   #    # #####   ##   #      #    #  ####   |', &
        &   /,8x,'|     #     #  #  ##   #   #    #  #  #      #    # #       |', &
        &   /,8x,'|     #    #    # # #  #   #   #    # #      #    #  ####   |', &
        &   /,8x,'|     #    ###### #  # #   #   ###### #      #    #      #  |', &
        &   /,8x,'|     #    #    # #   ##   #   #    # #      #    # #    #  |', &
        &   /,8x,'|     #    #    # #    #   #   #    # ######  ####   ####   |', &
        &   /,8x,'|                                                           |', &
        &   /,8x,'|                                                           |')"

    character(len=1000), parameter :: famheader = "(                              &
        &     8x,' ____________________________________________________________', &
        &   /,8x,'|                                                           |'  &
        &   /,8x,'|                                                           |', &
        &   /,8x,'|    #####  ##   #     # #####   ##   #      #    #  ####   |', &
        &   /,8x,'|    #     #  #  ##   ##   #    #  #  #      #    # #       |', &
        &   /,8x,'|    #### #    # # # # #   #   #    # #      #    #  ####   |', &
        &   /,8x,'|    #    ###### #  #  #   #   ###### #      #    #      #  |', &
        &   /,8x,'|    #    #    # #     #   #   #    # #      #    # #    #  |', &
        &   /,8x,'|    #    #    # #     #   #   #    # ######  ####   ####   |', &
        &   /,8x,'|                                                           |', &
        &   /,8x,'|                                                           |')"

    character(len=1000), parameter :: versioninfo = "(                            &
    &        8x,'|--------------- Version Information -----------------------|', &
    &      /,8x,'| Version tag = ', a44, '|',                                    &
    &      /,8x,'| ', a58, '|',                                                  &
    &      /,8x,'| ', a58, '|',                                                  &
    &      /,8x,'| ', a58, '|',                                                  &
    &      /,8x,'|                                                           |')"

    character(len=1000), parameter :: syminfo = "(                                &
    &         8x,'|-------------- Symmetry Information -----------------------|', &
    &       /,8x,'| S.p. generators        = ', a26, 7x, '|',                     &
    &       /,8x,'| Axis reduction  X Y Z  = ', 3i2, 27x, '|',                    &
    &       /,8x,'| SYM_CODE               = ', a26, 7x, '|',                     &
    &       /,8x,'| TRANS_CODE             = ', a26, 7x, '|')"

    character(len=1000), parameter :: compilationchoices = "(                     &
    &         8x,'|-------------- Compilation choices ------------------------|',&
    &       /,8x,'| Calculation type    = ', a36, '|',                           &
    &       /,8x,'| Boundary conditions = ', a36, '|',                           &
    &       /,8x,'| Derivatives of densities via ', a29, '|',                    &
    &       /,8x,'| ', a58, '|' )"

    character(len=200), parameter :: envinfo = "(                                 &
    &          8x,'|-------------- Environment Information --------------------|',&
    &        /,8x,'|  Number of MPI_ranks   = ', i6, 27x, '|')"

    character(len=1000), parameter :: compinfo = "(                               &
    &          8x,'|-------------- Compilation Information --------------------|',&
    &        /,8x,'| Compiled with:                                            |',&
    &        /,8x,'| ', a58, '|'                                                 ,&
    &        /,8x,'| Compilation flags reported:                               |',&
    &        /,8x,'| ', a58, '|'                                                 ,&
    &        /,8x,'| Optimisation flags reported:                              |',&
    &        /,8x,'| ', a58, '|',                                                 &
    &        /,8x,'|___________________________________________________________|')"

    ! intermediate character definitions
    character(len=26)                   :: symprint

    ! compilation choices determined by precompiler directives
#if(PASTA > 0)
    character(len=36), parameter :: calctype = 'PASTA '
#else
    character(len=36), parameter :: calctype = 'NUCLEI'
#endif
#if(USE_Periodic > 0)
    character(len=36), parameter  :: boundary_conditions= 'periodic'
#else
    character(len=36), parameter  :: boundary_conditions= 'anti-periodic'
#endif
#if(DENSUM == 1)
    character(len=29), parameter  :: den_deriv= 'density summation'
#else
    character(len=29), parameter  :: den_deriv= 'derivative routines'
#endif
#if(USE_MPI > 0)
    character(len=58), parameter  :: mpi_enabled= 'MPI enabled'
#else
    character(len=58), parameter  :: mpi_enabled= 'MPI disabled'
#endif


    if(MPI_RANK .eq. 0) then
    print *
    if(fam) then
        write(*, fmt=famheader)
    else
        write(*, fmt=header)
    endif
    write(*, fmt=versioninfo) versiontag, version1, version2, version3
    !----------------------------------------------------------------------------
    ! Information about symmetry choices
    symprint = adjustl(SYMSTRING)
    write(*, fmt=syminfo) symprint, reduX, reduY, reduZ, SYM_CODE, TRANS_CODE
    !----------------------------------------------------------------------------
    ! Other information about compile-time choices
    write(*,fmt=compilationchoices) &
    &  adjustl(calctype), adjustl(boundary_conditions), adjustl(den_deriv), &
    &  adjustl(mpi_enabled)
    !----------------------------------------------------------------------------
    ! Environment information
    write(*,fmt=envinfo) NPROCS
    !----------------------------------------------------------------------------
    ! Technical details about compilation
    write(*,fmt=compinfo) compiler, cflags, optflags
    endif

end subroutine print_header

end module version 