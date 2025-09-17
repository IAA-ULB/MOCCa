module IO_wf
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
  !
  !============================================================================== 
  ! Module governing the writing and reading of wave function files; 
  !   a submodule of the IO.f90 file.
  !------------------------------------------------------------------------------
  !
  ! Hephaestos keywords:
  !
  !     SYM_CODE   : $SYM_CODE  
  !        String encoding the symmetry choices of this particular version of 
  !        Tantalus. Written to .wf files created by this version.
  !
  !     TRANS_CODE : $TRANS_CODE
  !        String encoding the symmetry choices that this version of Tantalus
  !        can READ (in addition to its own type of files).
  ! 
  ! I.e. when compiled, the code can read files characterized by either 
  ! SYM_CODE or TRANS_CODE (employing an additional transformation in the second 
  ! case). The code will however ALWAYS write SYM_CODE .wf files. 
  !==============================================================================

  use compilation
  use GenInfo,       only: nx, ny, nz, dp, NPROCS, MPI_RANK, stp
  use wavefunctions, only: HFBLOCKS
  use functional,    only: ini_name_param, pairingtype, BCSGaps, HFBGaps, FermiEnergy
  use transform,     only: sym_transfo_needed

  implicit none

  !-----------------------------------------------------------------------------
  ! Version number of the .wf file written by this version of the code. 
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Some history:
  !   version 1 : initial version of the .wf files (fit of GSk1-2)
  !               (2018 - Nov. 2020)
  !   version 2 : implementation of symmetry encoding string
  !               (Nov. 2020 - March 2021)
  !   version 3 : inclusion of 
  !               * the full Bogoliubov transformation 
  !               * configuration matrix
  !               * HF transformation for calculations with diagsphamil=.false.
  !               (March 2021 - March 2021)
  !   version 4 : inclusion of 
  !               * blocking information 
  !               * single-particle hamiltonian
  !               (March 2021 - November 2021)
  !   version 5 :  inclusion of 
  !               * cranking frequencies omega_x/y/z
  !               (until September 2023)
  !   version 6 : separation of the HFPsi array from one read/write
  !               to "nwt" x read/writes for easier MPI useage
  !               (from September 2023)
  !-----------------------------------------------------------------------------
  integer, parameter  :: version_number = 6
  integer             :: file_version = 0

  !-----------------------------------------------------------------------------
  ! Encodings of the symmetry choices imposed by Hephaestos
  character(len=26), parameter :: SYM_CODE   = "$SYM_CODE"
  character(len=26), parameter :: TRANS_CODE = "$TRANS_CODE"

  !-----------------------------------------------------------------------------
  ! Characteristics of the calculation stored on the .wf file
  integer              :: filenx, fileny, filenz, filemv
  integer              :: filenwn, filenwp, filepairing
  integer              :: filenwt
  integer              :: fileblocks_global(8), fileblocks(8)
  integer, allocatable :: file_spwf_map(:),file_rank_map(:),file_spwf_inverse(:)
  real(KIND=dp)        :: filedx, fileneutrons, fileprotons
  !-----------------------------------------------------------------------------
  ! Did we succeed in reading a HFB configuration from file? 
  logical       :: readHFBinfofile= .false.
  ! Blocking information from file
  integer       :: fileblocktype=0, fileblocknumber = 0
  integer, allocatable          :: fileblockindices(:)
  character(len=2), allocatable :: fileBlockLowest(:)
  integer                       :: file_HFB_blocks(8)
  !-----------------------------------------------------------------------------
  logical             :: Allowtransform = .false.
  integer             :: extraspwfs(8) = 0


 contains

  subroutine read_tantalus_wf(chan, ifn)
    !---------------------------------------------------------------------------
    ! Reading all information from a previous Tantalus run stored in a .wf file.
    ! It does not (yet) exploit MPI I/O; reading is essentially done by rank 0
    ! and then broadcasted to the rest of the ranks.
    !
    ! This routine also performs a few sanity checks. 
    ! Currently:
    !   
    !   *) equality of (nx,ny,nz) between data and file
    !   *) equality of (nwn,nwp) between data and file
    !   *) the symmetry encoding matches either SYM_CODE or TRANS_CODE
    !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !   chan : integer, channel to open the file ifn on
    !   ifn  : character, name of the input file. 
    !          The code will first check for its existence.
    !---------------------------------------------------------------------------
    ! 
    ! Things read from file. (Not yet implemented ones are indicated by *)
    !
    ! Version
    ! Convergence information: E, dE                         (*)
    ! nx,ny,nz,dx,dt
    ! Symmetry information                                   (1)
    ! neutrons,protons
    ! nwn, nwp, Number of wavefunctions in every block
    ! spenergies, dispersions
    ! diagsphamil
    ! HFtransfo 
    ! (nwt) Wavefunctions
    ! Forcename
    ! single-particle hamiltonian
    ! Pairing information
    !    - Pairingtype
    !    - Rho_can = occupation factors 
    !      (HF)  nothing
    !      (BCS) Fermi level
    !        |   BCSGaps  
    !      (HFB) 
    !        |   blocktype, blocknumber
    !        |   block-sizes for HFB solver
    !        |   blocklowest/blockindices
    !        |   Fermi level
    !        |   rho_pairing    
    !        |   kappa_pairing
    !        |   can_transfo
    !        |   HFBgaps      
    !        |   Bogoliubov transformation 
    !        |   Configuration matrix      
    ! CrankingInfo
    ! Potentials                                             (2)
    ! Multipole Moments
    !     | The code writes the data on ALL the multipole moments.
    !     | For the format of the lines, see the Moments module.
    !
    !---------------------------------------------------------------------------    
    use functional
    use moments
    use cranking
    
    integer, intent(in)          :: chan
    character(len=*), intent(in) :: ifn
    character(len=20)            :: func_name_check
    character(len=26)            :: SYM_CODE_CHECK
    integer                      :: io,i, wave, wfcounter, targetrank
    logical                      :: exists
    real(KIND=dp)                :: Omega_file(3)
    real(KIND=dp), allocatable   :: filegaps(:,:), temp(:,:)
    logical                      :: filediagsphamil 
    logical                      :: check_x, check_y, check_z
    logical                      :: check_nwn, check_nwp

#if(USE_MPI>0)
    integer :: mpi_err
#endif   
    
    
    1 format ('Number of mesh points does not correspond to file.', / &
    &         'On file: nx= ', i3, ' ny= ', i3, ' nz= ',i3,            / &
    &         'In data: nx= ', i3, ' ny= ', i3, ' nz= ',i3)
    2 format ('Number of wavefunctions does not correspond to file.', / &
    &         'On file: nwn= ', i3, ' nwp=', i3,                      / &
    &         'In data: nwn= ', i3, ' nwp=', i3)

    3 format (' The symmetry choices  on file cannot be handled.')
    4 format (' SYM_CODE   = ', a26)
    5 format (' TRANS_CODE = ', a26)
    6 format (' ON FILE    = ', a26)
 
    if(MPI_RANK.eq.0) then
      !-------------------------------------------------------------------------
      ! First check if the file exists.
      inquire(file=ifn, exist=exists)
      if(.not.exists) then
        call stp('Input file specified does not exist!')
      endif
 
      open (chan,form='unformatted',file=ifn)

      read(chan, iostat=io) file_version
      if(file_version .gt. version_number) then
        call stp('Unsupported version number of the .wf file.')
      endif

      ! Convergence information                                (NOT IMPLEMENTED)
      read(chan,iostat=io) 
      !Parameters of the mesh
      read(Chan,iostat=io) filenx,fileny,filenz, filedx
      filemv = filenx*fileny*filenz

      ! Symmetry information       
      if(file_version .eq. 1) then 
        ! No symmetry information in version 1, only EV8-style calculations  
        read(Chan,iostat=io) 
      else
        read(Chan,iostat=io) SYM_CODE_CHECK

        if(SYM_CODE_CHECK .eq. SYM_CODE) then
          sym_transfo_needed = .false.
        elseif(SYM_CODE_CHECK .eq. TRANS_CODE) then
          sym_transfo_needed = .true.
        else
          print 3
          print 4, SYM_CODE 
          print 5, TRANS_CODE
          print 6, SYM_CODE_CHECK
          call stp('')
        endif
      endif

      !Number of protons and neutrons
        read(Chan,iostat=io) fileneutrons, fileprotons
      ! HFBLocks information 
      read(Chan,iostat=io) filenwn, filenwp, fileblocks_global
      filenwt = filenwn + filenwp
    endif

    !---------------------------------------------------------------------------
    ! Rank 0 now has a ton of information read from file, including the 
    ! dimensions of the symmetry blocks on the file.
#if(USE_MPI > 0)
    ! First, we broadcast this information
    call MPI_BCAST(filenx, 1, MPI_integer, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_BCAST(fileny, 1, MPI_integer, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_BCAST(filenz, 1, MPI_integer, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_BCAST(filemv, 1, MPI_integer, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_BCAST(filedx, 1, MPI_REAL8  , 0, MPI_COMM_WORLD, mpi_err)

    call MPI_BCAST(fileprotons , 1, MPI_integer, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_BCAST(fileneutrons, 1, MPI_integer, 0, MPI_COMM_WORLD, mpi_err)

    call MPI_BCAST(filenwn           ,1, MPI_integer,0, MPI_COMM_WORLD, mpi_err)
    call MPI_BCAST(filenwp           ,1, MPI_integer,0, MPI_COMM_WORLD, mpi_err)
    call MPI_BCAST(filenwt           ,1, MPI_integer,0, MPI_COMM_WORLD, mpi_err)
    call MPI_BCAST(fileblocks_global ,8, MPI_integer,0, MPI_COMM_WORLD, mpi_err)

    ! Seemingly useless to BCAST fileversion, but this is necessary for further
    ! logic further down in this routine
    call MPI_BCAST(file_version , 1, MPI_integer, 0, MPI_COMM_WORLD, mpi_err)

    call MPI_BCAST(sym_transfo_needed,1,MPI_LOGICAL, 0, MPI_COMM_WORLD, mpi_err)
#endif    
    ! .. now we have each rank decide what spwfs to take from file
    call loadbalance(fileblocks_global,                              & ! inputs
    &       fileblocks, file_spwf_map, file_rank_map,file_spwf_inverse)! outputs

    ! Arrays like these are stored on all ranks, hence "filenwt"
    allocate(spenergies (filenwt))
    allocate(dispersions(filenwt))
    allocate(HFtransfo  (filenwt,filenwt)) ; HFtransfo   = 0.0d0
    allocate(sphamil(filenwt,filenwt)) ; sphamil = 0.0d0

    if (allocated(rho_can)) deallocate(rho_can)
    allocate(rho_can(filenwt))

    if(MPI_RANK.eq.0) then
      read(chan,iostat=io) spenergies, dispersions
      if(file_version .ge. 3) then
        read(chan,iostat=io) filediagsphamil
        read(chan,iostat=io) HFtransfo
      endif
    endif
#if(USE_MPI>0)
    call MPI_BCAST(spenergies ,filenwt   , MPI_REAL8,0, MPI_COMM_WORLD, mpi_err)
    call MPI_BCAST(dispersions,filenwt   , MPI_REAL8,0, MPI_COMM_WORLD, mpi_err)
    call MPI_BCAST(dispersions,filenwt   , MPI_REAL8,0, MPI_COMM_WORLD, mpi_err)
    call MPI_BCAST(HFTRANSFO  ,filenwt**2, MPI_REAL8,0, MPI_COMM_WORLD, mpi_err)
#endif
    !---------------------------------------------------------------------------
    ! Reading the spwfs from file
    !---------------------------------------------------------------------------
    ! First allocate the needed space
    ! wavefunctions are distributed across ranks ....
    allocate(HFPsi(filenx*fileny*filenz,4, sum(fileblocks)))
    if(file_version .gt. 5) then
      if(MPI_RANK.eq.0) allocate(temp(filenx*fileny*filenz,4)) ! create space
      wfcounter = 0 ! this is the counter tracking how much spwfs each 
                    ! INDIVIDUAL rank has stored so far

      do wave=1, filenwt
        ! Read the spwf into dummy storage
        if(MPI_RANK.eq.0) read(chan,iostat=io) temp
        targetrank = file_rank_map(wave) ! rank to communicate the spwf to
        if(targetrank .eq. 0 .and. MPI_RANK.eq.0) then
            ! No communication is necessary for the spwfs stored on rank 0
            ! This case is ALWAYS executed for serial calculations.
            wfcounter = wfcounter + 1
            if(MPI_rank.eq.0) HFPsi(:,:,wfcounter) = temp
#if(USE_MPI > 0)
        else
          if(MPI_RANK.eq.0) then
            ! rank 0 sends the spwf to targetrank
            call MPI_SEND(temp, 4*filemv, MPI_REAL8, targetrank, 2,            &
            &                                           MPI_COMM_WORLD, mpi_err)
          else if(MPI_RANK .eq. targetrank) then
            ! ... which receives and stores in HFPSI
            wfcounter = wfcounter + 1
            call MPI_RECV(HFpsi(:,:,wfcounter), 4*filemv, MPI_REAL8, 0, 2,     &
              &                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_err)
          endif
#endif
        endif
      enddo
    else
      ! Originally, the .wf files contained the HFPsi array as one unformatted
      ! record. This is kind of unpractical for MPI applications.
      if(NPROCS .gt. 1) call stp('Old .wf files cannot be read with MPI runs.')

      ! We can safely read this in one go; a single rank is present
      read(chan,iostat=io) HFPsi
    endif
#if(USE_MPI > 0)    
    ! Possibly a superfluous barrier call, but good for my peace of mind
    call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)
#endif
    ! End of the most complicated part of this routine; back to easy sequential
    ! reads and some MPI broadcasting
    !---------------------------------------------------------------------------

    !--------------------------------------------------------------------------- 
    ! Name of the force and functional and full s.p. hamiltonian matrix
    if(MPI_RANK.eq.0) then
      read(chan, iostat=io) ini_name_param, func_name_check
      ! Single-particle hamiltonian
      if(file_version.ge.4) then
        read(chan, iostat=io) sphamil
      endif
    endif
#if(USE_MPI > 0)
    call MPI_BCAST(ini_name_param , len(ini_name_param) , MPI_CHARACTER,0,     &
    &                                                   MPI_COMM_WORLD, mpi_err)
    call MPI_BCAST(func_name_check, len(func_name_check), MPI_CHARACTER,0,     &
    &                                                   MPI_COMM_WORLD, mpi_err)

    call MPI_BCAST(sphamil, filenwt**2, MPI_REAL8,0,MPI_COMM_WORLD, mpi_err)
#endif

    !---------------------------------------------------------------------------
    ! Pairing information
    if(MPI_RANK.EQ.0) then
      read(chan, iostat=io) filepairing
      ! Write the occupation factors in all cases
      read(chan, iostat=io) rho_can
    endif

#if(USE_MPI > 0) 
    call MPI_BCAST(filepairing,      1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_BCAST(rho_can    ,filenwt, MPI_REAL8  , 0, MPI_COMM_WORLD, mpi_err)
#endif

    select case (filepairing)
    case(0) !-------------------------------------------------------------------
        ! HF: nothing to read
    case(1) !-------------------------------------------------------------------
        ! BCS calculation: read the gaps
        allocate(filegaps(filenwt,1))
        if(MPI_RANK .eq. 0) then
          read(chan, iostat=io) FermiEnergy       ! Lambda
          read(chan, iostat=io) filegaps
        endif

#if(USE_MPI > 0)
        call MPI_BCAST(Fermienergy,      2, MPI_REAL8,0, MPI_COMM_WORLD,mpi_err)
        call MPI_BCAST(filegaps  ,filenwt, MPI_REAL8,0, MPI_COMM_WORLD, mpi_err)
#endif

        ! Simply copy the gaps for now
        select case(pairingtype)
        case(0)
          ! Do nothing
        case(1)
          allocate(BCSGaps(filenwt)) ;  BCSGaps = filegaps(:,1)
        case(2)
          allocate(HFBgaps(filenwt, filenwt)) ; HFBgaps = 0
          do i=1, filenwt
              HFBgaps(i,i) = filegaps(i,1)
          enddo
        end select
    case(2) !-------------------------------------------------------------------
        ! HFB calculation
        if (allocated(filegaps)) deallocate(filegaps)

        allocate(filegaps(filenwt, filenwt)) 
        allocate(kappa_pairing(filenwt, filenwt)) 
        allocate(rho_pairing(filenwt, filenwt)) 

        if(file_version .gt. 3 ) then
          if(MPI_RANK .eq. 0) then
            read(chan, iostat=io) fileblocktype, fileblocknumber
            read(chan, iostat=io) file_HFB_blocks
          endif
#if(USE_MPI > 0)
          call MPI_BCAST(fileblocktype  , 1, MPI_INTEGER, 0, MPI_COMM_WORLD,   &
          &                                                             mpi_err)
          call MPI_BCAST(fileblocknumber, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,   &
          &                                                             mpi_err)
          call MPI_BCAST(file_HFB_blocks, 8, MPI_INTEGER, 0, MPI_COMM_WORLD,   &
          &                                                            mpi_err)
#endif
          ! Depending on the type of blocking, we read and BCAST different
          select case(fileblocktype) 
          case(0)
            if(MPI_RANK .eq. 0) read(chan, iostat = io)
          case(1,3,5)
            allocate(fileblockindices(fileblocknumber))
            if(MPI_RANK.eq.0) read(chan, iostat = io) fileblockindices
#if(USE_MPI > 0)
            call MPI_BCAST(fileblockindices,fileblocknumber, MPI_INTEGER, 0,   &
            &                                            MPI_COMM_WORLD,mpi_err)
#endif
          case(2,4,6)
            allocate(fileblocklowest(fileblocknumber))
            if(MPI_RANK .eq. 0) read(chan, iostat = io) fileblocklowest
#if(USE_MPI > 0)
            call MPI_BCAST(fileblocklowest,fileblocknumber, MPI_INTEGER, 0,    &
            &                                            MPI_COMM_WORLD,mpi_err)
#endif
          end select
        endif

        if(MPI_RANK .eq. 0) then
          read(chan, iostat=io) FermiEnergy       ! Lambda
          read(chan, iostat=io) rho_pairing
          read(chan, iostat=io) kappa_pairing     ! kappa
          read(chan, iostat=io) ! Canonical transformation
          read(chan, iostat=io) filegaps
        endif

#if(USE_MPI>0)
        call MPI_BCAST(FermiEnergy  ,          2, MPI_REAL8, 0,MPI_COMM_WORLD, &
        &                                                               mpi_err)
        call MPI_BCAST(rho_pairing  , filenwt**2, MPI_REAL8, 0,MPI_COMM_WORLD, &
        &                                                               mpi_err)
        call MPI_BCAST(kappa_pairing, filenwt**2, MPI_REAL8, 0,MPI_COMM_WORLD, & 
        &                                                               mpi_err)
        call MPI_BCAST(filegaps     , filenwt**2, MPI_REAL8, 0,MPI_COMM_WORLD, & 
        &                                                               mpi_err)
#endif
        ! For late-enough versions, we also read the full Bogoliubov 
        ! transformation and the associated configuration matrix.
        if(file_version .ge. 3) then
          readHFBinfofile = .true.
          allocate(Bogoliubov(2*filenwt, 2*filenwt)) 
          allocate(configmatrix(2*filenwt)) 

          ! We simply read these arrays here. If a transformation is needed,
          ! we will deal with it elsewhere.
          if(MPI_RANK .eq. 0) then
            read(chan, iostat=io) Bogoliubov
            if (io.ne.0) then
              call stp('ERROR: reading Bogoliubov transformation from file.')
            endif
            read(chan, iostat=io) configmatrix
            if (io.ne.0) then
              call stp('ERROR: reading the configuration matrix from file.')
            endif
          endif
#if(USE_MPI > 0)
          call MPI_BCAST(Bogoliubov, 4*filenwt**2, MPI_REAL8,0, MPI_COMM_WORLD,&
          &                                                             mpi_err)
          call MPI_BCAST(configmatrix, 2*filenwt , MPI_REAL8,0, MPI_COMM_WORLD,&
          &                                                             mpi_err)
#endif
        endif

        select case(pairingtype)
        case(0) !---------------------------------------------------------------
          ! This calculation is HF: Do nothing
        case(1) !---------------------------------------------------------------
          ! This calculation is HFB: use the diagonal matrix elements of 
          ! filegaps as BCSgaps. 
          allocate(BCSgaps(filenwt)) 
          do i=1, filenwt
            BCSgaps(i) = filegaps(i,i)
          enddo
        case(2) !---------------------------------------------------------------
          ! This calculation is HFB: simply copy the gaps for now
          if (allocated(HFBGaps)) then   
            deallocate(HFBGaps)        
          end if                       
          allocate(HFBGaps(filenwt, filenwt)) 
          HFBGaps = filegaps(1:filenwt, 1:filenwt)  
        end select
    case DEFAULT
      call stp('Something is seriously wrong with the .wf file.')
    end select 

    !---------------------------------------------------------------------------
    ! Cranking information
    if(MPI_RANK.eq.0) then
      if(file_version .gt. 4 ) then
        read(chan, iostat=io) omega_file
      else
        ! File-versions < 4 do not have this line
        read(chan, iostat=io) 
        omega_file = 0.0d0
      endif
      if(io.ne.0) then
        call stp('ERROR in reading cranking line of the wf file.')
      endif
    endif
#if(USE_MPI > 0)
    call MPI_BCAST(omega_file, 3, MPI_REAL8, 0, MPI_COMM_WORLD, mpi_err)
#endif
    ! Using the cranking frequencies read from file
    if(continueCrank) omega = omega_file

    !---------------------------------------------------------------------------
    ! Potentials: note that readpotentials handles all MPI affairs itself
    potentials_read = readpotentials(chan, filenx,fileny,filenz, sym_transfo_needed)
    !-------------------------------------------------------------------------
    ! Multipole moment information
    ! Note: ReadMoment handles all MPI affairs itself
    io = 0
    do while(io.eq.0)
      call ReadMoment(chan,io)
    enddo

    !-------------------------------------------------------------------------
    ! End of reading
    if(MPI_RANK.EQ. 0) close(chan)
    !-------------------------------------------------------------------------
    ! Sanity checks if transformation is not allowed
    if(.not.  AllowTransform) then
      if((filenx.ne.nx).or. (fileny.ne.ny) .or. (filenz.ne.nz)) then
          print 1, filenx, fileny, filenz, nx,ny,nz
          call stp('')
      endif
      if(filenwn.ne.nwn .or. filenwp.ne.nwp) then
          print 2, filenwn, filenwp, nwn, nwp
          call stp('')
      endif
    else
      ! We do not allow modification of the mesh, s.p. wavefunctions and 
      ! symmetry transformations at the same time. 
      check_x = (nx .ne. filenx) .and. (nx .ne. 2*filenx)
      check_y = (ny .ne. fileny) .and. (ny .ne. 2*fileny)
      check_z = (nz .ne. filenz) .and. (nz .ne. 2*filenz)

      check_nwn = (nwn .ne. filenwn) .and. (nwn .ne. 2*filenwn)
      check_nwp = (nwn .ne. filenwn) .and. (nwn .ne. 2*filenwn)

      if(sym_transfo_needed) then
       if(check_x .or. check_y .or. check_z) then 
        call stp("Please don't combine symmetry transformations and mesh modifications.")
       endif  
       if(check_nwn .or. check_nwp ) then 
        call stp("Please don't combine symmetry transformations and adding wavefunctions.")
       endif  
      endif
    endif
  end subroutine read_tantalus_wf
 
#if( $FAM == 0 )
  subroutine write_tantalus_wf(chan, ofn)
    !---------------------------------------------------------------------------
    ! Subroutine that dumps all information to a .wf file for future runs.
    ! Note: this does not (yet) use any MPI I/O operations, it simply relies on
    !       transferring all data to rank 0 and having it do the writing.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !   chan : integer, channel to open file ofn on
    !   ofn  : character, filename to write to. 
    !          If it does not exist, will get created.
    !---------------------------------------------------------------------------
    ! Things written to file. (Not yet implemented ones are indicated by (*) )
    !
    ! Version
    ! Convergence information: E, dE                         (*)
    ! nx,ny,nz,dx,dt                    
    ! Symmetry information                                   (1)
    ! neutrons,protons
    ! nwn, nwp, Number of wavefunctions in every block
    ! spenergies, dispersions
    ! diagsphamil
    ! HFtransfo 
    ! (nwt) Wavefunctions
    ! Forcename
    ! Single-particle hamiltonian
    ! Pairing information 
    !    - Pairingtype
    !    - Rho_can = occupation factors 
    !      (HF)  
    !        |   (nothing)
    !      (BCS) 
    !        |   Fermi level
    !        |   BCSGaps  
    !      (HFB)
    !        |   blocktype, blocknumber
    !        |   block-sizes for HFB solver
    !        |   blocklowest/blockindices
    !        |   Fermi level
    !        |   rho_pairing    
    !        |   kappa_pairing
    !        |   can_transfo
    !        |   HFBgaps      
    !        |   Bogoliubov transformation 
    !        |   Configuration matrix   
    ! CrankingInfo
    ! Potentials                                             (2)
    ! Multipole Moments
    !     | The code writes the data on ALL the multipole moments.
    !     | For the format of the lines, see the Moments module.
    !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Some remarks:
    !  (1) The symmetry information is encoded in a single string, the SYM_CODE.
    !  (2) Potentials are written on multiple lines, see the functional module.
    !---------------------------------------------------------------------------

    use functional
    use moments
    use Cranking

    integer, intent(in)          :: chan
    character(len=*), intent(in) :: ofn
    integer                      :: io, wave, wave_local, rank
#if(USE_MPI > 0)
    integer                      :: mpi_err
    real(KIND=dp), allocatable   :: tempwf(:,:)
#endif
     type(moment), pointer        :: mom

    call start_timer(T_wfoutput)

    open (chan,form='unformatted',file=ofn)
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Purely sequential part of the writing
    if(MPI_RANK .EQ. 0) then
      write(chan, iostat=io) version_number
      ! Convergence information                                (NOT IMPLEMENTED)
      write(chan,iostat=io) 
      !Parameters of the mesh
      write(Chan,iostat=io) nx,ny,nz, dx
      ! Symmetry information                                   
      write(Chan,iostat=io) SYM_CODE
      !Number of protons and neutrons
      write(Chan,iostat=io) neutrons,protons
      ! HFBLocks information (NOTE: this should be the GLOBAL information)
      write(Chan,iostat=io) nwn, nwp, hfblocks_global 
      ! Wavefunctions  
      write(chan,iostat=io) spenergies, dispersions
      ! information on the HF transformation
      write(chan, iostat=io) diagsphamil
      write(chan, iostat=io) HFtransfo
    endif
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Parallel part of the writing
    ! First, make sure the team is complete before proceeding
#if(USE_MPI > 0)
    call MPI_BARRIER(MPI_COMM_WORLD, mpi_err) 
    ! b) making space to receive spwfs from the other ranks
    if(MPI_RANK .eq.0) allocate(tempwf(mv,4))
#endif

    do wave = 1, nwt
        rank = rank_map(wave)           ! spwf wave is stored on which MPI rank?
        wave_local = spwf_inverse(wave) ! ... and has which local index? 

        if(rank.eq.0) then
          ! ------ Rank 0 writes its own wavefunctions ----------------
          if(MPI_RANK.eq.0) write(chan,iostat=io) HFPsi(:,:,wave_local)
#if(USE_MPI > 0)
        else
          ! ...  otherwise there is communication involved ....
          if(MPI_RANK .eq. 0) then
            ! -------rank 0 receives and writes -----------------------
            call MPI_RECV(        tempwf, 4*mv, MPI_REAL8, rank, 2, &
            &                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_err)
            write(chan,iostat=io) tempwf
          elseif(MPI_RANK .eq. rank) then
            ! -------rank "rank" sends ------- -----------------------
            call MPI_SEND(hfpsi(:,:,wave_local), 4*mv, MPI_REAL8, 0, 2,        &
            &                                           MPI_COMM_WORLD, mpi_err)
          endif
#endif
        endif
    enddo

    ! e) freeing up the space
#if(USE_MPI > 0)
    if(MPI_RANK .eq.0) deallocate(tempwf)
#endif
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Back to the sequential part of the writing
    if(MPI_RANK .EQ. 0) then
      ! Name of the force.
      write(chan, iostat=io) name_param, func_name
      ! Single-particle hamiltonian
      write(chan, iostat=io) sphamil
      !-------------------------------------------------------------------------
      ! Pairing information 
      write(chan, iostat=io) PairingType

      ! Write the occupation factors in all cases
      write(chan, iostat=io) rho_can

      select case (PairingType)
      case(0)
          ! HF: nothing to write
      case(1)
          ! BCS
          write(chan, iostat=io) FermiEnergy
          write(chan, iostat=io) BCSGaps 
      case(2)
          ! HFB
          write(chan, iostat=io) blocktype, blocknumber
          if(pairingscheme.eq.1) then
            write(chan, iostat=io) grad_blocks
          else
            write(chan, iostat=io) HFBlocks_global
          endif

          select case( blocktype)
          case(0)
            write(chan, iostat=io) 
          case(2,4,6)
            write(chan, iostat=io) blocklowest
          case(1,3,5)
            write(chan, iostat=io) blockindices
          end select

          write(chan, iostat=io) FermiEnergy       ! Lambda
          write(chan, iostat=io) rho_pairing       ! rho
          write(chan, iostat=io) kappa_pairing     ! kappa
          write(chan, iostat=io) Cantransfo        ! Canonical transformation
          write(chan, iostat=io) HFBgaps           ! Full matrix of gaps
          write(chan, iostat=io) Bogoliubov        ! Bogoliubov transformation
          write(chan, iostat=io) configmatrix      ! Configuration matrix
      end select
      ! Cranking information: frequencies in all Cartesian directions 
      write(chan, iostat=io) Omega(1:3)
      !-------------------------------------------------------------------------
      ! Potentials on file
      call writepotentials(chan,potentials)
      !-------------------------------------------------------------------------
      ! Multipole moment information                             
      !
      ! The Cray compilers on LUCIA want to inline the WriteMoment function while
      ! also flattening the linked list of multipole moments when optimisation 
      ! options -O2 or above are used. For reasons I do not understand, this 
      ! makes the executable segfault. Since this routine has absolutely no impact
      ! on execution time, I simply forbid the CRAY compiler to inline this function. 
      ! This magically solves the issue (which does not exist for ifort or gnu compilers) 
      ! 
      ! Cray version on LUCIA at the time of writing:
      ! Cray Fortran : Version 14.0.3
      ! 
      ! Note the double dollar-sign, to make sure Hephaestos does not replace these
      ! compiler directives. 
      ! 
      mom => root
      do while(associated(mom%next))
        mom => mom%next
        !DIR$$ NOINLINE
        call Writemoment(mom,chan)
        !DIR$$ INLINE
      enddo
    endif
    close(chan)

    call stop_timer(T_wfoutput)
  end subroutine write_tantalus_wf
#endif

#if (USE_HDF5 > 0 && $FAM == 0)
  subroutine write_tantalus_hdf5(ofn)
    !------------------------------------------------------------------------------------------
    ! Subroutine that writes a .hdf5 file to warmstart future runs or to do more
    ! analysis.
    !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    !  Groups
    !   /                   : general information
    !   |
    !   | > /wavefunctions/ : single-particle wavefunctions and their properties 
    !   | | > /wavefunctions/hfbasis/  : spwf-related quantities in the Hartree-Fock basis 
    !   | | > /wavefunctions/compbasis/: spwf-related quantities in the computational basis 
    !   | | > /wavefunctions/canbasis/ : spwf-related quantities in the canonical basis
    !   | 
    !   | > /fields/        : functions on the mesh
    !   | | > /fields/densities/   : mean-field densities on the mesh
    !   | | > /fields/potentials/  : mean-field potentials on the mesh
    !   | 
    !   | > /multipoles/     : quantities related to multipole moments    
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Technical note: the spwfs can be written to file in parallel when MPI
    ! is enabled, all other stuff is written serially by rank #0.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !   ofn  : character, filename to write to. 
    !          If the file does not exist, will get created.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! TODO
    ! - [ ] add planned values
    ! - [ ] add MPI BCAST calls
    ! - [ ] add MPI_RANK == 0 selection
    !------------------------------------------------------------------------------------------
    use HDF5
    use timing    , only : start_timer, stop_timer, T_wfoutput
    use functional, only : potentials 
    use densities , only : Density

    ! Outputfilename 
    character(len=*), intent(in)  :: ofn
    ! HDF5 integer identifiers
    integer(HID_T)                :: file_id, group_id, subgroup_id, dset_id, plist_id, space_id 
    ! error handling
    integer                       :: h5ferr
    integer(hsize_t),dimension(3) :: dims,data_dims

    call start_timer(T_wfoutput)
    ! - - - - - - - - - - - - - - -
    !  Administrative steps 
    ! - - - - - - - - - - - - - - -
    ! Initialize hdf5 interface 
    call h5open_f(h5ferr)
    ! Create new HDF5 file
    call h5fcreate_f(ofn, H5F_ACC_TRUNC_F, file_id, h5ferr)
    if(h5ferr.ne.0) call stp('Error creating HDF5 file.')

    call h5gcreate_f(file_id, 'fields'                 , group_id, h5ferr)    
    call h5gcreate_f(file_id, 'fields/densities'       , subgroup_id, h5ferr) ; call h5gclose_f(subgroup_id, h5ferr)
    call h5gcreate_f(file_id, 'fields/potentials'      , subgroup_id, h5ferr) ; call h5gclose_f(subgroup_id, h5ferr)
    call h5gclose_f(group_id, h5ferr)

    call h5gcreate_f(file_id, 'multipoles'             , group_id, h5ferr)   ; call h5gclose_f(group_id, h5ferr)
    ! - - - - - - - - - - - - - - -
    ! Start actual writing, organised per group in subroutines
    call write_hdf5_attributes(file_id)
    call write_hdf5_wavefunctions(file_id)
    call write_hdf5_fields(file_id, Density, Potentials)
    !call write_hdf5_multipoles(file_id)
    ! - - - - - - - - - - - - - - - -
    ! More technical closing steps
    call h5fclose_f(file_id,h5ferr)
    ! Close FORTRAN interface
    call h5close_f(h5ferr)

    call stop_timer(T_wfoutput)

  end subroutine write_tantalus_hdf5

  subroutine write_hdf5_attributes(file_id)
    !-----------------------------------------------------------------------------
    ! Subroutine that writes the attributes of the root group of the HDF5 file; 
    ! see the documentation for a full list of attributes written.
    !
    ! Input:
    !   file_id : integer(HID_T), identifier of the root group of the HDF5 file
    !-----------------------------------------------------------------------------
    use HDF5 
    use HDF5_auxiliary, only : hdf5_write_attr_integer, hdf5_write_attr_integer_1d
    use HDF5_auxiliary, only : hdf5_write_attr_double, hdf5_write_attr_double_1d
    use HDF5_auxiliary, only : hdf5_write_attr_char 

    use geninfo,        only : dx, nx, ny, nz, protons, neutrons
    use functional,     only : func_name, name_param
    use pairing,        only : PairingType, FermiEnergy
    use pairing,        only : BlockLowest, BlockType, BlockIndices
    use wavefunctions,  only : nwn, nwp, nwt, HFBlocks_global
    use cranking,       only : omega

    integer(HID_T), INTENT(IN) :: file_id 

    ! Versioning information
    call hdf5_write_attr_integer(file_id, 'version_number', version_number)
    ! Description; empty for now
    call hdf5_write_attr_char(file_id, 'description', ' ')
    ! Number of protons and neutrons
    call hdf5_write_attr_double(file_id, 'neutrons', neutrons)
    call hdf5_write_attr_double(file_id, 'protons',  protons)

    ! Type of Skyrme EDF 
    call hdf5_write_attr_char(file_id, 'func_name', trim(func_name))
    ! Parameterization name
    call hdf5_write_attr_char(file_id, 'name_param', trim(name_param))  
    ! Pairing ansatz
    call hdf5_write_attr_integer(file_id, 'PairingType', PairingType)  

    ! Parameters of the mesh: nx,ny,nz,dx,dy,dz
    call hdf5_write_attr_integer(file_id, 'nx', nx)
    call hdf5_write_attr_integer(file_id, 'ny', ny)
    call hdf5_write_attr_integer(file_id, 'nz', nz)
    ! A future generalization to non-cubic meshes would change these
    ! numbers; it is already in the .HDF5 files to be compatible with
    ! other mesh codes.
    call hdf5_write_attr_double( file_id, 'dx', dx)
    call hdf5_write_attr_double( file_id, 'dy', dx)
    call hdf5_write_attr_double( file_id, 'dz', dx)

    ! Type of boundary conditions/mesh
#if(USE_Periodic == 0)
    call hdf5_write_attr_char(file_id, 'BCtype'   , 'anti-periodic')
    call hdf5_write_attr_char(file_id, 'mesh_type', 'half-integer')
#else
    call hdf5_write_attr_char(file_id, 'BCtype'   , 'periodic')
    call hdf5_write_attr_char(file_id, 'mesh_type', 'integer')
#endif

#if(PASTA > 0)
    call hdf5_write_attr_char(file_id, 'calctype', 'PASTA')
#else 
    call hdf5_write_attr_char(file_id, 'calctype', 'NUCLEI')
#endif

    ! Number of wavefunctions for neutrons and protons
    call hdf5_write_attr_integer(file_id, 'nwn', nwn)
    call hdf5_write_attr_integer(file_id, 'nwp', nwp)
    call hdf5_write_attr_integer(file_id, 'nwt', nwt)
    call hdf5_write_attr_integer_1d(file_id, 'HFBlocks', HFBlocks_global,8)

    ! Symmetry information                                   
    call hdf5_write_attr_char(file_id, 'SYM_CODE', SYM_CODE)

    ! Pairing information
    call hdf5_write_attr_double_1d(file_id, "FermiEnergy", FermiEnergy, 2)
    ! TODO: include blocking information here!

    ! Cranking information
    call hdf5_write_attr_double_1d(file_id, "Omega", omega, 3)
  end subroutine write_hdf5_attributes

  subroutine write_hdf5_wavefunctions(file_id)
    !------------------------------------------------------------------------------------
    ! Subroutine that writes all information on the single-particle wavefunctions to 
    ! the /wavefunctions/ group and its subgroups on a HDF5 file.
    !
    ! Input:
    !   file_id     : integer(HID_T), identifier of the root group of the HDF5 file
    !
    ! Output:
    !   None
    !------------------------------------------------------------------------------------
    use HDF5
    use HDF5_auxiliary, only : hdf5_write_attr_double_1d,  hdf5_write_dataset_1d
    use wavefunctions,  only : HFPsi, nwt, dispersions, spenergies
    use pairing,        only : rho_can
    use BCS,            only : BCSgaps
    use HFB,            only : HFBgaps

    integer(HID_T), INTENT(IN) :: file_id
    integer(HID_T)             :: space_id, plist_id, dset_id, wf_id, hf_id, comp_id, can_id
    integer(hsize_t)           :: dims(3), dims_1d(1) !HDF5 requires a specific type of integer
    integer                    :: h5ferr

    ! Create groups for wavefunctions, fields and multipoles with associated subgroups 
    call h5gcreate_f(file_id, 'wavefunctions', wf_id, h5ferr)               
    call h5gcreate_f(file_id, 'wavefunctions/hfbasis'  , hf_id, h5ferr)  ; call h5gclose_f(hf_id, h5ferr)
    call h5gcreate_f(file_id, 'wavefunctions/compbasis', comp_id, h5ferr); call h5gclose_f(comp_id, h5ferr)
    call h5gcreate_f(file_id, 'wavefunctions/canbasis' , can_id, h5ferr) ; call h5gclose_f(can_id, h5ferr)
    call h5gclose_f(wf_id, h5ferr)

    !------------------------------------------------------------------------------------
    ! Writing the actual wavefunctions
#if(USE_MPI > 0)
    call stp('Parallel IO for wavefunctions not implemented yet.')
#else
    dims=(/nx*ny*nz,4,nwt/)
    ! Create dataspace for data_set 
    call h5screate_simple_f(3, dims, space_id, h5ferr) 
    if(h5ferr .ne. 0) call stp('Error creating dataspace for wavefunctions.')
    call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, h5ferr)
    if(h5ferr .ne. 0) call stp('Error creating property list for wavefunctions.')
    ! Create dataset with default properties "dset_id" is returned
    call h5dcreate_f(file_id,'/wavefunctions/states',H5T_NATIVE_DOUBLE,space_id,dset_id,h5ferr,plist_id)
    if(h5ferr .ne. 0) call stp('Error creating dataset for wavefunctions.')
    ! close access to plist
    call h5pclose_f(plist_id, h5ferr)
    ! close access to dataspace
    call h5sclose_f(space_id, h5ferr)
    ! Write dataset sequentally 
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, HFpsi, dims, h5ferr)
    if(h5ferr .ne. 0) call stp('Error writing wavefunctions to HDF5 file.')
    call h5dclose_f(dset_id, h5ferr)
#endif

    !------------------------------------------------------------------------------------
    ! Other information
    if(MPI_Rank.eq.0) then
      dims_1d = (/nwt/)
      ! Information on the computational basis
      select case(PairingType)
      case(0) ! HF
        ! Nothing to write for now
      case(1) ! BCS
        call hdf5_write_dataset_1d(file_id,'BCSgaps' , BCSgaps,nwt,'wavefunctions/compbasis/')
      case(2) ! HFB 
        ! Nothing to write for now
      end select
      ! Information on the HF basis
      call hdf5_write_dataset_1d(file_id,'spenergies' , spenergies,nwt,'wavefunctions/hfbasis/')
      call hdf5_write_dataset_1d(file_id,'dispersions',dispersions,nwt,'wavefunctions/hfbasis/')
      ! Information on the canonical basis
      call hdf5_write_dataset_1d(file_id,'rho_can', rho_can,nwt,'wavefunctions/canbasis/')
      !call h5screate_simple_f(1, dims_1d, space_id, h5ferr) 
      !call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, h5ferr)
      !call h5dcreate_f(file_id,'/wavefunctions/canbasis/rho_can',H5T_NATIVE_DOUBLE,space_id,dset_id,h5ferr,plist_id)
      !call h5pclose_f(plist_id, h5ferr)
      !call h5sclose_f(space_id, h5ferr)
      !call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, rho_can, dims_1d, h5ferr)
      !call h5dclose_f(dset_id, h5ferr)
    ! Information on the canonical basis
      ! TODO
    endif

  end subroutine write_hdf5_wavefunctions 

  subroutine write_hdf5_fields(file_id,R,F)
    !------------------------------------------------------------------------------------
    ! Subroutine that writes all information on densities and potentials to the 
    ! /fields/ group and its subgroups on the HDF5 file.
    ! 
    ! Input:
    !   file_id : integer(HID_T), identifier of the root group of the HDF5 file
    !   R       : densityvector, densities to be written
    ! Output:
    !   None
    !
    ! TODO
    ! - [ ] add MPI BCAST calls
    ! - [ ] add MPI_RANK == 0 selection
    ! - [ ] add 'linguistic' links 
    !------------------------------------------------------------------------------------
    use HDF5
    use densities, only : write_hdf5_densities, DensityVector
    use functional, only : write_hdf5_potentials, PotentialVector

    integer(HID_T), INTENT(IN)        :: file_id 
    type(DensityVector), INTENT(IN)   :: R
    type(PotentialVector), INTENT(IN) :: F

    call write_hdf5_densities(file_id,  R)
    call write_hdf5_potentials(file_id, F)

  end subroutine write_hdf5_fields

  subroutine write_hdf5_multipoles(file_id)
    !------------------------------------------------------------------------------------
    ! Subroutine that writes all information on the multipole operators to 
    ! /multipoles/ group on the HDF5 file.
    !
    ! Input:
    !   file_id : integer(HID_T), identifier of the root group of the HDF5 file
    !
    ! Output:
    !   None
    !------------------------------------------------------------------------------------
    use HDF5
    use moments

  integer(HID_T), INTENT(IN) :: file_id 

  end subroutine write_hdf5_multipoles
#endif 

#if( USE_HDF5 > 0)
  subroutine read_tantalus_hdf5(ifn, sym_transfo_needed)
    !------------------------------------------------------------------------------------------
    ! Subroutine that reads a .hdf5 file to warmstart future runs or to do more analysis.
    !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !   ifn  : character, filename to read from. (short for InputFileName)
    !
    ! Output:
    !   sym_transfo_needed : logical, whether symmetry transformations are needed
    !                        to adapt the wavefunctions on file to the current
    !                        symmetry settings.
    !------------------------------------------------------------------------------------------
    use HDF5
    use timing    , only : start_timer, stop_timer, T_wfinput
    use functional, only : potentials_read
    
    character(len=*), intent(in)  :: ifn
    logical, intent(out)          :: sym_transfo_needed
    integer                       :: h5ferr, i
    integer(HID_T)                :: file_id, root_id, dset_id, plist_id, space_id
    logical                       :: exists

    1 format ('Number of mesh points does not correspond to file.',    / &
    &         'On file: nx= ', i3, ' ny= ', i3, ' nz= ',i3,            / &
    &         'In data: nx= ', i3, ' ny= ', i3, ' nz= ',i3)
    2 format ('Number of wavefunctions does not correspond to file.',  / &
    &         'On file: nwn= ', i3, ' nwp=', i3,                       / &
    &         'In data: nwn= ', i3, ' nwp=', i3)

    3 format (' The symmetry choices on file cannot be handled by this executable.')
    4 format (' SYM_CODE   = ', a26)
    5 format (' TRANS_CODE = ', a26)
    6 format (' ON FILE    = ', a26)

    call start_timer(T_wfinput)

    !Initialize hdf5 interface 
    call h5open_f(h5ferr)

    ! First check if the file exists.
    inquire(file=ifn, exist=exists)
    if(.not.exists) then
      call stp('The input file you asked for does not exist! \n' // 'Input: ' // trim(ifn))
    endif
    !open the file
    call h5fopen_f(ifn,H5F_ACC_RDONLY_F,file_id,h5ferr)
    if(h5ferr.ne.0) call stp('Error opening HDF5 file.')

    ! Start actual reading
    call read_hdf5_attributes(file_id, sym_transfo_needed)
    call read_hdf5_fields(file_id, potentials_read)
    call read_hdf5_wavefunctions(file_id)
    ! TODO: add multipoles
    !call write_hdf5_multipoles(file_id)

    ! Close file and FORTRAN interface
    call h5fclose_f(file_id, h5ferr)
    call h5close_f(h5ferr)

    call stop_timer(T_wfinput)
  end subroutine read_tantalus_hdf5

  subroutine read_hdf5_attributes(file_id, sym_transfo_needed)
    !-----------------------------------------------------------------------------
    !
    ! Attention: surely not all information on file is consistently acted upon.
    !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !   file_id : integer(HID_T), identifier of the root group of the HDF5 file
    !
    ! Output:
    !   sym_transfo_needed : logical, whether symmetry transformations are needed
    !                        to adapt the wavefunctions on file to the current
    !                        symmetry settings.
    !-----------------------------------------------------------------------------
    use HDF5
    use HDF5_auxiliary, only : hdf5_read_attr_integer, hdf5_read_attr_integer_1d
    use HDF5_auxiliary, only : hdf5_read_attr_double , hdf5_read_attr_double_1d
    use HDF5_auxiliary, only : hdf5_read_attr_char 
    
    integer(HID_T), INTENT(IN)   :: file_id 
    character(len=26)            :: SYM_CODE_CHECK
    character(len=20)            :: func_name_check
    logical, INTENT(OUT)         :: sym_transfo_needed

    1 format ('Number of mesh points does not correspond to file.', / &
    &         'On file: nx= ', i3, ' ny= ', i3, ' nz= ',i3,            / &
    &         'In data: nx= ', i3, ' ny= ', i3, ' nz= ',i3)
    2 format ('Number of wavefunctions does not correspond to file.', / &
    &         'On file: nwn= ', i3, ' nwp=', i3,                      / &
    &         'In data: nwn= ', i3, ' nwp=', i3)

    3 format (' The symmetry choices  on file cannot be handled.')
    4 format (' SYM_CODE   = ', a26)
    5 format (' TRANS_CODE = ', a26)
    6 format (' ON FILE    = ', a26)

    ! Versioning information
    call hdf5_read_attr_integer(file_id, 'version_number', file_version)

    ! Number of protons and neutrons
    call hdf5_read_attr_double(file_id, 'neutrons', fileneutrons)
    call hdf5_read_attr_double(file_id, 'protons',  fileprotons)

    ! Type of Skyrme EDF 
    call hdf5_read_attr_char(file_id, 'func_name', func_name_check,len(func_name_check,kind=8))
    ! Parameterization name
    call hdf5_read_attr_char(file_id, 'name_param', ini_name_param,len(ini_name_param,kind=8))  
    ! Pairing ansatz
    call hdf5_read_attr_integer(file_id, 'PairingType', filepairing)  

    ! Parameters of the mesh: nx,ny,nz,dx,dy,dz
    call hdf5_read_attr_integer(file_id, 'nx', filenx)
    call hdf5_read_attr_integer(file_id, 'ny', fileny)
    call hdf5_read_attr_integer(file_id, 'nz', filenz)

    ! Read dx from file; note that information on dy and dz is discarded  
    call hdf5_read_attr_double( file_id, 'dx', filedx)

    ! Number of wavefunctions for neutrons and protons
    call hdf5_read_attr_integer(file_id, 'nwn', filenwn)
    call hdf5_read_attr_integer(file_id, 'nwp', filenwp)
    call hdf5_read_attr_integer(file_id, 'nwt', filenwt)
    call hdf5_read_attr_integer_1d(file_id, 'HFBlocks', fileblocks_global,8)

    ! Symmetry information                                   
    call hdf5_read_attr_char(file_id, 'SYM_CODE', SYM_CODE_CHECK, len(SYM_CODE_CHECK,kind=8))

    if(SYM_CODE_CHECK .eq. SYM_CODE) then
        sym_transfo_needed = .false.
    elseif(SYM_CODE_CHECK .eq. TRANS_CODE) then
        sym_transfo_needed = .true.
    else
      print 3
      print 4, SYM_CODE 
      print 5, TRANS_CODE
      print 6, SYM_CODE_CHECK
      call stp('')
    endif

    ! Pairing information
    call hdf5_read_attr_double_1d(file_id, "FermiEnergy", FermiEnergy, 2)
    ! TODO: include blocking information here!

    ! Cranking information
    !call hdf5_write_attr_double_1d(file_id, "Omega", omega, 3)
  end subroutine read_hdf5_attributes

  subroutine read_hdf5_fields(file_id, F)
    !------------------------------------------------------------------------------------
    ! Subroutine that reads information on densities and potentials from the
    ! /fields/ group and its subgroups on a HDF5 file.
    !
    ! Note: not all information is used, case in point being the mean-field densities.
    !
    ! Input:
    !   file_id : integer(HID_T), identifier of the root group of the HDF5 file
    !
    ! Output:
    !   F       : PotentialVector that will contain the potentials read from file
    !------------------------------------------------------------------------------------
    use HDF5
    use functional, only: PotentialVector, read_potentials_hdf5

    integer(HID_T), INTENT(IN)         :: file_id
    type(PotentialVector), INTENT(out) :: F

    F = read_potentials_hdf5(file_id, filenx, fileny, filenz, sym_transfo_needed)

    ! TODO: add MPI BCASTS here
  end subroutine read_hdf5_fields

  subroutine read_hdf5_wavefunctions(file_id) 
    !------------------------------------------------------------------------------------
    ! Subroutine that reads all information on the single-particle wavefunctions from
    ! the /wavefunctions/ group and its subgroups on a HDF5 file.
    !
    ! TODO: 
    ! - [ ] add MPI BCAST calls
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !   file_id     : integer(HID_T), identifier of the root group of the HDF5 file
    !   
    ! Output:
    ! 
    !------------------------------------------------------------------------------------
    use geninfo, only        : dv
    use HDF5
    use HDF5_auxiliary, only : hdf5_read_dataset_1d
    use wavefunctions,  only : HFPsi, dispersions, spenergies, loadbalance
    use pairing,        only : rho_can

    integer(HID_T), INTENT(IN) :: file_id
    integer(HID_T)             :: dset_id, group_id
    integer(hsize_t)           :: dims(3), dims_1d(1) 
    integer                    :: h5ferr, i

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! First, make a call to loadbalance in order to set all relevant arrays
    call loadbalance(fileblocks_global,                               &! inputs
    &       fileblocks, file_spwf_map, file_rank_map,file_spwf_inverse)! outputs

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! wavefunctions/compbasis/
    call h5gopen_f(file_id, '/wavefunctions/compbasis', group_id, h5ferr)
    select case(filepairing)
    case(0) ! HF
      ! No pairing gaps to read
    case(1) ! BCS
      allocate(BCSgaps(filenwt))
      call hdf5_read_dataset_1d(group_id, 'BCSgaps'     , BCSgaps, filenwt)
    case(2) ! HFB
    ! TODO
    end select
    call h5gclose_f(group_id, h5ferr)
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! wavefunctions/canbasis/
    if(allocated(rho_can)) deallocate(rho_can)
    allocate(rho_can(filenwt))
    call h5gopen_f(file_id, '/wavefunctions/canbasis', group_id, h5ferr)
    call hdf5_read_dataset_1d(group_id, 'rho_can'    ,     rho_can, filenwt)
    call h5gclose_f(group_id, h5ferr)

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! wavefunctions/hfbasis/
    allocate(spenergies (filenwt), dispersions(filenwt))
    call h5gopen_f(file_id, '/wavefunctions/hfbasis', group_id, h5ferr)
    call hdf5_read_dataset_1d(group_id, 'spenergies' ,  spenergies, filenwt)
    call hdf5_read_dataset_1d(group_id, 'dispersions', dispersions, filenwt)
    call h5gclose_f(group_id, h5ferr)
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Reading of the actual single-particle wavefunctions is kept for last
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    call h5dopen_f(file_id, '/wavefunctions/states', dset_id, h5ferr)
    if (h5ferr.ne.0) then
      call stp('ERROR: opening /wavefunctions/states dataset in hdf5 format')
    endif

#if(USE_MPI > 0)
    call stp('Parallel IO for wavefunctions not implemented yet.')
#else
    allocate(HFPsi(filenx*fileny*filenz,4, sum(fileblocks)))
    dims(1)=filenx*fileny*filenz
    dims(2)=4
    dims(3)=filenwt

    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, HFpsi, dims, h5ferr)

    print *, 'READING', fileblocks
    do i =1, sum(fileblocks)
      print *, 'NORM', i, sqrt(sum(HFPsi(:,:,i)**2)*dv)
    enddo

#endif
    call h5dclose_f(dset_id, h5ferr)
    if (h5ferr.ne.0) then
      call stp('ERROR: reading wafefunctions in hdf5 format')
    endif
  end subroutine read_hdf5_wavefunctions
#endif

  function check_blocking_structure() result(passed)
    !---------------------------------------------------------------------------
    !
    ! Some tests on the blocking structure on file vs. that demanded by the user. 
    !
    ! Sanity check
    !  #1 : no blocktypes that are not 0/2/4.
    !  #2 : blocklowest on file == blocklowest input by the user 
    !       * modulo permutations
    !  #3 : the number of blocked states in the Bogoliubov transformation in 
    !       each block matches the input and file
    !---------------------------------------------------------------------------
    use pairing, only: BlockLowest, blocktype

    integer :: i,j , NB, check_blocks(8), B
    integer, allocatable       :: check(:)
    logical                    :: identical, passed

    passed = .true.
    !---------------------------------------------------------------------------
    ! # 1 : no blocktypes that are not 0/2/4.
    if(blocktype.ne.0 .and. blocktype.ne.2 .and. blocktype.ne. 4) then
      call stp('Subroutine check_blocking_structure cannot deal (yet) with &
             &  blocktypes that are not 0/2/4.')
    endif    
    !---------------------------------------------------------------------------
    ! # 2 :  blocklowest on file == blocklowest input by the user 
    !        modulo permutations 
    if(allocated(fileblocklowest) .and. (.not. allocated(blocklowest))) then
      call stp('Blocklowest not allocated, while fileblocklowest is.')
    endif

    if(allocated(blocklowest) .and. (.not. allocated(fileblocklowest))) then
      call stp('Blocklowest allocated, while fileblocklowest is not.')
    endif
    
    if(size(fileblocklowest).ne.size(blocklowest)) then
      print *, ' Size of blocklowest on file:  ', size(fileblocklowest)
      print *, ' Size of blocklowest in input: ', size(blocklowest)
      call stp('')
    endif 
    
    NB = size(fileblocklowest)
    allocate(check(NB)) ; check = 0
    do i=1,NB
      do j=1, NB
        if(blocklowest(j) .eq. fileblocklowest(i)) then
          check(i) = check(i) + 1
        endif
      enddo
    enddo

    identical = .true.
    do i=1, NB
       if(check(i).ne.1) then
          identical = .false.
       endif
    enddo
    
    if(.not. identical) then
      print *, 'Blocklowest on file : ', fileblocklowest
      print *, 'Blocklowest on input: ', blocklowest
      print *, 'These are not identical.'
      call stp('')
    endif
    !---------------------------------------------------------------------------
    ! # 3: Check if the blocking structure on file actually matches the 
    !      structure asked for
    
    !---------------------------------------------------------------------------
    ! This method has turned out to NOT be a reliable indicator.
!    blocked_blocks =  figure_out_blocking_structure_agnostic(                  &
!    &                             sphamil, HFBgaps, FermiEnergy, Bogoliubov)
!  
!    check_blocks = 0
!    do i=1,NB
!      select case(blocklowest(i))
!      case('n+')
!        B = 1       
!      case('n-')
!        B = 3       
!      case('p+')
!        B = 5       
!      case('p-')
!        B = 7       
!      end select 
!      check_blocks(B) = check_blocks(B) + 1 
!    enddo
!  
!    do B=1,8
!      if(check_blocks(B).ne.blocked_blocks(B)) then
!        print *, 'Blocking structure of the Bogoliubov transformation on file'
!        print *, 'does not match that reported by the file.'
!        print *, ' Blocking structure of Bogoliubov matrix: ', blocked_blocks      
!        print *, ' Blocking structure asked for           : ', check_blocks      
!        stop
!      endif
!    enddo

    if(.not. sym_transfo_needed) then
      !---------------------------------------------------------------------------
      ! Instead, we check the "effective" block sizes. 
      ! A reference unblocked calculation will have HFBlocks = grad_blocks
      
      ! 04/01/2022: NOTE, that we can only do this for cases where we need NO
      !             symmetry transformation. If you transform symmetries, you
      !             are on your own! 
      
      check_blocks = HFBlocks
      
      do i=1,NB
        select case(blocklowest(i))
        case('n+')
          B = 1       
        case('n-')
          B = 3       
        case('p+')
          B = 5       
        case('p-')
          B = 7       
        end select 
        check_blocks(B)   = check_blocks(B)   - 1 
        check_blocks(B+1) = check_blocks(B+1) + 1 
      enddo

      do B=1,8
        if(file_HFB_blocks(B)+extraspwfs(B).ne.check_blocks(B)) then
          print *, 'Blocking structure of the Bogoliubov transformation on file'
          print *, 'does not match that reported by the file.'
          print *, ' Block structure of Bogoliubov matrix: ', file_HFB_blocks      
          print *, ' Block structure asked for           : ', check_blocks      
          call stp('')
        endif
      enddo
    endif
  end function check_blocking_structure

end module IO_wf
