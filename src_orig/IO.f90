module IO
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
 ! High-level module governing all aspects of in- and output. 
 ! The lower-level functions that this module relies on are grouped into
 !   - hdf5_auxiliary.f90 : auxiliary routines for hdf5 reading/writing
 !   - IO_wf.f90          : reading and writing of wavefunction files 
 !   - IO_aux.f90         : reading and writing of any other files
 !------------------------------------------------------------------------------
 !
 ! Hephaestos keywords:
 !
 !
 ! Parameters to pass into iniwavefunctions
 ! ININX  : $ININX
 ! ININY  : $ININY
 ! ININZ  : $ININZ
 !
 ! ININWN : $ININWN
 ! ININWP : $ININWP
 ! ININWT : $ININWT
 ! 
 !  TR    : $TR
 ! NTR    : $NTR
 !
 ! FAM    : $FAM
 !==============================================================================

use geninfo
use wavefunctions
use pairing
use functional
use momentsofinertia
use moments
use Coulombmod
use transform
use fission_MOI

use IO_wf, only: SYM_CODE, TRANS_CODE, allowtransform, extraspwfs
use IO_wf, only: version_number, file_version

#if (USE_HDF5 > 0)
use HDF5
#endif

implicit none


  !-----------------------------------------------------------------------------
  ! Filenames for in- and output of the code with respect to spwfs.
  character(len=100)  :: inputfilename, outputfilename
  ! Flag governing the reading of potentials from file
  ! If .true.  => attempt to read the potentials from file and use them
  !               to start iterating
  logical             :: potentials_from_file = .true.
  ! Signal the code to write extra output.
  character(len=100)   :: BXLFIT='', COMBI='', denfile='', potfile=''
  character(len=80)   :: sphffile='', spcanfile='', tofile='', blockfile=''
  character(len=80)   :: inertfile='', famfile='', xyfile='', xyinfile=''
  character(len=80)   :: finfile=''
  ! Signal the code to write the wavefunctions periodically to disk
  integer             :: checkpointiter = 0  

  logical                       :: passed_block_test = .true.

contains

  subroutine ReadInput(file_number, input_file)
    !---------------------------------------------------------------------------
    ! Subroutine to read all the data from the specified file (via the
    ! specified channel) or from STDIN if the variables are not present.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !   filenumber: (optional) integer, channel number 
    !   input_file: (opional) character, filename to look for input on
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! MPI represents a bit of a bookkeeping problem: namelist reading should
    ! be done by only one of the MPI ranks with results broadcasted to the rest.
    ! My philosophy here: 
    ! - To make things easier when adding/removing variables in the future, 
    !   I decided to do the MPI bookkeeping in each separate routine.
    ! - I broadcast ALL INPUT VARIABLES to ALL ranks, even if many variables
    !   will be acted upon by just one single rank. 
    !---------------------------------------------------------------------------
    use geninfo,       only : ReadGenInfo
    use evolution,     only : ReadEvolution
    use wavefunctions, only : ReadWFdata
    use scfiteration,  only : readscfiteration
    use moments,       only : readmomentdata
    use functional,    only : readfunctional
    use pairing,       only : initpairing
#if( $FAM == 1)
    use fam,           only : readfam
#endif
    implicit none

    ! These inputs control where the code will look for its input. Leaving them 
    ! empty will have the code rely on STDIN for input.
    integer(dp), intent(in), optional   :: file_number   
    character(26), intent(in), optional :: input_file 

    logical :: exists
#if(USE_MPI>0)
    integer :: mpi_err
#endif
    
    if(present(file_number)) then
      inquire(file=input_file, exist=exists)
      if(.not. exists) then
        print *, 'Specified input file does not exist!'
        call stp('')
      endif
      open(unit=file_number, file=input_file) 
    endif

    call ReadGenInfo(file_number)
    call readfunctional(file_number)
    call initpairing(file_number)
    call ReadEvolution(file_number)
    call ReadSCFIteration(file_number)
    call ReadWFdata(file_number)
    call ReadIOInput(file_number)
    if(N_inertia .gt. 4) then
      call read_inertia(file_number=file_number)
    else 
      inertia_l = inertia_l_hardcoded
      inertia_m = inertia_m_hardcoded
    endif
    call readmomentdata(file_number)
    call readcranking(file_number)

#if($FAM == 1)
    call readfam(file_number)

    if(xyfile .ne. '' .and. xyfile == xyinfile) then 
      print * ,"ERROR : XYfile and XYinfile carry the same name. Stopping ..."
      stop
    endif
#endif

    if(present(file_number)) then
      close(unit=file_number)
    endif

#if(USE_MPI > 0) 
  ! No MPI ranks can quit this routine before having received all information! 
  call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)
#endif

  end subroutine ReadInput

  subroutine ReadIOInput(file_number)
    !---------------------------------------------------------------------------
    ! Subroutine to read all the data on IO operations
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !   file_number : optional integer. If present, read from (open) channel
    !                 with this number. If absent, read from STDIN.
    !---------------------------------------------------------------------------

    integer(dp), intent(in), optional   :: file_number   
#if(USE_MPI>0)
    integer                             :: mpi_err
#endif

    NameList /IO/ InputFileName,OutputFileName, BXLFIT, COMBI, denfile,potfile,& 
    &           sphffile, spcanfile,checkpointiter, AllowTransform, extraspwfs,&
    &           tofile, blockfile, inertfile,  famfile, xyfile, xyinfile,      &
    &           finfile, N_inertia, potentials_from_file

    ! Only the first MPI RANK reads input
    if(MPI_RANK .eq. 0) then
      if(present(file_number)) then
        read (unit=file_number, nml=IO)
      else
        read (unit=*          , nml=IO)
      endif
    endif

    ! ... and then broadcasts information
    !      ( I am aware that these variables are likely to be useful only to 
    !        rank 0 core, but this might avoid future errors )
#if(USE_MPI > 0)
    call MPI_Bcast(InputFileName , len(InputFileName) , MPI_CHARACTER, 0, &
    &                                                   MPI_COMM_WORLD, mpi_err)
    call MPI_Bcast(OutputFileName, len(OutputFileName), MPI_CHARACTER, 0, &
    &                                                   MPI_COMM_WORLD, mpi_err)
    call MPI_Bcast(BXLFIT        , len(BXLFIT)        , MPI_CHARACTER, 0, &
    &                                                   MPI_COMM_WORLD, mpi_err)
    call MPI_Bcast(COMBI         , len(COMBI)         , MPI_CHARACTER, 0, &
    &                                                   MPI_COMM_WORLD, mpi_err)
    call MPI_Bcast(denfile       , len(denfile)       , MPI_CHARACTER, 0, &
    &                                                   MPI_COMM_WORLD, mpi_err)
    call MPI_Bcast(potfile       , len(potfile)       , MPI_CHARACTER, 0, &
    &                                                   MPI_COMM_WORLD, mpi_err)
    call MPI_Bcast(sphffile      , len(sphffile)      , MPI_CHARACTER, 0, &
    &                                                   MPI_COMM_WORLD, mpi_err)
    call MPI_Bcast(spcanfile     , len(spcanfile)     , MPI_CHARACTER, 0, &
    &                                                   MPI_COMM_WORLD, mpi_err)
    call MPI_Bcast(inertfile     , len(inertfile)     , MPI_CHARACTER, 0, &
    &                                                   MPI_COMM_WORLD, mpi_err)
    call MPI_Bcast(tofile        , len(tofile)        , MPI_CHARACTER, 0, &
    &                                                   MPI_COMM_WORLD, mpi_err)
    call MPI_Bcast(blockfile     , len(blockfile)     , MPI_CHARACTER, 0, &
    &                                                   MPI_COMM_WORLD, mpi_err)

    call MPI_Bcast(checkpointiter, 1                  , MPI_INTEGER, 0, &
    &                                                   MPI_COMM_WORLD, mpi_err)
    call MPI_Bcast(extraspwfs    , 8                  , MPI_INTEGER, 0, &
    &                                                   MPI_COMM_WORLD, mpi_err)

    call MPI_Bcast(allowtransform, 1                  , MPI_LOGICAL, 0, &
    &                                                   MPI_COMM_WORLD, mpi_err)

    call MPI_Bcast(N_inertia     , 1                  , MPI_INTEGER, 0, &
    &                                                   MPI_COMM_WORLD, mpi_err)
    call MPI_Bcast(potentials_from_file, 1            , MPI_INTEGER, 0, &
    &                                                   MPI_COMM_WORLD, mpi_err)
#endif
  
#if(USE_MPI == 0) 
  if(N_inertia .lt. 4) then
    call stp('N_inertia must be at least 4.')
  endif
#else 
  if(N_inertia .ne. 0) then
    call stp('N_inertia has to be 0 for MPI calculations.')
  endif
#endif

  end subroutine ReadIOInput

  subroutine PrintInput(file_number, input_file)
  !-----------------------------------------------------------------------------
  ! This subroutine prints all relevant information of the input, both from the
  ! user and from the wavefunction file.
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Input: 
  !  file_number : integer, only used for printing
  !  input_file  : character, only used for printing 
  !-----------------------------------------------------------------------------
   
    use wavefunctions
    use evolution
    use scfiteration
#if( $FAM == 1 )
    use fam, only : printfam
#endif
    use IO_wf, only : fileblockindices, fileBlockLowest, fileblocknumber, &
    &                 fileblocktype, file_version, readHFBinfofile

    integer*8, intent(in), optional     :: file_number
    character(11), intent(in), optional :: input_file 
#if(USE_MPI>0)
    integer                             :: mpi_err
#endif

    1 format ( 30('-'), 'General Information ', 30('-'))
    2 format ( ' Mesh parameters' )
    3 format ( '   nx = ', i5 , ' ny = ' , i5 , ' nz = ' , i5, ' mv = ' , i5)
    4 format ( '   dx = ', f20.10,' (fm  ) ')
    5 format ( '   dv = ', f5.2,' (fm^3) ')

    6 format ( ' Nucleus')
    7 format ( '    N = ', f10.5  ,'  Z = ', f10.5)
    8 format ( ' Wavefunctions')
    9 format ( '  nwt = ', i5, / &
    &          '  nwn = ', i5, / &
    &          '  nwp = ', i5 )
   98 format ( '  Derivatives stored explicitly: ', a3)
  990 format ( '  Nilsson initialization with hom = (', 3(f7.3) ,')')
  991 format ( '  Initialization with random spwfs')
  992 format ( '  HFBlocks = ', 8i5)
   10 format ( ' IO information', / &
    &          '  inputfilename  =', a32, / &
    &          '  outputfilename =', a32)
    
  101 format ( ' Information obtained from file ')  
  102 format ( '      - version number            : ', i5)
 1021 format ( '      - param. used on file       : ', 20a)
  103 format ( '      - Bogoliubov transfo read?  : ', l5)
 1031 format ( '      - Bogoliubov transfo used?  : ', l5)
  104 format ( '      - Blocking type             : ', i5)
  105 format ( '      - Blocknumber               : ', i5)
  106 format ( '      - Block indices             : ', 10i4)
  107 format ( '      - Block lowest              : ', 10a2)
  108 format ( '      - Passed blocking test      : ', l5)
  109 format ( '      - Potentials read from file : ', l5)

   11 format ( ' Filename for other output (not written if empty): ', /     &
             & '    BXL output     = ', a80, / &
             & '    DEN file       = ', a80, / &
             & '    POT file       = ', a80, / &
             & '    SPHF file      = ', a80, / &
             & '    SPCAN file     = ', a80, / &
             & '    TO file        = ', a80, / & 
             & '    BLOCK file     = ', a80, / &
             & '    INERT file     = ', a80, / &
             & '    FAM file       = ', a80, / &
             & '    XY file        = ', a80) 
  111 format ( '    Input data     = ', a26, / &
               '     on unit ', i10)
 1111 format ( ' Filename for FAM input (not used if empty): ', /     &
             & '    XY init file   = ', a80, / &
             & '    Ext.field file = ', a80)
  112 format ( ' Checkpointiter =', i10)
  113 format (' Printing spwf details during iterations: ', l5)
   12 format ( ' Convergence required', / &
    &          '  Energy convergence           < ', es8.1, / & 
    &          '  Multipole moment convergence < ', es8.1, / &
    &          '  Dispersion convergence       < ', es8.1, / &
    &          '  S.p. gradient convergenc     < ', es8.1, / &
    &          '  Fermi energy convergence     < ', es8.1, / &
    &          '  Angular momentum convergence < ', es8.1)
   13 format ( ' Inverse temperature Beta = ', f14.9)

    if(MPI_rank .eq. 0) then 
      ! Only one MPI rank needs to print information
      print *
      print 1
      print 2

      print 3 , nx, ny, nz, mv
      print 4 , dx
      print 5 , dv
      print 6
      print 7 , neutrons, protons
      print 8
      print 9 , nwt,nwn,nwp
      print 992, HFBlocks_global
      if(store_derivatives) then
        print 98, 'YES'
      else
        print 98, ' NO'
      endif
      if(trim(to_upper(inputfilename)).eq.'INIT' ) then !TODO: make this print too when starting from potentials
        if(adjustl(ini_strategy) .eq. 'NILSSON') then
          print 990, osc_freq
        elseif(adjustl(ini_strategy) .eq. 'RANDOM') then
          print 991
        endif
      endif
      print 13, inversetemp
      print 10, inputfilename, outputfilename
      if(trim(to_upper(inputfilename)).ne.'INIT') then
        print 101
        print 102, file_version
        print 1021, ini_name_param
        print 103, readHFBinfofile
        print 1031, Bogofromfile
        
        print 104, fileblocktype
        select case (fileblocktype)
          case(0)
          case(1,3,5)
            print 105, fileblocknumber
            print 106, fileblockindices
          case(2,4,6)
            print 105, fileblocknumber
            print 107, fileblocklowest
        end select
        print 108, passed_block_test

        print 109, potentials_from_file
      endif 

      print 112, checkpointiter
      print 113, print_adv_spwf_properties

      print 11, BXLFIT, DENFILE, POTFILE, SPHFFILE, SPCANFILE, TOFILE, BLOCKFILE, INERTFILE, FAMFILE, XYFILE
      if(present(file_number)) then
        print 111,  adjustl(trim(input_file)), file_number
      endif
#if( $FAM == 1 )
        print 1111,  XYINFILE, FINFILE
#endif    
      print 12, energy_prec, moment_prec, disp_prec, gradient_prec, fermi_prec,  &
      &         angmom_prec
      
      call printevolution
      call printscfiteration
      call printpairing_init
      call printmoment_init
      call printcranking_init
#if( $FAM == 1 )
      call printfam
#endif
      call printfunctional 
    endif    

  end subroutine PrintInput

  subroutine Readwavefunction()
    !---------------------------------------------------------------------------
    ! High-level routine to determine the starting point of a calculation. 
    !
    ! There are two main starting options, both of which has two suboptions
    !
    ! 1) Initialize in an EV8-style box with Nilsson orbitals
    !    a - start self-consistency cycles immediately
    !    b - read a set of potentials from file to start the calculations
    ! 
    ! 2) Read a set of spwfs from: 
    !    a - *.hdf5 file 
    !    b - *.wf file 
    ! 
    ! Which option is chosen based on the InputFileName keyword: 
    !  - INIT (case insensitive) : option 1a, no reading of any file
    !  - *.pot                   : option 1b, reading of a potential file
    !  - *.hdf5                  : option 2a, reading of a hdf5 wavefunction file
    !  - [any other filename]    : option 2b, reading of an unformatted fortran 
    !                                                      wavefunction file .wf
    !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! These inputs can then be amended by either
    !
    ! a) Breaking a symmetry, as coded by Hephaestos
    ! b) By adding points in the box and/or adding spwfs
    !    This option is currently limited to EV8-style calculations.
    !
    ! Options a) and b) cannot be combined in a single run, but can ofcourse 
    ! be achieved by running the code twice. 
    !
    ! None of a) or b) is allowed if the user does not set the AllowTransform
    ! flag to .true. This is coded like that as a general safeguard.
    !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Important note: MPI load balancing for the single-particle wavefunctions
    ! is rather complicated with respect to applying symmetry transformations 
    ! etc and is for this reason no performed "centrally". Rather, it is done
    ! within each kind of subroutine. 
    !    - iniwavefunctions => load balance based on the EV8 symmetries
    !    - readtantalus     => load balance based on the symmetries on file
    !    - transformspwfs   => load balance based on the actual symmetries
    !                          of the calculation.
    ! This kind of approach incurs some communication overheads that can 
    ! possibly be eliminated with a whole bunch of coding work, but this seems
    ! like an inefficient use of human time since it concerns only the set-up
    ! of a given calculation.
    !---------------------------------------------------------------------------
    use IO_wf, only : read_tantalus_wf, file_HFB_blocks
#if(USE_HDF5 == 1)
    use IO_wf, only : read_tantalus_hdf5
#endif
    use IO_wf, only : file_rank_map, file_spwf_map, file_spwf_inverse 
    use IO_wf, only : filenx, fileny, filenz, filenwn, filenwp,filedx, filemv
    use IO_wf, only : fileblocklowest, fileblocks, fileblocks_global,  &
    &                 fileblocktype, fileblocknumber, fileblockindices,&
    &                 readHFBinfofile, check_blocking_structure

    integer                       :: i, inputoption, lenchar
    character(len=:), allocatable :: standardized_input  
#if(USE_MPI>0)
    integer :: mpi_err
#endif

    call start_timer(T_wfini)
    
    standardized_input = trim(to_upper(inputfilename))
    lenchar=len(standardized_input)

    if(standardized_input.eq.'INIT') then
      ! option 1a : starting completely from scratch
      inputoption = 0 
    else if( standardized_input(lenchar-3:lenchar) .eq. '.POT') then
      ! option 1b : reading potentials
      inputoption = 1 
    else if( standardized_input(lenchar-4:lenchar) .eq. '.HDF5') then
      ! option 2a : reading complete .hdf5 file
      inputoption = 2 
    else
      ! option 2b : reading complete .wf file
      inputoption = 3        
    endif   
    !---------------------------------------------------------------------------
    ! Input options 
    if( inputoption.eq.0  .or. inputoption .eq. 1) then
      ! Generate starting point with Nilsson wavefunctions.
      call iniwavefunctions($ININX, $ININY, $ININZ, $ININWN, $ININWP)

      guessgaps         = .true.
      fileblocks        = HFBlocks
      
      if(inputoption.eq.1) then
        Potentials_read = readpotentials_separate(12, inputfilename)
        Coulomb_read_from_file = .true. 
        ! Signalling that we have direct and Exchange potentials read
      endif

      if( SYM_CODE .ne. "0 1 001 000 10 000 010 111" ) then
        ! Initialisation with nil8 wavefunctions is always EV8-style
        ! Thus we signal that a symmetry transformation is needed
        sym_transfo_needed = .true.
        if( TRANS_CODE .ne. "0 1 001 000 10 000 010 111") then
          print *, "---------------------------------------------------"
          print *, "| Calculations cannot be initialized from scratch |"
          print *, "| for this particular symmetry option.            |"
          print *, "---------------------------------------------------"
          call stp('')
        endif
      endif

      ! ... and then make the code think these things were read from file
      ! by setting all file_X quantities to those initialized
      guessgaps         = .true.

      filenx = $ININX ; fileny = $ININY ; filenz = $ININZ ; filedx = dx
      filemv = filenx*fileny*filenz
      fileblocks_global = HFBlocks_global
      fileblocks        = HFBlocks
      file_spwf_map     = spwf_map
      file_rank_map     = rank_map
      file_spwf_inverse = spwf_inverse
    else if(inputoption.eq.2) then
#if (USE_HDF5 > 0)
      ! Option 2a) start from a previous calculation with hdf5 input file
      call read_tantalus_hdf5(inputfilename, sym_transfo_needed)
#else
      call stp('HDF5 support was not enabled at compilation.')
#endif
    else
      ! Option 2b) start from a previous calculation with .wf input file
      call read_tantalus_wf(12, inputfilename)
      ! No need to guess gaps every time (unless the user asked for it)
    endif
    !---------------------------------------------------------------------------
    ! Transformation options
    if(allowtransform ) then
      if(  sym_transfo_needed ) then 
          ! Option a): break a symmetry and transform the spwfs appropriately
          call Transformspwfs( HFPsi, filenx, fileny, filenz,fileblocks_global,&
          &                    fileblocks, file_rank_map, file_spwf_inverse)
          ! ----> this features a call to load_balance and hence sets the 
          !       correct spwfs mappings everywhere
      else
          ! Option b): add points and/or add spwfs
          call  TransformInput(filenx,fileny,filenz,filenwn,filenwp,filedx,    & 
          &                    fileblocks,file_HFB_blocks, file_spwf_map,      &
          &                    file_rank_map, file_spwf_inverse, extraspwfs)
          ! ----> this features a call to load_balance and hence sets the 
          !       correct spwfs mappings everywhere
      endif
    else  
      ! Sanity check the input
      if(sym_transfo_needed) then
        call stp('Symmetry transformation needed, but not allowed by user.')
      endif
      ! We still need to set the information regarding spwf mapping
      HFblocks     = fileblocks
      spwf_map     = file_spwf_map
      rank_map     = file_rank_map
      spwf_inverse = file_spwf_inverse
    endif
    !---------------------------------------------------------------------------
    ! The following information needs to be transferred in every case
    nwt_local = sum(HFBlocks)
#if(USE_MPI>0)
      HFblocks_global = 0
      call MPI_ALLREDUCE(HFblocks,HFBlocks_global,8, MPI_INTEGER,MPI_SUM,      &
      &                                                 MPI_COMM_WORLD, mpi_err)
#else
      HFBlocks_global = HFBlocks
#endif
    !---------------------------------------------------------------------------
    ! with everything safely in memory, we add in an orthonormalisation to 
    ! guarantee we can start calculating stuff.
#if(USE_MPI > 0)
    ! Copy the 1D wavefunctions to the 2D layout, since that is how we 
    ! orthonormalize ...
    call transfer_1D_to_2D(HFPsi, HFPsi_2D)
#endif
    call orthonormalize
#if(USE_MPI > 0)
    ! ... and make sure the results get back to the original layout
    call transfer_2D_to_1D(HFPsi_2D, HFPsi)
#endif
    !---------------------------------------------------------------------------
    ! Failsafe for the HF transformation
    if(.not.allocated(HFTransfo)) then
        allocate(HFTransfo(nwt_local,nwt_local))
        HFtransfo = 0.0d0
        do i=1,nwt_local
            HFtransfo(i,i) = 1.0d0
        enddo
    endif
    !---------------------------------------------------------------------------
    call set_spwf_symmetries(sx, sy, sz, HFblocks)
    !---------------------------------------------------------------------------
    if(guessgaps) then
      ! Guess some pairing gaps if asked for (always if starting from INIT)
      call initializeGaps(gapvalue)
    endif
    !---------------------------------------------------------------------------
    ! Checking the blocking options: making sure things on the file are in line 
    ! with what the user asked for
    if(Bogofromfile .and. readHFBinfofile .and. pairingscheme.eq.1) then
      if(.not.allocated(fileblocklowest)) then
          ! This is the one case which we will accept: no blocking on the file, 
          ! but blocking in the input. In this case, we need to do an 
          ! explicit diagonalization from the start.
          !Bogofromfile = .false.
      else
          ! In any other case, we check all things we can check.
          passed_block_test =  check_blocking_structure()      
      endif
    endif
    call stop_timer(T_wfini)
  end subroutine ReadWaveFunction

subroutine write_header(iochannel)
    !---------------------------------------------------------------------------
    ! This routine writes a header to file that contains a bunch of information 
    ! on the calculation. It looks like:
    !
    !     #   N = i3, Z = i3, A = i3
    !     #   nwn = i3, nwp = i3
    !     #   Name of the parameterization
    !     #   type of functional
    !     #   Fermi energies of both nucleon species
    !     #   Quadrupole deformation in terms of Q20 and Q22
    !     #   Quadrupole deformation in terms of B20 and B22
    !     #   Quadrupole deformation in terms of B2 and gamma
    !     #   BI 1: Blocktype, Blocknumber
    !     #   BI 2: BlockIndices
    !     #   BI 3: Blocklowest
    !     #   [EMPTY CURRENTLY, RESERVED FOR SYMMETRY INFORMATION]
    !
    !---------------------------------------------------------------------------
    integer, intent(in) :: iochannel
!     type(moment), pointer :: Q20, Q22
   
    1 format("# N = ", i8, ' Z = ', i8, ' A = ', i8)
    2 format("# nwn = ", i9, ", nwp = ", i9)
    3 format("# (nx,ny,nz) = (", 3i5, "), dx = ", f10.8, ' fm')
    4 format("# Parameterisation    : ", a40)
    5 format("# Functional type     : ", a40)
    6 format("# Fermi energies      : ", 2f15.4)
!     7 format("# Quadrupole   Q20,Q22: ", 2f15.4)
!     8 format("# Quadrupole   B20,B22: ", 2f15.4)
!     9 format("# Quadrupole    Q, gam: ", 2f15.4)

!    10 format("# BI 1: ", 2i3)
!    11 format("# BI 2: ", 99i4)
!    12 format("# BI 3: ", 99a3)

    write(iochannel, fmt=1)  int(neutrons), int(protons),int(neutrons+protons)
    write(iochannel, fmt=2)  nwn, nwp
    write(iochannel, fmt=3)  nx, ny, nz, dx
    write(iochannel, fmt=4)  name_param
    write(iochannel, fmt=5)  func_name
    write(iochannel, fmt=6)  FermiEnergy
  
    !Q20 =>FindMoment(2,0,.false.     )
    !Q22 =>FindMoment(2,2,.false., Q20)
    !write(iochannel, fmt=7) sum(Q20%value), sum(Q22%value)
    !write(iochannel, fmt=8)    Q20%beta(4), Q22%beta(4)
    !write(iochannel, fmt=9)    Q(3), G(3)
    
   ! write(iochannel, fmt=10)  blocktype, blocknumber
   ! if(blocknumber .gt. 0) then
   !   write(iochannel, fmt=11) Blockindices
   !   write(iochannel, fmt=12) Blocklowest
   ! else
   !   write(iochannel, fmt=11)
   !   write(iochannel, fmt=12)
   ! endif

    write(iochannel, fmt='(a1)') '#'

  end subroutine write_header

#if( $FAM == 0 )
  subroutine WriteWaveFunction(chan, ofn)
    !--------------------------------------------------------------------------------------------
    ! Simple routine to call the appropriate subroutine depending on type of output
    ! asked for by the user.
    !
    ! Input:
    !------------
    !   chan : integer, channel to open file ofn on
    !   ofn  : character, filename to write to.
    !          If it does not exist, will get created.
    !          If this ends in "HDF5" (case-insensitive), then the code will write an HDF5 file.
    !          If not, then a simple fortran unformatted file will be written.
    !--------------------------------------------------------------------------------------------
    use IO_wf, only: write_tantalus_wf
#if(USE_HDF5 == 1)
    use IO_wf, only : write_tantalus_hdf5
#endif

    integer, intent(in)          :: chan
    character(len=*), intent(in) :: ofn

    if(trim(to_upper(ofn(len_trim(ofn)-3:))).eq.'HDF5') then
#if(USE_HDF5>0)
      call write_tantalus_hdf5(ofn) !new hdf5 format
#else
      call stp('HDF5 support was not enabled at compilation.')
#endif
    else
      call write_tantalus_wf(chan, ofn) ! old style in .wf file
    endif

  end subroutine WriteWaveFunction
#endif 

#if ($FAM == 0)
  subroutine write_advanced_output(iter, iomsg)
    !--------------------------------------------------------------------------- 
    ! Collection routine for writing advanced output. 
    ! Current options:
    !     BXLFIT   : one-line file incorporating essential info for a fit by
    !                the Brussels group
    !     POTFILE  : file containing 3D information on various potentials
    !     DENFILE  : file containing 3D information on the densities
    !     SPHFFILE : information on the single-particle wavefunctions in the 
    !                Hartree-Fock basis.
    !     SPCANFILE: information on the single-particle wavefunctions in the 
    !                Canonical basis.
    !---------------------------------------------------------------------------
    integer          :: iter
    character(len=*) :: iomsg

    ! Bonus file for quick feedback into the fit
    if(BXLFIT .ne. '') then
        call Brussels_output(iter, iomsg)
    endif  
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Write the neutron, proton and charge density to a file for postprocessing 
    if(DENFILE .ne. '') then
      call write_densities(Density, DENFILE)
    endif
    if(TOFILE .ne. '') then
$TR   call stp('Time-odd densities do not figure in a calculation that assumes time-reversal.')
$NTR  call write_timeodd_densities(Density, TOFILE)
    endif
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Write the relevant potentials to a file for postprocessing
    if(POTFILE .ne. '') then
      call write_potentialfile(potentials, POTFILE)
    endif
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Single-particle wave function information 
    ! a) in the HF-basis
    if(SPHFFILE .ne. '') then
      call write_sp_info(SPHFFILE)
    endif 
    ! b) in the canonical basis
    if(pairingtype.eq.2 .and. SPCANFILE .ne. '') then
      call write_sp_info_can(SPCANFILE)
    endif
    ! c) for the combinatorial level density code
    if(COMBI .ne. '') then
      call combi_output(COMBI)
    endif
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! 
    if(pairingtype.eq.2 .and. BLOCKFILE .ne. '') then
      call write_blocked_sps(BLOCKFILE)
    endif        
    
    if(inertfile .ne. '') then
      call write_inertias(inertfile)
    endif

  end subroutine write_advanced_output

  subroutine Brussels_output(iter, iomsg)
    !---------------------------------------------------------------------------
    ! Write an extra file for use in the Brussels fitting protocol to
    !     zXXXnYYY.out"
    ! where XXX = proton number
    !       YYY = neutron number
    !
    ! It contains on a single line
    !
    !      N, Z, Total energy, Enocor, Erot+Evib,  
    ! &    b20, b22, b30, b32, b40                                 &
    ! &    sqrt(<r^2_p>/Z), B(1:3), J2(1:3),                       &
    ! &    avgap_uv(n), avgap_uv(p), Epairn, Epairp                &
    ! &    DEhistory, iter, iomsg
    !
    ! Notes:
    ! *  sqrt(<r^2_p>/Z) is calculated as in the moments module, i.e. it is 
    !    calculated from the charge density, which is not necessarily the proton
    !    density.
    ! * B is the Belyaev moment of inertia, along every axis
    ! * J2 is the expectation value of J^2, along every axis
    !
    ! * Erot is the sum of the rotational corrections
    ! * Evib is the sum of the vibrational corrections
    ! 
    ! * avgap_uv(p/n)    is the average gap for protons or neutrons as
    !                     calculated for the selected pairing approximation
    !                     by averaging the gaps with the anomalous density 
    !                     kappa (uv). This is of course zero on the HF level.
    !
    ! * DEhistory is the maximal deviation in the last 5 iterations
    !
    ! * io is a character that indicates if problems have been detected.
    !   Currently:
    !      * 'CONVERGED'     =>  The calculation exited when it was judged 
    !                            complete.
    !      * 'MAXITER'       => The calculation reached the maximum specified 
    !                           number of iterations without thinking it was 
    !                           complete. Does not necessarily mean something 
    !                           is wrong.
    !---------------------------------------------------------------------------

    use Moments    
    use functional
    character(len=*), intent(in) :: iomsg

    type(Moment), pointer :: Q20, Q22, r2, Q30, Q32, Q40
    real(KIND=dp)         :: E, Enocor,Erot_vib, rms,B(3), DEhistory
    real(KIND=dp)         :: b20, b22, b30, b32, b40
    integer, intent(in)   :: iter
    integer               :: iterh

    character(len=len(BXLFIT)+12) :: filedone
    
    Q20 =>FindMoment( 2,0,.false.)
    Q22 =>FindMoment( 2,2,.false., Q20)
    r2  =>FindMoment(-2,0,.false., Q20) ! The rms radius is associated with l=-2

    Q30 => FindMoment(3,0,.false., Q20)
    Q32 => FindMoment(3,2,.false., Q20) ! We start searching from Q20, as that is guaranteed to exist
    Q40 => FindMoment(4,0,.false., Q20)

    write(filedone,'(a,"z",i3.3,"n",i3.3,".out")')        &  
     &     trim(adjustl(BXLFIT)),int(protons),int(neutrons) 
  
    open(unit=10,file=filedone)

    E = TotalE 
    Enocor = totalE - sum(rotcorrection)      &
    &                 - sum(COMcorrection(2,:)) &
    &                 - vibcorrection
    Erot_vib = sum(rotcorrection)+ vibcorrection

    b20 = Q20%beta(4)
    b22 = Q22%beta(4)    
    rms     =     r2%chargevalue

    if(associated(Q30)) then
      b30 = Q30%beta(4)
    else
      b30 = 0.0d0
    endif

    if(associated(Q32)) then
      b32 = Q32%beta(4)
    else
      b32 = 0.0d0
    endif
    
    b40 = Q40%beta(4)

    select case(pairingtype)
    case(0,1)
      B = Belyaev(:,3)
    case(2)
      if(inversetemp.lt.0) then
        B = Bely_coll(:,3)
      else
        B = Belyaev(:,3)
      endif
    end select
  
    iterh=min(iter,5)
    DEhistory=maxval(Ehistory(1:iterh))-minval(Ehistory(1:iterh))
    
    
    write(10,'(2i4,20(1x,f20.10), i6)', advance='NO')               &
    &     int(protons), int(neutrons), E, Enocor, Erot_vib,         &
    &     b20, b22, b30, b32, b40,                                  &
    &     sqrt(rms/protons),  B(1:3), J2_coll(1:3,3),               &
    &     average_gap(2,1), average_gap(2,2),  PairingEnergy(1:2),  &
    &     DEhistory, iter
  
    write(10, '(2x, a99)') adjustl(iomsg)
    close(10)

  end subroutine Brussels_output
  
  subroutine write_densities(R, fname)
    !---------------------------------------------------------------------------
    ! Write the following densities to a file named "fname"
    !    rho(neutron), rho(proton), rho(charge), 
    !    tau(neutron), tau(proton),
    !    \tilde{rho}(neutron), \tilde{rho}(proton)
    !---------------------------------------------------------------------------
    ! The file contains a header written by the subroutine write_header,
    ! supplemented by a dedicated line explaining the content of each column.
    ! The format of the body of said file is
    ! 
    !   x , y , z,  rho_n, rho_p, rho_c, tau_n, tau_p, DP_I_I_n, DP_I_I_p
    !
    ! where the first three numbers are the Cartesian coordinates in fm, withµ
    ! the densities all in their natural units. The mesh points are traverse in 
    ! column-major order ('Fortran order'), which might not be how your favorite 
    ! plotting tool prefers it. Note that the densities are written "as-is" to 
    ! file, i.e. only in part of the box that is actually represented 
    ! numerically. It is up to postprocessing to construct the densities in the 
    ! simulation volume.
    !---------------------------------------------------------------------------
    type(DensityVector), intent(in), target :: R
    character(len=*), intent(in)            :: fname
    integer                                 :: io, i,j,k, mi

    1 format('#', 6x, 'X[fm]',20x,'Y[fm]', 20x,'Z[fm]', 20x,   &
      &               'rho_n', 20x, 'rho_p', 20x,'rho_c', 20x, &
      &               'tau_n', 20x, 'tau_p', 20x,              &
      &               'tilde{rho}_n', 13x, 'tilde{rho}_p')

    open(1,file=fname, iostat=io)
    if(io.ne.0) then    
      print *, 'Something went wrong with writing a density to file.'
      print *, 'filename = ', fname
      call stp('')
    endif
    
    call write_header(1)
    write(1, fmt=1) 
    do k=1,nz
      do j=1,ny
        do i=1,nx
          write(1, fmt='(3es25.12)', advance='no') &
          &          meshx(i), meshy(j), meshz(k)
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          ! the contributions above are indexed according to (x,y,z) but 
          ! we do not have this luxury for most of the densities
          mi = meshindex(i,j,k)
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          ! The ordinary and charge density; always defined
          write(1, fmt='(2es25.12)', advance='no') R%D_I_I(mi,1),R%D_I_I(mi,2)
          write(1, fmt='( es25.12)', advance='no') R%chargedensity(i,j,k)
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          ! The kinetic density; its definition depends on the type of EDF used
!$TAUSCALAR write(1, fmt='(2es25.12)', advance='no') &
!$TAUSCALAR &         R%D_Nm_Nm(mi,1), R%D_Nm_Nm(mi,2)
!$TAUTENSOR write(1, fmt='(2es25.12)', advance='no') &
!$TAUTENSOR &         R%D_N_N(mi,1,1,1) + R%D_N_N(mi,2,2,1) + R%D_N_N(mi,3,3,1),&
!$TAUTENSOR &         R%D_N_N(mi,1,1,2) + R%D_N_N(mi,2,2,2) + R%D_N_N(mi,3,3,2)
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          ! The pairing fields FP_I_I
          write(1, fmt='(2es25.12)',advance='no') R%DP_I_I(mi,1), R%DP_I_I(mi,2)
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          ! We are done writing this line in the output
          write(1, fmt='()') !  newline character
        enddo
      enddo
    enddo

    close(1)
  end subroutine write_densities
  
 subroutine write_potentialfile(F,fname)
    !---------------------------------------------------------------------------
    ! Write the mean-field potentials to a file named "fname":
    !  central       coulomb   coulomb   kinetic  pairing   spin-orbit  
    !                 direct   exchange
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Remarks:
    ! 
    !   *) Note that G_I_NS has 9 components for each isospin, corresponding 
    !      to the gradient and spin indices. These are arranged in lexographical
    !      order.  
    !   *) F_I_I should be the exact potential corresponding to the derivative 
    !      of the Skyrme energy with respect to D_I_I. This means that it should
    !      not include
    !       - the contributions from constraints
    !       - the contribution of the Coulomb interaction
    !   *) To reproduce the complete state of the code, it is important that 
    !      the Coulomb potentials written to file are the potentials 
    !      corresponding to the charge density; for the direct potential this
    !      is the U that satisfies 
    ! 
    !                 Delta U = 4 pi rho_charge
    !
    !      If the finite extent of the nucleons charge density is taken into 
    !      account selfconsistently, this means that U is NOT the potential 
    !      which should be added to F_I_I. 
    !
    !---------------------------------------------------------------------------
    !
    ! The file contains a header written by the subroutine write_header,
    ! supplemented by
    ! #   X[fm] Y[fm] Z[fm] V_nuc(n) V_nuc(p) V_c(n) V_c(p) 
    !               V_so(xx,n), V_so(xy,n), ..., V_so(zz,p)
    ! 
    ! where the # are included so that Numpy (or other plotting tools) can 
    ! ignore these lines when naively plotting stuff. Note that the fourth
    ! line is currently empty, but is reserved for future additions concerning
    ! symmetry options of the current run.
    !
    ! where the first three numbers are the Cartesian coordinates (units of fm).
    ! The points are written down in column-major order ('Fortran order'), 
    ! which might not be how your favorite plotting tool prefers it.
    !---------------------------------------------------------------------------
    use Coulombmod ! module explicitly 'used' in order to be able to place the 
                   ! values of the direct and exchange Coulomb potentials 
                   ! correctly on the mesh


    type(PotentialVector), intent(in) :: F
    character(len=*), intent(in)      :: fname
    real(KIND=dp), pointer            :: Vnucp(:,:,:), Vnucn(:,:,:)
    real(KIND=dp), allocatable        :: Coulp(:,:,:), Excp(:,:,:)

    real(KIND=dp), allocatable, target   :: temp(:,:)
    integer                              :: io, i,j,k, mu, nu, ox, oy, oz, mi
    character(len=1) :: directions(3) 

    1 format('#', 6x, 'X[fm]',20x,'Y[fm]', 20x,'Z[fm]', 20x, 'V_nuc(n)', 17x, 'V_nuc(p)', 17x, &
      &      'V_cd', 21x, 'V_ce', 21x, 'V_kin(n)', 17x, 'V_kin(p)', 17x, 'FP_n', 21x, 'FP_p',21x) 
    2 format('W_', 2a1,'(n)', 18x, 'W_', 2a1,'(p)', 18x )

    open(1,file=fname, iostat=io)
    if(io.ne.0) then    
      print *, 'Something went wrong with writing a potential to file.'
      print *, 'filename = ', fname
      call stp('')
    endif

    call write_header(1)
    write(1, fmt=1, advance='no') 
  
    directions = (/'x', 'y', 'z'/)
    do mu=1,3
      do nu=1,3
        write(1, fmt=2, advance='no') directions(mu), directions(nu), &
        &                             directions(mu), directions(nu)
      enddo
    enddo
    write(1, fmt='()')
    
    ! The central nuclear potential is the potential associated with D_I_I, but 
    ! it should not include the contribution of the constraints, nor the 
    ! contribution of the direct and exchange Coulomb potentials
    allocate(temp(nx*ny*nz,2), coulp(nx,ny,nz), excp(nx,ny,nz))
    
    temp = F%F_I_I(:,1:2) !- constraint_I_I

    Vnucn(1:nx,1:ny,1:nz)  => temp(:,1)
    Vnucp(1:nx,1:ny,1:nz)  => temp(:,2)

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Subtracting the coulomb potential depends on our treatment of the 
    ! proton and neutron finite size effect
    ox = coul_offset_x ; oy = coul_offset_y ; oz = coul_offset_z
    if((all(protonsize.eq.0.0) .and. all(neutronsize.eq.0.0)) .or.         &
    &                             (.not. nucleonsize_selfconsistent)) then
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      ! No finite size effect
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      ! The index juggling is ugly, but necessary, because the Coulomb 
      ! potential has a different size than the Lagrange mesh.
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      do k=1,nz
        do j=1,ny
          do i=1,nx
            Vnucp(i,j,k)          =   Vnucp(i,j,k) &
            &                     - F%CoulombPotential(i+ox,j+oy,k+oz)    &
            &                     - F%ExchangePotential(i,j,k)
          enddo
        enddo
      enddo
    else
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! Finite size effects taken into account
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! Note that there is no index juggling since these matrices are 
      ! conveniently defined on the ordinary mesh.
      Vnucn = Vnucn - F%FoldedCoul(:,:,:,1) &
      &             - F%FoldedExchange(:,:,:,1)
      Vnucp = Vnucp - F%FoldedCoul(:,:,:,2) &
      &             - F%FoldedExchange(:,:,:,2)
    endif
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! The potentials related to the charge density
    Coulp = F%CoulombPotential (ox+1:ox+nx,oy+1:oy+ny,oz+1:oz+nz)
    Excp  = F%ExchangePotential(ox+1:ox+nx,oy+1:oy+ny,oz+1:oz+nz)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          ! Mesh coordinates and F_I_I and coulomb contribution to it.
          write(1, fmt='(7es25.12)', advance='no') &
          &          meshx(i), meshy(j), meshz(k),        &
          &            Vnucn(i,j,k), Vnucp(i,j,k),        & 
          &            Coulp(i,j,k), Excp(i,j,k) 
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          ! the contributions above are indexed according to (x,y,z) but 
          ! we do not have this luxury for the following potentials
          mi = meshindex(i,j,k)
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          ! The kinetic potential is the potential F_Nm_Nm if D_Nm_Nm is used 
          ! in the calculation. If instead the kinetic density is decontracted, 
          ! i.e. D_N_N is used, then we write the scalar component of the tensor
!$TAUSCALAR write(1, fmt='(2es25.12)', advance='no') &
!$TAUSCALAR &         F%F_Nm_Nm(mi,1), F%F_Nm_Nm(mi,2)
!$TAUTENSOR write(1, fmt='(2es25.12)', advance='no') &
!$TAUTENSOR &         F%F_N_N(mi,1,1,1) + F%F_N_N(mi,2,2,1) + F%F_N_N(mi,3,3,1),&
!$TAUTENSOR &         F%F_N_N(mi,1,1,2) + F%F_N_N(mi,2,2,2) + F%F_N_N(mi,3,3,2)
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          ! The pairing fields FP_I_I
          write(1, fmt='(2es25.12)',advance='no') F%FP_I_I(mi,1), F%FP_I_I(mi,2)
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          ! The spin-orbit potential is the potential G_I_NS
          !do mu=1,3
          !  do nu=1,3
          !    write(1, fmt='(2es25.12)', advance='no') &
          !    &            F%G_I_NS(mi,mu,nu,1), F%G_I_NS(mi,mu,nu,2)
          !  enddo
          !enddo
          write(1, fmt='()') !  newline character
        enddo
      enddo
    enddo

    close(1)
  end subroutine write_potentialfile

    !!NS_t0t3: fname is deleted from write nabla argument
  subroutine write_nablaJ()
    !---------------------------------------------------------------------------
    ! Debugging routine that can be used to write both
    !    (1)  the vector component of Jmunu 
    !    (2)  the distinct components of the divergence of this vector
    ! to file.
    !---------------------------------------------------------------------------
    !NS_t0t3:
    ! real(KIND=dp), pointer           :: dxn(:,:,:), dxp(:,:,:), lapdn(:,:,:)
    ! real(KIND=dp), pointer           :: dyn(:,:,:), dyp(:,:,:), lapdp(:,:,:)
    ! real(KIND=dp), pointer           :: dzn(:,:,:), dzp(:,:,:)
    ! character(len=*), intent(in)     :: fname
    ! integer                          :: io, i,j,k, it
    ! real(KIND=dp), target             :: divJ(nx*ny*nz,3,2), lapd(nx*ny*nz,2)
    ! real(KIND=dp), target             :: Jmn(nx*ny*nz,3,2)

    ! 1 format('#  X[fm]   Y[fm]   Z[fm]       ')
    ! open(1,file=fname, iostat=io)
    ! if(io.ne.0) then    
    !   print *, 'Something went wrong with writing a density to file.'
    !   print *, 'filename = ', fname
    !   call stp('')
    ! endif
    
    ! do it=1,2
    !   Jmn(:,1,it) = C_I_NS(:,2,3,it) - C_I_NS(:,3,2,it) 
    !   Jmn(:,2,it) = C_I_NS(:,3,1,it) - C_I_NS(:,1,3,it) 
    !   Jmn(:,3,it) = C_I_NS(:,1,2,it) - C_I_NS(:,2,1,it) 
    ! enddo

    ! do it=1,2
    !   call Derive_X(Jmn(:,1,it), -1, divJ(:,1,it)) 
    !   call Derive_Y(Jmn(:,2,it), -1, divJ(:,2,it)) 
    !   call Derive_Z(Jmn(:,3,it), -1, divJ(:,3,it)) 
    ! enddo

    ! dxn(1:nx,1:ny,1:nz)  => Jmn(:,1,1)
    ! dxp(1:nx,1:ny,1:nz)  => divJ(:,1,1)
    ! dyn(1:nx,1:ny,1:nz)  => Jmn(:,2,1)
    ! dyp(1:nx,1:ny,1:nz)  => divJ(:,2,1)
    ! dzn(1:nx,1:ny,1:nz)  => Jmn(:,3,1)
    ! dzp(1:nx,1:ny,1:nz)  => divJ(:,3,1)
    
    ! lapd  = LAP_D_I_I
    ! lapdn(1:nx,1:ny,1:nz) => LAP_D_I_I(:,1)
    ! lapdp(1:nx,1:ny,1:nz) => LAP_D_I_I(:,2)
    
    ! call write_header(1)
    ! write(1, fmt=1) 
    ! do k=1,nz
    !   do j=1,ny
    !     do i=1,nx
    !       write(1, fmt='(3f8.3, 8es25.12E3)') meshx(i), meshx(j), meshz(k),      &
    !       &                           dxn(i,j,k), dxp(i,j,k), &
    !       &                           dyn(i,j,k), dyp(i,j,k), &
    !       &                           dzn(i,j,k), dzp(i,j,k), &
    !       &                           lapdn(i,j,k), lapdp(i,j,k)
    !     enddo
    !   enddo
    ! enddo

    ! close(1)
  end subroutine write_nablaJ

  subroutine write_timeodd_densities(R, fname)
    !---------------------------------------------------------------------------
    ! Write the following densities to a file named "fname"
    !    D_I_S, C_I_N
    ! All are total densities, i.e. neutron + proton.
    !---------------------------------------------------------------------------
    ! The file contains a header written by the subroutine write_header,
    ! supplemented by
    !
    !     #   X[fm] Y[fm] Z[fm]
    ! 
    ! where the # are included so that Numpy (or other plotting tools) can 
    ! ignore these lines when naively plotting stuff. Note that the fourth
    ! line is currently empty, but is reserved for future additions concerning
    ! symmetry options of the current run.
    !
    ! The format of the body of said file is
    ! 
    !   x , y , z, D_I_Sx/y/z (n),D_I_Sx/y/z (p), C_I_Nx/y/z (n),C_I_Nx/y/z (p) 
    !
    ! where the first three numbers are the Cartesian coordinates (units of fm).
    ! The points are written down in column-major order ('Fortran order'), 
    ! which might not be how your favorite plotting tool prefers it.
    !---------------------------------------------------------------------------
    ! Note that the densities are written "as they are" to file, i.e. only in
    ! part of the box that is actually represented numerically. It is up to
    ! postprocessing to actually construct the densities in the entire box.
    !---------------------------------------------------------------------------

!$NTR    real(KIND=dp), pointer           :: Sxn(:,:,:), Sxp(:,:,:)
!$NTR    real(KIND=dp), pointer           :: Jxn(:,:,:), Jxp(:,:,:)
!$NTR    real(KIND=dp), pointer           :: Syn(:,:,:), Syp(:,:,:)
!$NTR    real(KIND=dp), pointer           :: Jyn(:,:,:), Jyp(:,:,:)
!$NTR    real(KIND=dp), pointer           :: Szn(:,:,:), Szp(:,:,:)
!$NTR    real(KIND=dp), pointer           :: Jzn(:,:,:), Jzp(:,:,:)
!$NTR    real(KIND=dp), pointer           :: Txn(:,:,:), Txp(:,:,:)
!$NTR    real(KIND=dp), pointer           :: Tyn(:,:,:), Typ(:,:,:)
!$NTR    real(KIND=dp), pointer           :: Tzn(:,:,:), Tzp(:,:,:)

    type(DensityVector), intent(in), target :: R
    character(len=*), intent(in)            :: fname!
!
!    real(KIND=dp), allocatable, target      :: totalangmom(:,:,:)
!    integer                                 :: io, i,j,k
!$NTR integer                                :: it
!$TR  real(KIND=dp)                          :: trash!

!    1 format('#  X[fm]   Y[fm]   Z[fm] ', &
!    &        '   Sxn     Syn     Szn   ', &
!    &        '   Sxp     Syp     Szp   ', &
!    &        '   jxn     jyn     jzn   ', &
!    &        '   jxp     jyp     jzp   ', &
!    &        '   Jxn     Jyn     Jzn   ', &
!    &        '   Jxp     Jyp     Jzp   ')
!    2 format('#  0       1       2     ', &
!    &        '   3       4       5     ', &
!    &        '   6       7       8     ', &
!    &        '   9      10      11     ', &
!    &        '  12      13      14     ', &
!    &        '  15      16      17     ', &
!    &        '  18      19      20     ')

!$TR trash = R%D_I_I(1,1) ! to stop compiler complaints when TR is conserved

!    open(1,file=fname, iostat=io)
!    if(io.ne.0) then    
!      print *, 'Something went wrong with writing a density to file.'
!      print *, 'filename = ', fname
!      call stp('')
!    endif!

!$NTR    Sxn(1:nx,1:ny,1:nz)  => R%D_I_S(:,1,1) ; Sxp(1:nx,1:ny,1:nz)  => R%D_I_S(:,1,2)
!$NTR    Syn(1:nx,1:ny,1:nz)  => R%D_I_S(:,2,1) ; Syp(1:nx,1:ny,1:nz)  => R%D_I_S(:,2,2)
!$NTR    Szn(1:nx,1:ny,1:nz)  => R%D_I_S(:,3,1) ; Szp(1:nx,1:ny,1:nz)  => R%D_I_S(:,3,2)

!$NTR    Jxn(1:nx,1:ny,1:nz)  => R%C_I_N(:,1,1) ; Jxp(1:nx,1:ny,1:nz)  => R%C_I_N(:,1,2)
!$NTR    Jyn(1:nx,1:ny,1:nz)  => R%C_I_N(:,2,1) ; Jyp(1:nx,1:ny,1:nz)  => R%C_I_N(:,2,2)
!$NTR    Jzn(1:nx,1:ny,1:nz)  => R%C_I_N(:,3,1) ; Jzp(1:nx,1:ny,1:nz)  => R%C_I_N(:,3,2)


    ! Calculate the total angular momentum density
!$NTR    allocate(totalangmom(nx*ny*nz,3,2)) ; totalangmom = 0.0d0
!$NTR    do it=1,2
!$NTR      totalangmom(:,1,it) = 0.5 * R%D_I_S(:,1,it) ! spin part
!$NTR      totalangmom(:,2,it) = 0.5 * R%D_I_S(:,2,it) ! spin part
!$NTR      totalangmom(:,3,it) = 0.5 * R%D_I_S(:,3,it) ! spin part
!$NTR      do i=1, nx*ny*nz
!$NTR        ! X component : J_x ~ y j_z - z j_y
!$NTR        TotalAngMom(i,1,it) = TotalAngMom(i,1,it) &
!$NTR        & + meshgrid(i,2) * R%C_I_N(i,3,it) - meshgrid(i,3) * R%C_I_N(i,2,it)
!$NTR
!$NTR        ! Y component : J_y ~ z j_x - x j_z
!$NTR        TotalAngMom(i,2,it) = TotalAngMom(i,2,it) &
!$NTR        & + meshgrid(i,3) * R%C_I_N(i,1,it) - meshgrid(i,1) * R%C_I_N(i,3,it)
!$NTR
!$NTR        ! Z component : J_z ~ x j_y - y j_x
!$NTR        TotalAngMom(i,3,it) = TotalAngMom(i,3,it) &
!$NTR        & + meshgrid(i,1) * R%C_I_N(i,2,it) - meshgrid(i,2) * R%C_I_N(i,1,it)
!$NTR      enddo
!$NTR    enddo

!   
!$NTR    Txn(1:nx,1:ny,1:nz)  => TotalAngMom(:,1,1) 
!$NTR    Txp(1:nx,1:ny,1:nz)  => TotalAngMom(:,1,2)
!$NTR    Tyn(1:nx,1:ny,1:nz)  => TotalAngMom(:,2,1) 
!$NTR    Typ(1:nx,1:ny,1:nz)  => TotalAngMom(:,2,2)
!$NTR    Tzn(1:nx,1:ny,1:nz)  => TotalAngMom(:,3,1) 
!$NTR    Tzp(1:nx,1:ny,1:nz)  => TotalAngMom(:,3,2)!!

!    call write_header(1)
!    write(1, fmt=2) 
!    write(1, fmt=1) 
!    do k=1,nz
!      do j=1,ny
!        do i=1,nx
!$NTR          write(1, fmt='(3f8.3, 18es25.12E3)') meshx(i), meshx(j), meshz(k),   &
!$NTR          &                                Sxn(i,j,k), Syn(i,j,k), Szn(i,j,k), & 
!$NTR          &                                Sxp(i,j,k), Syp(i,j,k), Szp(i,j,k), & 
!$NTR          &                                Jxn(i,j,k), Jyn(i,j,k), Jzn(i,j,k), & 
!$NTR          &                                Jxp(i,j,k), Jyp(i,j,k), Jzp(i,j,k), &
!$NTR          &                                Txn(i,j,k), Tyn(i,j,k), Tzn(i,j,k), &
!$NTR          &                                Txp(i,j,k), Typ(i,j,k), Tzp(i,j,k)

!$TR          write(1, fmt='(3f8.3, 18es25.12E3)') meshx(i), meshx(j), meshz(k),   &
!$TR          &                                0.0d0,0.0d0,0.0d0, &
!$TR          &                                0.0d0,0.0d0,0.0d0, &
!$TR          &                                0.0d0,0.0d0,0.0d0, &
!$TR          &                                0.0d0,0.0d0,0.0d0, &
!$TR          &                                0.0d0,0.0d0,0.0d0, &
!$TR          &                                0.0d0,0.0d0,0.0d0
!        enddo
!      enddo
!    enddo
!
!    deallocate(TotalAngMom)
!    close(1)

  end subroutine write_timeodd_densities

  subroutine write_sp_info(fname)
    !---------------------------------------------------------------------------
    ! Write detailed information on the single-particle spectrum to file.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Note that this is all information in the Hartree-Fock basis. Note that 
    ! time-reversal symmetry is hardcoded for the moment, hence every line
    ! represents a single-particle wavefunction and its partner.
    !---------------------------------------------------------------------------
    ! The file contains a header written by the subroutine write_header, 
    ! supplemented by
    !     #   Information in the HartreeFock basis
    !     #   i  iso  P  occ E  JxT JyT Jz J SxT SyT Sz
    !
    ! In the body of the file, it contains the following information  
    !      wave, isospin, parity, rho, spenergy, JX, JY, JZ, JJ, SxT, SyT, Sz
    !
    !   wave     : numbering 
    !   isospin  : -1 for neutrons, +1 for protons
    !   parity   : parity quantum number
    !   rho      : diagonal element of rho(wave,wave), with an extra factor
    !              of two for time-reversal
    !   spenergy : expectation value of the sp hamiltonian for this state
    !              (at convergence, these are eigenstates of h)
    !   JX       : matrix element of Jx T for this state.
    !              Note extra time-reversal reversal!
    !   JY       : matrix element of Jy T for this state.
    !              Note extra time-reversal reversal!
    !   JZ       : matrix element of Jz T for this state.
    !   JJ       : J quantum number (real number) that corresponds to this state
    !              such that 
    !                    (JJ+1) JJ = <Jx^2> + <Jy^2> + <Jz^2> 
    !
    !   SxT      : matrix element of S_x T (real part)
    !   SyT      : matrix element of S_y T (imaginary part)
    !   Sz       : matrix element of S_z 
    !---------------------------------------------------------------------------
    use wavefunctions

    character(len=*), intent(in) :: fname
    integer                      :: io, i,  wave
    integer                      :: ProtonOrder(nwp), NeutronOrder(nwn)
    real(KIND=dp)                :: Jx, Jy, Jz, JJ, Spinx, Spiny, Spinz, P
 
    1 format(2i5, 10f10.4)
    2 format("# Neutron spwfs")
    3 format("# Proton spwfs")
    4 format("# Information in the Hartree-Fock basis")

    60 format ("#",3x,'i',3x,'iso',3x,'P',4x,'occ',7x,'<h>',7x,  &
     &        'JxT',7x,'JyT', 7x ,'Jz', 8x, 'J', 9x, 'SxT', 7x,'SyT',7x,'Sz')    

    open(1,file=fname, iostat=io)
    if(io.ne.0) then    
      print *, 'Something went wrong with the sp. info to file.'
      print *, 'filename = ', fname
      call stp('')
    endif
    call write_header(1)
    write(1, fmt=4)
    write(1, fmt='(a1)') '#'
    write(1, fmt=60) 
    write(1, fmt=2) 

    neutronorder = OrderSpwfsISO(-1)
    protonorder  = OrderSpwfsISO(+1)
    !---------------------------------------------------------------------------
    ! First do the neutron wavefunctions    
    do i=1,nwn
      wave = neutronorder(i)
      P = P_hf(wave)

      Jx = HF_JTR(1,wave) ; Spinx = HF_STR (1,wave)
      Jy = HF_JTI(2,wave) ; Spiny = HF_STI (2,wave)
      Jz = HF_J  (3,wave) ; Spinz = HF_spin(3,wave)
      JJ = HF_JJ(wave)

      write(1, fmt=1) wave, -1, p, rho_HF(wave), spenergies(wave),   & 
      &               Jx, Jy,Jz, JJ, Spinx, Spiny, Spinz
    enddo      
    write(1, fmt=3) 
    !---------------------------------------------------------------------------
    ! Then do the proton wavefunctions    
    do i=1,nwp
      wave =  protonorder(i)
      P = P_hf(wave)

      Jx = HF_JTR(1,wave) ; Spinx = HF_STR (1,wave)
      Jy = HF_JTI(2,wave) ; Spiny = HF_STI (2,wave)
      Jz = HF_J  (3,wave) ; Spinz = HF_spin(3,wave)
      JJ = HF_JJ(wave)

      write(1, fmt=1) wave, +1, p, rho_HF(wave), spenergies(wave), & 
      &                Jx, Jy,Jz,JJ,Spinx,Spiny,Spinz
    enddo
    close(1)
  end subroutine write_sp_info

  subroutine write_sp_info_can(fname)
    !---------------------------------------------------------------------------
    ! Write detailed information on the single-particle spectrum to file.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Note that this is all information in the basis where rho is diagonal. 
    ! (This may or may not be the canonical basis, depending on the options 
    !  of the calculation.) 
    ! Note that time-reversal symmetry is hardcoded for the moment, hence every 
    ! line represents a single-particle wavefunction and its partner.
    !---------------------------------------------------------------------------
    ! The file contains a header written by the subroutine write_header, 
    ! supplemented by
    !     #   Information in the basis that diagonalizes RHO
    !     #   i  iso  P  occ E  JxT JyT Jz J
    ! 
    ! In the body of the file, it contains the following information  
    !      wave, isospin,parity,rho_can, canenergy, JX, JY, JZ, JJ, SxT, SyT, Sz
    !
    !   wave     : numbering 
    !   isospin  : -1 for neutrons, +1 for protons
    !   parity   : parity quantum number
    !   rho_can  : diagonal element of rho, with an extra factor of two for
    !              time-reversal
    !   canenergy: expectation value of the sp hamiltonian for this state
    !              NOTE: these states are not eigenstates of the sp hamiltonian!
    !   JX       : matrix element of Jx T for this state.
    !              Note extra time-reversal reversal!
    !   JY       : matrix element of Jy T for this state.
    !              Note extra time-reversal reversal!
    !   JZ       : matrix element of Jz T for this state.
    !   JJ       : J quantum number (real number) that corresponds to this state
    !              such that 
    !                    (JJ+1) JJ = <Jx^2> + <Jy^2> + <Jz^2> 
    !   SxT      : matrix element of S_x T (real part)
    !   SyT      : matrix element of S_y T (imaginary part)
    !   Sz       : matrix element of S_z 
    !---------------------------------------------------------------------------
    use wavefunctions

    character(len=*), intent(in) :: fname
    integer                      :: io, i, wave
    integer                      :: ProtonOrder(nwp), NeutronOrder(nwn)
    real(KIND=dp)                :: Jx, Jy, Jz, JJ, Spinx, Spiny, Spinz, P
 
    1 format(2i5, 10f10.4)
    2 format("# Neutron spwfs")
    3 format("# Proton spwfs")
    4 format("# Information in the basis that diagonalizes RHO")

    60 format ("#",3x,'i',3x,'iso',3x,'P',4x,'occ',7x,'<h>',7x,  &
     &        'JxT',7x,'JyT', 7x ,'Jz', 8x, 'J', 9x, 'SxT', 7x,'SyT',7x,'Sz') 

    open(1,file=fname, iostat=io)
    if(io.ne.0) then    
      print *, 'Something went wrong with the sp. info to file.'
      print *, 'filename = ', fname
      call stp('')
    endif
    
    call write_header(1)
    write(1, fmt=4)
    write(1, fmt='(a1)') '#'
    write(1, fmt=60) 
    write(1, fmt=2) 

    neutronorder = OrderSpwfsISO(-1,.true.)
    protonorder  = OrderSpwfsISO(+1,.true.)
    !---------------------------------------------------------------------------
    ! First do the neutron wavefunctions    
    do i=1,nwn
      wave = neutronorder(i)    
      P = P_can(wave)

      Jx = can_JTR(1,wave) ; Spinx = can_STR (1,wave)
      Jy = can_JTI(2,wave) ; Spiny = can_STI (2,wave)
      Jz = can_J  (3,wave) ; Spinz = can_spin(3,wave)
      JJ = can_JJ(wave)

      write(1, fmt=1) wave, -1, p, rho_can(wave), canenergies(wave), Jx, Jy,Jz,&
      &               JJ, Spinx, Spiny, Spinz
    enddo      
    write(1, fmt=3) 
    !---------------------------------------------------------------------------
    ! Then do the proton wavefunctions    
    do i=1,nwp
      wave =  protonorder(i)
      P = P_can(wave)

      Jx = can_JTR(1,wave) ; Spinx = can_STR (1,wave)
      Jy = can_JTI(2,wave) ; Spiny = can_STI (2,wave)
      Jz = can_J  (3,wave) ; Spinz = can_spin(3,wave)
      JJ = can_JJ(wave)

      write(1, fmt=1) wave, +1, p, rho_can(wave), canenergies(wave), Jx, Jy,Jz,&
      &               JJ, Spinx, Spiny, Spinz
    enddo
    close(1)
  end subroutine write_sp_info_can
  
  subroutine write_blocked_sps(fname)
    !---------------------------------------------------------------------------
    ! Write the density of the 'blocked' single-particle states to file, i.e. 
    ! the single-particle states in the canonical basis that come out with 
    ! rho_can = 1 due to blocking.
    !
    ! The file contains a header written by the subroutine write_header, 
    ! supplemented by the following 
    !
    ! # X [fm] Y [fm] Z [fm]   Psi_1  Psi_2  Psi_3 ....
    !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    !
    !---------------------------------------------------------------------------
    character(len=*), intent(in) :: fname
    integer                      :: io
    integer                      :: i,j,k, wave
    real(KIND=dp), allocatable   :: psis(:,:,:,:) 
    real(KIND=dp), pointer       :: tempwf_one(:,:,:), tempwf_two(:,:,:)
    real(KIND=dp), pointer       :: tempwf_three(:,:,:), tempwf_four(:,:,:)

    1 format('#  X[fm]   Y[fm]   Z[fm]')
    2 format(7x, ' |Psi_', i1, '|^2' , 10x)

    allocate(psis(nx,ny,nz,blocknumber)) ; psis = 0

    if(.not.allocated(blocked_sps)) then
      print *, 'Cannot write single-particle wavefunctions to file.'
      deallocate(psis)
      return
    endif

    open(1,file=fname, iostat=io)
    if(io.ne.0) then    
      print *, 'Something went wrong with writing blocked states to file.'
      print *, 'filename = ', fname
      call stp('')
    endif
    
    call write_header(1)
    write(1, fmt=1, advance='no')
    do k=1, blocknumber
      write(1, fmt=2, advance ='no') k
    enddo
    write(1, fmt=*)
    
    do wave=1,blocknumber
      tempwf_one(1:nx,1:ny,1:nz)   => canpsi(1:nx*ny*nz,1,wave)
      tempwf_two(1:nx,1:ny,1:nz)   => canpsi(1:nx*ny*nz,2,wave)
      tempwf_three(1:nx,1:ny,1:nz) => canpsi(1:nx*ny*nz,3,wave)
      tempwf_four(1:nx,1:ny,1:nz)  => canpsi(1:nx*ny*nz,4,wave)
      do k=1, nz
        do j=1,ny
          do i=1,nx
            psis(i,j,k,wave) =     tempwf_one(i,j,k)**2   &
            &                     +tempwf_two(i,j,k)**2   &
            &                     +tempwf_three(i,j,k)**2 &
            &                     +tempwf_four(i,j,k)**2 
          enddo
        enddo
      enddo
    enddo
    
    do k=1,nz
      do j=1,ny
        do i=1,nx
          write(1, fmt='(3f8.3)', advance='no')  meshx(i), meshx(j), meshz(k)
          do wave=1,blocknumber
            write(1, fmt='(1x,es25.12,1x)', advance='no') psis(i,j,k,wave) 
          enddo
          write(1, fmt=*)
        enddo
      enddo
    enddo
      
  end subroutine write_blocked_sps

  subroutine write_inertias(fname)
    !---------------------------------------------------------------------------
    ! Write information on collective inertia's to file:
    !
    !    a) the total collective inertia tensor for 
    !       1. neutrons
    !       2. protons
    !       3. total
    !    b) the M1 intermediate matrix
    !       1. neutrons
    !       2. protons
    !    c) the M3 intermediate matrix
    !       1. neutrons
    !       2. protons
    !
    ! Input :
    !   fname : filename to write information to.
    !
    !---------------------------------------------------------------------------
    use fission_MOI
    use moments
      
    character(len=*), intent(in) :: fname
    integer                      :: io
    integer                      :: k,it, l, m
    type(Moment), pointer        :: current 

    1 format ('#', 18x)
    2 format (' B_{ ', i2, i2, '}     ')
   21 format (14x)
    3 format (f15.3)
    4 format ('# Neutrons')
    5 format ('# Protons')
    6 format ('# Total')
    7 format (' Q_{ ', i2, 1x, i2, '} | ', 99es15.5 )
   
    9 format ('#', 50('-'), 'Collective inertia tensor', 50('-'))
   91 format ('#', 50('-'), '    Intermediate M^1     ', 50('-'))
   92 format ('#', 50('-'), '    Intermediate M^3     ', 50('-'))
   
    open(1,file=fname, iostat=io)
    if(io.ne.0) then    
      print *, 'Something went wrong with writing collective inertias to file.'
      print *, 'filename = ', fname
      call stp('')
    endif
    
    call write_header(1)
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Values of all relevant multipole moments
    write(1, fmt=1, advance ='no') 
    do k=1, N_inertia
       l = inertia_l(k)
       m = inertia_m(k)
       write(1, fmt=2, advance ='no') l, m
    enddo
    write(1, fmt=*)
    write(1, fmt=21, advance ='no')
    do k=1, N_inertia
       l = inertia_l(k)
       m = inertia_m(k)
       current => findmoment(l, m, .false.)      
       if(associated(current)) then
         ! We only write the value to file if the multipole moment can be 
         ! found, i.e. when it is not restricted by symmetry
         write(1, fmt=3, advance ='no') current%beta(4)
       else
         ! This multipole moment was not found, we just write 0.0 to file
         write(1, fmt=3, advance ='no') 0.0d0
       endif
    enddo
    write(1, *)
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Collective inertia
    write(1, fmt=1)
    write(1, fmt=9)
    write(1, fmt=6)
    do k=1, N_inertia
        l = inertia_l(k)
        m = inertia_m(k)
        write(1, fmt=7) l,m,collective_inertia(k,1:N_inertia)
    enddo

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Intermediate matrix M1
    write(1, fmt=1)
    write(1, fmt=91)
    do it=1,2
      select case(it)
      ! FORTRAN does not seem to allow for calculated fmt = it + 3 statements,
      ! so hardcoding it is.
      case(1)
        write(1, fmt=4)
      case(2)
        write(1, fmt=5)
      case(3)
        write(1, fmt=6)
      end select      
      do k=1, N_inertia
        l = inertia_l(k)
        m = inertia_m(k)
        write(1, fmt=7) l,m,M1(k,1:N_inertia ,it)
      enddo    
    enddo

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Intermediate matrix M3
    write(1, fmt=1)
    write(1, fmt=92)
    do it=1,2
      select case(it)
      ! FORTRAN does not seem to allow for calculated fmt = it + 3 statements,
      ! so hardcoding it is.
      case(1)
        write(1, fmt=4)
      case(2)
        write(1, fmt=5)
      case(3)
        write(1, fmt=6)
      end select      
      do k=1, N_inertia
        l = inertia_l(k)
        m = inertia_m(k)
        write(1, fmt=7) l,m,M3(k,1:N_inertia ,it)
      enddo    
    enddo
  
  end subroutine write_inertias
#endif

#if( $FAM == 1)
   subroutine init_fam_file(l, m, eff_charge_n, eff_charge_p, fname)
    !---------------------------------------------------------------------------
    ! Create file to write strength function S(omega, F) obtained from FAMtalus
    !---------------------------------------------------------------------------
    ! The file contains a header written by the subroutine write_header,
    ! supplemented by a dedicated line explaining the content of each column.
    ! The format of the body of said file is
    !   omega[MeV]    S_free[...]     S[...]   iter
    !                         '-> unit depends on the external field                  
    ! 
    ! The actual strength is written to this file by subroutine append_fam_file()
    ! called each time a frequency is converged. 
    !---------------------------------------------------------------------------
    integer, intent(in)               :: l, m
    real(kind=DP), intent(in)         :: eff_charge_n, eff_charge_p
    character(len=*), intent(in)      :: fname
    integer                           :: io

    print *, ' writing strength function to file :  ', fname

    1 format ( '# external field:   ', /, &
    &          '#    F = Q_', i1, i1,/, &
    &          '#    neutron eff charge = ', f10.3, ' e', /, &
    &          '#    proton eff charge  = ', f10.3, ' e')
    2 format('#', 2x, 'omega',12x,'S_free', 21x, 'S', 17x, 'iter')


    open(1,file=fname, iostat=io)
    if(io.ne.0) then    
      print *, 'filename = ', fname
      call stp('')
    endif
    
    call write_header(1) ! write general header info

    write(1, fmt=1) l, m, eff_charge_n, eff_charge_p ! write info of extrenal field 
    write(1, fmt=2)      ! write column names

    close(1)

  end subroutine init_fam_file

  subroutine init_fam_file_new(fname)
    !---------------------------------------------------------------------------
    ! Create file to write strength function S(omega, F) obtained from FAMtalus
    !---------------------------------------------------------------------------
    ! The file contains a header written by the subroutine write_header,
    ! supplemented by a dedicated line explaining the content of each column.
    ! The format of the body of said file is
    ! 
    !  omega[MeV] gamma[MeV] iter S_complex_re S_complex_im S_n+ S_n- S_p+ S_p- S_tot
    ! 
    ! The units of the strength depend on the external field F. 
    !
    ! The actual strength is written to this file by subroutine append_fam_file()
    ! called each time a frequency is converged. 
    !---------------------------------------------------------------------------
    use fam
    character(len=*), intent(in)      :: fname
    integer                           :: io

    print *, ' writing strength function to file :  ', fname

    1 format ( '# external field:   ', /, &
    &          '#    F = Q_', i1, i1,/, &
    &          '#    neutron eff charge = ', f10.3, ' e', /, &
    &          '#    proton eff charge  = ', f10.3, ' e')

    2 format ( '# sum rules: ', / , '#   m1 = ', es20.8)
    3 format('#', 4x, 'omega', 5x, 'gamma',4x, 'iter',  8x,'S_complex_re', 13x, &
    &  'S_complex_im', 18x, 'S_n+', 21x, 'S_n-', 21x, 'S_p+', 21x, 'S_p-', 20x, &
    &  'S_tot') 


    open(1,file=fname, iostat=io)
    if(io.ne.0) then    
      print *, 'filename = ', fname
      call stp('')
    endif
    
    call write_header(1) ! write general header info

    write(1, fmt=1) l, m, eff_charge_n, eff_charge_p ! write info of extrenal field 
    write(1, fmt=2) ewsr ! write sum rules  
    write(1, fmt=3)      ! write column names

    close(1)

  end subroutine init_fam_file_new

  subroutine append_fam_file(omega, S, iter, S_free, fname)
    real(kind=dp), intent(in)         :: omega, S, S_free
    integer, intent(in)               :: iter
    character(len=*), intent(in)      :: fname
    integer                           :: io

    print *, ' append fam file :  ', fname

    open(1, file=fname, status='old', position='append', iostat=io)
    if(io.ne.0) then    
      print *, 'filename = ', fname
      call stp('')
    endif
    
    write(1, fmt='(f8.3, es25.12E3, es25.12E3, i10)') omega, S_free, S, iter
      
    close(1)

  end subroutine append_fam_file


  subroutine append_fam_file_new(S_decomp, iter, fname)
    use fam
    real(KIND=dp), intent(in)         :: S_decomp(8)
    integer, intent(in)               :: iter
    character(len=*), intent(in)      :: fname
    integer                           :: io

    1 format (f10.3, f10.3, i7, es25.12E3, es25.12E3, es25.12E3, es25.12E3, &
    & es25.12E3, es25.12E3, es25.12E3) 

    print *, ' append fam file :  ', fname

    open(1, file=fname, status='old', position='append', iostat=io)
    if(io.ne.0) then    
      print *, 'filename = ', fname
      call stp('')
    endif
    
    write(1, 1) omega_fam, smear, iter, strength_complex%re, strength_complex%im, &
    & sum(S_decomp(1:2)), sum(S_decomp(3:4)),sum(S_decomp(5:6)), sum(S_decomp(7:8)), sum(S_decomp(:))

    close(1)

  end subroutine append_fam_file_new

  subroutine init_xy_file(fname)
    !---------------------------------------------------------------------------
    ! Create file to write X and Y amplitudes for each freq obtained from FAMtalus
    ! as well as the perturbing operator F. 
    !---------------------------------------------------------------------------
    ! The file contains a header written by the subroutine write_header,
    ! The complex matrices X_mn Y_mn are written in a sparse format as 
    !    m    n    X_mn%re   X_mn%im     Y_mn%re   Y_mn%im
    ! The file is appended for each FAM frequency, different blocks seperated by 
    ! a single line:
    !   & omega = [omega]  [smear]
    ! At the top of the file, the perturbing operator F is written. 
    !---------------------------------------------------------------------------
    use fam
    character(len=*), intent(in)      :: fname
    integer                           :: io

    print *, ' writing XY to file :  ', fname

    1 format ( '# external field:   ', /, &
    &          '#    F = Q_', i1, i1,/, &
    &          '#    neutron eff charge = ', f10.3, ' e', /, &
    &          '#    proton eff charge  = ', f10.3, ' e')

    2 format ( '# sum rules: ', / , '#   m1 = ', es20.8)
    3 format('#', 5x, 'm', 6x, 'n',14x, 'X(/F)_mn_re', 14x, 'X(/F)_mn_im', 14x, 'Y(/F)_mn_re', 14x, 'Y(/F)_mn_im') 


    open(1,file=fname, iostat=io)
    if(io.ne.0) then    
      print *, 'filename = ', fname
      call stp('')
    endif
    
    call write_header(1) ! write general header info

    write(1, fmt=1) l, m, eff_charge_n, eff_charge_p ! write info of extrenal field 
    write(1, fmt=2) ewsr ! write sum rules  
    write(1, fmt=3)      ! write column names

    close(1)


  end subroutine init_xy_file


  subroutine append_xy_file(fname, O20, O02)
    use fam
    character(len=*), intent(in)           :: fname
    complex(KIND=dp), intent(in), optional :: O20(:,:), O02(:,:)
    integer                                :: io, mu, nu

    1 format ( '& omega = ', f10.3, f10.3) 
    2 format (i7, i7, es25.12E3, es25.12E3, es25.12E3, es25.12E3) 

    print *, ' append fam file :  ', fname

    open(1, file=fname, status='old', position='append', iostat=io)
    if(io.ne.0) then    
      print *, 'filename = ', fname
      call stp('')
    endif
    

    if (present(O20) .and. present(O02)) then
      do nu = 1, nwt
        do mu = 1, nwt
          if(abs(O20(mu,nu)) > 1e-10 .or. abs(O02(mu,nu)) > 1e-10) then
            write(1, fmt=2) mu, nu, O20(mu,nu)%re, O20(mu,nu)%im, O02(mu,nu)%re, O02(mu,nu)%im
          end if
        enddo
      enddo

    else
      write(1, fmt=1) omega_fam, smear

      do nu = 1, nwt
        do mu = 1, nwt
          if(abs(X(mu,nu)) > 1e-10 .or. abs(Y(mu,nu)) > 1e-10) then
            write(1, fmt=2) mu, nu, X(mu,nu)%re, X(mu,nu)%im, Y(mu,nu)%re, Y(mu,nu)%im
          end if
        enddo
      enddo
    endif

    close(1)

  end subroutine append_xy_file


subroutine init_perturbed_denfile(fname)
    !---------------------------------------------------------------------------
    ! Create file named "fname"to write the following perturbed densities 
    !    drho(neutron), drho(proton), drho(charge)
    !---------------------------------------------------------------------------
    ! The file contains a header written by the subroutine write_header,
    ! supplemented by a dedicated line explaining the content of each column.
    ! The format of the body of said file is
    ! 
    !   x, y, z, drho_n_sym%re, drho_n_sym%im, drho_n_asym%re, drho_n_asym%im, 
    !     drho_p_sym%re, drho_p_sym%im, drho_p_asym%re, drho_p_asym%im, 
    !     drho_c_sym%re, drho_c_sym%im, drho_c_asym%re, drho_c_asym%im
    !
    ! where the first three numbers are the Cartesian coordinates in fm, withµ
    ! the densities all in their natural units. The mesh points are traverse in 
    ! column-major order ('Fortran order'), which might not be how your favorite 
    ! plotting tool prefers it. Note that the densities are written "as-is" to 
    ! file, i.e. only in part of the box that is actually represented 
    ! numerically. It is up to postprocessing to construct the densities in the 
    ! simulation volume.
    !---------------------------------------------------------------------------
    use fam
    character(len=*), intent(in)            :: fname
    integer                                 :: io

    print *, ' writing perturbed densities to file :  ', fname

    1 format ( '# external field:   ', /, &
    &          '#    F = Q_', i1, i1,/, &
    &          '#    neutron eff charge = ', f10.3, ' e', /, &
    &          '#    proton eff charge  = ', f10.3, ' e')

    2 format ( '# sum rules: ', / , '#   m1 = ', es20.8)
    3 format('#', 19x, 'X[fm]',20x,'Y[fm]', 20x,'Z[fm]',  &
      & 12x, 'drho_n_sym_re', 12x , 'drho_n_sym_im', 11x, 'drho_n_asym_re', 11x , 'drho_n_asym_im', &
      & 12x, 'drho_p_sym_re', 12x , 'drho_p_sym_im', 11x, 'drho_p_asym_re', 11x , 'drho_p_asym_im', &
      & 12x, 'drho_c_sym_re', 12x , 'drho_c_sym_im', 11x, 'drho_c_asym_re', 11x , 'drho_c_asym_im')


    open(1,file=fname, iostat=io)
    if(io.ne.0) then    
      print *, 'filename = ', fname
      call stp('')
    endif
    
    call write_header(1) ! write general header info

    write(1, fmt=1) l, m, eff_charge_n, eff_charge_p ! write info of extrenal field 
    write(1, fmt=2) ewsr ! write sum rules  
    write(1, fmt=3)      ! write column names

    close(1)

  end subroutine init_perturbed_denfile

subroutine append_perturbed_denfile(Rs, Ra, fname)
    !---------------------------------------------------------------------------
    ! Append the the file named "fname" the following perturbed densities 
    !    drho(neutron), drho(proton), drho(charge)
    !---------------------------------------------------------------------------
    ! Different blocks corresponding to different fam frequencies omega are 
    ! separated by a single line 
    !     & omega = [omega%re] [omega%im]
    !---------------------------------------------------------------------------
    use fam
    type(DensityVector), intent(in), target :: Rs, Ra
    character(len=*), intent(in)            :: fname
    integer                                 :: io, i,j,k, mi


    1 format ( '& omega = ', f10.3, f10.3) 

    print *, ' append fam file :  ', fname

    open(1, file=fname, status='old', position='append', iostat=io)
    if(io.ne.0) then 
      print *, 'Something went wrong with writing a density to file.'
      print *, 'filename = ', fname
      call stp('')
    endif
    
    write(1, fmt=1) omega_fam, smear

    do k=1,nz
      do j=1,ny
        do i=1,nx
          write(1, fmt='(3es25.12)', advance='no') &
          &          meshx(i), meshy(j), meshz(k)
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          ! the contributions above are indexed according to (x,y,z) but 
          ! we do not have this luxury for most of the densities
          mi = meshindex(i,j,k)
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          ! The ordinary and charge density; always defined
          write(1, fmt='(2es25.12)', advance='no') Rs%D_I_I(mi,1)%re,Rs%D_I_I(mi,1)%im
          write(1, fmt='(2es25.12)', advance='no') Ra%D_I_I(mi,1)%re,Ra%D_I_I(mi,1)%im
          write(1, fmt='(2es25.12)', advance='no') Rs%D_I_I(mi,2)%re,Rs%D_I_I(mi,2)%im
          write(1, fmt='(2es25.12)', advance='no') Ra%D_I_I(mi,2)%re,Ra%D_I_I(mi,2)%im
          write(1, fmt='(2es25.12)', advance='no') Rs%chargedensity(i,j,k)%re,Rs%chargedensity(i,j,k)%im
          write(1, fmt='(2es25.12)', advance='no') Ra%chargedensity(i,j,k)%re,Ra%chargedensity(i,j,k)%im
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          ! We are done writing this line in the output
          write(1, fmt='()') !  newline character
        enddo
      enddo
    enddo

    close(1)
  end subroutine append_perturbed_denfile

#endif

  function force_halfinteger(j) result(jforced)
      !-------------------------------------------------------------------------
      ! Small function to round a real number to a half integer number, 
      ! useful for angular momenta.
      !
      !-------------------------------------------------------------------------
      real(KIND=dp) :: j, jforced
      integer       :: i

      i = 1
      do while(i .lt. 2*abs(j))     
          i = i + 2 
      enddo

      if(abs(2*abs(j) - i) .lt. abs(2*abs(j)-i+2)) then
          jforced = i/2.0_dp
      else
          jforced = i/2.0_dp - 1
      endif
  end function force_halfinteger

#if ($FAM == 0)
  subroutine combi_output(COMBI)
    !---------------------------------------------------------------------------
    ! Write an extra file to serve as input to a combinatorial calculation 
    ! of the nuclear level density.
    !
    ! ATTENTION: this output assumes an axial nucleus with a symmetry axis 
    !            along the z-axis. If the single-particle states are not  
    !            (at least approximately) eigenstates of J_z, then this output
    !            will effectively be nonsense.
    !
    ! Additional note: the MOI that are written are the "COLLECTIVE" Belyaev 
    !                  values, i.e. those without the contributions from any
    !                  blocked qps.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! TODO: Document the precise output of this routine.
    !
    !---------------------------------------------------------------------------
    use Moments, only: multipolefactor, calculatetotalQl
    
    character(len=*), intent(in) :: combi
    integer, allocatable :: indices(:)
    integer              :: i,ii, p1, p2,jj
    real(KIND=dp)        :: A, mstate1, mstate2
    real(KIND=dp)        :: b2, b3,b4, q2(4), q3(4), q4(4)
    real(KIND=dp), allocatable :: tempgaps(:,:)

    1 format (a1, 3i4)
    2 format (2(f5.1,i2,3f8.3))
    3 format ( 2i4,2(x,f8.4),19(x,f7.3),2f12.3)

    open(unit=6, file=COMBI)
    
    select case(pairingtype)
    case(0)
      ! HF: we don't do nothing for the moment
    case(1)
      ! BCS case: things are straightforward
      !-------------------------------------------------------------------------
      ! a) single-particle neutron states
$TR      write(6, fmt=1) '*', int(protons), int(neutrons+protons),  +nwn/2
$NTR     write(6, fmt=1) '#', int(protons), int(neutrons+protons),  +nwn/2
      indices = OrderSpwfsISO(-1)
      
      do i=1,nwn/2
          ii    = indices(i)
          jj    = indices(i+nwn/2)

          if(ii .lt. (HFBlocks(1))) p1 =  0
          if(ii .gt. (HFBlocks(1))) p1 =  1

          if(jj .lt. (HFBlocks(1))) p2 =  0
          if(jj .gt. (HFBlocks(1))) p2 =  1

          mstate1 = angmom_z_real(HFPsi(:,:,ii),HFPsi(:,:,ii),HFdPsi(:,:,:,ii))
          mstate2 = angmom_z_real(HFPsi(:,:,jj),HFPsi(:,:,jj),HFdPsi(:,:,:,jj))

          mstate1 = force_halfinteger(mstate1)
          mstate2 = force_halfinteger(mstate2)

          write(6,fmt=2) mstate1,p1,spenergies(ii),rho_can(ii)/2,BCSgaps(ii),& 
          &              mstate2,p2,spenergies(jj),rho_can(jj)/2,BCSgaps(jj) 
      enddo
      !-------------------------------------------------------------------------
      ! b) single-particle proton states
$TR      write(6, fmt=1) ' ', int(protons), int(neutrons+protons),+nwp/2
$NTR     write(6, fmt=1) ' ', int(protons), int(neutrons+protons),+nwp/2
      indices = OrderSpwfsISO(+1)
      do i=1,nwp/2
          ii     = indices(i)
          jj     = indices(i+nwp/2)

          if(ii .lt. sum(HFBlocks(1:3))) p1 =  0
          if(ii .gt. sum(HFBlocks(1:3))) p1 =  1

          if(jj .lt. sum(HFBlocks(1:3))) p2 =  0
          if(jj .gt. sum(HFBlocks(1:3))) p2 =  1

          mstate1 = angmom_z_real(HFPsi(:,:,ii),HFPsi(:,:,ii),HFdPsi(:,:,:,ii))
          mstate2 = angmom_z_real(HFPsi(:,:,jj),HFPsi(:,:,jj),HFdPsi(:,:,:,jj))

          mstate1 = force_halfinteger(mstate1)
          mstate2 = force_halfinteger(mstate2)
          
          write(6,fmt=2) mstate1,p1,spenergies(ii),rho_can(ii)/2,BCSgaps(ii),& 
          &              mstate2,p2,spenergies(jj),rho_can(jj)/2,BCSgaps(jj) 
      enddo
    !---------------------------------------------------------------------------
    case(2)
      ! HFB case: things are less straightforward
      tempgaps = HFBGaps

      if((.not. diagsphamil) .and. allocated(HFTransfo)) then
        call calc_gaps_HF(tempgaps, HFTransfo)
      endif
      !-------------------------------------------------------------------------
      ! a) single-particle neutron states
$TR      write(6, fmt=1) '*', int(protons), int(neutrons+protons),+nwn/2
$NTR     write(6, fmt=1) '#', int(protons), int(neutrons+protons),+nwn/2
      indices = OrderSpwfsISO(-1)
      
       do i=1,nwn/2
          ii    = indices(i)
          jj    = indices(i+nwn/2)

          if(ii .le. sum(HFBlocks(1:2))) p1 =  0
          if(ii .gt. sum(HFBlocks(1:2))) p1 =  1

          if(jj .le. sum(HFBlocks(1:2))) p2 =  0
          if(jj .gt. sum(HFBlocks(1:2))) p2 =  1
          
          mstate1 = angmom_z_real(HFPsi(:,:,ii),HFPsi(:,:,ii),HFdPsi(:,:,:,ii))
          mstate2 = angmom_z_real(HFPsi(:,:,jj),HFPsi(:,:,jj),HFdPsi(:,:,:,jj))

          mstate1 = force_halfinteger(mstate1)
          mstate2 = force_halfinteger(mstate2)
          
$TR          write(6,fmt=2) mstate1,p1,spenergies(ii),rho_HF(ii),tempgaps(ii,ii),& 
$TR          &              mstate2,p2,spenergies(jj),rho_HF(jj),tempgaps(jj,jj)

$NTR         write(6,fmt=2) mstate1,p1,spenergies(ii),rho_HF(ii),maxval(abs(tempgaps(ii,:))),& 
$NTR         &              mstate2,p2,spenergies(jj),rho_HF(jj),maxval(abs(tempgaps(jj,:))) 
      enddo
      !-------------------------------------------------------------------------
      ! b) single-particle proton states
$TR       write(6, fmt=1) ' ', int(protons), int(neutrons+protons),+nwp/2
$NTR      write(6, fmt=1) ' ', int(protons), int(neutrons+protons),+nwp/2
      indices = OrderSpwfsISO(+1)
      do i=1,nwp/2
          ii     = indices(i)
          jj     = indices(i+nwp/2)

          if(ii .le. sum(HFBlocks(1:6))) p1 =  0 
          if(ii .gt. sum(HFBlocks(1:6))) p1 =  1

          if(jj .le. sum(HFBlocks(1:6))) p2 =  0
          if(jj .gt. sum(HFBlocks(1:6))) p2 =  1

          mstate1 = angmom_z_real(HFPsi(:,:,ii),HFPsi(:,:,ii),HFdPsi(:,:,:,ii))
          mstate2 = angmom_z_real(HFPsi(:,:,jj),HFPsi(:,:,jj),HFdPsi(:,:,:,jj))

          mstate1 = force_halfinteger(mstate1)
          mstate2 = force_halfinteger(mstate2)

$TR      write(6,fmt=2) mstate1,p1,spenergies(ii),rho_HF(ii),tempgaps(ii,ii),& 
$TR      &              mstate2,p2,spenergies(jj),rho_HF(jj),tempgaps(jj,jj) 
          
$NTR      write(6,fmt=2) mstate1,p1,spenergies(ii),rho_HF(ii),maxval(abs(tempgaps(ii,:))),& 
$NTR      &              mstate2,p2,spenergies(jj),rho_HF(jj),maxval(abs(tempgaps(jj,:)))  
    enddo
    !---------------------------------------------------------------------------
    end select
    !  The final line is composed of various informations read by the 
    !  gennew code. The FORTRAN read statement is
    !
    !            READ(15,*,err=5) IZ,IA,BETA,HFBET4,HFGN,HFGP,
    ! &       HFDN,HFDP,HFDDN,HFDDP,HFEN,HFEP,HFUN,HFUP,HFLN,
    ! &       HFLP,HFINX,HFIPX,HFINY,HFIPY,HFINZ,HFIPZ,HFJ2,HFE1,HFE2    
    ! 
    ! NOTE: the level density code assumes that, if the nucleus is axially
    !       symmetric, the symmetry axis is the z-axis. It is the single-particle
    !       expectation values of Jz that are written on file.
    !---------------------------------------------------------------------------

    ! Calculate the total deformations Q_l
    q2 = CalculateTotalQl(2) 
    q3 = CalculateTotalQl(3)
    q4 = CalculateTotalQl(4)
    
    ! Rescale to dimensionless quantities
    A = neutrons + protons
    b2 = q2(4) * MultipoleFactor(A,A,2)
    b3 = q3(4) * MultipoleFactor(A,A,3)
    b4 = q4(4) * MultipoleFactor(A,A,4)

    ! Note: items marked with (*) are written as zero and, to the best of
    ! my (=W.R.) knowledge, not used by the level density code.  
    !                        IZ          IA    BETA     B4
    write(unit=6, fmt=3)  int(protons),int(A),  b2, b4,  &
    !                      HGN    HFGP  HFDN             HFDP
    &                      b3,    0.0,  average_gap(2,1),average_gap(2,2), &
    !                      HFDDN,HFDDP,HFEN,HFEP,HFUN,HFUP
    &                      0.0,   0.0,  0.0, 0.0, 0.0, 0.0, & 
    !                      HFLN, HFLP, 
    &                      FermiEnergy(1), FermiEnergy(2),  &
    !                      HFINX /     ,    HFIPX
    &                      Bely_coll(1,1), Bely_coll(1,2),       &
    !                      HFINY /     ,    HFIPY,      
    &                      Bely_coll(2,1), Bely_coll(2,2),       &
    !                      HFINZ /     ,    HFIPZ,      
    &                      Bely_coll(3,1), Bely_coll(3,2),       &
    !                      HFJ2               HFE1  , HE2
    &                      J2_coll(3,3),  totalE, 0.0

    close(unit=6)
  end subroutine combi_output
#endif    
end module IO
