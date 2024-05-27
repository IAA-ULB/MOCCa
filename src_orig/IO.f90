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
 ! Module governing the in- and output of Tantalus. 
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
 !==============================================================================

use geninfo
use wavefunctions
use pairing
use functional
use momentsofinertia
use moments
use Coulombmod
use transform

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
  ! Filenames for in- and output of the code with respect to spwfs.
  character(len=100)  :: inputfilename, outputfilename
  ! Signal the code to write extra output.
  character(len=80)   :: BXLFIT='', COMBI='', denfile='', potfile=''
  character(len=80)   :: sphffile='', spcanfile='', tofile='', blockfile=''
  character(len=80)   :: inertfile=''
  ! Signal the code to write the wavefunctions periodically to disk
  integer             :: checkpointiter = 0  
  !-----------------------------------------------------------------------------
  logical             :: Allowtransform = .false.
  integer             :: extraspwfs(8) = 0
  !-----------------------------------------------------------------------------
  ! Encodings of the symmetry choices imposed by Hephaestos
  character(len=26), parameter :: SYM_CODE   = "$SYM_CODE"
  character(len=26), parameter :: TRANS_CODE = "$TRANS_CODE"
  !-----------------------------------------------------------------------------
  ! Characteristics of the calculation stored on the .wf file
  integer              :: filenx, fileny, filenz, filemv
  integer              :: filenwn, filenwp, filepairing
  integer              :: filenwt, fileneutrons, fileprotons
  integer              :: fileblocks_global(8), fileblocks(8)
  integer, allocatable :: file_spwf_map(:),file_rank_map(:),file_spwf_inverse(:)
  real(KIND=dp) :: filedx
  !-----------------------------------------------------------------------------
  ! Did we succeed in reading a HFB configuration from file? 
  logical       :: readHFBinfofile= .false.
  ! Blocking information from file
  integer       :: fileblocktype=0, fileblocknumber = 0
  integer, allocatable          :: fileblockindices(:)
  character(len=2), allocatable :: fileBlockLowest(:)
  integer                       :: file_HFB_blocks(8)
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
    use fission_moi,   only : read_inertia
  
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
    call read_inertia(file_number)
    call readmomentdata(file_number)
    call readcranking(file_number)

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
    use fission_moi,   only : N_inertia

    integer(dp), intent(in), optional   :: file_number   
#if(USE_MPI>0)
    integer                             :: mpi_err
#endif

    NameList /IO/ InputFileName,OutputFileName, BXLFIT, COMBI, denfile,potfile,& 
    &           sphffile, spcanfile,checkpointiter, AllowTransform, extraspwfs,&
    &           tofile, blockfile, inertfile, N_inertia

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
    call MPI_Bcast(N_inertia     , 1                  , MPI_INTEGER, 0, &
    &                                                   MPI_COMM_WORLD, mpi_err)

    call MPI_Bcast(allowtransform, 1                  , MPI_LOGICAL, 0, &
    &                                                   MPI_COMM_WORLD, mpi_err)
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
   99 format ( '  Nilsson initialization with hom = (', 3(f7.3) ,')')
   10 format ( ' IO information', / &
    &          '  inputfilename  =', a32, / &
    &          '  outputfilename =', a32)
    
  101 format ( ' Information obtained from file ')  
  102 format ( '      - version number          : ', i5)  
 1021 format ( '      - param. used on file     : ', 20a)
  103 format ( '      - Bogoliubov transfo read?: ', l5)  
 1031 format ( '      - Bogoliubov transfo used?: ', l5)  
  104 format ( '      - Blocking type           : ', i5)
  105 format ( '      - Blocknumber             : ', i5)
  106 format ( '      - Block indices           : ', 10i4)
  107 format ( '      - Block lowest            : ', 10a2)
  108 format ( '      - Passed blocking test    : ', l5)
    
   11 format ( ' Filename for other output (not written if empty): ', /     &
             & '    BXL output     = ', a80, / &
             & '    DEN file       = ', a80, / &
             & '    POT file       = ', a80, / &
             & '    SPHF file      = ', a80, / &
             & '    SPCAN file     = ', a80, / &
             & '    TO file        = ', a80, / & 
             & '    BLOCK file     = ', a80, / &
             & '    INERT file     = ', a80) 
 1111 format ( '    Input data     = ', a26, / &
               '     on unit ', i10)
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
      if(trim(to_upper(inputfilename)).eq.'INIT') print 99, osc_freq

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
      endif 

      print 112, checkpointiter
      print 113, print_adv_spwf_properties

      print 11, BXLFIT, DENFILE, POTFILE, SPHFFILE, SPCANFILE, TOFILE, BLOCKFILE, INERTFILE
      if(present(file_number)) then
        print 1111,  adjustl(trim(input_file)), file_number
      endif
      print 12, energy_prec, moment_prec, disp_prec, gradient_prec, fermi_prec,  &
      &         angmom_prec
      
      call printevolution
      call printscfiteration
      call printpairing_init
      call printmoment_init
      call printcranking_init
      call printfunctional  
    endif    

  end subroutine PrintInput

  subroutine Readwavefunction()
    !---------------------------------------------------------------------------
    ! High-level routine to determine the starting point of a calculation. 
    !
    ! There are two main starting options, one of which has two suboptions
    !
    ! 1) Initialize in an EV8-style box with Nilsson orbitals
    !    a - start self-consistency cycles immediately
    !    b - read a set of potentials from file to start the calculations
    ! 
    ! 2) Read a set of spwfs from file 
    ! 
    ! Which option is chosen based on the InputFileName keyword: 
    !  - INIT (case insensitive) : option 1a, no reading of any file
    !  - *.pot                   : option 1b, reading of a potential file
    !  - [any other filename]    : option 2, reading of a wavefunction file
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
    else
      ! option 2 : reading complete .wf file
      inputoption = 2        
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
        symtransfo_needed = .true.
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
    else
      ! Option 2) start from a previous calculation.
      call ReadTantalus(12, inputfilename)
      ! No need to guess gaps every time (unless the user asked for it)
    endif
    !---------------------------------------------------------------------------
    ! Transformation options
    if(allowtransform ) then
      if(  symtransfo_needed ) then 
          ! Option a): break a symmetry and transform the spwfs appropriately
          call Transformspwfs( HFPsi, filenx, fileny, filenz,fileblocks_global,&
          &                    fileblocks, file_rank_map, file_spwf_inverse)
      else
          ! Option b): add points and/or add spwfs
          call  TransformInput(filenx,fileny,filenz,filenwn,filenwp,filedx,    & 
          &                               fileblocks,file_HFB_blocks,extraspwfs)
      endif
    else  
      ! Sanity check
      if(symtransfo_needed) then
        call stp('Symmetry transformation needed, but not allowed by user.')
      endif
      ! We still need to set this particular information
      HFblocks  = fileblocks
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
    ! ... and these if (and only if) transformspwfs was not called above
    ! If transformspwfs was called, this assignment was taken care of inside 
    ! that routine.
    if(.not. symtransfo_needed) then
      spwf_map     = file_spwf_map
      rank_map     = file_rank_map
      spwf_inverse = file_spwf_inverse
    endif

    !---------------------------------------------------------------------------
    ! with everything safely in memory, we add in an orthonormalisation to 
    ! guarantee we can start calculating stuff.
#if(USE_MPI > 0)
    ! Copy the 1D wavefunctions to the 2D layout, since that is how we 
    ! orthonormalize ...
    call transfer_1D_to_2D(HFPsi, HFPsi_2D)
#endif
    call  orthonormalize
#if(USE_MPI > 0)
    ! ... and make sure the results get back to the original layout
    call transfer_2D_to_1D(HFPsi_2D, HFPsi)
#endif
    !---------------------------------------------------------------------------
    ! Failsafe for the HF transformation
    if(.not.allocated(HFTransfo)) then
        allocate(HFTransfo(nwt,nwt)) 
        HFtransfo = 0.0d0
        do i=1,nwt
            HFtransfo(i,i) = 1.0d0
        enddo
    endif
    !---------------------------------------------------------------------------
    call set_spwf_symmetries(sx, sy, sz, HFblocks)
    call update_spwf_symmetries(.true.)
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

  subroutine ReadTantalus(chan, ifn)
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
      inquire(file=inputfilename, exist=exists)
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
          symtransfo_needed = .false.
        elseif(SYM_CODE_CHECK .eq. TRANS_CODE) then
          symtransfo_needed = .true.
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

    call MPI_BCAST(symtransfo_needed,1,MPI_LOGICAL, 0, MPI_COMM_WORLD, mpi_err)
#endif    
    ! .. now we have each rank decide what spwfs to take from file
    call loadbalance(fileblocks_global,balancing_strategy, &           ! inputs
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
    potentials_read = readpotentials(chan, filenx,fileny,filenz, symtransfo_needed)
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

      if(symtransfo_needed) then
       if(check_x .or. check_y .or. check_z) then 
        call stp("Please don't combine symmetry transformations and mesh modifications.")
       endif  
       if(check_nwn .or. check_nwp ) then 
        call stp("Please don't combine symmetry transformations and adding wavefunctions.")
       endif  
      endif
    endif
  end subroutine ReadTantalus

  subroutine WriteTantalus(chan, ofn)
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

  end subroutine WriteTantalus

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
      call write_timeodd_densities(TOFILE)
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

    write(filedone,'(a,"z",i3.3,"n",i3.3".out")')        &  
     &     trim(adjustl(BXLFIT)),int(protons),int(neutrons) 
  
    open(unit=10,file=filedone)

    E = TotalE 
    if(rotcorr.ne.0) then
        Enocor = totalE - sum(rotcorrection)      &
        &                 - sum(COMcorrection(2,:)) & 
        &                 - sum(vibcorrection)
        Erot_vib = sum(rotcorrection)+ sum(vibcorrection)
    else
        Enocor = totalE - sum(COMcorrection(2,:)) 
        Erot_vib = 0.
    endif

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

    if(.not. symtransfo_needed) then
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
  
  subroutine massage_Bogoliubov()
    !---------------------------------------------------------------------------
    ! Subroutine to transform the Bogoliubov transformation found on file 
    ! towards one of the type demanded by the user. 
    ! 
    ! Currently only works for .wf that were created with
    !    blocktype   = 2 
    !    blocklowest = a combination of n+,n-,p+,p  (WITH NO REPEATS!)
    !
    ! and that break time-reversal symmetry.
    !    
    ! The routine checks 
    !   (i)   checks in what blocks excitations are         (array UNDO)
    !   (ii)  checks in what blocks excitations should be   (array DODO)
    !   (iii) flips quasiparticles of the lowest qp energy in every 
    !         block where an either 
    !          (a) excitation should be and isn't ; or
    !          (b) excitation shoud not be and is
    !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! 13/04/21: Some undesired behaviour because of this routine in the 
    !           the systematic calculations by G. Scamps. For this reason, this
    !           subroutine is currently "dead code", i.e. never called again.
    !---------------------------------------------------------------------------
    integer :: i, B, sb, N, N2, T
    integer, allocatable       :: undo(:), dodo(:)
    real(KIND=dp), allocatable :: temp(:)
    logical  :: flip

    !---------------------------------------------------------------------------    
    ! First, we see in what blocks we need to undo some qp excitations
    if(fileblocknumber .gt. 0) then
      allocate(undo(fileblocknumber))
      do i=1, fileblocknumber
        select case(fileblocklowest(i))
        case('n+')
          B = 1
        case('n-')
          B = 3
        case('p+')
          B = 5
        case('p-')
          B = 7
        case('n0', 'p0')
          call stp('Tantalus cannot handle FILEFROMBOGO=.true. with n0 or p0')
        end select
        
        undo(i) = B     
      enddo
    endif 
    !---------------------------------------------------------------------------
    ! Then, we check in with what the user asked
    if(blocknumber .gt. 0) then
      allocate(dodo(blocknumber))
      do i=1, blocknumber
        select case(blocklowest(i))
        case('n+')
          B = 1
        case('n-')
          B = 3
        case('p+')
          B = 5
        case('p-')
          B = 7
        case('n0', 'p0')
          call stp('Tantalus cannot handle FILEFROMBOGO=.true. with n0 or p0')
        end select
      
        dodo(i) = B     
      enddo 
    endif
    !---------------------------------------------------------------------------
    ! Then we flip some quasiparticle excitations
    sb = 0 
    do B=1,8,2
      N = file_HFB_blocks(B)   ; if(N.eq.0) cycle
      N2= file_HFB_blocks(B+1)
      T = N + N2

      flip = .false.
      if(allocated(undo)) then
        do i=1,fileblocknumber
           if(undo(i) .eq. B) flip = .true.
        enddo
      endif
      
      if(flip) then
        temp = Bogoliubov(sb+1:sb+2*T, sb+T+N+1)

        Bogoliubov(sb  +1:sb+  T, sb+T+N+1) = temp(T+1:2*T)
        Bogoliubov(sb+T+1:sb+2*T, sb+T+N+1) = temp(  1:  T)
      endif
      
      flip = .false.
      if(allocated(dodo)) then
        do i=1,blocknumber
           if(dodo(i) .eq. B) then
            flip = .true.
           endif
        enddo
      endif
      
      if(flip) then
        configmatrix(sb+T+N+1) = 0
        configmatrix(sb+    1) = 1
      endif

      sb = sb + 2 *T 
    enddo
    !---------------------------------------------------------------------------

  end subroutine massage_Bogoliubov

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
    type(moment), pointer :: Q20, Q22
   
    1 format("# N = ", i8, ' Z = ', i8, ' A = ', i8)
    2 format("# nwn = ", i9, ", nwp = ", i9)
    3 format("# (nx,ny,nz) = (", 3i5, "), dx = ", f10.8, ' fm')
    4 format("# Parameterisation    : ", a40)
    5 format("# Functional type     : ", a40)
    6 format("# Fermi energies      : ", 2f15.4)
    7 format("# Quadrupole   Q20,Q22: ", 2f15.4)
    8 format("# Quadrupole   B20,B22: ", 2f15.4)
    9 format("# Quadrupole    Q, gam: ", 2f15.4)

   10 format("# BI 1: ", 2i3)
   11 format("# BI 2: ", 99i4)
   12 format("# BI 3: ", 99a3)

    write(iochannel, fmt=1)  int(neutrons), int(protons),int(neutrons+protons)
    write(iochannel, fmt=2)  nwn, nwp
    write(iochannel, fmt=3)  nx, ny, nz, dx
    write(iochannel, fmt=4)  name_param
    write(iochannel, fmt=5)  func_name
    write(iochannel, fmt=6)  FermiEnergy
  
    Q20 =>FindMoment(2,0,.false.     )
    Q22 =>FindMoment(2,2,.false., Q20)    
    write(iochannel, fmt=7) sum(Q20%value), sum(Q22%value)
    write(iochannel, fmt=8)    Q20%beta(4), Q22%beta(4)
    write(iochannel, fmt=9)    Q(3), G(3)
    
    write(iochannel, fmt=10)  blocktype, blocknumber
    if(blocknumber .gt. 0) then
      write(iochannel, fmt=11) Blockindices
      write(iochannel, fmt=12) Blocklowest
    else
      write(iochannel, fmt=11) 
      write(iochannel, fmt=12)
    endif

    write(iochannel, fmt='(a1)') '#'

  end subroutine write_header

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
$TAUSCALAR write(1, fmt='(2es25.12)', advance='no') &
$TAUSCALAR &         R%D_Nm_Nm(mi,1), R%D_Nm_Nm(mi,2)
$TAUTENSOR write(1, fmt='(2es25.12)', advance='no') &
$TAUTENSOR &         R%D_N_N(mi,1,1,1) + R%D_N_N(mi,2,2,1) + R%D_N_N(mi,3,3,1),&
$TAUTENSOR &         R%D_N_N(mi,1,1,2) + R%D_N_N(mi,2,2,2) + R%D_N_N(mi,3,3,2)
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
    
    temp = F%F_I_I(:,1:2) - constraint_I_I

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
$TAUSCALAR write(1, fmt='(2es25.12)', advance='no') &
$TAUSCALAR &         F%F_Nm_Nm(mi,1), F%F_Nm_Nm(mi,2)
$TAUTENSOR write(1, fmt='(2es25.12)', advance='no') &
$TAUTENSOR &         F%F_N_N(mi,1,1,1) + F%F_N_N(mi,2,2,1) + F%F_N_N(mi,3,3,1),&
$TAUTENSOR &         F%F_N_N(mi,1,1,2) + F%F_N_N(mi,2,2,2) + F%F_N_N(mi,3,3,2)
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          ! The pairing fields FP_I_I
          write(1, fmt='(2es25.12)',advance='no') F%FP_I_I(mi,1), F%FP_I_I(mi,2)
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          ! The spin-orbit potential is the potential G_I_NS
          do mu=1,3
            do nu=1,3
              write(1, fmt='(2es25.12)', advance='no') &
              &            F%G_I_NS(mi,mu,nu,1), F%G_I_NS(mi,mu,nu,2)
            enddo
          enddo
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

  subroutine write_timeodd_densities(fname)
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
    
$NTR    real(KIND=dp), pointer           :: Sxn(:,:,:), Sxp(:,:,:)
$NTR    real(KIND=dp), pointer           :: Jxn(:,:,:), Jxp(:,:,:)
$NTR    real(KIND=dp), pointer           :: Syn(:,:,:), Syp(:,:,:)
$NTR    real(KIND=dp), pointer           :: Jyn(:,:,:), Jyp(:,:,:)
$NTR    real(KIND=dp), pointer           :: Szn(:,:,:), Szp(:,:,:)
$NTR    real(KIND=dp), pointer           :: Jzn(:,:,:), Jzp(:,:,:)
$NTR    real(KIND=dp), pointer           :: Txn(:,:,:), Txp(:,:,:)
$NTR    real(KIND=dp), pointer           :: Tyn(:,:,:), Typ(:,:,:)
$NTR    real(KIND=dp), pointer           :: Tzn(:,:,:), Tzp(:,:,:)

    real(KIND=dp), allocatable, target  :: totalangmom(:,:,:)
  
    character(len=*), intent(in)     :: fname
    integer                          :: io, i,j,k
$NTR integer                         :: it

    1 format('#  X[fm]   Y[fm]   Z[fm] ', &
    &        '   Sxn     Syn     Szn   ', &
    &        '   Sxp     Syp     Szp   ', &
    &        '   jxn     jyn     jzn   ', &
    &        '   jxp     jyp     jzp   ', &
    &        '   Jxn     Jyn     Jzn   ', &
    &        '   Jxp     Jyp     Jzp   ')
    2 format('#  0       1       2     ', &
    &        '   3       4       5     ', &
    &        '   6       7       8     ', &
    &        '   9      10      11     ', &
    &        '  12      13      14     ', &
    &        '  15      16      17     ', &
    &        '  18      19      20     ')

    open(1,file=fname, iostat=io)
    if(io.ne.0) then    
      print *, 'Something went wrong with writing a density to file.'
      print *, 'filename = ', fname
      call stp('')
    endif

$NTR    Sxn(1:nx,1:ny,1:nz)  => D_I_S(:,1,1) ; Sxp(1:nx,1:ny,1:nz)  => D_I_S(:,1,2)
$NTR    Syn(1:nx,1:ny,1:nz)  => D_I_S(:,2,1) ; Syp(1:nx,1:ny,1:nz)  => D_I_S(:,2,2)
$NTR    Szn(1:nx,1:ny,1:nz)  => D_I_S(:,3,1) ; Szp(1:nx,1:ny,1:nz)  => D_I_S(:,3,2)

$NTR    Jxn(1:nx,1:ny,1:nz)  => C_I_N(:,1,1) ; Jxp(1:nx,1:ny,1:nz)  => C_I_N(:,1,2)
$NTR    Jyn(1:nx,1:ny,1:nz)  => C_I_N(:,2,1) ; Jyp(1:nx,1:ny,1:nz)  => C_I_N(:,2,2)
$NTR    Jzn(1:nx,1:ny,1:nz)  => C_I_N(:,3,1) ; Jzp(1:nx,1:ny,1:nz)  => C_I_N(:,3,2)


    ! Calculate the total angular momentum density
$NTR    allocate(totalangmom(nx*ny*nz,3,2)) ; totalangmom = 0.0d0
$NTR    do it=1,2
$NTR      totalangmom(:,1,it) = 0.5 * D_I_S(:,1,it) ! spin part
$NTR      totalangmom(:,2,it) = 0.5 * D_I_S(:,2,it) ! spin part
$NTR      totalangmom(:,3,it) = 0.5 * D_I_S(:,3,it) ! spin part
$NTR      do i=1, nx*ny*nz
$NTR        ! X component : J_x ~ y j_z - z j_y
$NTR        TotalAngMom(i,1,it) = TotalAngMom(i,1,it) &
$NTR        & + meshgrid(i,2) * C_I_N(i,3,it) - meshgrid(i,3) * C_I_N(i,2,it)
$NTR
$NTR        ! Y component : J_y ~ z j_x - x j_z
$NTR        TotalAngMom(i,2,it) = TotalAngMom(i,2,it) &
$NTR        & + meshgrid(i,3) * C_I_N(i,1,it) - meshgrid(i,1) * C_I_N(i,3,it)
$NTR
$NTR        ! Z component : J_z ~ x j_y - y j_x
$NTR        TotalAngMom(i,3,it) = TotalAngMom(i,3,it) &
$NTR        & + meshgrid(i,1) * C_I_N(i,2,it) - meshgrid(i,2) * C_I_N(i,1,it)
$NTR      enddo
$NTR    enddo

   
$NTR    Txn(1:nx,1:ny,1:nz)  => TotalAngMom(:,1,1) 
$NTR    Txp(1:nx,1:ny,1:nz)  => TotalAngMom(:,1,2)
$NTR    Tyn(1:nx,1:ny,1:nz)  => TotalAngMom(:,2,1) 
$NTR    Typ(1:nx,1:ny,1:nz)  => TotalAngMom(:,2,2)
$NTR    Tzn(1:nx,1:ny,1:nz)  => TotalAngMom(:,3,1) 
$NTR    Tzp(1:nx,1:ny,1:nz)  => TotalAngMom(:,3,2)

    call write_header(1)
    write(1, fmt=2) 
    write(1, fmt=1) 
    do k=1,nz
      do j=1,ny
        do i=1,nx

$NTR          write(1, fmt='(3f8.3, 18es25.12E3)') meshx(i), meshx(j), meshz(k),   &
$NTR          &                                Sxn(i,j,k), Syn(i,j,k), Szn(i,j,k), & 
$NTR          &                                Sxp(i,j,k), Syp(i,j,k), Szp(i,j,k), & 
$NTR          &                                Jxn(i,j,k), Jyn(i,j,k), Jzn(i,j,k), & 
$NTR          &                                Jxp(i,j,k), Jyp(i,j,k), Jzp(i,j,k), &
$NTR          &                                Txn(i,j,k), Tyn(i,j,k), Tzn(i,j,k), &
$NTR          &                                Txp(i,j,k), Typ(i,j,k), Tzp(i,j,k)

$TR          write(1, fmt='(3f8.3, 18es25.12E3)') meshx(i), meshx(j), meshz(k),   &
$TR          &                                0.0d0,0.0d0,0.0d0, &
$TR          &                                0.0d0,0.0d0,0.0d0, &
$TR          &                                0.0d0,0.0d0,0.0d0, &
$TR          &                                0.0d0,0.0d0,0.0d0, &
$TR          &                                0.0d0,0.0d0,0.0d0, &
$TR          &                                0.0d0,0.0d0,0.0d0

        enddo
      enddo
    enddo

    deallocate(TotalAngMom)
    close(1)

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
    do it=1,3
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
        write(1, fmt=7) l,m,collective_inertia(k,1:N_inertia ,it)
      enddo    
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
    !
    !---------------------------------------------------------------------------
    type(moment), pointer        :: quadrupole
    character(len=*), intent(in) :: combi
    integer, allocatable :: indices(:)
    integer              :: i,ii, p1, p2,jj
    real(KIND=dp)        :: A, mstate1, mstate2
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

    quadrupole => FindMoment(2,0,.false.)
    A = protons + neutrons

    ! Note: items marked with (*) are written as zero and, to the best of
    ! my (=W.R.) knowledge, not used by the level density code.  
    !                        IZ          IA    BETA     B4
    write(unit=6, fmt=3)  int(protons),int(A),quadrupole%beta(4), 0.0,  &
    !                      HGN    HFGP  HFDN             HFDP
    &                      0.0,   0.0,  average_gap(2,1),average_gap(2,2), & 
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
    
end module IO
