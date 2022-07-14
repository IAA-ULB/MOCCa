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
  !-----------------------------------------------------------------------------
  integer, parameter  :: version_number = 5
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
  integer       :: filenx, fileny, filenz, filenwn, filenwp, filepairing
  integer       :: filenwt, fileneutrons, fileprotons, fileblocks(8)
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
    !---------------------------------------------------------------------------

    use geninfo,       only : ReadGenInfo
    use evolution,     only : ReadEvolution
    use wavefunctions, only : ReadWFdata
    use scfiteration,  only : readscfiteration
    use moments,       only : readmomentdata
    use functional,    only : readfunctional
    use pairing,       only : initpairing
    use fission_moi,   only : N_inertia, read_inertia
  
    implicit none

    ! These inputs control where the code will look for its input. Leaving them 
    ! empty will have the code rely on STDIN for input.
    integer(dp), intent(in), optional   :: file_number   
    character(26), intent(in), optional :: input_file 

    logical :: exists

    NameList /IO/ InputFileName,OutputFileName, BXLFIT, COMBI, denfile,potfile,& 
    &           sphffile, spcanfile,checkpointiter, AllowTransform, extraspwfs,&
    &           Counter, run, tofile, blockfile, inertfile, N_inertia
    
    if(present(file_number)) then
      inquire(file=input_file, exist=exists)
      if(.not. exists) then
        print *, 'Specified input file does not exist!'
        stop
      endif
      open(unit=file_number, file=input_file) 
    endif

    call ReadGenInfo(file_number)
    call readfunctional(file_number)
    call initpairing(file_number)
    call ReadEvolution(file_number)
    call ReadSCFIteration(file_number)
    call ReadWFdata(file_number)
    
    if(present(file_number)) then
      read (unit=file_number, nml=IO)
    else
      read (unit=*, nml=IO)
    endif
    
    if(N_inertia .gt. 0) then
      call read_inertia(file_number)
    else
      if(N_inertia .lt. 0) then
        print *, 'Wrong value for N_inertia.'
        stop
      endif
    endif
    
    call readmomentdata(file_number)
    call readcranking(file_number)

    if(present(file_number)) then
      close(unit=file_number)
    endif

  end subroutine ReadInput

  subroutine PrintInput(file_number, input_file)
  !-----------------------------------------------------------------------------
  ! This subroutine prints all relevant information of the input, both from the
  ! user and from the wavefunction file.
  !-----------------------------------------------------------------------------
   
    use wavefunctions
    use evolution
    use scfiteration

    integer*8, intent(in), optional     :: file_number
    character(11), intent(in), optional :: input_file 
   
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
   12 format ( ' Convergence required', / &
    &          '  Energy convergence           < ', es8.1, / & 
    &          '  Multipole moment convergence < ', es8.1, / &
    &          '  Dispersion convergence       < ', es8.1, / &
    &          '  S.p. gradient convergenc     < ', es8.1, / &
    &          '  Fermi energy convergence     < ', es8.1, / &
    &          '  Angular momentum convergence < ', es8.1)
   13 format ( ' Inverse temperature Beta = ', f14.9)

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
    
  end subroutine PrintInput
  
  subroutine Readwavefunction()
    !---------------------------------------------------------------------------
    ! High-level routine to determine the starting point of a calculation. 
    !
    ! There are two starting options
    !
    ! 1) Initialize in an EV8-style box with Nilsson orbitals
    ! 2) Read a set of spwfs from file 
    ! 
    ! Which option is chosen based on the InputFileName keyword: it is is
    !  'INIT' (case insensitive) then the code performs option 1). Otherwise
    ! it will attempt to read said file. 
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !
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
    ! flag to .true. This behavior is coded like that as a general safeguard.
    !---------------------------------------------------------------------------
    integer :: i
    !---------------------------------------------------------------------------
    ! Input options 
    if(trim(to_upper(inputfilename)).eq.'INIT') then
      ! Option 1) generate starting point with Nilsson wavefunctions.
      call iniwavefunctions($ININX, $ININY, $ININZ, $ININWN, $ININWP)
      guessgaps         = .true.
      fileblocks        = HFBlocks

      if( SYM_CODE .ne. "0 1 001 000 10 000 010 111" ) then
        ! Initialisation with nil8 wavefunctions is always EV8-style
        ! Thus we signal that a symmetry transformation is needed
        symtransfo_needed = .true.
        if( TRANS_CODE .ne. "0 1 001 000 10 000 010 111") then
          print *, "---------------------------------------------------"
          print *, "| Calculations cannot be initialized from scratch |"
          print *, "| for this particular symmetry option.            |"
          print *, "---------------------------------------------------"
          stop        
        endif
      endif
      filenx = $ININX ; fileny = $ININY ; filenz = $ININZ ; filedx = dx
    else
      ! Option 2) start from a previous calculation.
      call ReadTantalus(12, inputfilename)
      ! No need to guess gaps by default (unless the user asked for it)
    endif
    !---------------------------------------------------------------------------  
    ! Transformation options
    if(allowtransform ) then
      if(  symtransfo_needed ) then 
          ! Option a): break a symmetry and transform the spwfs appropriately
          call  Transformspwfs( HFPsi, fileblocks, filenx, fileny, filenz)
          call  GramSchmidt  
      else
          ! Option b): add points and/or add spwfs
          call  TransformInput(filenx,fileny,filenz,filenwn,filenwp,filedx,    & 
          &                               fileblocks,file_HFB_blocks,extraspwfs)
          call  GramSchmidt  
          ! The added spwfs are added somewhat randomly, hence we add an extra
          ! orthonormalisation in the mix.
      endif
    else  
      ! We still need to set this particular information
      HFblocks = fileblocks
      if(symtransfo_needed) then
        print *, 'Symmetry transformation needed, but not allowed by user.'
        stop
      endif
    endif
    
    ! Failsafe for the HF transformation
    if(.not.allocated(HFTransfo)) then
        allocate(HFTransfo(nwt,nwt)) 
        HFtransfo = 0.0d0
        do i=1,nwt
            HFtransfo(i,i) = 1.0d0
        enddo
    endif
    
    call set_spwf_symmetries(sx, sy, sz, HFblocks)
    call update_spwf_symmetries()
    !---------------------------------------------------------------------------
    if(guessgaps) then
      ! Guess some pairing gaps if asked for (always if starting from INIT)
      call initializeGaps()
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
    
  end subroutine ReadWaveFunction

  subroutine ReadTantalus(chan, ifn)
    !---------------------------------------------------------------------------
    ! Reading all information from a previous Tantalus run. 
    ! Heavily based on the MOCCa input routine.
    !
    ! Also performs a few sanity checks. 
    ! Currently:
    !   
    !   *) equality of (nx,ny,nz) between data and file
    !   *) equality of (nwn,nwp) between data and file
    !   *) the symmetry encoding matches either SYM_CODE or TRANS_CODE
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
    integer                      :: io,i
    logical                      :: exists
    real(KIND=dp)                :: Omega_file(3)
    real(KIND=dp), allocatable   :: filegaps(:,:), temp(:,:)
    logical                      :: filediagsphamil 
    logical                      :: check_x, check_y, check_z
    logical                      :: check_nwn, check_nwp
    
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
    !---------------------------------------------------------------------------
    ! First check if the file exists.
    inquire(file=inputfilename, exist=exists)
    if(.not.exists) then
      print *, 'Input file specified does not exist!'
      stop
    endif
    !---------------------------------------------------------------------------
    open (chan,form='unformatted',file=ifn)
    
    read(chan, iostat=io) file_version
    if(file_version .gt. version_number) then
      print *, 'Unsupported version number of the .wf file.'
      print *, 'Maximum current version: ', version_number
      stop
    endif

    ! Convergence information                                  (NOT IMPLEMENTED)
    read(chan,iostat=io) 
    !Parameters of the mesh
    read(Chan,iostat=io) filenx,fileny,filenz, filedx

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
        stop  
      endif
    endif

    !Number of protons and neutrons
    read(Chan,iostat=io) fileneutrons, fileprotons
    ! HFBLocks information 
    read(Chan,iostat=io) filenwn, filenwp, fileblocks
    filenwt = filenwn + filenwp
    ! Wavefunctions
    !- - - - - - - - - - - - - - - -
    ! First allocate the needed space
    allocate(HFPsi(filenx*fileny*filenz,4, filenwn+filenwp))
    allocate(spenergies(filenwn+filenwp))
    allocate(dispersions(filenwn+filenwp))

    if (allocated(rho_can)) then 
      deallocate(rho_can)        
    end if                       

    allocate(rho_can(filenwn + filenwp))
    
    read(chan,iostat=io) spenergies, dispersions    
    
    if(file_version .ge. 3) then
      read(chan,iostat=io) filediagsphamil
      allocate(HFtransfo(filenwt,filenwt)) 
      read(chan,iostat=io) HFtransfo
    endif
    
    read(chan,iostat=io) HFPsi    
    ! Name of the force and functional
    read(chan, iostat=io) ini_name_param, func_name_check
    ! Single-particle hamiltonian
    if(file_version.ge.4) then
      allocate(current_sph(filenwt,filenwt))
      read(chan, iostat=io) current_sph
    endif
    !---------------------------------------------------------------------------
    ! Pairing information                                      
    read(chan, iostat=io) filepairing
    ! Write the occupation factors in all cases
    read(chan, iostat=io) rho_can

    select case (filepairing)
    case(0)
        ! HF: nothing to read
    case(1)
        ! BCS: read the gaps
        allocate(filegaps(filenwn+filenwp,1))
        read(chan, iostat=io) FermiEnergy       ! Lambda
        read(chan, iostat=io) filegaps
        
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
    case(2)
        ! HFB
        if (allocated(filegaps)) then 
          deallocate(filegaps)        
        end if                       

        allocate(filegaps(filenwt, filenwt)) 
        allocate(kappa_pairing(filenwt, filenwt)) 
        allocate(rho_pairing(filenwt, filenwt)) 
        
        if(file_version .gt. 3 ) then
          read(chan, iostat=io) fileblocktype, fileblocknumber
          read(chan, iostat=io) file_HFB_blocks
          select case(fileblocktype) 
            case(0)
              read(chan, iostat = io)
            case(1,3,5)
              allocate(fileblockindices(fileblocknumber))
              read(chan, iostat = io) fileblockindices
            case(2,4,6)
              allocate(fileblocklowest(fileblocknumber))
              read(chan, iostat = io) fileblocklowest
          end select
        endif

        read(chan, iostat=io) FermiEnergy       ! Lambda
        read(chan, iostat=io) rho_pairing
        read(chan, iostat=io) kappa_pairing     ! kappa
        read(chan, iostat=io) ! Canonical transformation

        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        ! Do some gymnastics to read the gaps
        allocate(temp(filenwt, filenwt))
        io = 0
        read(chan, iostat=io) temp ! HFBGaps
        filegaps = temp(1:filenwt, 1:filenwt)

        !-----------------------------------------------------------------------        
        ! We no longer do these gymnastics, which were only necessary to support
        ! old .wf files, none of which still exist (I think/hope).
!        if(io.ne.0) then
!          rewind(chan)
!          
!          rewindc=17
!          if(file_version.lt.3) rewindc=15
!          do c=1,rewindc
!                read(chan, iostat=io)
!          enddo
!          deallocate(temp) ; allocate(temp(filenwt, filenwt))
!          read(chan, iostat=io) temp
!        endif    
        !-----------------------------------------------------------------------        

        if (io.ne.0) then
          print *, 'ERROR in reading the gaps from file.'
          stop
        endif
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        ! For late-enough versions, we can  read also the full Bogoliubov 
        ! transformation and the associated configuration matrix.
        if(file_version .ge. 3) then
        
          readHFBinfofile = .true.
          allocate(Bogoliubov(2*filenwt, 2*filenwt)) 
          allocate(configmatrix(2*filenwt)) 

          ! We simply read these arrays here. If a transformation is needed,
          ! we will deal with it elsewhere.
          read(chan, iostat=io) Bogoliubov
          if (io.ne.0) then
            print *, 'ERROR in reading the Bogoliubov transformation from file.'
            stop
          endif
          read(chan, iostat=io) configmatrix
          if (io.ne.0) then
            print *, 'ERROR in reading the configuration matrix from file.'
            stop
          endif
        endif

        select case(pairingtype)
        case(0)
          ! Do nothing
        case(1)
          ! Use the diagonal HFBgaps for the BCSgaps
          ! If a transformation is needed we will deal with it elsewhere
          allocate(BCSgaps(filenwt)) 
          do i=1, filenwt
            BCSgaps(i) = filegaps(i,i)
          enddo
        case(2)
          ! Simply copy the gaps for now
          if (allocated(HFBGaps)) then   
            deallocate(HFBGaps)        
          end if                       
          allocate(HFBGaps(filenwt, filenwt)) 
          HFBGaps = filegaps(1:filenwt, 1:filenwt)  
        end select
    case DEFAULT
      print *, 'Something is seriously wrong with the .wf file.'
      stop
    end select   
    ! Cranking information       
    if(file_version .gt. 4 ) then                              
      read(chan, iostat=io) omega_file(1:3)
    else
      ! File-versions < 4 do not have this line
      read(chan, iostat=io) 
      omega_file = 0.0d0
    endif
    if(io.ne.0) then
        print *, 'ERROR in reading cranking line of the wf file.'
        stop
    endif
    
    if(continueCrank) then
      ! Using the cranking frequencies read from file
      omega = omega_file
    endif
    
    !---------------------------------------------------------------------------
    ! Potentials                                               
    call readpotentials(chan, filenx,fileny,filenz, symtransfo_needed)
    !---------------------------------------------------------------------------
    ! Multipole moment information                             
    io = 0
    do while(io.eq.0)
      call ReadMoment(chan,io)
    enddo
    
    ! End of reading
    close(chan)
    !---------------------------------------------------------------------------
    if(.not.  AllowTransform) then
      !-------------------------------------------------------------------------
      ! Sanity checks if transformation is not allowed
      if((filenx.ne.nx).or. (fileny.ne.ny) .or. (filenz.ne.nz)) then
          print 1, filenx, fileny, filenz, nx,ny,nz
          stop
      endif
      if(filenwn.ne.nwn .or. filenwp.ne.nwp) then
          print 2, filenwn, filenwp, nwn, nwp
          stop
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
          print *, "Please don't combine symmetry transformations and mesh modifications."
          stop
         endif  
         if(check_nwn .or. check_nwp ) then 
          print *, "Please don't combine symmetry transformations and adding wavefunctions."
          stop
         endif  
      endif
    endif

  end subroutine ReadTantalus

  subroutine WriteTantalus(chan, ofn)
    !---------------------------------------------------------------------------
    ! Subroutine that dumps all information to a .wf file for future runs.
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
    integer                      :: io
    type(moment), pointer        :: mom

    open (chan,form='unformatted',file=ofn)

    write(chan, iostat=io) version_number
    ! Convergence information                                  (NOT IMPLEMENTED)
    write(chan,iostat=io) 
    !Parameters of the mesh
    write(Chan,iostat=io) nx,ny,nz, dx
    ! Symmetry information                                   
    write(Chan,iostat=io) SYM_CODE
    !Number of protons and neutrons
    write(Chan,iostat=io) neutrons,protons
    ! HFBLocks information 
    write(Chan,iostat=io) nwn, nwp, hfblocks
    ! Wavefunctions  
    write(chan,iostat=io) spenergies, dispersions
    ! information on the HF transformation
    write(chan, iostat=io) diagsphamil
    write(chan, iostat=io) HFtransfo
    
    write(chan,iostat=io) HFPsi                              
    ! Name of the force.
    write(chan, iostat=io) name_param, func_name
    ! Single-particle hamiltonian
    write(chan, iostat=io) current_sph
    !---------------------------------------------------------------------------
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
          write(chan, iostat=io) HFBlocks
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
    !---------------------------------------------------------------------------
    ! Potentials on file
    call writepotentials(chan)
    !---------------------------------------------------------------------------
    ! Multipole moment information                             
    mom => root
    do while(associated(mom%next))
      mom => mom%next
      call Writemoment(mom,chan)
    enddo
    close(chan)

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
      call write_densities(DENFILE)
    endif
    if(TOFILE .ne. '') then
$TR   print *, 'Time-odd densities do not figure in a calculation that assumes time-reversal.'
$TR   stop
      call write_timeodd_densities(TOFILE)
    endif
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Write the relevant potentials to a file for postprocessing
    if(POTFILE .ne. '') then
      call write_potentials(POTFILE)
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
      b32 = Q32%beta(4)
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
      print *, 'Subroutine check_blocking_structure cannot deal (yet) with'
      print *, 'blocktypes that are not 0/2/4.'
      stop
    endif    
    !---------------------------------------------------------------------------
    ! # 2 :  blocklowest on file == blocklowest input by the user 
    !        modulo permutations 
    if(allocated(fileblocklowest) .and. (.not. allocated(blocklowest))) then
      print *, 'Blocklowest not allocated, while fileblocklowest is.'
      stop      
    endif

    if(allocated(blocklowest) .and. (.not. allocated(fileblocklowest))) then
      print *, 'Blocklowest allocated, while fileblocklowest is not.'
      stop      
    endif
    
    if(size(fileblocklowest).ne.size(blocklowest)) then
      print *, ' Size of blocklowest on file:  ', size(fileblocklowest)
      print *, ' Size of blocklowest in input: ', size(blocklowest)
      stop
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
      stop
    endif
    !---------------------------------------------------------------------------
    ! # 3: Check if the blocking structure on file actually matches the 
    !      structure asked for
    
    !---------------------------------------------------------------------------
    ! This method has turned out to NOT be a reliable indicator.
!    blocked_blocks =  figure_out_blocking_structure_agnostic(                  &
!    &                             current_sph, HFBgaps, FermiEnergy, Bogoliubov)
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
          stop
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
          print *, 'Tantalus cannot handle FILEFROMBOGO=.true. with n0 or p0'
          stop
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
          print *, 'Tantalus cannot handle FILEFROMBOGO=.true. with n0 or p0'
          stop
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
   
    1 format("# N = ", i3, ' Z = ', i3, ' A = ', i3)
    2 format("# nwn = ", i3, ", nwp = ", i3)
    3 format("# (nx,ny,nz) = (", 3i3, "), dx = ", f8.6, ' fm')
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

  subroutine write_densities(fname)
    !---------------------------------------------------------------------------
    ! Write the following densities to a file named "fname"
    !    rho(neutron), rho(proton), rho(charge)
    !---------------------------------------------------------------------------
    ! The file contains a header written by the subroutine write_header,
    ! supplemented by
    !
    !     #   X[fm] Y[fm] Z[fm] rho_n[fm^{-3}] rho_p[fm^{-3}] rho_c[fm^{-3}]
    ! 
    ! where the # are included so that Numpy (or other plotting tools) can 
    ! ignore these lines when naively plotting stuff. Note that the fourth
    ! line is currently empty, but is reserved for future additions concerning
    ! symmetry options of the current run.
    !
    ! The format of the body of said file is
    ! 
    !        x , y , z,  rho_n, rho_p, rho_c
    !
    ! where the first three numbers are the Cartesian coordinates (units of fm)
    ! and the densities are all in units of fm^{-3}. 
    ! The points are written down in column-major order ('Fortran order'), 
    ! which might not be how your favorite plotting tool prefers it.
    !---------------------------------------------------------------------------
    ! Note that the densities are written "as they are" to file, i.e. only in
    ! part of the box that is actually represented numerically. It is up to
    ! postprocessing to actually construct the densities in the entire box.
    !---------------------------------------------------------------------------
    real(KIND=dp), pointer           :: dn(:,:,:), dp(:,:,:)
    character(len=*), intent(in)     :: fname
    integer                          :: io, i,j,k

    1 format('#  X[fm]   Y[fm]   Z[fm]       rho_n[fm^{-3}]           rho_p[fm^{-3}]           rho_c[fm^{-3}]')
    open(1,file=fname, iostat=io)
    if(io.ne.0) then    
      print *, 'Something went wrong with writing a density to file.'
      print *, 'filename = ', fname
      stop
    endif

    dn(1:nx,1:ny,1:nz)  => D_I_I(:,1)
    dp(1:nx,1:ny,1:nz)  => D_I_I(:,2)

    call write_header(1)
    write(1, fmt=1) 
    do k=1,nz
      do j=1,ny
        do i=1,nx
          write(1, fmt='(3f8.3, 3es25.12E3)') meshx(i), meshx(j), meshz(k),      &
          &                           dn(i,j,k), dp(i,j,k), chargedensity(i,j,k) 
        enddo
      enddo
    enddo

    close(1)
  end subroutine write_densities

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
    
    real(KIND=dp), pointer           :: Sxn(:,:,:), Sxp(:,:,:)
    real(KIND=dp), pointer           :: Jxn(:,:,:), Jxp(:,:,:)
    real(KIND=dp), pointer           :: Syn(:,:,:), Syp(:,:,:)
    real(KIND=dp), pointer           :: Jyn(:,:,:), Jyp(:,:,:)
    real(KIND=dp), pointer           :: Szn(:,:,:), Szp(:,:,:)
    real(KIND=dp), pointer           :: Jzn(:,:,:), Jzp(:,:,:)
    character(len=*), intent(in)     :: fname
    integer                          :: io, i,j,k

    1 format('#  X[fm]   Y[fm]   Z[fm]   ')
    open(1,file=fname, iostat=io)
    if(io.ne.0) then    
      print *, 'Something went wrong with writing a density to file.'
      print *, 'filename = ', fname
      stop
    endif

$NTR    Sxn(1:nx,1:ny,1:nz)  => D_I_S(:,1,1) ; Sxp(1:nx,1:ny,1:nz)  => D_I_S(:,1,2)
$NTR    Syn(1:nx,1:ny,1:nz)  => D_I_S(:,2,1) ; Syp(1:nx,1:ny,1:nz)  => D_I_S(:,2,2)
$NTR    Szn(1:nx,1:ny,1:nz)  => D_I_S(:,3,1) ; Szp(1:nx,1:ny,1:nz)  => D_I_S(:,3,2)

$NTR    Jxn(1:nx,1:ny,1:nz)  => C_I_N(:,1,1) ; Jxp(1:nx,1:ny,1:nz)  => C_I_N(:,1,2)
$NTR    Jyn(1:nx,1:ny,1:nz)  => C_I_N(:,2,1) ; Jyp(1:nx,1:ny,1:nz)  => C_I_N(:,2,2)
$NTR    Jzn(1:nx,1:ny,1:nz)  => C_I_N(:,3,1) ; Jzp(1:nx,1:ny,1:nz)  => C_I_N(:,3,2)

    call write_header(1)
    write(1, fmt=1) 
    do k=1,nz
      do j=1,ny
        do i=1,nx
          write(1, fmt='(3f8.3, 12es25.12E3)') meshx(i), meshx(j), meshz(k),   &
          &                                Sxn(i,j,k), Syn(i,j,k), Szn(i,j,k), & 
          &                                Sxp(i,j,k), Syp(i,j,k), Szp(i,j,k), & 
          &                                Jxn(i,j,k), Jyn(i,j,k), Jzn(i,j,k), & 
          &                                Jxp(i,j,k), Jyp(i,j,k), Jzp(i,j,k)

        enddo
      enddo
    enddo

    close(1)

  end subroutine write_timeodd_densities

  subroutine write_potentials(fname)
    !---------------------------------------------------------------------------
    ! Write the following potentials to a file named "fname"
    !  F_I_I(n/p),  F_c(n/p),  V_so(n/p),   V_pair(n/p)
    !  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  central       coulomb   spin-orbit   pairing  
    !
    ! Remarks:
    !   *) no definition yet of V_so or V_pair
    !   *) F_c is the potential of the DIRECT Coulomb energy, directly obtained
    !      from the charge density. It is in general NOT this potential that 
    !      the protons (and neutrons) feel. 
    !---------------------------------------------------------------------------
    !
    ! The file contains a header written by the subroutine write_header,
    ! supplemented by
    ! #   X[fm] Y[fm] Z[fm] V_nuc(n) V_nuc(p) V_c(n) V_c(p) 
    !                                             V_so(n) V_so(p) V_p(n) V_p(p)
    ! 
    ! where the # are included so that Numpy (or other plotting tools) can 
    ! ignore these lines when naively plotting stuff. Note that the fourth
    ! line is currently empty, but is reserved for future additions concerning
    ! symmetry options of the current run.
    !
    ! The format of the body of said file is
    ! 
    !     x,y,z,F_I_I(n),F_I_I(p), F_c,V_so(n), V_so(p),V_pair(n),V_pair(p)
    !
    ! where the first three numbers are the Cartesian coordinates (units of fm).
    ! The points are written down in column-major order ('Fortran order'), 
    ! which might not be how your favorite plotting tool prefers it.
    !---------------------------------------------------------------------------
    ! IMPORTANT:
    !  While this routine now claims to write V_so and V_pair to file, right now
    !  it just writes zeros in those columns.
    !---------------------------------------------------------------------------
    character(len=*), intent(in) :: fname
    real(KIND=dp), pointer       :: Vnucp(:,:,:), Vnucn(:,:,:)
    real(KIND=dp), allocatable   :: Couln(:,:,:), Coulp(:,:,:)

    real(KIND=dp), allocatable, target   :: temp(:,:)
    integer                              :: io, i,j,k

    1 format('#  X[fm]   Y[fm]   Z[fm]  V_nuc(n) V_nuc(p) V_c(n) V_c(p) V_so(n) V_so(p) V_p(n) V_p(p)')

    open(1,file=fname, iostat=io)
    if(io.ne.0) then    
      print *, 'Something went wrong with writing a potential to file.'
      print *, 'filename = ', fname
      stop
    endif

    call write_header(1)
    write(1, fmt=1) 
  
    ! The central nuclear potential is the field associated with D_I_I, but it
    ! should not include the constraints, nor the contribution of the 
    ! Coulomb potential
    allocate(temp(nx*ny*nz,2), couln(nx,ny,nz), coulp(nx,ny,nz))
    temp = F_I_I(:,1:2) - constraint_I_I

    Vnucn(1:nx,1:ny,1:nz)  => temp(:,1)
    Vnucp(1:nx,1:ny,1:nz)  => temp(:,2)

    ! Subtracting the coulomb potential depends on our treatment of the 
    ! proton and neutron finite size effect
    if((all(protonsize.eq.0.0) .and. all(neutronsize.eq.0.0)) .or.         &
    &                             (.not. nucleonsize_selfconsistent)) then
      ! No finite size effect for either protons or neutrons
      Vnucp = Vnucp - CoulombPotential  - ExchangePotential

      Couln = 0.0d0
      Coulp = CoulombPotential
    else
      ! Finite size effects taken into account
      Vnucn = Vnucn - FoldedCoul(:,:,:,1)       - FoldedExchange(:,:,:,1)
      Vnucp = Vnucp - FoldedCoul(:,:,:,2)       - FoldedExchange(:,:,:,2)

      Couln = FoldedCoul(:,:,:,1)
      Coulp = FoldedCoul(:,:,:,2)
    endif

    do k=1,nz
      do j=1,ny
        do i=1,nx
          write(1, fmt='(3f8.3, 8es25.12)') meshx(i), meshy(j), meshz(k),      &
          &          Vnucn(i,j,k), Vnucp(i,j,k), Couln(i,j,k), Coulp(i,j,k),   &
          &          0.0, 0.0 ,0.0, 0.0 
        enddo
      enddo
    enddo

    deallocate(couln, coulp)
    close(1)
  end subroutine write_potentials

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
      stop
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

      write(1, fmt=1) wave, -1, p, 2*rho_HF(wave), spenergies(wave),   & 
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

      write(1, fmt=1) wave, +1, p, 2*rho_HF(wave), spenergies(wave), & 
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
      stop
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

    if(.not.allocated(blocked_sps)) then
      print *, 'Cannot write single-particle wavefunctions to file.'
      return
    endif

    open(1,file=fname, iostat=io)
    if(io.ne.0) then    
      print *, 'Something went wrong with writing blocked states to file.'
      print *, 'filename = ', fname
      stop
    endif
    
    call write_header(1)
    write(1, fmt=1, advance='no')
    do k=1, blocknumber
      write(1, fmt=2, advance ='no') k
    enddo
    write(1, fmt=*)
    
    allocate(psis(nx,ny,nz,blocknumber)) ; psis = 0
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
      stop
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
    ! Write an extra file for input of the combinatorial level density code.
    !
    ! ATTENTION: this output assumes an axial nucleus with a symmetry axis 
    !            along the z-axis. If the single-particle states are not  
    !            (at least approximately) eigenstates of J_z, then this output
    !            will effectively be nonsense.
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

          if(ii .lt. (HFBlocks(1))) p1 =  0
          if(ii .gt. (HFBlocks(1))) p1 =  1

          if(jj .lt. (HFBlocks(1))) p2 =  0
          if(jj .gt. (HFBlocks(1))) p2 =  1

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

          if(ii .lt. sum(HFBlocks(1:3))) p1 =  0
          if(ii .gt. sum(HFBlocks(1:3))) p1 =  1

          if(jj .lt. sum(HFBlocks(1:3))) p2 =  0
          if(jj .gt. sum(HFBlocks(1:3))) p2 =  1

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
    &                      Belyaev(1,1), Belyaev(1,2),       &
    !                      HFINY /     ,    HFIPY,      
    &                      Belyaev(2,1), Belyaev(2,2),       &
    !                      HFINZ /     ,    HFIPZ,      
    &                      Belyaev(3,1), Belyaev(3,2),       &
    !                      HFJ2               HFE1  , HE2
    &                       J2(2,3),  totalE, 0.0

    close(unit=6)
  end subroutine combi_output
    
end module IO
