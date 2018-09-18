module IO
 !==============================================================================
 !  #######   ##   #    # #####   ##   #      #    #  ####
 !     #     #  #  ##   #   #    #  #  #      #    # #
 !     #    #    # # #  #   #   #    # #      #    #  ####
 !     #    ###### #  # #   #   ###### #      #    #      #
 !     #    #    # #   ##   #   #    # #      #    # #    #
 !     #    #    # #    #   #   #    # ######  ####   ####
 !
 !  Copyright W. Ryssens & M. Bender
 !
 !============================================================================== 
 ! Module governing the in- and output of Tantalus. 
 !
 !==============================================================================

use geninfo
use wavefunctions
use pairing

implicit none

  !-----------------------------------------------------------------------------
  ! Filenames for in- and output of the code with respect to spwfs.
  character(len=100) :: inputfilename, outputfilename

contains

  subroutine ReadInput
    !---------------------------------------------------------------------------
    ! Subroutine to read all the data from STDIN.
    !
    !---------------------------------------------------------------------------

    use geninfo,       only : ReadGenInfo
    use evolution,     only : ReadEvolution
    use wavefunctions, only : ReadWFdata
    use scfiteration,  only : readscfiteration
    use moments,       only : readmomentdata
    use functional,    only : readfunctional
    use pairing,       only : initpairing
  
    implicit none

    NameList /IO/ InputFileName,OutputFileName
    
    call ReadGenInfo
    call readfunctional
    call initpairing
    call ReadEvolution
    call ReadSCFIteration
    call ReadWFdata
    
    read (unit=*, nml=IO)
    
    call readmomentdata
    
  end subroutine ReadInput

  subroutine PrintInput
  !-----------------------------------------------------------------------------
  ! This subroutine prints all relevant information of the input, both from the
  ! user and from the wavefunction file.
  !-----------------------------------------------------------------------------
   
    use wavefunctions
    use evolution
    use scfiteration
   
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
   10 format ( ' IO information', / &
    &          '  inputfilename  =', a20, / &
    &          '  outputfilename =', a20)
   11 format ( ' Convergence required', / &
    &          '  Energy convergence           < ', e8.1, / & 
    &          '  Multipole moment convergence < ', e8.1, / &
    &          '  Dispersion convergence       < ', e8.1 )

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
    print 10, inputfilename, outputfilename
    print 11, energy_prec, moment_prec, disp_prec

    call printevolution
    call printscfiteration
    call printfunctional  
    
  end subroutine PrintInput
  
  subroutine Readwavefunction
    !---------------------------------------------------------------------------
    ! Read information from a finished Tantalus calculation.
    ! This is the part that should decide on how to read from different inputs.
    !---------------------------------------------------------------------------
    
    if(trim(to_upper(inputfilename)).eq.'INIT') then  
      ! Generate starting point with Nilsson wavefunctions.
      call iniwavefunctions()
      ! Guess some pairing gaps
      call GuessGaps()
    else
      ! If not, start from a previous calculation.
      call ReadTantalus(12, inputfilename)
      ! No need to guess gaps, they should read from file. 
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
    !
    !---------------------------------------------------------------------------
    ! 
    ! Things read from file. (Not yet implemented ones are indicated by *)
    !
    ! Version
    ! Convergence information: E, dE                         (*)
    ! nx,ny,nz,dx,dt                    
    ! Symmetry information                                   (*)
    ! neutrons,protons
    ! Number of wavefunctions in every block
    ! (nwt) Wavefunctions                                    
    ! Forcename                                              
    ! Pairing information                                    
    !    - Pairingtype
    !    - Rho_can = occupation factors 
    !      (HF)  nothing
    !      (BCS) BCSGaps
    !      (HFB) rho_pairing    
    !        |   kappa_pairing
    !        |   can_transfo
    !        |   HFBgaps            
    ! Densities                                              (*)
    ! CrankingInfo                                           (*)                      
    ! Moments                                                (*)
    !
    !---------------------------------------------------------------------------
    
    use functional
    
    integer, intent(in)          :: chan
    character(len=*), intent(in) :: ifn
    character(len=20)            :: func_name_check
    integer                      :: io, version
    logical                      :: exists
    
    integer       :: filenx, fileny, filenz, filenwn, filenwp,i, filepairing
    integer       :: filenwt 
    real(KIND=dp) :: filedx
    real(KIND=dp), allocatable :: filegaps(:,:)
    
    1 format ('Number of mesh points does not correspond to file.', / &
    &         'On file: nx= ', i3, ' ny= ', i3, ' nz= ',i3,            / &
    &         'In data: nx= ', i3, ' ny= ', i3, ' nz= ',i3)
    2 format ('Number of wavefunctions does not correspond to file.', / &
    &         'On file: nwn= ', i3, ' nwp=', i3,                      / &
    &         'In data: nwn= ', i3, ' nwp=', i3)
    
    !---------------------------------------------------------------------------
    ! First check if the file exists.
    inquire(file=inputfilename, exist=exists)
    if(.not.exists) then
      print *, 'Input file specified does not exist!'
      stop
    endif
    !---------------------------------------------------------------------------
    open (chan,form='unformatted',file=ifn)
    
    read(chan, iostat=io) version
    ! Convergence information                                  (NOT IMPLEMENTED)
    read(chan,iostat=io) 
    !Parameters of the mesh
    read(Chan,iostat=io) filenx,fileny,filenz, filedx
    ! Symmetry information                                     (NOT IMPLEMENTED)
    read(Chan,iostat=io) 
    !Number of protons and neutrons
    read(Chan,iostat=io) neutrons,protons
    ! HFBLocks information 
    read(Chan,iostat=io) filenwn, filenwp, hfblocks
    filenwt = filenwn + filenwp
    ! Wavefunctions
    !- - - - - - - - - - - - - - - -
    ! First allocate the needed space
    allocate(HFPsi(filenx*fileny*filenz,4, filenwn+filenwp))
    allocate(spenergies(filenwn+filenwp))
    allocate(dispersions(filenwn+filenwp))
    allocate(rho_can(filenwn + filenwp))
    
    read(chan,iostat=io) spenergies, dispersions    
    read(chan,iostat=io) HFPsi                              
    ! Name of the force and functional
    read(chan, iostat=io) name_param, func_name_check

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
        allocate(BCSGaps(filenwt)) ;  BCSGaps = filegaps(:,1)        
    case(2)
        ! HFB

        allocate(filegaps(2*filenwt, 2*filenwt)) 
        read(chan, iostat=io) FermiEnergy       ! Lambda
        read(chan, iostat=io) ! rho
        read(chan, iostat=io) ! kappa
        read(chan, iostat=io) ! Canonical transformation
        read(chan, iostat=io) filegaps     ! Full matrix of gaps

        ! Simply copy the gaps for now
        allocate(HFBGaps(2*filenwt, 2*filenwt)) ; HFBGaps = filegaps  
    end select   
    !---------------------------------------------------------------------------
    ! Densities                                                (NOT IMPLEMENTED)
    ! No idea yet on how to implement this, as the nature of the densities
    ! calculated every calculation can be very different.      
    read(chan, iostat=io)

    ! Cranking information                                     (NOT IMPLEMENTED)
    read(chan, iostat=io)
    ! Multipole moment information                             (NOT IMPLEMENTED)
    read(chan, iostat=io)
    
    !---------------------------------------------------------------------------
    ! Sanity checks
    if((filenx.ne.nx).or. (fileny.ne.ny) .or. (filenz.ne.nz)) then
        print 1, filenx, fileny, filenz, nx,ny,nz
        stop
    endif
    if(filenwn.ne.nwn .or. filenwp.ne.nwp) then
        print 2, filenwn, filenwp, nwn, nwp
        stop
    endif
    !---------------------------------------------------------------------------
    ! Assign correct reflection symmetries for the derivative routines. 
    ! Should be handled by HEPHAESTOS in the future though.
    allocate(sx(4,nwt), sy(4,nwt), sz(4,nwt))
    do i=1, HFBlocks(1)
        sx(1,i) =  1 ; sy(1,i) = +1 ; sz(1,i) = +1
        sx(2,i) = -1 ; sy(2,i) = -1 ; sz(2,i) = +1 
        sx(3,i) = -1 ; sy(3,i) = +1 ; sz(3,i) = -1
        sx(4,i) =  1 ; sy(4,i) = -1 ; sz(4,i) = -1
    enddo
    do i=HFBlocks(1) + 1,HFBlocks(1) + HFBlocks(3)
        sx(1,i) =  1 ; sy(1,i) = +1 ; sz(1,i) = -1
        sx(2,i) = -1 ; sy(2,i) = -1 ; sz(2,i) = -1 
        sx(3,i) = -1 ; sy(3,i) = +1 ; sz(3,i) = +1
        sx(4,i) =  1 ; sy(4,i) = -1 ; sz(4,i) = +1
    enddo
    do i=HFBlocks(1) + HFBlocks(3)+1,HFBlocks(1) + HFBlocks(3) +HFBlocks(5)
        sx(1,i) =  1 ; sy(1,i) = +1 ; sz(1,i) = +1
        sx(2,i) = -1 ; sy(2,i) = -1 ; sz(2,i) = +1 
        sx(3,i) = -1 ; sy(3,i) = +1 ; sz(3,i) = -1
        sx(4,i) =  1 ; sy(4,i) = -1 ; sz(4,i) = -1
    enddo
    do i=HFBlocks(1)+HFBlocks(3)+HFBlocks(5) + 1,                      &
    &       HFBlocks(1)+HFBlocks(3)+HFBlocks(5) + HFBLocks(7)
        sx(1,i) =  1 ; sy(1,i) = +1 ; sz(1,i) = -1
        sx(2,i) = -1 ; sy(2,i) = -1 ; sz(2,i) = -1 
        sx(3,i) = -1 ; sy(3,i) = +1 ; sz(3,i) = +1
        sx(4,i) =  1 ; sy(4,i) = -1 ; sz(4,i) = +1
    enddo
    close(chan)
  end subroutine ReadTantalus

  subroutine WriteTantalus(chan, ofn)
    !---------------------------------------------------------------------------
    ! Subroutine that dumps all information to file for future runs.
    ! Heavily based on the MOCCa output routine.
    !---------------------------------------------------------------------------
    ! 
    ! Things written to file. (Not yet implemented ones are indicated by *)
    !
    ! Version
    ! Convergence information: E, dE                         (*)
    ! nx,ny,nz,dx,dt                    
    ! Symmetry information                                   (*)
    ! neutrons,protons
    ! Number of wavefunctions in every block
    ! (nwt) Wavefunctions                                    
    ! Forcename                                              
    ! Pairing information                                    
    !    - Pairingtype
    !    - Rho_can = occupation factors 
    !      (HF)  nothing
    !      (BCS) Fermi level
    !        |   BCSGaps  
    !      (HFB) Fermi level
    !        |   rho_pairing    
    !        |   kappa_pairing
    !        |   can_transfo
    !        |   HFBgaps            
    ! Densities                                              (*)
    ! CrankingInfo                                           (*)                      
    ! Moments                                                (*)
    !
    !---------------------------------------------------------------------------

    use functional

    integer, intent(in)          :: chan
    character(len=*), intent(in) :: ofn
    integer                      :: io
    
    open (chan,form='unformatted',file=ofn)

    write(chan, iostat=io) 1
    ! Convergence information                                  (NOT IMPLEMENTED)
    write(chan,iostat=io) 
    !Parameters of the mesh
    write(Chan,iostat=io) nx,ny,nz, dx
    ! Symmetry information                                     (NOT IMPLEMENTED)
    write(Chan,iostat=io) 
    !Number of protons and neutrons
    write(Chan,iostat=io) neutrons,protons
    ! HFBLocks information 
    write(Chan,iostat=io) nwn, nwp, hfblocks
    ! Wavefunctions  
    write(chan,iostat=io) spenergies, dispersions
    write(chan,iostat=io) HFPsi                              
    ! Name of the force.
    write(chan, iostat=io) name_param, func_name
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
        write(chan, iostat=io) FermiEnergy       ! Lambda
        write(chan, iostat=io) rho_pairing       ! rho
        write(chan, iostat=io) kappa_pairing     ! kappa
        write(chan, iostat=io) Cantransfo        ! Canonical transformation
        write(chan, iostat=io) HFBgaps           ! Full matrix of gaps
    end select
    !---------------------------------------------------------------------------

    ! Densities                                                (NOT IMPLEMENTED)
    ! No idea yet on how to implement this, as the nature of the densities
    ! calculated every calculation can be very different.      
    write(chan, iostat=io)   
     
    ! Cranking information                                     (NOT IMPLEMENTED)
    write(chan, iostat=io)
    ! Multipole moment information                             (NOT IMPLEMENTED)
    write(chan, iostat=io)
    
  end subroutine WriteTantalus
end module IO
