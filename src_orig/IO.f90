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
use functional
use momentsofinertia
use moments

implicit none

  !-----------------------------------------------------------------------------
  ! Filenames for in- and output of the code with respect to spwfs.
  character(len=100) :: inputfilename, outputfilename

  ! Signal the code to write extra output.
  character(len=40) :: BXLFIT = '', COMBI='', denfile=''

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
  
    implicit none

    ! These inputs control where the code will look for its input. Leaving them 
    ! empty will have the code rely on STDIN for input.
    integer(dp), intent(in), optional   :: file_number   
    character(11), intent(in), optional :: input_file 

    logical :: exists

    NameList /IO/ InputFileName,OutputFileName, BXLFIT, COMBI, denfile
    
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
    
    call readmomentdata(file_number)

    if(present(file_number)) then
      close(unit=file_number)
    endif

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
   99 format ( '  Nilsson initialization with hom = (', 3(f7.3) ,')')
   10 format ( ' IO information', / &
    &          '  inputfilename  =', a20, / &
    &          '  outputfilename =', a20)
   11 format ( '  BXL output     =', a20)
  111 format ( '  Combi ouput    =', a20)
   12 format ( ' Convergence required', / &
    &          '  Energy convergence           < ', e8.1, / & 
    &          '  Multipole moment convergence < ', e8.1, / &
    &          '  Dispersion convergence       < ', e8.1 )
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
    if(BXLFIT .ne. '') then
        print 11, BXLFIT 
    endif
    if(COMBI .ne. '') then
        print 111, COMBI
    endif
    print 12, energy_prec, moment_prec, disp_prec
    
    call printevolution
    call printscfiteration
    call printpairing_init
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
    ! CrankingInfo                                           (*)                      
    ! Potentials                                             
    ! Moments                                                (*)
    !
    !---------------------------------------------------------------------------
    
    use functional
    use moments
    
    integer, intent(in)          :: chan
    character(len=*), intent(in) :: ifn
    character(len=20)            :: func_name_check
    integer                      :: io, version
    logical                      :: exists
    
    integer       :: filenx, fileny, filenz, filenwn, filenwp,i, filepairing
    integer       :: filenwt, fileneutrons, fileprotons
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
    read(Chan,iostat=io) fileneutrons, fileprotons
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
    ! Cranking information                                     (NOT IMPLEMENTED)
    read(chan, iostat=io)
    !---------------------------------------------------------------------------
    ! Potentials                                               
    call readpotentials(chan)
    !---------------------------------------------------------------------------
    ! Multipole moment information                             
    io = 0
    do while(io.eq.0)
      call ReadMoment(chan,io)
    enddo
    
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

  subroutine WriteTantalus(chan, ofn, iter, iomsg)
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
    ! Multipole Moments                                                 
    !     | The code writes the data on ALL the multipole moments.
    !     | For the format of the lines, see the Moments module.
    !---------------------------------------------------------------------------

    use functional
    use moments

    integer, intent(in)          :: chan
    character(len=*), intent(in) :: ofn
    integer                      :: io
    type(moment), pointer        :: mom
    integer, intent(in)          :: iter
    character(len=*), intent(in) :: iomsg

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
    ! Cranking information                                     (NOT IMPLEMENTED)
    write(chan, iostat=io)
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
    
    ! Bonus file for quick feedback into the fit
    if(BXLFIT .ne. '') then
        call Brussels_output(iter, iomsg)
    endif  

    ! Output for the level density code
    if(COMBI .ne. '') then
        call Combi_output()
    endif

    ! Write the density to a file for postprocessing 
    if(DENFILE .ne. '') then
      call writedensity(D_I_I, DENFILE)
    endif
  end subroutine WriteTantalus

  subroutine Brussels_output(iter, iomsg)
    !---------------------------------------------------------------------------
    ! Write an extra file for use in the Brussels fitting protocol to
    !     Out/zXXXnXXX.out   
    !
    ! It contains on a single line
    !
    !      N, Z, Total energy, Q20(t), Q22(t), Q(t), Q(t),  &   
    ! &    Gamma(n), Gamma(p),  <r^2_p>, Rotcorrection, iter, io
    !
    ! Notes:
    ! *  <r^2_p> is calculated as in the moments module, i.e. it is calculated  
    !    from the charge density, which is not necessarily the proton density.
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

    type(Moment), pointer :: Q20, Q22, r2
    real(KIND=dp)         :: E, quad(2), rms, q2(3)
    integer, intent(in)   :: iter

    character(len=len(BXLFIT)+12) :: filedone
    
    Q20 =>FindMoment( 2,0,.false.)
    Q22 =>FindMoment( 2,2,.false., Q20)
    r2  =>FindMoment(-2,0,.false., Q22) ! The rms radius is associated with l=-2

    ! Write the filename
    write(filedone,'(a,"z",i3.3,"n",i3.3,".out")'), trim(adjustl(BXLFIT)),     &
    & int(protons),int(neutrons)
  
    open(unit=10,file=filedone)

    E = TotalE
    quad(1) = sum(Q20%value)
    quad(2) = sum(Q22%value)    
    rms     =     r2%value(2)
    q2      = CalculateTotalQl(2)

    write(10,'(2i4,7f15.6, i6)', advance='NO')  &
    &     int(protons),int(neutrons),E,quad, q2(3), G(3)*180/pi, sqrt(rms/protons),&
    &     Rotcorrection, iter
  
    write(10, '(2x, a99)') adjustl(iomsg)
    close(10)

  end subroutine Brussels_output

  subroutine combi_output
    !---------------------------------------------------------------------------
    ! Write an extra file for input of the combinatorial level density code.
    !
    !
    ! ATTENTION: this output assumes an axial nucleus with a symmetry axis 
    !            along the z-axis. If the single-particle states are not  
    !            (at least approximately) eigenstates of J_z, then this output
    !            will effectively be nonsense.
    !---------------------------------------------------------------------------

    integer, allocatable :: indices(:)
    integer              :: i,ii, p1, p2,jj
    real(KIND=dp)        :: R0, A, fac, mstate1, mstate2

    1 format (a1, 3i4)
    2 format (2(f5.1,i2,3f8.3))
    3 format ( 2i4,17(x,f7.4),2f9.2)

    open(unit=6, file=COMBI)

    ! a) single-particle neutron states
    write(6, fmt=1) '*', int(protons), int(neutrons+protons), nwn/2
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

        write(6,fmt=2),mstate1,p1,spenergies(ii),rho_can(ii)/2,BCSgaps(ii),& 
        &              mstate2,p2,spenergies(jj),rho_can(jj)/2,BCSgaps(jj) 
    enddo

    ! b) single-particle neutron states
    write(6, fmt=1) ' ', int(protons), int(neutrons+protons), nwp/2
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
        
        write(6,fmt=2),mstate1,p1,spenergies(ii),rho_can(ii)/2,BCSgaps(ii),& 
        &              mstate2,p2,spenergies(jj),rho_can(jj)/2,BCSgaps(jj) 
    enddo

    !  The final line is composed of various informations
    !  Z, A, beta2, beta4, Gn, Gp, Deltan, Deltap, ddmn, ddmp,
    !       econdn,econdp,eshcorn,eshcorp,lambdan,lambdap,ainer,rigid,
    !       etott,etable

    R0  = 1.2
    A   = neutrons + protons
    fac = 4. * pi/(3. * (R0 *(A))**2 * A) * Q(3) * sqrt(5/(16*pi))

    !                                              Q40   Gn   Gp  Deltan Deltap  
    write(unit=6, fmt=3), int(protons),int(A),fac*Q(3), 0.0, 0.0, 0.0, 0.0,0.0,&
    !                     ddmn, ddmp, econdn, econdp, eshcorn, eshcorp 
    &                      0.0,  0.0,    0.0,    0.0,     0.0,     0.0,        & 
     !                     lambdan, lambdap, ANONYMOUOS NUMBER,  ainer, rigid, 
    &                      FermiEnergy(1),  FermiEnergy(2), 0.0, Belyaev(3,3), &
    &                      Rigid(3,3), totalE, 0.0

    close(unit=6)
  end subroutine combi_output
    
  function force_halfinteger(j) result(jforced)
      
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
end module IO
