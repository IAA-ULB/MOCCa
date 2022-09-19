module pairing
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
 !
 ! High-level module delegating the solving of the pairing subproblem, as well
 ! as some routines of general interest.
 !
 !==============================================================================
 ! Hephaestos keywords:
 !
 ! FORBIDBCS : $FORBIDBCS 
 ! TR        : $TR
 ! NTR       : $NTR
 ! PBROKEN   : $PBROKEN
 !==============================================================================

 use compilation
 use wavefunctions
 use hartreefock
 use BCS
 use HFB
 use HFB_gradient
 use pairingcutoffs
 use timing
 use parameterization  

 implicit none
 
 !------------------------------------------------------------------------------
 ! Pairing density matrix and anomalous density matrix in the "as is" basis.
 real(KIND=dp), allocatable :: rho_pairing(:,:),  kappa_pairing(:,:)
 ! ... and in the physical HF basis, but only the DIAGONAL matrix elements
 real(KIND=dp), allocatable :: rho_HF(:)
 ! ... and in the canonical basis ...
 real(KIND=dp), allocatable :: rho_can(:), kappa_can(:)
 ! Note: the object kappa_can only has an effect on the calculation in the BCS 
 !       case, when kappa is actually off-diagonal. In the HFB case, no 
 !       guarantees are given as to the canonical form of kappa, and the 
 !       code does not rely on it being this way.
 ! ...  and the "configuration matrix" ...
 real(KIND=dp), allocatable :: configmatrix(:)
 ! ... and finally, the Bogoliubov transformation.
 real(KIND=dp), allocatable, target :: Bogoliubov(:,:)
 ! Do we start with the Bogoliubov transformation from file? 
 ! This is important for the gradient solver, though not as much for the 
 ! HFB_direct solver. 
 logical             :: BogoFromFile   = .true.
 !------------------------------------------------------------------------------
 ! Quasiparticle excitation energies, either HF, BCS or HFB.
 real(KIND=dp), allocatable :: QPenergies(:)
 !------------------------------------------------------------------------------
 ! Transformation from the HFBasis into the canonical basis
 real(KIND=dp), allocatable :: CanTransfo(:,:)
 ! Transformation from the HFBasis into the basis where Kappa (with cutoffs)
 ! is canonical. (NOT CALCULATED AT THE MOMENT)
 real(KIND=dp), allocatable :: CanCutTransfo(:,:)
 !------------------------------------------------------------------------------
 ! Fermi energy for neutrons and protons.
 real(KIND=dp) :: FermiEnergy(2) , FermiHistory(2)
 !------------------------------------------------------------------------------
 ! Particle number dispersion
 real(KIND=dp) :: Dispersion(2)
 !------------------------------------------------------------------------------
 ! Type of pairing to employ.
 ! (0): Hartree-Fock
 ! (1): Hartree-Fock + BCS
 ! (2): Hartree-Fock-Bogoliubov
 integer :: PairingType = 0
 !------------------------------------------------------------------------------
 ! Decide which Fermisolver to use for direct HFB solution strategies
 ! "Brent"  => use a modified bisection solver
 ! "Secant" => use a secant routine
 character(len=99) :: FermiSolver='Brent'
 !------------------------------------------------------------------------------
 ! Decide which module gets to calculate the pairing gaps.
 procedure(calcBCSgaps), pointer :: CalcGaps
 !------------------------------------------------------------------------------
 ! Type of blocking we want. 
 ! (0) no blocking
 ! (1) ordinary blocking, based on indices
 ! (2) ordinary blocking, asking for lowest energy configurations
 ! (3) EFA blocking, based on indices.
 ! (4) EFA blocking, asking for lowest energy configurations.
 ! (5) ordinary blocking, index selection through overlap
 ! (6) EFA blocking, index selection through overlap
 !  
 ! If this is nonzero, the code will look for a new namelist "Indices"
 !
 !------------------------------------------------------------------------------
 integer :: BlockType  = 0
 integer :: BlockNumber= 0 
 !------------------------------------------------------------------------------
 ! Indices of the levels to block. 
 integer, allocatable :: BlockIndices(:) 
 !------------------------------------------------------------------------------
 ! Indices of the quasi-particles that ended up blocked, as well as their 
 ! closest partners under time-reversal. This information is  only used to 
 ! regularize the calculation of the rotational energy.
 integer, allocatable       :: blocked_qps(:), partner_qps(:)
 real(KIND=dp), allocatable :: partner_overlaps(:)
 !------------------------------------------------------------------------------
 ! Indices (in the CANONICAL basis) of the single-particle states that obtained
 ! v^2 = 1 through blocking
 integer, allocatable :: blocked_sps(:)
 !------------------------------------------------------------------------------ 
 ! EFA blocking for the lowest qp. 
 ! Indicated by either
 ! 'n+' : positive parity neutron
 ! 'n-' : negative parity neutron
 ! 'n0' : (any parity)    neutron
 ! 'p+' : positive parity proton
 ! 'p-' : negative parity proton
 ! 'p0' : (any parity)    proton
 ! '  ' : nothing
 character(len=2), allocatable :: BlockLowest(:)
 !------------------------------------------------------------------------------
 ! Entropy of the statistical mixture in the case of finite temperature
 real(KIND=dp) :: entropy(2) = 0
 !------------------------------------------------------------------------------
 ! Values for the average pairing gap
 ! First index is for 
 ! (1) f_k v_k ^2 weighted
 ! (2) f_k u_k v_k weighted
 !  The second index is for isospin.
 real(KIND=dp) :: average_gap(2,2)
 !------------------------------------------------------------------------------
 ! Integer tracking the way the contribution from a free gas is counted for the 
 ! readjustment of the Fermi energy.
 ! 0) No correction: the number of particles is the ordinary counting.
 ! 1) Subtraction  : the nucleus is modelled as being immersed in a gas of 
 !                   free particles. The number of particles is the calculated
 !                   one MINUS the number of particles in a free gas at the 
 !                   same chemical potential.
 ! 2) Nogas        : the code is only allowed to occupy the bound states.
 integer(sp)   :: particles_in_gas = 0
 real(KIND=dp) :: ngas(2)
 !------------------------------------------------------------------------------
 ! Whether or not to guess some pairing gaps when starting the code.
 logical :: guessgaps = .false.
 ! ... and with which values to initialize them. Negative values signal
 ! the code to use the default value, 1.5 MeV.
 real(KIND=dp) :: gapvalue(2) = -2
 !-----------------------------------------------------------------------------
 ! Determine the scheme used for solving the HFB problem in the HF basis
 !  (0) => Direct solution, i.e. construction and diagonalisation of the 
 !         HFB Hamiltonian
 !  (1) => Gradient solution, i.e. following the manifold of HFB solutions
 integer :: pairingscheme = 0

contains

  subroutine initpairing(file_number)
    !---------------------------------------------------------------------------
    ! Read and initialize pairing options from the namelists
    !   /Pairing/
    !   /Indices/
    !---------------------------------------------------------------------------
    character(len=20)                   :: Type = 'HF'
    integer(dp), intent(in), optional   :: file_number   
    integer                             :: i
    
    NameList /Pairing/ Type, Constantgap,                                      &
    &                  BlockType, BlockNumber, particles_in_gas, maxhfbiter,   & 
    &                  FermiSolver, guessgaps,  pairingscheme,                 &
    &                  gradient_precon, bogofromfile, gapvalue

    NameList /Indices/ BlockIndices, blocklowest, blockfname

    if(present(file_number)) then
      read(unit=file_number, NML=Pairing)
    else
      read(unit=*, NML=Pairing)
    endif  

    Type        = to_upper(Type)
    if('HF' .eq.adjustl(type)) then
      pairingtype = 0
    elseif('BCS' .eq. adjustl(type)) then
      pairingtype = 1
    elseif('HFB' .eq. adjustl(type)) then
      pairingtype = 2
    elseif('' .eq. adjustl(type)) then
      pairingtype = 0
    else
      print *, 'This type of pairing is not implemented yet.'
      stop
    endif

$FORBIDBCS if( pairingtype .eq. 1) then
$FORBIDBCS    print *, "BCS pairing treatment not allowed."
$FORBIDBCS    stop
$FORBIDBCS endif

    FermiSolver = to_upper(FermiSolver)    
    if(adjustl(FermiSolver).eq.'SECANT') then
      FindFermi => FindFermi_secant
    elseif(adjustl(FermiSolver).eq.'BRENT') then
      FindFermi => FindFermi_brent
    else
      print *, 'Unknown FermiSolver', FermiSolver, ' selected.'
      stop
    endif
    
    if(Blocktype.lt.0 .or. BlockType.gt.6) then
        print *, 'This value of BlockType is not accepted.'
        stop
    endif
    
    if(particles_in_gas .lt. 0 .or. particles_in_gas .gt. 2) then
      print *, 'This value for particles_in_gas is not accepted.'
      stop
    endif

    if(particles_in_gas .ne. 0 .and. inversetemp .eq. -1) then
       print *, 'Particles_in_gas should be zero for T=0 calculations.'
      stop
    endif
    !---------------------------------------------------------------------------
    ! Reading information on the blocking if needed.
    if(BlockNumber.ne.0) then
        ! Sanity check: only allow for blocking in HFB mode
        if(pairingtype.ne.2) then 
          print *, 'Blocking only allowed when doing HFB calculations.'
          stop
        endif

        allocate(BlockIndices(BlockNumber)) ; BlockIndices = 0
        allocate(BlockLowest(BlockNumber))  ; BlockLowest  = ' ' 
        read(unit=*, nml=Indices)

        ! Sanity checks on the BlockLowest array: 
        ! (a) do not proceed with empty list
        ! (b) do not allow for selection on parity of the blocked state if
        !  parity is broken
        if(blocktype.eq.2 .or. blocktype .eq. 4) then
          do i=1, blocknumber
            select case(blocklowest(i))
            case('n+', 'n-')
$PBROKEN              print *, 'Cannot block a neutron qp with definite parity.'
$PBROKEN              stop
            case('p+', 'p-')
$PBROKEN              print *, 'Cannot block a proton qp with definite parity.'
$PBROKEN              stop
            case('n0', 'p0')
              ! allowed
            case DEFAULT
              ! something else went wrong
              print *, 'Did not read all elements in blocklowest correctly.'
              stop
            end select
          enddo
        endif
        
        ! Sanity checks on the BlockIndices array:
        ! (a) check if everything was read
        if(blocktype.eq.1 .or. blocktype.eq.3) then
          do i=1,blocknumber
            if(blockindices(i) .eq.0) then
              print *, 'Did not read all elements in blockindices correctly.'
              stop
            endif
          enddo
        endif

        ! Sanity check on the useage of the gradient solver for EFA
!        if((blocktype.eq.3 .or. blocktype.eq.4) .and. pairingscheme .eq.1) then
!          print *,' The gradient solver cannot yet handle EFA blocking.'
!          stop
!        endif 
        
        ! Sanity check: cannot do full blocking if time-reversal is not broken
$TR        if(blocktype.eq.1 .or. Blocktype.eq.2) then
$TR         print *, 'Cannot do true blocking when time-reversal is conserved.'
$TR         stop
$TR        endif

        ! Reading model spwf to block
        if(blockfname .ne. "") then
           call read_modelwf(blockfname)
           if(blocknumber.gt. 1) then
              print *, 'Cannot block more than one modelspwf.'
              stop
           endif
        endif
    endif

    !---------------------------------------------------------------------------
    ! Cutoff decision
    select case(CutType)
    case(1)
       PairingCutoff => SymmetricFermi
    case(2)
       PairingCutoff => FermiAbove
    case(3)
       PairingCutoff => CosineCut
    case(4)
       PairingCutoff => SymmetricFermizero
    case DEFAULT
       print *, 'Unknown cutoff type CutType. Valid options are 1-4.'
       stop
    end select
    pairingcut(1) = cutneutron
    pairingcut(2) = cutproton
    !---------------------------------------------------------------------------
    select case(PairingType)
    case(0)
      CalcGaps => calcHFgaps
    case(1)
      CalcGaps => calcBCSgaps
    case(2)
      CalcGaps => calcHFBgaps
    end select
    !---------------------------------------------------------------------------
    ! 
    if((pairingscheme .ne. 0) .and. (pairingscheme.ne.1)) then
      print *, 'Invalid pairingscheme value.'
      stop
    endif
  end subroutine initpairing

  subroutine printpairing_init
    !---------------------------------------------------------------------------
    ! Print information on the treatment of pairing at the start of the run.
    !---------------------------------------------------------------------------
    1 format(80('-'))
    2 format(' Pairing treatment: ', a60)
   21 format('   Fermi-solver: ', a99 )
  211 format('   Pairing strategy:', a60)
 2111 format('     -> Fermi tolerance: ', es10.3)
    6 format('   Cutoff parameters  = ', a20)
    7 format('     dE (n,p) = ', 2f5.2, ' MeV ')
    8 format('     mu (n,p) = ', 2f5.2, ' MeV ')
   81 format('   Stabilisation active')
   82 format('    Estab(p,n)= ', 2f4.1, ' MeV')
   83 format('   Initialisation:')
  831 format('   -> guessed initial gaps Delta')
 8311 format('       Gn = ', f4.1, ' MeV ', ' Gp = ', f4.1, ' MeV.' )
  832 format('   -> started with HFB transformation from file')
  833 format('   -> started with explicit vacuum construction')

   13 format('   Gas-treatment:  Normal'            )    
   14 format('   Gas-treatment:  Subtraction method')    
   15 format('   Gas-treatment:  Bound states only')    

  100 format('   Fixed Fermi= ', 2f12.4)

   90 format ('  Blocking Options')
   91 format ('    Blocking type: ', i1)
   92 format ('    Ordinary Blocking' )
   93 format ('    Equal Filling    ' )
   
   10 format ('    Blocknumber  = ', i2 )
   11 format ('    Blocklowest  = ', 20(1x, a2))
   12 format ('    BlockIndices = ', 20i3)
   16 format ('    Block through overlap')
   17 format ('    Blockfile    = ', 40a)

    character(len=60) :: ptreat, pscheme

    print 1

    select case (pairingtype)
    case(0)
        ptreat = 'Hartree-Fock (HF)'
        print 2, adjustl(ptreat)
    case(1)
        ptreat = 'Bardeen-Cooper-Schrieffer (HF+BCS)'
        print 2, adjustl(ptreat)
        print 2111, pairing_prec
    case(2)
        ptreat = 'Hartree-Fock-Bogoliubov (HFB)'
        print 2, ptreat
 
        select case(pairingscheme)
        case(0)
          pscheme = 'Direct diagonalisation'
          print 211, adjustl(pscheme)
          print 21, adjustl(FermiSolver)
        case(1)
          pscheme = 'Gradient solver'
          print 211,  adjustl(pscheme)
        end select 
        print 2111, pairing_prec
    end select

    select case(CutType)
    case(1)
       print 6, 'Symmetric Fermi'
    case(2)
       print 6, 'Fermi above'
    case(3)
       print 6, 'Cosine'
    case(4)
       print 6, 'Sym. Fermi + Heaviside'
    case DEFAULT
       print *, 'Unrecognized type of pairing cutoff.'
       print *, 'Accepted values are 1-4.'
       stop
    end select

    print 7, pairingcut
    print 8, PairingMu

    if(fixfermi) then
        print 100, mun, mup
    endif
  
    if(abs(Estabp).gt.1d-10 .or. abs(Estabn).gt.1d-10) then
      print 81
      print 82, Estabp, Estabn
    endif

    print 83
    if(pairingtype.eq.2) then
      if(bogofromfile) then
        print 832
      else
        print 833
      endif 
    endif
    if(guessgaps) then
      print 831
      print 8311, gapvalue
    endif

    if(particles_in_gas .eq.1) then
      print 14
    elseif(particles_in_gas .eq.2) then
      print 15
    else
      print 13
    endif

    if(Blocktype .ne. 0) then
        print 90
        print 91, Blocktype
        select case(Blocktype)
        case(1)
            print 92
            print 10, Blocknumber
            print 12, Blockindices
        case(2)
            print 92
            print 10, Blocknumber   
            print 11, Blocklowest
        case(3)
            print 93
            print 10, Blocknumber
            print 12, BlockIndices
        case(4)
            print 93
            print 10, Blocknumber
            print 11, Blocklowest
        case(5)
            print 92
        case(6)
            print 93
            print 16
            print 17, adjustl(blockfname)
        end select
    endif

  end subroutine printpairing_init
  
  subroutine initializeGaps( gapvalue )
    !---------------------------------------------------------------------------
    ! Allocate the pairing gaps (BCSgaps or HFBgaps) and assign them a 
    ! value according to the users wishes. 
    !
    ! Note that 
    !
    ! a) this routine respects symmetries, i.e. gaps forbidden by symmetry
    !    are initialized to zero, independently from gapvalue.
    !
    ! b) This routine is not necessary (but will function) when starting the 
    !    code from a .wf file. In that case, the gaps will be overwritten.
    !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    !
    ! Input:
    !   gapvalue :  values to assign to the pairing gaps for neutrons and protons
    !               in MeV. If gapvalue < 0, use the default value, 1.5 MeV.
    !---------------------------------------------------------------------------
    integer                   :: wave, wave2, si, B, N, s, N2, it
    real(KIND=dp), intent(in) :: gapvalue(2) 
    real(KIND=dp)             :: fill(2)
    
    do it=1,2
      if(gapvalue(it) .lt. 0) then
        fill(it) = 1.5
      else
        fill(it) = gapvalue(it)
      endif
    enddo
    
    

    select case (PairingType)
    case(0)
      !-------------------------------------------------------------------------
      ! Nothing to do for HF calculations
      return
    case(1)
      !-------------------------------------------------------------------------
      !BCS Calculation
      if(.not.allocated(BCSGaps)) then
        allocate(BCSGaps(nwt)) ; BCSGaps = 0.0
      endif
      do wave=1,nwt
        BCSgaps(wave) = 0.1
      enddo
    case(2)
      !-------------------------------------------------------------------------
      ! HFB calculations
      if(.not.allocated(HFBGaps)) then
        allocate(HFBGaps(nwt,nwt)) ; HFBGaps = 0.0
      endif
  
      HFBGaps = 0.0
      si = 0 
      do B= 1,8,2 ! Loop over only half of the blocks
        N = HFBlocks(B) ; if(N.eq.0) cycle
        N2= HFblocks(B+1)
        
        if(B .ge. 5) then
          it = 2
        else
          it = 1
        endif
        
        do wave=si+1,si+N
             !------------------------------------------------------------------
             ! If there is a conserved time-like symmetry, then we store only
             ! half of the gaps, see HFB.f90
$TR             do wave2=si+1,si+N
$NTR          do wave2=si+N+1,si+N+N2


            ! We check if the spwfs are (close to) having opposite parity
            if(abs(P_hf(wave) * P_hf(wave2) + 1 ).lt. 0.1) cycle
            !------------------------------------------------------------------
            if(allocated(kappa_pairing)) then
              ! We've found a kappa on file and can use it to guess better 
              ! signs and sizes
              if(abs(kappa_pairing(wave, wave2)).gt.1d-8) then
                s = int(kappa_pairing(wave, wave2)/abs(kappa_pairing(wave, wave2)))
              else
                s = 1
              endif
              ! Guess something with the sign of kappa beween 0.3 and 1.5
              HFBgaps( wave, wave2) = &
              &     s*max(min(10*abs(kappa_pairing(wave, wave2)),1.5),0.3)
            else
              ! Either we don't have a kappa in storage, or the pairing has 
              ! collapsed in this subblock.
              HFBgaps( wave, wave2) = fill(it)
            endif
$NTR        HFBgaps(wave2, wave) = -HFBgaps(wave, wave2)
          enddo 
        enddo
        si = si + N + N2
      enddo
    end select  
  end subroutine initializeGaps
  
  subroutine SolvePairing(scheme,ifail)
    !---------------------------------------------------------------------------
    ! Master routine for the solving of the pairing equations.
    !
    ! Input:  
    !   Scheme : determines the type of solution method 
    !
    !        ( 0) => Direct diagonalisation of the HFB Hamiltonian, followed by
    !                explicit construction of a Bogoliubov vacuum state
    !        (+1) => Perform a heavy-ball step in the limited subspace
    !        (-1) => Do everything similar to a heavy-ball step, but don't 
    !                actually update the Bogoliubov transformation
    ! Output:
    !   ifail  : if non-zero, something went wrong with a diagonalization 
    !---------------------------------------------------------------------------
    use parameterization, only : hbm

    integer, intent(in)        :: scheme
    integer, intent(out)       :: ifail
    real(KIND=dp), allocatable :: sphamil(:,:)

    call start_timer(T_pairing)

    if(.not.allocated(rho_can)) then
      allocate(rho_can(nwt))              ; rho_can    = 0.0
    endif
    if(.not.allocated(kappa_can)) then
      allocate(kappa_can(nwt))            ; kappa_can  = 0.0 
    endif
    if(.not.allocated(qpenergies)) then
        if(PairingType .le. 1)   allocate(QPenergies(nwt))         
        if(PairingType .eq. 2)   allocate(QPenergies(2*nwt))         
        qpenergies = 0.0
    endif

    ! Always allocate the configuration matrix C
    if(.not.allocated(configmatrix)) then
       allocate(configmatrix(2*nwt))      ; configmatrix = 0.0
    endif
        
    select case (Pairingtype)
    case(0)
        if(inversetemp .eq. -1) then
            call NaiveFill(rho_can)
        else
            call FiniteTemperatureHF(rho_can,FermiEnergy, particles_in_gas)
        endif
        ! Make sure the program does not stop because this is uninitialized

        ifail = 0
    case(1)
      !-------------------------------------------------------------------------
      ! BCS-type pairing
      ! The diagonal elements of rho and kappa are only set. 
      call solvepairing_BCS(FermiEnergy, rho_can, kappa_can, qpenergies,       &
      &                     particles_in_gas, BlockType, Blockindices,         & 
      &                     blocklowest, blocked_qps)
      ! Make sure the program does not stop because this is uninitialized
      ifail = 0
    case(2)
      !-------------------------------------------------------------------------
      ! HFB-type pairing
      if(.not.allocated(CanTransfo)) then
        ! Allocate the full matrices
        allocate(CanTransfo(nwt, nwt))     ; CanTransfo    = 0.0
      endif
      if(.not.allocated(CanCutTransfo)) then
        allocate(CanCutTransfo(nwt, nwt))  ; CanCutTransfo = 0.0
      endif
      if(.not. allocated(rho_pairing)) then
        allocate(rho_pairing(nwt,nwt))     ; rho_pairing   = 0.0
      endif
      if(.not.allocated(kappa_pairing)) then
        allocate(kappa_pairing(nwt,nwt))   ; kappa_pairing = 0.0
      endif     
      if(.not.allocated(Bogoliubov)) then
        allocate(Bogoliubov(2*nwt,2*nwt))  ; Bogoliubov    = 0.0
      endif
      if(.not.allocated(qpdispersions)) then
        allocate(qpdispersions(2*nwt)) ; qpdispersions = 0.0d0
      endif
      
      ! Depending on the algorithm in use, we build a different single-particle
      ! hamiltonian matrix.
      sphamil = build_sph(scheme, efficientHFB)

      !-------------------------------------------------------------------------
      ! Find the Fermi energy
      select case(scheme)
      case( 0)
        call solvepairing_HFB_direct(  &
        &   sphamil,HFBgaps,FermiEnergy,Bogoliubov,rho_pairing,kappa_pairing,  &
        &   configmatrix, qpenergies,BlockType, Blockindices, blocklowest,     &
        &   blocked_qps, partner_qps, partner_overlaps, ifail)
      case(+1)
        call solvepairing_HFB_gradient( &
        &   sphamil,HFBgaps,FermiEnergy,Bogoliubov,rho_pairing,                &
        &   kappa_pairing, configmatrix, qpenergies,BlockType, Blockindices,   &
        &   blocklowest, blocked_qps, partner_qps, partner_overlaps,           &
        &   .true. , 1, ifail)
      case(-1)
        call solvepairing_HFB_gradient( &
        &   sphamil,HFBgaps,FermiEnergy,Bogoliubov,rho_pairing,                &
        &   kappa_pairing, configmatrix, qpenergies,BlockType, Blockindices,   &
        &   blocklowest, blocked_qps, partner_qps, partner_overlaps,           &
        &  .false., 1, ifail)
      end select
    end select
    ! Construct the density in the Hartree-Fock basis 
    ! (which is generally only used for printing)
    if(pairingtype.eq.2) then
      rho_hf =  construct_rho_HF(rho_pairing, HFtransfo)
    else
      rho_hf = rho_can/2.0d0
    endif  
    !---------------------------------------------------------------------------
    ! If beta != infty, we calculate the number of particles in the gas
    if(inversetemp.ne.-1) then
      ngas(1) = gasoccupations(fermienergy(1), hbm(1))
      ngas(2) = gasoccupations(fermienergy(2), hbm(2))
    endif
    !---------------------------------------------------------------------------
    ! Compute the cutoffs
    call ComputePairingCutoffs(fermienergy)
    call stop_timer(T_pairing)
   
    if(blocknumber.ne.0) then
      blocked_sps = identify_blocked_particle(Bogoliubov) 
    endif   
  end subroutine SolvePairing

  subroutine calc_avg_gap()
    !---------------------------------------------------------------------------
    ! Calculate the average gap
    !---------------------------------------------------------------------------
    select case (PairingType)
    case(0)
      average_gap = 0
    case(1)
      average_gap = average_gap_BCS()
    case(2)
      average_gap = average_gap_HFB()
    end select
  end subroutine calc_avg_gap

  function build_sph(pscheme, efficient) result(sph)
    !---------------------------------------------------------------------------
    !
    !   
    !---------------------------------------------------------------------------
    real(KIND=dp), allocatable :: sph(:,:)
    integer, intent(in)        :: pscheme
    logical, intent(in)        :: efficient 
    integer                    :: i

    allocate(sph(nwt,nwt)) ; sph = 0.0d0

    if((pscheme.eq. 0 .and. (.not. efficient)) &
    &   .or. (.not. allocated(current_sph))) then
      ! Diagonal part
      do i=1, nwt
        sph(i,i) = spenergies(i)
      enddo
    else
      ! Full matrix
      sph = current_sph
    endif
 
  end function build_sph

  subroutine printpairing(stabfactor)
    !---------------------------------------------------------------------------
    !
    !---------------------------------------------------------------------------

    real*8, intent(in) :: stabfactor(2)

    1 format (26('-'), ' Pairing ', 25('-'))
    2 format (25x, ' N ',7x, ' P ')
    3 format (' Fermi Level (MeV) ',2x,f13.8,2x,f13.8)
    4 format (' Particles         ',2x,f13.8,2x,f13.8)
    5 format (' Dispersion        ',2x,f13.8,2x,f13.8)
    6 format (' Average gap   v^2 ',2x, f13.8, 2x, f13.8)
   61 format (' Average gap   uv  ',2x, f13.8, 2x, f13.8)
    7 format (60('-'))

    8 format ('  gas-like          ', 2x, f13.8, 2x, f13.8)
    9 format ('  nucleus           ', 2x, f13.8, 2x, f13.8)

   10 format (' Stab. factor       ', 2x, f13.8, 2x, f13.8)
   11 format (' Overlap with model ', 2x, f13.8)
 
!   12 format ('                               ++  +-  -+  --')
!   13 format (' Number parity    n:', 2x, 4i3)
!   14 format (' Number parity    p:', 2x, 4i3)

    select case(PairingType)
    case (0)
        if(inversetemp .eq. -1) return
        
        print 1
        print 2
        print 3, FermiEnergy
        print 4, sum(rho_can(1:nwn)), sum(rho_can(nwn+1:nwt))

        if(inversetemp .ne. -1)then
          print 8, ngas
          print 9, sum(rho_can(1:nwn))-ngas(1), sum(rho_can(nwn+1:nwt))-ngas(2)
        endif
  
        print 5,  HFdispersion
    case (1,2)
        ! BCS and HFB
        print 1    
        print 2
        print 3, FermiEnergy
        print 4, sum(rho_can(1:nwn)), sum(rho_can(nwn+1:nwt))
        select case(PairingType)
        case(1)
            print 5, BCSdispersion
        case(2)
            print 5, HFBdispersion
        end select
        print 6,  average_gap(1,:)
        print 61, average_gap(2,:) 

        if(abs(Estabp).gt.1d-10 .or. abs(Estabn).gt.1d-10) then
          print 10, stabfactor
        endif
        if(blocktype.ge.5) then
            print 11, blockoverlap
        endif
!        
!$NTR    nb = number_parity_throughU(Bogoliubov, configmatrix, grad_blocks)        
!$NTR    print 12
!$NTR    print 13, nb(1:4)
!$NTR    print 14, nb(5:8)

        if(pairingtype.eq.2) then
          call PrintHFBConvergence(rho_pairing, kappa_pairing, Bogoliubov)
        endif
    end select
    print 7
  end subroutine PrintPairing
  
  function calcpairingenergy() result(E)
    !---------------------------------------------------------------------------
    ! Calculate the pairing energy.
    !-----------------------------------------------------f----------------------
    integer       :: wave, it, wave2
    real(KIND=dp) :: E(2)
    
    E = 0.0
    
    select case(PairingType) 
    case(0)
      ! HF case
    case(1)
      ! BCS case
      do wave = 1,nwt
          it = 1
          if(wave .gt. nwn) it = 2
          E(it) = E(it) - BCSgaps(wave)*Kappa_can(wave)
      enddo
    case(2)
      ! HFB case
      do wave = 1,nwt
        do wave2 = 1,nwt
          it = 1
          if(wave .gt. nwn) it = 2
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
          ! The pairing energy, full-blown is 
          ! E_pair = 1/2  Tr (Delta kappa^{\dagger})
          !        = 1/2 sum_ij Delta_{ij} \kappa^{\dagger}_{ji}
          !        = 1/2 sum_ij Delta_{ij} \kappa^*_{ij}
          ! 
          ! When we have no antilinear, antihermitian symmetry then we sum 
          ! exactly this formula. Otherwise, both kappa and delta have the 
          ! following structure
          !
          !  ( 0      Delta)     (  0      kappa )
          !  (-Delta     0 )     ( -kappa    0   )
          !
          ! but only the part on the upper right is actually stored in memory
          ! 

$TR         E(it) = E(it) - Kappa_pairing(wave,wave2)*HFBgaps(wave,wave2)
$NTR         E(it) = E(it) + 0.5 * Kappa_pairing(wave,wave2)*HFBgaps(wave,wave2)
        enddo
      enddo
    end select
  end function calcpairingenergy

  subroutine CalcEntropy()
    !---------------------------------------------------------------------------
    ! Calculate the entropy associated with a given statistical mixture at
    ! finite temperature.
    !---------------------------------------------------------------------------

    integer       :: i, it
    real(KIND=dp) :: fac

    if(inversetemp == -1) then  
        entropy  = 0
        return
    endif

    select case(PairingType)
    case(0)
        !-----------------------------------------------------------------------
        entropy = 0
        do i=1,nwt
            if(i .gt. nwn) then
                it = 2
            else
                it = 1
            endif
            if(rho_can(i) .gt. 1d-10) then
                ! Notice the factor 2 for time-reversal. The rho_cans are
                ! double the true occupations!                 
                entropy(it) = entropy(it) &
                &          -       rho_can(i)/2.0  * dlog(     rho_can(i)/2.0)
                if( rho_can(i) .lt. 2.0d0) then 
                    entropy(it) = entropy(it) &
                    &      -  (1 - rho_can(i)/2.0) * dlog( 1 - rho_can(i)/2.0)
                endif
            endif
        enddo
        ! Restoring the factor 2.
        entropy = 2 * entropy
    case (1)
        !-----------------------------------------------------------------------
        entropy = 0
        do i=1,nwt
            if(i .gt. nwn) then
                it = 2
            else
                it = 1
            endif

            ! We use the BCSf array, and not a recalculated version of the 
            ! occupation factors: depending on other options they might have 
            ! changed. 
            fac = BCSf(i)        
            if(fac .gt. 1d-10) then
              entropy(it) = entropy(it) - fac * log(fac) 
              if( (1 - fac) .gt. 1d-10) then            
                entropy(it) = entropy(it)- (1-fac) * log(1-fac)
              endif
            endif
        enddo
        entropy = 2 * entropy
    case (2)
        entropy = 0
        !-----------------------------------------------------------------------
        ! Notice that this sum is over ALL quasiparticles, not just the chosen
        ! ones!
        ! S = - sum_i f_i ln(f_i)
        do i=1,2*nwn
            ! neutron  quasiparticles
            if(configmatrix(i).gt.0d0) then
             entropy(1) = entropy(1) - configmatrix(i) * log(configmatrix(i))
            endif
        enddo
        do i=2*nwn+1, 2*nwt
            ! proton quasiparticles
            if(configmatrix(i).gt.0d0) then
              entropy(2) = entropy(2) - configmatrix(i) * log(configmatrix(i))
            endif
        enddo
        ! And a factor of two for time-reversal  
        entropy = 2* entropy
    end select
    
  end subroutine CalcEntropy

  function average_gap_HFB() result(gap)
     !--------------------------------------------------------------------------
     ! Calculation of two types of "average gap", based on 
     !
     ! M. Bender et al., EPJA 8, 59-75 (2000).
     !
     ! Note that this routine is located here, instead of in the HFB module, 
     ! because we need access to the full rho and kappa matrices.
     !--------------------------------------------------------------------------
     ! The original BCS formulation is given by:
     !
     ! <v2 Delta > = sum_k f_k v^2_k   Delta_k / sum f_k**2 v^2_k
     ! <uv Delta > = sum_k f_k u_k v_k Delta_k / sum f_k**2 u_k v_k
     !
     ! Note that the f_k in the original reference is the square of our cutoff!
     !
     ! As in the BCS case, we calculate the average gap here WITHOUT cutoffs, 
     ! as our definition of the gaps already includes all of the cutoff factors
     ! already, in contrast to the EPJA paper. 
     !
     ! We do the summation in the canonical basis, where rho and kappa take 
     ! a simple form. We transform the gaps to this basis, but I would like 
     ! to remark that (due to the presence of cutoffs) they need not be 
     ! canonical in that basis. Kappa however, picks out only the canonical 
     ! part.
     !--------------------------------------------------------------------------

    real(KIND=dp) :: gap(2,2), norm(2,2), v2, uv
    real(KIND=dp), allocatable :: gaps_can(:,:)
    integer       :: it1, wave,i
$NTR integer      :: wavebar

    gap = 0 ; norm = 0
    if(.not.allocated(HFBgaps)) return

    allocate(gaps_can(nwt,nwt)) ; gaps_can = 0.0
    gaps_can = matmul(transpose(cantransfo), HFBgaps)
    gaps_can = matmul(gaps_can, cantransfo)
  
    do wave    =1, nwt
      it1 = 1 ; if(wave .gt.nwn) it1 = 2

$NTR  wavebar = conjugp(wave)

      v2  = rho_can(wave)
$TR      gap(1,it1) = gap(1,it1)  +  v2 * abs(gaps_can(wave,wave))                
$NTR     if(wavebar .ne. 0) then ! The conjugate partner has not necessarily
$NTR                             ! been found, in which case gap = 0 anyway
$NTR       gap(1,it1) = gap(1,it1)  +  v2 * abs(gaps_can(wave,wavebar))                
$NTR     endif      
      norm(1,it1)= norm(1,it1) +  v2                

      uv  = kappa_can(wave)
$TR      gap(2,it1) = gap(2,it1)  +  abs(uv * gaps_can(wave,wave))                 
$NTR      if(wavebar .ne. 0) then
$NTR        gap(2,it1) = gap(2,it1)  +  abs(uv * gaps_can(wave,wavebar))     
$NTR      endif            
      norm(2,it1)= norm(2,it1) +  abs(uv)                                       
    enddo

    do i=1,2
      do it1 =1,2
        gap(i,it1) = gap(i,it1)/norm(i,it1)
        ! numerical safeguard
        if(gap(i,it1) .gt. 1d2) gap(i,it1) = 0.0d0
      enddo
    enddo

    deallocate(gaps_can)
  end function average_gap_HFB
  
  pure function construct_rho_HF(rho, HF_T) result(r_HF)
    !---------------------------------------------------------------------------
    ! Construct the diagonal matrix elements of the density matrix in the 
    ! Hartree-Fock basis.
    !
    ! Input:
    !   rho : density matrix in the "as is" basis
    !   HF_T: Hartree-Fock transformation 
    !
    ! Output:
    !   r_hf: Diagonal matrix elements of rho in the Hartree-Fock basis
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in)  :: rho(:,:), HF_T(:,:)
    real(KIND=dp), allocatable :: r_HF(:)

    integer :: k,l,j
    
    allocate(r_HF(nwt))
    
    do k=1,nwt  
      r_HF(k) = 0.0d0
      do l=1,nwt
        do j=1,nwt
          r_HF(k) = r_HF(k) +  HF_T(l,k) * rho(l,j) * HF_T(j,k)
        enddo
      enddo
    enddo
    
  end function construct_rho_HF
  
  function identify_blocked_particle(Bogo) result (indices)
    !---------------------------------------------------------------------------
    ! Finds the indices in the canonical basis of the levels that obtain v^2 
    ! as a cause of blocking. 
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Input:
    !   Bogo   : Bogoliubov transformation
    ! Output:
    !   indices: indices of the individual levels identified as 'blocked'
    !            in the canonical basis.
    !---------------------------------------------------------------------------
    real(KIND=dp), allocatable, intent(in) :: Bogo(:,:)
    integer, allocatable :: indices(:)
    integer              :: B, si, sb, N, N2, k, ind, j, qp, sp
    real(KIND=dp), allocatable :: U(:), V(:)
    real(KIND=dp)        :: over
    
    allocate(indices(blocknumber)) ; indices= 0
      
    si  = 0
    sb  = 0
    ind = 0
    do B=1,8,2
      N = HFBlocks(B) ; if(N.eq.0) cycle
      N2= HFBlocks(B+1)
      
      do k=1,blocknumber
        qp = blocked_qps(k)
        
        if( qp .gt. si .and. qp .le. si+N+N2) then
          ! This particular blocked qp lives in this symmetry block(s)
          U = Bogo(sb     +1:sb+  N+  N2, sb + N + N2 - qp + si + 1)
          V = Bogo(sb+N+N2+1:sb+2*N+2*N2, sb + N + N2 - qp + si + 1)
          ! These are matrix elements in the HF basis, so we transform to 
          ! the canonical basis
          U = matmul(transpose(cantransfo(si+1:si+N+N2, si+1:si+N+N2)), U)
          V = matmul(transpose(cantransfo(si+1:si+N+N2, si+1:si+N+N2)), V)
            
          over = 0          
          sp   = 0
          do j=1,N+N2
            ! We check for the sp state with the largest V element
            if( V(j)**2 .gt. over .and. abs(rho_can(si+j) - 1.0d0).lt.0.1d0) then
              over = V(j)**2
              sp   = j
            endif 
          enddo
          !
          ind = ind +1 
          indices(ind) = si + sp         
        endif
      enddo    
      si = si +   N +   N2
      sb = sb + 2*N + 2*N2
    enddo
   end function identify_blocked_particle

  subroutine read_modelwf(fname)
      !-------------------------------------------------------------------------
      !
      !-------------------------------------------------------------------------
      logical                       :: exists = .true.
      character(len=40), intent(in) :: fname 
      integer                       :: io, filenx,fileny,filenz,fileit,filepar
      integer                       :: i,j,k,l !, sxh(4), syh(4), szh(4)
      real(KIND=dp)                 :: filedx
      real(KIND=dp), pointer        :: model3d(:,:,:) 
      !real(KIND=dp), allocatable    :: dmodel3d(:,:,:), ddmodel3d(:,:,:)

      1 format (3i3, f8.3, 2i3)
      2 format (99f18.15)


      inquire(file=fname, EXIST = exists)

      if( .not. exists) then
        print *, 'File for model spwf does not exist.'
        stop
      else
        allocate(modelspwf(nx*ny*nz,4,2)) ; modelspwf = 0
        open(unit = 12, file=fname, iostat=io)
        !-----------------------------------------------------------------------
        ! Read the header:
        ! nx ny nz dx it parity 
        read(unit=12, fmt=1) filenx, fileny, filenz,filedx, fileit, filepar
        ! Sanity checks
        if((filenx .ne. nx) .or. &
        &  (fileny .ne. ny) .or. & 
        &  (filenz .ne. nz) .or. &
        &  (filedx .ne. dx)) then
          print *, 'Mesh of the model spwf does not match the calculation.'
          stop
        endif

        !-----------------------------------------------------------------------
        ! Read U(r)
        do l=1,4
          model3d(1:nx, 1:ny, 1:nz) => modelspwf(1:nx*ny*nz,l,1)
          do k=1,nz
            do j=1,ny
              do i=1,nx
               read(unit=12,fmt=2) model3d(i,j,k)
              enddo
            enddo
          enddo
        enddo
        ! Read V(r)
        do l=1,4
          model3d(1:nx, 1:ny, 1:nz) => modelspwf(1:nx*ny*nz,l,2)
          do k=1,nz
            do j=1,ny
              do i=1,nx
               read(unit=12,fmt=2) model3d(i,j,k) 
              enddo
            enddo
          enddo
        enddo

        !-----------------------------------------------------------------------
        ! Some lines of code for checking the correct construction of the 
        ! model spwfs on the mesh.
        !
        !-----------------------------------------------------------------------
        
!       sxh(1) =  1 ; syh(1) = +1 ; szh(1) = -1
!        sxh(2) = -1 ; syh(2) = -1 ; szh(2) = -1 
!        sxh(3) = -1 ; syh(3) = +1 ; szh(3) = +1
!        sxh(4) =  1 ; syh(4) = -1 ; szh(4) = +1

!        allocate(dmodel3d(nx*ny*nz,3,4)) ; dmodel3d = 0.0
!        allocate(ddmodel3d(nx*ny*nz,6,4)) ; ddmodel3d = 0.0

!        call inilag
!        do l=1, 4
!          call derive_tot_1D(modelspwf(:,l,1),sxh(l), syh(l), szh(l), &
!                                             dmodel3d(:,:,l),ddmodel3d(:,:,l))
!        enddo
!        print *
!        print *, 'Jz', &
!                angmom_z_real(modelspwf(:,:,1),modelspwf(:,:,1), dmodel3d) & 
!                &                                /(sum(modelspwf(:,:,1)**2)*dv)
!        do l=1, 4
!          call derive_tot_1D(modelspwf(:,l,2),-sxh(l), syh(l), -szh(l), &
!                &                             dmodel3d(:,:,l),ddmodel3d(:,:,l))
!        enddo
!        print *
!        print *, 'Jz',  & 
!               &  angmom_z_real(modelspwf(:,:,2),modelspwf(:,:,2), dmodel3d) & 
!               &  /(sum(modelspwf(:,:,2)**2)*dv)
!        stop
        !-----------------------------------------------------------------------
        ! Assigning the right blocking blocks
        if(fileit .eq. 1) then
            if (filepar.gt.0) then
              modelblock = 1
            else
              modelblock = 3
            endif            
        else
            if (filepar.gt.0) then
              modelblock = 5
            else
              modelblock = 7
            endif            
        endif

      endif 
  end subroutine read_modelwf

  subroutine clean_pairing()
    !---------------------------------------------------------------------------
    ! Clean up various arrays that might have been allocated to start fresh.
    !---------------------------------------------------------------------------
    if(allocated(rho_pairing))   deallocate(rho_pairing)
    if(allocated(kappa_pairing)) deallocate(kappa_pairing)
    if(allocated(configmatrix))  deallocate(configmatrix)
    if(allocated(Bogoliubov))    deallocate(Bogoliubov)
    if(allocated(QPenergies))    deallocate(QPenergies)
    if(allocated(QPdispersions)) deallocate(QPdispersions)
    if(allocated(CanTransfo))    deallocate(CanTransfo)
    if(allocated(CanCutTransfo)) deallocate(CanCutTransfo)
    if(allocated(BlockIndices))  deallocate(BlockIndices)
    if(allocated(Blocklowest))   deallocate(BlockLowest)

  end subroutine clean_pairing

end module
