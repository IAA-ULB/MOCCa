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

 use compilation
 use wavefunctions
 use hartreefock
 use BCS
 use HFB
 use pairingcutoffs
  
 implicit none
 
 !------------------------------------------------------------------------------
 ! Pairing density matrix and anomalous density matrix in the HF basis. 
 real(KIND=dp), allocatable :: rho_pairing(:,:),  kappa_pairing(:,:)
 ! ... and in the canonical basis ...
 real(KIND=dp), allocatable :: rho_can(:), kappa_can(:)
 ! Note that the object kappa_can only has an effect on the calculation 
 ! in the BCS case, when kappa is actually off-diagonal. In the HFB case, no
 ! guarantees are given as to the canonical form of kappa, and the code does not
 ! rely on it being this way.
 ! ...  and the "configuration matrix" ...
 real(KIND=dp), allocatable :: configmatrix(:)
 ! ... and finally, the Bogliubov transformation.
 real(KIND=dp), allocatable :: Bogoliubov(:,:)
 !------------------------------------------------------------------------------
 ! Quasiparticle excitation energies, either HF, BCS or HFB.
 real(KIND=dp), allocatable :: QPenergies(:)
 !------------------------------------------------------------------------------
 ! Transformation from the HFBasis into the canonical basis
 real(KIND=dp), allocatable :: CanTransfo(:,:)
 ! Transformation from the HFBasis into the basis where Kappa (with cutoffs)
 ! is canonical.
 real(KIND=dp), allocatable :: CanCutTransfo(:,:)
 !------------------------------------------------------------------------------
 ! Fermi energy for neutrons and protons.
 real(KIND=dp) :: FermiEnergy(2) , FermiHistory(2)
 !------------------------------------------------------------------------------
 ! Particle number dispersion
 real(KIND=dp) :: Dispersion(2)
 !------------------------------------------------------------------------------
 ! Cutoff functions
 integer :: CutType = 1
 !------------------------------------------------------------------------------
 ! Type of pairing to employ.
 ! (0), (1), (2)
 integer :: PairingType = 0
 !------------------------------------------------------------------------------
 ! Decide which Fermisolver to use. 
 ! "Brent"  => use a modified bisection solver
 ! "Secant" => use a secant routine
 character(len=99) :: FermiSolver='Brent'
 !------------------------------------------------------------------------------
 ! Decide which module gets to calculate the pairing gaps.
 procedure(calcBCSgaps), pointer :: CalcGaps
 !------------------------------------------------------------------------------
 ! Mixing parameter for the HFB equations
 real(KIND=dp) :: HFBMix = 1.0
 !------------------------------------------------------------------------------
 ! Mixtype HFB
 ! 0 => linear mixing of rho and kappa
 ! 1 => linear mixing of the generalized density
 integer       :: HFBmixtype = 0
 !------------------------------------------------------------------------------
 ! Type of blocking we want. 
 ! (0) no blocking
 ! (1) ordinary blocking, based on indices
 ! (2) ordinary blockgin, asking for lowest energy configurations
 ! (3) EFA blocking, based on indices.
 ! (4) EFA blocking, asking for lowest energy configurations.
 ! 
 ! If this is nonzero, the code will look for a new namelist "Indices"
 !
 ! This currently only works for HFB calculations!
 !------------------------------------------------------------------------------
 integer :: BlockType  = 0
 integer :: BlockNumber= 0 
 !------------------------------------------------------------------------------
 ! Indices of the levels to block. 
 integer, allocatable :: BlockIndices(:) 
 ! Indices of the quasi-particles that ended up blocked.
 integer, allocatable :: blocked_qps(:)
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
 integer       :: particles_in_gas = 0
 real(KIND=dp) :: ngas(2)


contains

  subroutine initpairing(file_number)
    !---------------------------------------------------------------------------
    ! Read and initialize pairing options. 
    !
    !---------------------------------------------------------------------------
    character(len=20) :: Type = 'HF'
    integer(dp), intent(in), optional   :: file_number   
    
    NameList /Pairing/ Type, CutType, Constantgap, hfbmix, hfbmixtype,         &
    &                  BlockType, BlockNumber, cutneutron, cutproton,          &
    &                  particles_in_gas, maxhfbiter, FermiSolver

    NameList /Indices/ BlockIndices, blocklowest

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

    FermiSolver = to_upper(FermiSolver)    
    if(adjustl(FermiSolver).eq.'SECANT') then
      FindFermi => FindFermi_secant
    elseif(adjustl(FermiSolver).eq.'BRENT') then
      FindFermi => FindFermi_brent
    else
      print *, 'Unknown FermiSolver', FermiSolver, ' selected.'
      stop
    endif
    
    if(Blocktype.lt.0 .or. BlockType.gt.4) then
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
        allocate(BlockIndices(BlockNumber)) ; BlockIndices = 0
        allocate(BlockLowest(BlockNumber))  ; BlockLowest  = ' ' 
        read(unit=*, nml=Indices)
    endif

    !---------------------------------------------------------------------------
    ! Cutoff decision
    select case(CutType)
    case(1)
       PairingCutoff => SymmetricFermi
    case(2)
       PairingCutoff => CosineCut
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
  end subroutine initpairing

  subroutine printpairing_init
    !---------------------------------------------------------------------------
    ! Print information on the treatment of pairing at the start of the run.
    !---------------------------------------------------------------------------
    1 format(80('-'))
    2 format(' Pairing treatment: ', a60)
   21 format('   Fermi-solver: ', a99 )
    3 format('   Linear mixing of (rho,kappa)')    
    4 format('   Linear mixing of eigenvalues of R')
    5 format('   HFBmix = ', f5.3)
    6 format('   Cutoff parameters  = ', a20)
    7 format('     dE (n,p) = ', 2f4.1, ' MeV ')
    8 format('     mu (n,p) = ', 2f4.1, ' MeV ')

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


    print 1

    select case (pairingtype)
    case(0)
        print 2, 'Hartree-Fock (HF)'
    case(1)
        print 2, 'Bardeen-Cooper-Schrieffer (HF+BCS)'
    case(2)
        print 2, 'Hartree-Fock-Bogoliubov (HFB)'
        print 21, adjustl(FermiSolver)
    end select

    if(pairingtype.eq.2) then
        if(HFBmixtype.eq.0) then
            print 3
        elseif(HFBmixtype.eq.1) then
            print 4
        endif 
        print 5, HFBmix
    endif    

    select case(CutType)
    case(1)
       print 6, 'Symmetric Fermi'
    case(2)
       print 6, 'Cosine'
    end select

    print 7, pairingcut
    print 8, PairingMu

    if(fixfermi) then
        print 100, mun, mup
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
        end select
    endif

  end subroutine printpairing_init
  
  subroutine GuessGaps()
    !---------------------------------------------------------------------------
    ! 
    !
    !
    !---------------------------------------------------------------------------
    integer :: wave, wave2, si, B, N

    
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
      ! Simply put 1.0 
      do wave=1,nwt
        BCSgaps(wave) = 1.0 
      enddo
    case(2)
      !-------------------------------------------------------------------------
      ! HFB calculations
      if(.not.allocated(HFBGaps)) then
        ! Factor of 2 through time-reversal
        allocate(HFBGaps(2*nwt,2*nwt)) ; HFBGaps = 0.0
      endif
  
      HFBGaps = 0.0
      si = 0 
      do B= 1,8
        N = HFBlocks(B)
        do wave=si+1,si+N
          do wave2=si+1,si+N
            HFBgaps( wave, wave2) = 1.0
          enddo 
        enddo
        si = si + N
      enddo
    end select  
  end subroutine GuessGaps
  
  subroutine SolvePairing
    !---------------------------------------------------------------------------
    ! Master routine for the solving of the pairing equations.
    !---------------------------------------------------------------------------
    use parameterization, only : hbm

    if(.not.allocated(rho_can)) then
      allocate(rho_can(nwt))              ; rho_can    = 0.0
    endif
    if(.not.allocated(kappa_can)) then
      allocate(kappa_can(nwt))            ; kappa_can  = 0.0 
    endif
    if(.not.allocated(qpenergies)) then
        allocate(QPenergies(nwt))         ; qpenergies = 0.0
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
    case(1)
      !-------------------------------------------------------------------------
      ! BCS-type pairing
      ! The diagonal elements of rho and kappa are only set. 
      call solvepairing_BCS(FermiEnergy, rho_can, kappa_can, qpenergies,       &
      &                     particles_in_gas)

      ! Calculate the average gap
      average_gap = average_gap_BCS()
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

      if(.not. allocated(HFBsizes)) call inithfb

      !-------------------------------------------------------------------------
      ! Find the Fermi energy
      call solvepairing_HFB(FermiEnergy, Bogoliubov,rho_pairing, kappa_pairing,&
      &                     configmatrix, qpenergies, HFBmix, HFBmixtype,      &
      &                     BlockType, Blockindices, blocklowest, blocked_qps)
    end select

    !---------------------------------------------------------------------------
    ! If beta != infty, we calculate the number of particles in the gas
    if(inversetemp.ne.-1) then
      ngas(1) = gasoccupations(fermienergy(1), hbm(1))
      ngas(2) = gasoccupations(fermienergy(2), hbm(2))
    endif

    !---------------------------------------------------------------------------
    ! Compute the cutoffs
    call ComputePairingCutoffs(fermienergy)
    !---------------------------------------------------------------------------

  end subroutine SolvePairing

  subroutine printpairing
    !---------------------------------------------------------------------------
    !
    !
    !
    !---------------------------------------------------------------------------

    1 format (26('-'), ' Pairing ', 25('-'))
    2 format (25x, ' N ',7x, ' P ')
    3 format (' Fermi Level (MeV) ',2x,f13.8,2x,f13.8)
    4 format (' Particles         ',2x,f13.8,2x,f13.8)
    5 format (' Dispersion        ',2x,f13.8,2x,f13.8)
!    6 format (' dN/da             ',2x,f13.8,2x,f13.8,/,                       & 
!    &         ' dZ/da             ',2x,f13.8,2x,f13.8 )

    6 format (' Average gap   v^2 ',2x, f13.8, 2x, f13.8,/,                    &
              ' Average gap   uv  ',2x, f13.8, 2x, f13.8)
    7 format (60('-'))

    8 format ('  gas-like         ', 2x, f13.8, 2x, f13.8)
    9 format ('  nucleus          ', 2x, f13.8, 2x, f13.8)

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
            print 6, average_gap
        case(2)
            print 5, HFBdispersion
            call PrintHFBConvergence(rho_pairing, kappa_pairing)
        end select
    end select

!    if(inversetemp.ne.-1) then
!        call EstimateDNDA()
!        print 6, dNda
!    endif
    
    print 7
  end subroutine PrintPairing
  
  function calcpairingenergy() result(E)
    !---------------------------------------------------------------------------
    ! Calculate the pairingenergy
    !
    !---------------------------------------------------------------------------
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
      do wave = 1,nwt
        do wave2 = 1,nwt
          it = 1
          if(wave .gt. nwn) it = 2
          E(it) = E(it) - Kappa_pairing(wave,wave2)*HFBgaps(wave,wave2)
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

  subroutine clean_pairing()

    if(allocated(rho_pairing))   deallocate(rho_pairing)
    if(allocated(kappa_pairing)) deallocate(kappa_pairing)
    if(allocated(configmatrix))  deallocate(configmatrix)
    if(allocated(Bogoliubov))    deallocate(Bogoliubov)
    if(allocated(QPenergies))    deallocate(QPenergies)
    if(allocated(CanTransfo))    deallocate(CanTransfo)
    if(allocated(CanCutTransfo)) deallocate(CanCutTransfo)
    if(allocated(BlockIndices))  deallocate(BlockIndices)
    if(allocated(Blocklowest))   deallocate(BlockLowest)

  end subroutine clean_pairing

end module
