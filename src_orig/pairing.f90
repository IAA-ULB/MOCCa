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
 !
 ! The overall strategy is to delegate as much work as possible to the specific
 ! different pairing solvers, leaving only the truely 'global' quantities here.
 ! These are: 
 !  a) The density matrix rho_pairing
 !  b) The anomalous density kappa_pairing
 !  c) The pairing cutoffs
 !  d) The Fermi energy
 !  e) The dispersion 
 ! All of the rest, including the U and V matrices and all other possible 
 ! quantities are not truly 'global'. It is only these quantities that determine
 ! the final many-body auxiliary state.
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
 real(KIND=dp), allocatable :: rho_pairing(:,:), kappa_pairing(:,:)
 real(KIND=dp), allocatable :: rho_can(:), kappa_can(:)
 real(KIND=dp), allocatable :: configmatrix(:)
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
 real(KIND=dp) :: FermiEnergy(2) 
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
 ! (1) ordinary blocking
 ! (2) EFA blocking. 
 ! Signals the code to look for a new namelist "Indices"
 !
 ! This currently only works for HFB calculations, but could be envisaged for 
 ! BCS calculations as well. 
 !------------------------------------------------------------------------------
 integer :: BlockType  = 0
 integer :: BlockNumber= 0 
 !------------------------------------------------------------------------------
 ! Indices of the levels to block. 
 integer, allocatable :: BlockIndices(:) 

 !------------------------------------------------------------------------------
 ! Entropy of the statistical mixture in the case of finite temperature
 real(KIND=dp) :: entropy(2) = 0
contains

  subroutine initpairing
    !---------------------------------------------------------------------------
    ! Read and initialize pairing options. 
    !
    !---------------------------------------------------------------------------
    character(len=20) :: Type = 'HF'
    
    NameList /Pairing/ Type, CutType, Constantgap, hfbmix, hfbmixtype,         &
    &                  BlockType, BlockNumber, cutneutron, cutproton
    NameList /Indices/ BlockIndices

    read(unit=*, NML=Pairing)
  
    Type = to_upper(Type)
    
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
    
    if(Blocktype.lt.0 .or. BlockType.gt.2) then
        print *, 'This value of BlockType is not accepted.'
        stop
    endif
    !---------------------------------------------------------------------------
    ! Reading information on the blocking if needed.
    if(BlockNumber.ne.0) then
        allocate(BlockIndices(BlockNumber)) ; BlockIndices = 0
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
    1 format(80('-'))
    2 format(' Pairing treatment: ', a30)
    3 format('   Linear mixing of (rho,kappa)')    
    4 format('   Linear mixing of eigenvalues of R')
    5 format('   HFBmix = ', f5.3)
    6 format('   Cutoff parameters  = ', a20)
    7 format('     dE (n,p) = ', 2f4.1, ' MeV ')
    8 format('     mu (n,p) = ', 2f4.1, ' MeV ')

    print 1

    select case (pairingtype)
    case(0)
        print 2, 'Hartree-Fock (HF)'
    case(1)
        print 2, 'Bardeen-Cooper-Schrieffer (HF+BCS)'
    case(2)
        print 2, 'Hartree-Fock-Bogoliubov (HFB)'
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
    integer :: wave
    
    if(.not.allocated(rho_can)) then
      allocate(rho_can(nwt))   ; rho_can    = 0.0
    endif
    if(.not.allocated(kappa_can)) then
      allocate(kappa_can(nwt)) ; kappa_can  = 0.0 
    endif
    if(.not.allocated(qpenergies)) then
        allocate(QPenergies(nwt)) ; qpenergies = 0.0
    endif
    ! Always allocate the configuration matrix C
    if(.not.allocated(configmatrix)) then
       allocate(configmatrix(2*nwt))      ;  configmatrix = 0.0
    endif
    
    select case (Pairingtype)
    case(0)
      call NaiveFill(rho_can)
    case(1)
      !-------------------------------------------------------------------------
      ! BCS-type pairing
      ! The diagonal elements of rho and kappa are only set. 
      call solvepairing_BCS(FermiEnergy, rho_can, kappa_can, qpenergies)

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

      if(all(HFBsizes.eq.0)) call inithfb

      !-------------------------------------------------------------------------
      ! Find the Fermi energy
      call solvepairing_HFB(FermiEnergy, rho_pairing, kappa_pairing,           &
      &    configmatrix, qpenergies,HFBmix, HFBmixtype, BlockType, Blockindices)
          
    end select
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
    3 format (' Fermi Level (MeV) ',2x,f10.5,2x,f10.5)
    4 format (' Particles         ',2x,f10.5,2x,f10.5)
    5 format (' Dispersion        ',2x,f10.5,2x,f10.5)
    7 format (60('-'))

    print 1

    select case(PairingType)
    case (0)
        return
    case (1,2)
        ! BCS and HFB
        print 2
        print 3, FermiEnergy
        print 4, sum(rho_can(1:nwn)), sum(rho_can(nwn+1:nwt))
        select case(PairingType)
        case(1)
            print 5, BCSdispersion
        case(2)
            print 5, HFBdispersion
        end select
    end select
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
      do wave=1,nwt
          it = 1
          if(wave .gt. nwn) it =2
          E(it) = E(it) - BCSgaps(wave)*Kappa_can(wave)
      enddo
    case(2)
      do wave=1,nwt
        do wave2=wave,nwt
          it = 1
          if(wave .gt. nwn) it =2
          ! Factor of two due to symmetricity
          E(it) = E(it) + Kappa_pairing(wave2,wave)*HFBgaps(wave,wave2)
        enddo
      enddo
    end select
  end function calcpairingenergy

  subroutine CalcEntropy()
    !---------------------------------------------------------------------------
    ! Calculate the entropy associated with a given statistical mixture at
    ! finite temperature.
    !---------------------------------------------------------------------------

    integer :: i

    select case(PairingType)
    case(0)
        print *, 'Entropy calculation not yet implemented for HF.' 
        stop
    case (1)
        print *, 'Entropy calculation not yet implemented for BCS.'
    case (2)
        entropy = 0
        ! S = - sum_i f_i ln(f_i)
        do i=1,2*nwn
            entropy(1) = entropy(1) - configmatrix(i) * log(configmatrix(i))
        enddo
        do i=2*nwn+1, 2*nwt
            entropy(2) = entropy(2) - configmatrix(i) * log(configmatrix(i))
        enddo
    end select
    
  end subroutine CalcEntropy


end module
