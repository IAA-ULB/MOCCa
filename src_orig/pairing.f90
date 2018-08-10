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
 ! Transformation from the HFBasis into the canonical basis
 real(KIND=dp), allocatable :: CanTransfo(:,:)
 !------------------------------------------------------------------------------
 ! Fermi energy for neutrons and protons.
 real(KIND=dp) :: FermiEnergy(2) 
 !------------------------------------------------------------------------------
 !
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
 
contains

  subroutine initpairing
    !---------------------------------------------------------------------------
    ! Read and initialize pairing options. 
    !
    !---------------------------------------------------------------------------
    character(len=20) :: Type
    
    NameList /Pairing/ Type, CutType, Constantgap 
    read(unit=*, NML=Pairing)
  
    Type = to_upper(Type)
    
    if('HF' .eq.adjustl(type)) then
      pairingtype = 0
    elseif('BCS' .eq. adjustl(type)) then
      pairingtype = 1
    elseif('HFB' .eq. adjustl(type)) then
      pairingtype = 2
    else
      print *, 'This type of pairing is not implemented yet.'
    endif
    
    !---------------------------------------------------------------------------
    ! Cutoff decision
    select case(CutType)
    case(1)
       PairingCutoff => SymmetricFermi
    case(2)
       PairingCutoff => CosineCut
    end select
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

      ! Determine the size of the HFB matrices
      call initHFB()  
      HFBGaps = 0.0
      si = 0 
      do B= 1,8
        N = HFBlocks(B)
        do wave=si+1,si+N
          do wave2=wave+1,si+N
            HFBgaps( wave, wave2) = 1.0
            HFBgaps(wave2, wave)  =-1.0
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
      allocate(kappa_can(nwt)) ; kappa_can  = 0.0 
    endif
    
    select case (Pairingtype)
    case(0)
      call NaiveFill(rho_can)
    case(1)
      !-------------------------------------------------------------------------
      ! BCS-type pairing
      ! The diagonal elements of rho and kappa are only set. 
      call solvepairing_BCS(FermiEnergy, rho_can, kappa_can)
      
    case(2)
      !-------------------------------------------------------------------------
      ! HFB-type pairing
      if(.not.allocated(CanTransfo)) then
        ! Allocate the full matrices
        allocate(CanTransfo(nwt, nwt))     ; CanTransfo    = 0.0
        allocate(rho_pairing(nwt,nwt))     ; rho_pairing   = 0.0
        allocate(kappa_pairing(nwt,nwt))   ; kappa_pairing = 0.0
      endif
      ! Find the Fermi energy
      call solvepairing_HFB(FermiEnergy, rho_pairing, kappa_pairing)
      
      ! Find the transformation to the canonical basis
      call Canonical(rho_pairing, kappa_pairing, rho_can, kappa_can, cantransfo)
      
      ! Apply this transformation
      call ConstructCanonicalBasis(cantransfo, rho_can)
    end select
    
    ! Compute the cutoffs
    call ComputePairingCutoffs(fermienergy)

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
    7 format (60('-'))

    print 1

    select case(PairingType)
    case (0)
        return
    case (1,2)
        ! BCS and HFB
        print 2
        print 3, FermiEnergy
!        print 4, sum(D_I_I(:,1))*dv, sum(D_I_I(:,2))*dv
    end select
    print 7
  end subroutine PrintPairing
  
  function calcpairingenergy() result(E)
    !---------------------------------------------------------------------------
    ! Calculate the pairingenergy
    !
    !---------------------------------------------------------------------------
    integer       :: wave, it
    real(KIND=dp) :: E(2)
    
    E = 0.0
    
    select case(PairingType) 
    case(0)
    
    case(1)
      ! BCS case
      do wave=1,nwt
          it = 1
          if(wave .gt. nwn) it =2
          E(it) = E(it) - BCSgaps(wave)*Kappa_can(wave)
      enddo
    case(2)
    
    end select
  end function calcpairingenergy

end module
