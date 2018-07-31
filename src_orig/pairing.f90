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
 use pairingcutoffs
  
 implicit none
 
 !------------------------------------------------------------------------------
 ! Pairing density matrix and anomalous density matrix in the HF basis. 
 real(KIND=dp), allocatable :: rho_pairing(:,:), kappa_pairing(:,:)
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

 
contains

  subroutine initpairing
    !---------------------------------------------------------------------------
    ! Read and initialize pairing options. 
    !
    !---------------------------------------------------------------------------
    character(len=20) :: Type
    
    NameList /Pairing/ Type, CutType
    read(unit=*, NML=Pairing)
  
    Type = to_upper(Type)
    
    if('HF' .eq.adjustl(type)) then
      pairingtype = 0
    elseif('BCS' .eq. adjustl(type)) then
      pairingtype = 1
    else
      print *, 'This type of pairing is not implemented yet.'
    endif
    
    select case(CutType)
    case(1)
       PairingCutoff => SymmetricFermi
    case(2)
       PairingCutoff => CosineCut
    end select
    
  end subroutine initpairing
  
  subroutine SolvePairing
    !---------------------------------------------------------------------------
    ! Master routine for the solving of the pairing equations.
    !---------------------------------------------------------------------------
    integer :: wave
    
    if(.not.allocated(rho_pairing)) then
      allocate(rho_pairing(nwt,nwt))     ; rho_pairing   = 0.0
      allocate(kappa_pairing(nwt,nwt))   ; kappa_pairing = 0.0
      allocate(occupations(nwt))         ; occupations   = 0.0
    endif
    
    select case (Pairingtype)
    case(0)
      call NaiveFill()
    case(1)
    
      !-------------------------------------------------------------------------
      ! BCS-type pairing
      ! The diagonal elements of rho and kappa are only set. 
      call solvepairing_BCS(FermiEnergy, rho_pairing, kappa_pairing)
      ! The density matrices need not be changed, HFBasis = canonical basis.
      do wave=1,nwt
        occupations(wave) = rho_pairing(wave,wave)
      enddo 
    end select
    
  end subroutine SolvePairing
  
  subroutine printpairing
    !---------------------------------------------------------------------------
    !
    !
    !
    !---------------------------------------------------------------------------
    
    use densities
    
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
        print 4, sum(D_I_I(:,1))*dv, sum(D_I_I(:,2))*dv
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
          E(it) = E(it) - 0.5*BCSgaps(wave)*Kappa_pairing(wave,wave)
      enddo
    case(2)
    
    end select
  end function calcpairingenergy

end module
