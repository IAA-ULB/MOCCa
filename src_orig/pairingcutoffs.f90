module pairingcutoffs
!===============================================================================
! Cutoff routines
!===============================================================================

  use wavefunctions
  
  abstract interface
    real(KIND=dp) function Cut(E, Lambda, it) result(Cutoff)
      import                    :: dp
      real(KIND=dp), intent(in) :: E, Lambda
      integer, intent(in)       :: it
    end function Cut
  end interface

 procedure(Cut), pointer :: PairingCutoff 

 !------------------------------------------------------------------------------
 ! Parameters for the symmetric Fermi function cutoff
 real(KIND=dp) :: cutneutron = 5, cutproton= 5.
 real(KIND=dp) :: PairingCut(2) = 0,  PairingMu(2) =0.5
 !------------------------------------------------------------------------------
 ! Storage for all cutoffs. 
 ! Currently these are ALWAYS in the HF basis.
 real(KIND=dp), allocatable :: PCutoffs(:)

contains

  subroutine ComputePairingCutoffs(Lambda)
  !-----------------------------------------------------------------------------
  ! Computes and stores all pairing cutoffs for further use, as a function of 
  ! the Fermi energy
  !-----------------------------------------------------------------------------
    integer :: wave, it
    real(KIND=dp), intent(in) :: Lambda(2)
    
    if(.not.allocated(PCutoffs)) then
        allocate(PCutoffs(nwt))
    endif
    
    do wave=1,nwt
        it = 1
        if(wave .gt. nwn) it = 2
        PCutoffs(wave) = PairingCutoff(spenergies(wave), Lambda(it), it)         
    enddo
  
  end subroutine ComputePairingCutoffs

  real(KIND=dp) function SymmetricFermi(E, Lambda, it) result(Cutoff)
    !---------------------------------------------------------------------------
    ! Calculates a cutoff that utilises two fermi functions, above and below the
    ! Fermi level.
    ! Cutoff is taken with a fermi function both above and below the fermi 
    !    energy.
    !        f^-2 = [1 + exp((  epsilon - lambda - DeltaE)/mu)]
    !             * [1 + exp((- epsilon + lambda - DeltaE)/mu)]
    !    with mu and DeltaE being read from input.
    !    This is described in S.J. Krieger et al., Nucl. Phys.A517 (1990) 275
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in) :: E, Lambda
    integer, intent(in)       :: it
    real(Kind=dp)             :: Up, Down
    
    Up   =     (E - Lambda - PairingCut(it))/PairingMu(it)
    Down =   - (E - Lambda + PairingCut(it))/PairingMu(it)
    Cutoff = sqrt(sqrt(1.0_dp/(1.0_dp + exp(Up))))
    Cutoff = Cutoff * sqrt(sqrt(1.0_dp/(1.0_dp + exp(Down))))

    return
  end function SymmetricFermi
  
  real(KIND=dp) function CosineCut(E, Lambda, it) result(Cutoff)
    !---------------------------------------------------------------------------
    ! Calculates a cutoff using a cosine function.
    ! 
    ! Define:
    ! xd = DeltaE - mu/2
    ! xu = DeltaE + mu/2
    !
    ! then
    ! f = 1                                     
    !   when abs(Epsilon - Lambda) <= xd
    ! f = 0.5 * cos ( (abs(Epsilon - lambda) -xd )*pi/2) + 0.5 
    !   when xd <= abs(Epsilon - Lambda) <= xu
    ! f = 0.0  
    !   when xu <= abs(Epsilon - Lambda)
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in) :: E, Lambda
    integer, intent(in)       :: it
    real(Kind=dp)             :: Up, Down
    
    Down = PairingCut(it) - PairingMu(it)/2.0_dp
    Up   = PairingCut(it) + PairingMu(it)/2.0_dp
  
    if( abs(E - Lambda).le. Down) then
      Cutoff = 1.0_dp      
    elseif( abs(E-Lambda) .le. Up) then
      Cutoff = 0.5_dp * cos((abs(E-Lambda) - Down)*pi/PairingMu(it)) + 0.5
    else
      Cutoff = 0.0_dp
    endif
  
  end function Cosinecut
  
end module pairingcutoffs
