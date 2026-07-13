!===============================================================================
!     __  __  ___   ____ ____
!    |  \/  |/ _ \ / ___/ ___|__ _
!    | |\/| | | | | |  | |   / _` |
!    | |  | | |_| | |__| |__| (_| |
!    |_|  |_|\___/ \____\____\__,_|
!
!    Copyright (C) 2026 W. Ryssens and M. Bender
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU Affero General Public License as published
!    by the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU Affero General Public License for more details.
!
!    You should have received a copy of the GNU Affero General Public License
!    along with this program.  If not, see <https://www.gnu.org/licenses/>.
!
!===============================================================================
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

 ! Procedure pointer to the type of pairing cutoff
 procedure(Cut), pointer :: PairingCutoff 
 !------------------------------------------------------------------------------
 ! Cutoff functions
 integer :: CutType = 1
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
    
    !---------------------------------------------------------------------------
    ! We first calculate the cutoffs in the HF basis, where the single-particle
    ! hamiltonian is diagonal
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
  
  real(KIND=dp) function FermiAbove(E, Lambda, it) result(Cutoff)
    !---------------------------------------------------------------------------
    ! Calculates a cutoff that utilises a single Fermi function above the 
    ! Fermi energy. 
    ! 
    !        f^-2 = [1 + exp((  epsilon - lambda - DeltaE)/mu)]
    !
    !    with mu and DeltaE being read from input.
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in) :: E, Lambda
    integer, intent(in)       :: it
    real(Kind=dp)             :: Up
    
    Up   =     (E - Lambda - PairingCut(it))/PairingMu(it)
    Cutoff = sqrt(sqrt(1.0_dp/(1.0_dp + exp(Up))))

    return
  end function FermiAbove
  
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
  
  real(KIND=dp) function SymmetricFermizero(E, Lambda, it) result(Cutoff)
    !---------------------------------------------------------------------------
    ! Calculates a cutoff that utilises two fermi functions, above and below the
    ! Fermi level. In addition, there is a cutoff that is (close to) a Heaviside  
    ! function theta(-E).
    !
    !        f^-2 = [1 + exp((  epsilon - lambda - DeltaE)/mu)]
    !             * [1 + exp((- epsilon + lambda - DeltaE)/mu)]
    !             * [1 + exp((- epsilon )/mu2)]
    ! 
    ! with mu and DeltaE being read from input. and mu2 = 1e-4
    ! The additional cut-off prevents the BCS-gas problem.
    !  To converge when lambda is close or above zero, the cut-off has to be 
    ! 2 MeV above the fermi energy.
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in) :: E, Lambda
    integer, intent(in)       :: it
    real(Kind=dp)             :: Up, Down, Up2, ECut2
    
    Up   =     (E - Lambda - PairingCut(it))/PairingMu(it)
    Ecut2=  max(Lambda+2,0.)
    Up2  =     (E-Ecut2)/1e-4
    Down =   - (E - Lambda + PairingCut(it))/PairingMu(it)
    Cutoff = sqrt(sqrt(1.0_dp/(1.0_dp + exp(Up))))
    Cutoff = Cutoff * sqrt(sqrt(1.0_dp/(1.0_dp + exp(Down))))
    Cutoff = Cutoff * sqrt(sqrt(1.0_dp/(1.0_dp + exp(Up2))))

    return
  end function SymmetricFermizero

  subroutine clean_pairingcutoffs
    if(allocated(Pcutoffs)) deallocate(Pcutoffs)
  end subroutine clean_pairingcutoffs
end module pairingcutoffs
