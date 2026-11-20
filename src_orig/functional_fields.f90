submodule (functional) functional_fields
!==============================================================================
!
! Implementation submodule of module functional: functional_fields
!  Interface block lives in src_orig/functional.f90.
!
!==============================================================================
implicit none
contains

function Add_potentialvector(F1, F2) result(F)
    !---------------------------------------------------------------------------
    ! Sum two potential-vectors F1 and F2 to get a new potential-vector F, so
    ! that F_a = F1_a + F2_a for every potential with latin index a.
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !     F1, F2 : potential-vectors
    ! Output:
    !     F      : new potential vector
    !---------------------------------------------------------------------------
    type(PotentialVector), intent(in) :: F1, F2
    type(PotentialVector)             :: F

    ! Automatically generated expressions for all Skyrme potentials
$INIPOTENTIALS
$ADD_POTENTIALS

    ! Additions for the Coulomb potentials are always executed
    F%CoulombPotential  = F1%CoulombPotential  + F2%CoulombPotential
    F%ExchangePotential = F1%ExchangePotential + F2%ExchangePotential
  
    if(allocated(F1%Foldedcoul)) then
        F%Foldedcoul      = F1%Foldedcoul     + F2%Foldedcoul
        F%Foldedexchange  = F1%Foldedexchange + F2%Foldedexchange
    endif
end function Add_potentialvector


function multiply_potentialvector(a, F1) result(F)
    !---------------------------------------------------------------------------
    ! Multiply a potential-vector F1 with a scalar a such that 
    ! F_b = a * F1_b for every potential with latin index b.
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !     a      : scalar 
    !     F1     : potential-vector
    ! Output:
    !     F      : new potential vector
    !---------------------------------------------------------------------------
    type(PotentialVector), intent(in) :: F1
    real(KIND=dp), intent(in)         :: a
    type(PotentialVector)             :: F

$INIPOTENTIALS
$MULTIPLY_POTENTIALS

    F%CoulombPotential  = a*F1%CoulombPotential 
    F%ExchangePotential = a*F1%ExchangePotential
 
    if(allocated(F1%Foldedcoul)) then
      F%Foldedcoul      = a*F1%Foldedcoul     
      F%Foldedexchange  = a*F1%Foldedexchange 
    endif
end function multiply_potentialvector

!  function calcRouth_onthefly(Rin, Fin) result(routh)
!     !---------------------------------------------------------------------------
!     ! This routine calculates the total Routhian for a given set of densities.
!     !
!     ! Note:
!     !   (*) It only includes contributions that are included in the calculation
!     !       in a variational way, i.e. any rotational correction or two-body
!     !       COM-correction is not included.
!     !
!     !   (*) This routine has NO side-effects, i.e. it calculates the total
!     !       routhian and ONLY the total Routhian. No module-level variables
!     !       are set.
!     !
!     ! ATTENTION: testing routine, results not guaranteed!
!     !
!     !
!     !---------------------------------------------------------------------------
!     use Coulombmod
!     use Densities
!
!     real(KIND=dp)                     :: routh
!     type(DensityVector), intent(in)   :: Rin
!     type(PotentialVector), intent(in) :: Fin
!
!     real(KIND=dp) :: K(2), TE, TO, PE(2), S, C, CE, COM(2,2)
!
!     routh = 0.0d0
!
!     ! Kinetic energy
!     K = CompKinetic_density(Rin)
!     ! Note that this estimate of the energy does NOT include the 2-body COM
!     ! To save CPU cycles, we override its calculation
! #if(PASTA == 0)
!     call compComCorrection(K, .true., COM)
! #else
!     ! A waste of CPU time for pasta calculations
!     COM = 0.0d0
! #endif
!     ! Skyrme energy
!     call compSkyrme(Rin, Fin, S, TE, TO, PE)
!     ! Coulomb energies
!     C  = CoulombEnergy_Direct(Rin,Fin)
!     CE = CoulombEnergy_Exchange(Rin)
!
!     ! The total energy is comprised of
!     !      Kinetic part + Skyrme part + corrections + Coulomb energy
!     routh =sum(S + K) + sum(COM) + C + CE
!     ! Note again: no 2-body COM or rotational correction in here
!
!     ! Addition of the various constraints
!     ! TODO: correct this!
!     !routh = routh  + sum(Constraint_I_I * Rin%D_I_I)*dv/2.0_dp
!     ! Missing = cranking contribution and correct multipole contribution
!
!
!  end function calcRouth_onthefly
 

#if( $FAM == 0)
  function PVectorInproduct(F1, F2) result(x)
    !---------------------------------------------------------------------------
    ! Define a basic inproduct on the space of the potential vectors.
    ! 
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Input: 
    !   F1, F2 : potentialvectors to calculate the inproduct of. 
    !       
    ! Output: 
    !      x: inproduct value, i.e. < F_1 | F_2 >
    !---------------------------------------------------------------------------
    type(PotentialVector), intent(in) :: F1,F2
    real(KIND=dp) :: x

    x = 0.0d0
$PVECTORINPRODUCT
    
 end function PVectorInproduct
#endif


end submodule functional_fields
