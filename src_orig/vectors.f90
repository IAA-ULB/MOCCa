module vectors
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
 !  Module containing type declarations for both density- and potential-vectors.
 !
 !  All of these declations are separated from the densities.f90 and 
 !  functional.f90 modules such that they can be employed by lower-level 
 !  modules without much fuss.
 !
 !==============================================================================
 ! Hephaestos keywords:
 !
 ! DECLARATION         : [WAY TOO LONG]
 ! DECLARATION_FIELDS  : [WAY TOO LONG]
 !==============================================================================

 use geninfo

 implicit none

  type DensityVector
    !-------------------------------------------------------------------------
    ! Custom type for containing all Skyrme densities
    !
    !   \mathcal{R} = (D_I_I, D_N_N, ....)
    !
    ! but also their derivatives, as well as the charge density!
$DECLARATION
    !-------------------------------------------------------------------------
    ! Charge density of the protons, possibly including the correction for the
    ! finite size of the proton. It is stored here, as both the moments module
    ! and the coulomb module need it, even though Coulomb depends on the
    ! moments module.
    real(KIND=dp), allocatable :: chargedensity(:,:,:)
  end type DensityVector

  type PotentialVector
    !-------------------------------------------------------------------------
    ! Custom type for containing all potentials
    !
    !   \mathcal{F} = (F_I_I, F_N_N, ....)
$DECLARATION_POTENTIALS
    !-------------------------------------------------------------------------
    ! Coulomb fields
    ! Attention: CoulombPotential is defined on an EXTENDED box, see the
    !            coulomb module
    real(KIND=dp), allocatable :: CoulombPotential(:,:,:)
    real(KIND=dp), allocatable :: ExchangePotential(:,:,:)
    ! Array with the folded Coulomb potential, necessary if we take the
    ! finite size of the nucleons into account
  real(KIND=dp), allocatable :: FoldedCoul(:,:,:,:), FoldedExchange(:,:,:,:)

  end type PotentialVector

end module vectors


