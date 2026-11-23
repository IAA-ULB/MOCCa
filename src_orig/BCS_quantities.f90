submodule (BCS) BCS_quantities
!==============================================================================
!
! Implementation submodule of BCS: quantities derived from the BCS solution
! (quasiparticle energies, occupations, dispersion, average gaps, ordering).
! The interfaces live in module BCS (src_orig/BCS.f90).
!
!==============================================================================

 use wavefunctions
 use pairingcutoffs

 implicit none

contains

 module subroutine BCSQPEnergies(Fermi)
   !----------------------------------------------------------------------------
   ! This subroutine calculates the BCS quasiparticle energies as a function of
   ! the Fermi energies and the pairing gaps.
   !
   ! eqp_k = sqrt[(epsilon_k - \lambda)**2 + Delta_{k, kbar}**2]
   !
   ! which is formula (6.72) on page 235 in Ring & Schuck.
   !----------------------------------------------------------------------------
    real(KIND=dp), intent(in) :: Fermi(2)
    integer                   :: wave, it
    real(KIND=dp)             :: epsilon, lambda

    if(.not.allocated(BCSqps)) then
      allocate(BCSqps(nwt)) ; BCSqps=0
    endif

    do wave=1,nwt
        it = 1
        if(wave .gt. nwn) it = 2

        epsilon = spenergies(wave)
        lambda  = Fermi(it)
        BCSqps(wave)  = sqrt((epsilon - lambda)**2 + abs(BCSgaps(wave))**2)
    enddo

  end subroutine BCSQPEnergies

  module subroutine calcBCSOccupations(Fermi)
   !----------------------------------------------------------------------------
   ! Find the occupation numbers of the HFBasis from the quasiparticle energies
   ! and the Fermi energies.
   !
   ! v^2_k = 0.5 * (1 - (Epsilon - Lambda)/(E_{qp}))
   !
   ! or formula (6.51) on page 231 in Ring & Schuck. Note that these are only
   ! entries of the density matrix in the case of zero-temperature calculations.
   !----------------------------------------------------------------------------
    real(KIND=dp), intent(in)   :: Fermi(2)
    integer                     :: wave,it
    real(KIND=dp)               :: eqp

    if(.not.allocated(BCSoccupations)) then
      allocate(BCSoccupations(nwt)) ; BCSoccupations = 0.0d0
    endif

    do wave=1,nwt
      it  = 1
      if(wave .gt. nwn) it = 2

      eqp = BCSqps(wave)

      ! Zero-temperature BCS
      BCSOccupations(wave) = 0.5d0*(1 - (spenergies(wave) - Fermi(it))/eqp)
    enddo
    ! Time reversal symmetry
    BCSOccupations = 2 * BCSOccupations

  end subroutine calcBCSOccupations

   module subroutine calcBCSdispersion(rho_can, kappa_can)
      !-------------------------------------------------------------------------
      ! Calculate the dispersion
      !  Tr rho - rho^2
      !-------------------------------------------------------------------------
      ! We calculate the dispersion of the particle number
      !
      !  <N^2> = sum_{ab} <a^{\dagger}_{a} a_{a} a^{\dagger}_{b} a_{b} >
      !        = sum_{ab} rho_{aa} rho_{bb}
      !                +  rho_{ab} ( 1 - rho^*_{ab})
      !                +  kappa_{ab}^* kappa_{ab}
      !
      ! So <N^2> - <N>^2 = Tr(rho ( 1 -rho)) +  Tr(kappa * kappa^{\dagger})
      !
      !-------------------------------------------------------------------------
      ! Note, that at T = 0, we have that (kappa * kappa^{\dagger}) = rho(1-rho).
      ! So in that case, we have
      !  < Delta N^2 > = < N^2 > - <N>^2 = 2 * Tr(rho(1-rho))
      ! which is the old formula from EV8, CR8, etc...
      !
      ! We implement however the formula above, since this is the one that
      ! correctly generalizes to T != 0.
      !-------------------------------------------------------------------------
      real(KIND=dp), intent(in):: rho_can(:), kappa_can(:)
      integer :: wave

      ! BCS dispersion
      BCSdispersion = 0.0d0
      do wave=1,nwn
        BCSdispersion(1) = BCSdispersion(1)                                    &
        &                             + 0.5d0*rho_can(wave)*(1-0.5d0*rho_can(wave))&
        &                             + kappa_can(wave)**2
      enddo

      do wave=nwp+1,nwt
        BCSdispersion(2) = BCSdispersion(2)                                    &
        &                             + 0.5d0*rho_can(wave)*(1-0.5d0*rho_can(wave))&
        &                             + kappa_can(wave)**2
      enddo
      BCSdispersion  = 2 * BCSdispersion

   end subroutine calcBCSdispersion

   module function average_gap_BCS() result(gap)
      !-------------------------------------------------------------------------
      ! Calculation of two types of "average gap", based on
      !
      ! M. Bender et al., EPJA 8, 59-75 (2000).
      !
      ! <v2 Delta > = sum f_k v^2_k   Delta_k / sum f_k v^2_k
      ! <uv Delta > = sum f_k u_k v_k Delta_k / sum f_k u_k v_k
      !
      ! Note that the f_k in the reference is the square of our cutoff!
      !
      ! I calculate the average gap here WITHOUT cutoffs, as our definition
      ! of the gaps already includes all of the cutoff factors already, in
      ! contrast to the EPJA paper.
      !-------------------------------------------------------------------------

      real(KIND=dp) :: gap(2,2), norm(2,2), v2, uv
      integer       :: it, wave

      gap = 0 ; norm = 0
      if(.not.allocated(Pcutoffs)) return
      do wave=1,nwt
        it = 1
        if(wave.gt.nwn) it = 2

        uv = 0.5d0 * BCSgaps(wave)/(BCSqps(wave))
        v2 =       BCSoccupations(wave)

        ! Note that the definition of the gaps include the cutoff factors.
        !  v^2 weighted
        gap(1,it) = gap(1,it)    + v2 * BCSgaps(wave) !*  Pcutoffs(wave)**2
        norm(1,it)= norm(1,it)   + v2                 !*  Pcutoffs(wave)**2
        ! uv weighted
        gap(2,it) = gap(2,it)    + uv * BCSgaps(wave) !*  Pcutoffs(wave)**2
        norm(2,it)= norm(2,it)   + uv                 !*  Pcutoffs(wave)**2
      enddo
      gap = gap/norm

   end function average_gap_BCS

   module subroutine clean_BCS()

     if(allocated(BCSgaps))        deallocate(BCSgaps)
     if(allocated(BCSqps))         deallocate(BCSqps)
     if(allocated(BCSoccupations)) deallocate(BCSoccupations)
     if(allocated(BCSf))           deallocate(BCSf)

   end subroutine clean_BCS

  module function Order(energies, nwf) result(Indices)
    !---------------------------------------------------------------------------
    ! Returns the indices for an ordered traversal of the input array.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !  energies :  real*8, a set of energies to be ordered
    !  nwf      :  size of the array to be ordered
    ! Output:
    !  indices  :  integer, the indices to get the energies in ascending order
    !---------------------------------------------------------------------------
    integer, intent(in)        :: nwf
    real(Kind=dp),intent(in)   :: Energies(nwf)
    real(Kind=dp)              :: Eswap(nwf)
    integer                    :: Indices(nwf)
    integer                    :: i, HolePos, ToInsertIndex
    real(Kind=dp)              :: ToInsert

    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !Filling Energies & Indices
    do i=1,nwf
       Indices(i) = i
    enddo

    Eswap = Energies
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !Sort the energies
    do i=2,nwf
      !Make a hole at index i
      ToInsert = Eswap(i)
      HolePos  = i
      ToInsertIndex = Indices(i)
      do while(ToInsert.lt.Eswap(HolePos-1))
        !Move the hole one place down
        Eswap(HolePos) = Eswap(HolePos-1)
        Indices(HolePos) = Indices(HolePos-1)
        HolePos = HolePos - 1
        if(HolePos.eq.1) exit
      enddo
      !Insert the energy at the correct place
      Eswap(HolePos)    = ToInsert
      Indices(HolePos)  = ToInsertIndex
    enddo

  end function Order

end submodule BCS_quantities
