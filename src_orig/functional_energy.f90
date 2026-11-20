submodule (functional) functional_energy
!==============================================================================
!
! Implementation submodule of module functional: functional_energy
!  Interface block lives in src_orig/functional.f90.
!
!==============================================================================
implicit none
contains

 subroutine CalcEnergy(Rin, Fin, calc_expensive)
    !---------------------------------------------------------------------------
    ! Calculate (i)   the energy
    !           (ii)  the Routhian
    !           (iii) the free energy (when T!= 0)
    ! for a given set of mean-field densities.
    !
    ! Input: 
    !    Rin           : density-vector
    !    Fin           : potential-vector
    !    calc_expensive: controls the calculation of the numerically expensive 
    !                    parts of the total energy. 
    !                    Right now these are:
    !                      (i) the 2-body centre-of-mass correction
    !                     (ii) the Belyaev moments of inertia and <J^2> for
    !                          the rotational and (optional) vibration corrections
    !
    ! This routine is used to calculate "printable" energies, i.e. it is used
    ! for complete calculations of the energy in preparation of a print-out. 
    ! In practice, this means that several module-level quantities are set by 
    ! this routine, all kinds of possible additions to the energy are 
    ! calculated, and the history of the evolution of the energy is saved.
    !
    ! Thus, this routine has a large amount of side-effects, in contrast to the
    ! function calcRouth_onthefly, which has no side-effects and returns ONLY
    ! the value of the Routhian for a given set of mean-field densities. 
    !---------------------------------------------------------------------------
    use momentsofinertia, only: calcJ2andBelyaev, calcrigid
    use fission_MOI, only: calc_collective_inertia
    use Coulombmod
    use densities

    logical, intent(in) :: calc_expensive
    type(DensityVector), intent(in)   :: Rin
    type(PotentialVector), intent(in) :: Fin

    call start_timer(T_energy)

    !---------------------------------------------------------------------------
    ! First we calculate all the individual terms/parts
    ! Kinetic energy
#if ( $TAUPRESENT == 1 )
    if(store_derivatives) then
      Kinetic = CompKinetic_spwfs()
    else
      Kinetic = CompKinetic_density(Rin)
    endif
#else 
    Kinetic = CompKinetic_spwfs() 
#endif

    ! COM correction
    ! (pass signal if we want to skip the calculation of the two-body part)
#if(PASTA == 0)
    call CompCOMCorrection(Kinetic, .not. calc_expensive, COMCorrection)
#else
    ! A waste of CPU time for pasta calculations
    COMCorrection = 0.0d0
#endif

    ! Skyrme functional
    call compSkyrme(Rin, Fin, Skyrme, tot_even, tot_odd, pairdenenergy)

    ! Pairing energy: can be used to check the validity of the calculation.
    ! It is summed by integrating Delta instead of the pairing densities.
    PairingEnergy = CalcPairingEnergy()

    ! If the stabilisation for the pairing is active, calculate the
    ! stabilisationfactor and rescale the pairing energies
    if(abs(Estabp).gt.1d-10 .or. abs(Estabn).gt.1d-10) then
      ! We use the pairing energy obtained  by integrating the pairing
      ! densities
      pairstabfactor = CompStabilisingFactor(PairDenEnergy)

      ! The energy as deduced from the densities
      PairDenE_stab = PairDenEnergy * ( 1 - pairstabfactor)

      ! The energy as deduced from the pairing gaps is no longer right, when
      ! stabilisation is active.
      !-------------------------------------------------------------------------
      ! Assuming a pairing functional that is linear in (kappa kappa*),
      ! the non-stabilised gaps and pair energy are related by
      !   Delta  = d E_pair / d kappa*
      !    E_pair = sum kappa* Delta
      ! The stabilised quantities are
      !   Delta^s  = [ 1 + StabilisingGapFactor ] Delta
      !   E_pair^s = [ 1 - StabilisingGapFactor ] E_pair
      !            = [ 1 - StabilisingGapFactor ] sum kappa* Delta
      ! where StabilisingGapFactor = PairingStabCut^2 / E_pair^2 is a global
      ! state-independent factor. Therefore
      !              [ 1 - StabilisingGapFactor ]
      !   E_pair^s = ---------------------------- sum kappa* Delta^s
      !              [ 1 + StabilisingGapFactor ]
      !-------------------------------------------------------------------------
      PairE_stab    = PairingEnergy  * (1 - pairstabfactor)/(1 + pairstabfactor)

      ! We correct the 'Skyrme' energy here, as the pairing energy was already
      ! summed in there. Hence, we subtract it and add the stabilised one.
      Skyrme = Skyrme - sum(PairDenEnergy) + sum(PairdenE_stab)
    endif

    ! Direct contribution of the Coulomb potential
    CoulombDirect   = coulomb_energy_direct(Rin, Fin)   
    ! Exchange contribution
    CoulombExchange = coulomb_energy_exchange(Rin) 

    call calcrigid(Rin)
#if(PASTA == 0 && $FAM == 0 )
    ! Perform the calculation of all relevant inertias
    if(calc_expensive) then
      call start_timer(T_MOI)
      call calcJ2andBelyaev()
      call stop_timer(T_MOI)
      call calc_collective_inertia()
    endif
    ! ... and only then calculate the contributions to the energy
    Rotcorrection    = calc_rotational_correction(rotcorr)
    Vibcorrection    = calc_vibrational_correction(vibcorr)
#else
    ! These collective corrections are not relevant to
    ! (1) not relevant to infinite systems (PASTA == 1)
    ! (2) expensive to calculate when the number of spwfs is high (FAM == 1)
    Vibcorrection = 0.0d0
    Rotcorrection = 0.0d0
    COMCorrection = 0.0d0
#endif

    ! Entropy calculation when temperature is finite
    call calcentropy()

    !NS: Calculate electrons energy
#if(PASTA==0)
    ElectronEnergyKin=0.0d0
    ElectronEnergyExch=0.0d0
#else
    call calcElectronEnergy()
#endif
    ! The total energy is comprised of 
    !      Kinetic part + Skyrme part + Coulomb energy
    TotalE = Skyrme + sum(Kinetic)
    TotalE = TotalE + CoulombDirect + CoulombExchange
    ! Plus schematic corrections for the collective energy
    TotalE = TotalE + sum(Rotcorrection) + Vibcorrection + sum(COMCorrection)

    ! Total energy from single-particle energies
    SpwfEnergy = calcspwfenergy()

    ! The free energy
    FreeEner = TotalE
    if(inversetemp .gt. 0.0d0) FreeEner = FreeEner - sum(entropy)/inversetemp

    ! Calculate the Routhian
    Routhian = TotalE                                       &
    !                              cranking contribution
    !                               -  omega_mu <J_mu>
    &                         - sum(crankenergy_cut)/2.0_dp &
    !                              multipole contribution
    !                               -  lambda_ml < Q_ml > 
    &                 + sum(Constraint_I_I(:,1:2) * Rin%D_I_I(:,1:2))*dv/2.0_dp

    call stop_timer(T_energy)

 end subroutine CalcEnergy


 subroutine CompSkyrme(R, F, S, T_even, T_odd, PE)
    !---------------------------------------------------------------------------
    ! Integrate the Skyrme energy density for the given set of densities.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !  R : a density-vector, i.e. a set of densities
    ! Output:
    !  S     : total Skyrme energy
    !  T_even: contribution of time-even terms to S
    !  T_odd : contribution of time-odd terms to S
    !  PE    : pairing energy, separated by isospin
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! This routine has a ton of side-effects:
    !  - sets of the values of the individual terms
    !  - sets the bilinear, trilinear, quadrilinear and density-dependent 
    !    partial symmations of the energy.
    !
    ! TODO: further 'functionalize' these routines, eliminate side-effects.
    !---------------------------------------------------------------------------
    real(KIND=dp)                     :: Edensity(mv)
    type(DensityVector), intent(in)   :: R
    type(PotentialVector), intent(in) :: F
    real(KIND=dp), intent(out)        :: S, T_even, T_odd, PE(2)
    
$CALCULATION    

    T_even = &
$TOTAL_EVEN

    T_odd  = &
$TOTAL_ODD

    bilinear = &
$TOTAL_BI

    trilinear = &
$TOTAL_TRI

    quadrilinear = &
$TOTAL_QUAD

    densitydependent = &
$TOTAL_DD

    S = T_even + T_odd
      
    PE(1) = &
$TOTALPAIR_NEUTRON

    PE(2) = &
$TOTALPAIR_PROTON

 end subroutine CompSkyrme


 subroutine PrintSkyrme()
    !---------------------------------------------------------------------------
    ! Print all contributions to the Skyrme energy, automatically generated by
    ! Hephaestos.
    !---------------------------------------------------------------------------
    
    1 format (91('-'))
    3 format (' Skyrme Energy',16x, 'Isospin  1 2 3 4', 33x, 'Energy [MeV]')
    !4 format (17x, 'Total Skyrme:', 3f20.6)
    5 format (17x, 'Total Skyrme:', 41x, f20.6)    

   51 format (17x, '   time-even:', 41x, f20.6)
   52 format (17x, '   time-odd :', 41x, f20.6)

   53 format (17x,     '   bilinear :', 41x, f20.6)
   54 format (17x,     '   trilinear:', 41x, f20.6)
   55 format (17x,     'quadrilinear:', 41x, f20.6)
   56 format (12x, 'density-dependent:', 41x, f20.6)
   
   97 format ( a38, 2a2, 29x, f20.6)
   98 format ( a38, 3a2, 27x, f20.6)
$QUADRI   99 format ( a38, 4a2, 25x, f20.6)

     print 1
     print 3
     print 1
$PRINT
     print 1
     print 5,   Skyrme
     print 1
     print 51,  tot_even
     print 52,  tot_odd
     print 53,  bilinear
     print 54,  trilinear
     print 55, quadrilinear
     print 56, densitydependent
     print 1
 end subroutine PrintSkyrme


 function CompKinetic_spwfs() result(kinetic)
    !---------------------------------------------------------------------------
    ! This subroutine computes the total kinetic energy,
    ! according to the following formula:
    !    E_k = -\hbar/2m \int d^3x \sum_{k} v_{k} \Psi_k^* \Delta \Psi_k
    !---------------------------------------------------------------------------
    ! Note that the 1-body c.o.m. correction is not taken into account here!
    !---------------------------------------------------------------------------
    use Constants

    integer          :: wave, it,k,i, wave_global
    real(KIND=dp)    :: Inproduct
    real(KIND=dp)    :: Kinetic(2)
#if(USE_MPI>0)
    integer          :: mpi_err
#endif
#if($FAM==1)
    integer          :: wave_2, B, N, si
#endif 

    Kinetic = 0.0_dp
#if($FAM == 0)
    ! Kinetic Energy
    do wave=1,nwt_local              ! local spwf index
        wave_global = spwf_map(wave) ! global spwf index

        ! Isospin is neutron in the first half of blocks, proton in the rest
        it = 2
        if(wave_global.le.nwn) it = 1

        Inproduct = 0.0_dp
        do k=1,4
                do i=1,mv
                       Inproduct = Inproduct + DenPsi(i,k,wave) *  &
                       &  ( DenddPsi(i,1,k,wave) + &
                       &    DenddPsi(i,4,k,wave) + &
                       &    DenddPsi(i,6,k,wave))
                enddo
        enddo
        Kinetic(it)= Kinetic(it) + rho_can(wave_global)*Inproduct
    enddo
#else
    ! Nondiagonal sum for FAM
    si = 0
    do B=1,8
      N = HFBLocks(B); if(N.eq.0) cycle 

      ! Isospin is neutron in the first half of blocks, proton in the rest
      it = 1
      if(B.ge.5) it = 2

      do wave=si+1,si+N 
        do wave_2=si+1,si+N
          Inproduct = 0.0_dp
          do k=1,4
            do i=1,mv
              Inproduct = Inproduct + HFPsi(i,k,wave_2) *  &
              &  ( HFddPsi(i,1,k,wave) + &
              &    HFddPsi(i,4,k,wave) + &
              &    HFddPsi(i,6,k,wave))
            enddo
          enddo          
          Kinetic(it)= Kinetic(it) + rho_pairing(wave, wave_2)*Inproduct
        enddo
      enddo
      si = si + N
    enddo
$TR  Kinetic = 2 * Kinetic ! Time-reversal factor 2
#endif

#if(USE_MPI > 0)
    ! Sum the contributions across all MPI ranks
    call MPI_ALLREDUCE(MPI_IN_PLACE, Kinetic, 2, MPI_REAL8, MPI_SUM,           &
    &                                                   MPI_COMM_WORLD, mpi_err)
#endif

    Kinetic=-Kinetic * hbm * dv
    return
  end function CompKinetic_spwfs

#if( $TAUPRESENT == 1 )

  function CompKinetic_density(Rin) result(kinetic)
    !---------------------------------------------------------------------------
    ! This subroutine computes the total kinetic energy from the kinetic density
    !    E_k = -\hbar/2m \int d^3x tau
    !---------------------------------------------------------------------------
    ! Note that the 1-body c.o.m. correction is not taken into account here!
    !---------------------------------------------------------------------------
    real(KIND=dp)                   :: Kinetic(2)
    type(DensityVector), intent(in) :: Rin
    integer                         :: it

    do it=1,2
$TAUSCALAR    Kinetic(it) = hbm(it) *dv * sum(Rin%D_Nm_Nm(:,it))
$TAUTENSOR    Kinetic(it) = hbm(it) *dv * sum(Rin%D_N_N(:,1,1,it)  &
$TAUTENSOR            &                     + Rin%D_N_N(:,2,2,it)  &
$TAUTENSOR            &                     + Rin%D_N_N(:,3,3,it),1)
    enddo
  end function CompKinetic_density
#endif

#if ($FAM == 0)

function calculate_epasta(nucE) result(Epasta)
  !------------------------------------------------------------------------
  ! Calculate the energy per particle for a pasta calculation. This is a
  ! rather simple expression, but it pays to not have it in multiple places.
  !
  ! For now, only nucE is input as all other parts are constant throughout
  ! any given calculation.
  !
  ! Input:
  !  nucE : total energy of the configuration of protons and neutrons as
  !         calculated.
  !
  ! Output:
  !  Epasta: E/A for the entire configuration that accounts for the
  !          contributions of the sea of electrons.
  !------------------------------------------------------------------------

  real(KIND=dp), intent(in) :: nucE
  real(KIND=dp)             :: Epasta

  Epasta = nucE + ElectronEnergyKin+ElectronEnergyExch-protons*Qnp
  Epasta = Epasta/dble(protons+neutrons)

end function calculate_epasta
#endif


  function calcspwfenergy() result(spwfenergy)
    !---------------------------------------------------------------------------
    ! Calculates the total energy from the single-particle energies.
    ! TODO: define inputs!
    !---------------------------------------------------------------------------

    use wavefunctions
    use moments

    integer       :: wave, wave2
    real(KIND=dp) :: spwfenergy, e_rear

    ! Start by summing the single-particle energies
    spwfenergy = 0
    if(pairingtype.lt.2) then
      do wave=1,nwt
        spwfenergy = spwfenergy + rho_can(wave) * spenergies(wave)
      enddo
    elseif(allocated(canenergies)) then
      do wave=1,nwt
        spwfenergy = spwfenergy + rho_can(wave) * canenergies(wave)
      enddo
    else 
      ! Safeguard case
      do wave=1,nwt 
        do wave2=1,nwt 
          spwfenergy = spwfenergy + rho_pairing(wave,wave2) * sphamil(wave2,wave)
        enddo
      enddo
$TR   spwfenergy = 2 * spwfenergy      
    endif
    !
    spwfenergy = 0.5d0 * spwfenergy

    ! Calculation of rearrangement energy (without Coulomb Exchange)
    e_rear = 0
    e_rear = e_rear - 0.5d0*  trilinear
    e_rear = e_rear -         quadrilinear
$EREAR


    spwfenergy = spwfenergy + e_rear
    ! Add kinetic and CoulombExchange contributions
    spwfenergy = spwfenergy + 0.5d0 * sum(kinetic) + CoulombExchange/3.d0

    ! Add the 1-body COMcorrection. 
    if(COM1body.eq.1) then
        ! Perturbative use: add it completely
        SpwfEnergy = SpwfEnergy  + sum(COMCorrection(1,:))
    else if(COM1body.eq.2) then
        ! self-consistent use; add only half
        SpwfEnergy = SpwfEnergy  + sum(COMCorrection(1,:))/2.0_dp
    endif

    if(COM2body.gt.0) then
        SpwfEnergy = SpwfEnergy  + sum(COMCorrection(2,:))
    endif

    ! Subtract contribution by multipole constraints
    Constraint_I_I = constraints_sph_elmult(.false.)
    SpwfEnergy = SpwfEnergy - &
    &          sum(Constraint_I_I(:,1:2) * DBLE(Density%D_I_I(:,1:2)))*dv/2.0_dp

    ! Subtract contribution by cranking constraints
    SpwfEnergy = SpwfEnergy - sum(crankenergy_cut)/2.0_dp

    ! Add the pairing energy (with the stabilisation)
    if(abs(Estabp).gt.1d-10 .or. abs(Estabn).gt.1d-10) then
      SpwfEnergy = SpwfEnergy + sum(PairdenE_stab)
    else
      SpwfEnergy = SpwfEnergy + sum(PairdenEnergy)
    endif
    ! Add the rotational correction
    Spwfenergy = Spwfenergy + sum(Rotcorrection)
    ! And the vibrational correction
    Spwfenergy = Spwfenergy + vibcorrection
  end function calcspwfenergy


  subroutine CompCOMCorrection(kin, override_2body, Comcorr)
    !---------------------------------------------------------------------------
    ! Calculation of the centre-of-mass correction, both one-body and two-body.
    ! General reference for the actual calculation of the entire correction
    !
    ! M. Bender et al., Eur. Phys. J. A 7, 467-478 (2000)
    !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    !  Input:
    ! -------- 
    !  kin           : real     Kinetic energy for the scaling of the 1-body COM
    !
    !  override_2body: logical. If true, skip the 2-body COM calculation, even 
    !                           if the functional requires it.
    !  Output:
    ! --------
    !  ComCorr : 2x2 array, 1/2-body centre-of-mass correction for p/n 
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
    ! REMARK FOR FUTURE GENERALISATIONS
    ! ----------------------------------
    ! For all current symmetry options available, the calculation as implemented
    ! is complete and correct to the best of my (W.R.) knowledge.
    ! HOWEVER, once time-reversal and parity both are broken, the expectation
    ! value of the total momentum of the nucleus is no longer restricted by
    ! symmetry, i.e.
    !
    !             < P > != 0
    !
    ! although individual components might still be restricted by remaining
    ! symmetries.
    !
    ! If that is the case, the correction calculated here should have an extra
    ! contribution that still needs to be implemented namely,
    !
    !        extra term =  - f < P >^2
    !
    ! which will, however, be computationally cheap as it is the square of a
    ! one-body expectation value.
    !
    ! TODO:
    ! -------
    !  (a) add <P>^2 term to the 2-body COM calculation for parity and
    !       time-reversal broken calculations.
    !  (b) allow computation for Hartree-Fock calculations
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !
    ! For future reference (all sums over the entire sp. space, unless
    ! explicitly mentioned)
    !
    !  E_cm = - f < P^2 >     ( f^{-1} = 2 m A)
    !
    !       = - f (sum_{ijkl} P_ij P_kl < a^dagger_i a_j a^dagger_k a_l >
    !
    !  The Wick theorem gives us
    !
    !   < a^dagger_i a_j a^dagger_k a_l >
    !     =    rho_{ji}      rho_{lk}            (a)
    !       -  kappa^*_{ik}  kappa_{lj}          (b)
    !       +  rho_{li} ( delta_{jk} - rho_{jk}) (c1) and (c2)
    !
    ! In the canonical basis, each of these gives rise to
    !
    ! (a) => -f  ( sum_i P_ii rho_ii )^2 = 0  as it is the square of <P>
    !
    ! The one-body component is given by (c1)
    !
    ! (c1) => -f  sum_ij P_ij P_ji rho_{ii}
    !       = -f sum_i P^2_ii rho_ii
    !       = -f (-i hbar)^2 sum_i Delta_ii rho_ii
    !       = +f hbar^2 sum_i Delta_ii rho_ii
    !
    ! The ph-part of the two-body component is given by (c2)
    !
    ! (c2) => +f sum_ij rho_ii rho_jj P_ij P_ji
    !       = +f (-i hbar)^2 sum_ij rho_ii rho_jj Nabla_ij Nabla_ji
    !       = -f hbar^2 sum_ij rho_ii rho_jj Nabla_ij Nabla_ji
    !
    !   (noting that Nabla_ji = - Nabla_ij^*)
    !       = +f hbar^2 sum_ij rho_ii rho_jj |Nabla_ij|^2
    !
    ! The pp-part of the two-body component is given by (b)
    !
    ! (b) => +f sum_ij P_ij P_{ibar jbar}  kappa^*_{i ibar} kappa_{jbar j}
    !      = +f (-i hbar)^2 sum_ij kappa^*_{i ibar} kappa_{jbar j}
    !                                       Nabla_{i j} \Nabla_{ibar jbar}
    !      = -f hbar^2 sum_ij kappa^*{i ibar} kappa^*_{i ibar} kappa_{jbar j}
    !                                       Nabla_{i j} \Nabla_{ibar jbar}
    !
    !     where ibar and jbar are the canonical partners of i and j.
    !
    ! Note
    ! (*) that no symmetries have been used yet at this point
    ! (*) I have not explicitly kept track of the vector nature of P.
    ! (*) The matrix elements of nabla are calculated in compnablamelements
    !     in the densities module.
    !---------------------------------------------------------------------------
    ! There is also a phenomenological way to include the two-body part as a
    ! rescaling of the one-body part, as documented in
    !
    !  M. Butler, D. Sprung and J. Martorell
    !  A improved approximate treatment of c.m. motion in DDHF calculations.
    !  Nucl. Phys. A422 157-166 (1984).
    !
    ! The one-body part is obtained as
    !
    !   hbar^2/2m => hbar^2/(2*m) * (1 - f(A)/A)
    !
    ! with
    !
    !   f(A) = 2/(t + 1/(3t)) with t = (1.5 * A)**(1/3).
    !
    ! It is activated by putting COM1Body = 3, COM2BODY = 0.
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(out)    :: comcorr(2,2)
    logical, intent(in)           :: override_2body
    real(KIND=dp), intent(in)     :: kin(2)

    integer       :: it, i,j
$NTR integer       :: B, ibar, jbar, ii, jj, N, N2, N3, N4, si
    real(KIND=dp) :: NablaMElements(3,2,nwt,nwt), fac
    real(KIND=dp) :: Butler_t, Butler_f, prefac(2)

    call start_timer(T_com)
    call start_timer(T_com1)

    !Only reset the one-body part
    COMCorr(1,:) = 0.0_dp

    select case(COM1Body)
    case(0)
      ! No contribution
    case(1,2)
      ! Deduce 1-body COM correction from the Kinetic Energy
      COMCorr(1,:) = - Kin(:) * nucleonmass/ &
      &                 (neutrons * nucleonmass(1) + protons * nucleonmass(2))
    case(3)
      ! Deduce 1-body COM correction from the Kinetic Energy with Butlers
      ! formula.

      Butler_t = (1.5d0 * (neutrons + protons))**(1.d0/3.d0)
      Butler_f = 2.d0/(Butler_t + 1.d0/(3*Butler_t))

      COMCorr(1,:) = - Kin(:) * nucleonmass * Butler_f/ &
      &                 (neutrons * nucleonmass(1) + protons * nucleonmass(2))
    end select

    call stop_timer(T_com1)

    ! The calculations is not yet implemented for Hartree-Fock calculations
$NTR    if(COM2body .ne. 0 .and. pairingtype .eq. 0) then
$NTR      call stp('Two-body COM not implemented yet for Hartree-Fock calculations with time-reversal breaking.')
$NTR    endif
    
    if((COM2body .eq. 1) .and. (.not. override_2body)) then
      !-------------------------------------------------------------------------
      ! The 2-body COM correction, calculated as discussed above
      !-------------------------------------------------------------------------
      call start_timer(T_com2)

      NablaMElements = compNablaMelements()
      call start_timer(T_com2_summation)

      COMCorrection(2,:) = 0.0
      COM2pp = 0.0 ; COM2ph = 0.0 ; COM2_pp_debug = 0.0 ; COM2_ph_debug=0.0
      do i=1,nwt
         ! We sum over all possible (i,j) pairs, the matrix elements are
         ! correctly calculated either way.
         it = 1
         if(i.gt.nwn) it = 2
         do j=1,nwt
            ! + sum_ij rho_ii rho_jj |Nabla_ij|^2
            ! v^2 v^2 part
            fac = rho_can(i)*rho_can(j)
$TR         fac = fac / 4.0d0 ! rho_can is twice too large if T is conserved
            COM2_ph_debug(1,it) = COM2_ph_debug(1,it) + fac*NablaMElements(1,1,i,j)**2
            COM2_ph_debug(2,it) = COM2_ph_debug(2,it) + fac*NablaMElements(2,2,i,j)**2
            COM2_ph_debug(3,it) = COM2_ph_debug(3,it) + fac*NablaMElements(3,1,i,j)**2
$TR         ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
$TR         ! - sum_ij kappa^*{i ibar} kappa^*_{i ibar} kappa_{jbar j}
$TR         !                                    Nabla_{i j} \Nabla_{ibar jbar}
$TR         ! uv uv part
$TR         ! (in the case of conserved time-reversal)
$TR         fac = kappa_can(i)*kappa_can(j)
$TR         COM2_pp_debug(1,it) = COM2_pp_debug(1,it) + fac*NablaMElements(1,1,i,j)**2
$TR         COM2_pp_debug(2,it) = COM2_pp_debug(2,it) + fac*NablaMElements(2,2,i,j)**2
$TR         COM2_pp_debug(3,it) = COM2_pp_debug(3,it) + fac*NablaMElements(3,1,i,j)**2
         enddo
      enddo

$NTR      si = 0
$NTR      do B=1,8,4 ! This is essentially an isospin loop now
$NTR        N  = HFblocks_global(B) ; if(N.eq.0) cycle
$NTR        N2 = HFblocks_global(B+1)
$NTR        N3 = HFblocks_global(B+2)
$NTR        N4 = HFblocks_global(B+3)
$NTR        it = 1  ; if(B .eq.5) it = 2
$NTR        ! We simply loop over all possible combinations of spwfs
$NTR        ! as the NablaMElements array is zero in the right places
$NTR        do i=1, N+N2+N3+N4
$NTR          ii   = si + i
$NTR          ibar = conjugp(ii) ; if(ibar .eq.0) cycle
$NTR          do j=1,  N+N2+N3+N4
$NTR            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
$NTR            ! - sum_ij kappa^*{i ibar} kappa^*_{i ibar} kappa_{jbar j}
$NTR            !                                Nabla_{i j} \Nabla_{ibar jbar}
$NTR            ! (in the case of broken time-reversal)
$NTR            jj   = si +  j
$NTR            jbar = conjugp(jj) ; if(jbar .eq.0) cycle
$NTR            fac = -  kappa_can(ii)*kappa_can(jbar)
$NTR            COM2_pp_debug(3,it) = COM2_pp_debug(3,it) &
$NTR            &                   + fac * NablaMElements(3,1,ii,jj)          &
$NTR            &                         * NablaMElements(3,1,ibar,jbar)
$NTR            COM2_pp_debug(1,it) = COM2_pp_debug(1,it) &
$NTR            &                   + fac * NablaMElements(1,1,ii,jj)          &
$NTR            &                         * NablaMElements(1,1,ibar,jbar)
$NTR            COM2_pp_debug(2,it) = COM2_pp_debug(2,it) &
$NTR            &                   - fac * NablaMElements(2,2,ii,jj)          &
$NTR            &                         * NablaMElements(2,2,ibar,jbar)
$NTR            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
$NTR          enddo
$NTR        enddo
$NTR        si = si + N + N2 + N3 + N4
$NTR      enddo
      call stop_timer(T_com2_summation)

      ! The prefactor f
      prefac=hbm*nucleonmass/(neutrons*nucleonmass(1) + protons*nucleonmass(2))
      do it=1,2
        COM2_ph_debug(:,it) = prefac(it) * COM2_ph_debug(:,it)
        COM2_pp_debug(:,it) = prefac(it) * COM2_pp_debug(:,it)
      enddo
      ! In the case of Time-reversal conservation, we summed over only half
      ! the states
$TR   COM2_ph_debug = 2*COM2_ph_debug
$TR   COM2_pp_debug = 2*COM2_pp_debug

      ! Summing the three directions
      do it=1,2
        COM2ph(it) = sum(COM2_ph_debug(:,it))
        COM2pp(it) = sum(COM2_pp_debug(:,it))
      enddo

      ! The final relevant energy is the sum of ph and pp contributions
      do it=1,2
        COMCorrection(2,it) = COM2ph(it) + COM2pp(it)
      enddo
      call stop_timer(T_com2)
     endif
     call stop_timer(T_com)

  end subroutine CompCOMCorrection


  function calc_rotational_correction(rot_option) result(Erot)
    !---------------------------------------------------------------------
    !
    ! Input:
    !   rot_option
    ! Output:
    !   Erot:
    !---------------------------------------------------------------------
    use momentsofinertia, only: select_J2_and_MOI

    integer, intent(in)        :: rot_option
    real(KIND=dp)              :: Erot(3)

    real(KIND=dp)              :: B(3),f_rot(3)
    real(KIND=dp)              :: J2(3), MOI(3), MOI_c

    Erot = 0.0d0
    select case(rot_option)
    case(0)
      !----------------------------------------------------------------------
      ! No rotational correction
      return
    case(1)
      !----------------------------------------------------------------------
      ! Calculate a phenomenological collective correction based on the
      ! perturbative Belyaev moment of inertia.
      !
      !    Erot = - \sum_{\mu} f^rot_mu <J_mu^2>/(2 * I_{\mu})
      !
      ! where
      !      * the sum is over all three Cartesian directions
      !      *  < J^2_{mu} > is the expectation value of the angular momentum
      !                      squared in a given direction
      !      *  I_mu is the Belyaev moment of inertia along a given direction
      !      * f^rot_mu is a cutoff function for the rotational correction
      !
      ! We take
      !
      !      f^rot_mu = b tanh( c B_mu )
      !
      ! where B_mu = I_mu / I_c is the ratio between the calculated Belyaev
      ! moment of inertia and (one third of) the classical moment of inertia,
      ! and the parameters b and c are encoded as rotcorrb and rotcorrc.
      call select_J2_and_MOI(J2,MOI)
      MOI_c = classical_MOI(neutrons,protons)
      B = MOI/MOI_c * 3.0d0
      f_rot         = rotcorrb * tanh(rotcorrc * B)
      Erot          = - f_rot * J2/(2*MOI)
    case DEFAULT
      call stp('Unknown value for rot_opion in calc_rotational_correction')
    end select

  end function calc_rotational_correction


  function calc_vibrational_correction(vib_option) result(Evib)
    !---------------------------------------------------------------------------
    ! Calculate the energy correction for spurious vibrational motion.
    !
    ! Input:
    !  vib_option: type of vibrational correction to calculate
    !          (0) : no correction
    !          (1) : simple formula based on the rotational one
    !          (2) : correction inspired by the ZPE of a collective Hamiltonian
    ! Output:
    !   Evib: the total vibrational correction
    !---------------------------------------------------------------------------
    use momentsofinertia, only :  select_J2_and_MOI
    use fission_MOI, only      :  M2, M3

    integer, intent(in)        :: vib_option
    real(KIND=dp)              :: Evib

    real(KIND=dp)              :: MOI_c, fvib(3), MOI(3), J2(3), B(3)
    real(KIND=dp)              :: M3_temp(2,2), M3_inv(2,2), temp(2,2)

    select case(vib_option)
    case (0)
      ! No vibrational correction include
      Evib = 0
      return
    case (1)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! Simple vibrational correction based on the rotational moment of inertia
      !
      ! Evib = - \sum_{\mu}  f^vib_mu  <J_mu^2>/(2 * I_{\mu}
      !
      ! where
      !      * the sum is over all three Cartesian directions
      !      *  < J^2_{mu} > is the expectation value of the angular momentum
      !                      squared in a given direction
      !      *  I_mu is the Belyaev moment of inertia along a given direction
      !      *  f^vib_mu = d B_mu exp ( -l  (B - b_vib)**2  )
      !
      ! The quantity  B_mu = I_mu / I_c is the ratio between the calculated
      ! Belyaev moment of inertia and one third of the classical moment of
      ! inertia, 1.0/3.0 * I_c as calculated in the function classical_MOI.
      !
      ! The free parameters are d, l and b_vib; these are specified by a
      !  parameterization and are called vibcorrd, vibcorrl and vibcorrb in
      !  parameterization.f90.
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      MOI_c = classical_MOI(neutrons, protons)
      call select_J2_and_MOI(J2, MOI)
      B     = MOI/MOI_c * 3.0d0

      fvib  = vibcorrd * B * exp( - vibcorrl * (B - vibcorrb)**2)
      Evib  = - sum(fvib * J2/(2*MOI))
    case (2)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! Vibrational correction from the zero-point energy of a collective Hamiltonian
      !
      !  Evib = - 1/4 * d * Trace M_3^{-1} M_2
      !
      ! where
      !  - M_{n} are the QRPA-like response matrices for the Q20 and Q22
      !     collective degrees of freedom; see the fission_moi.f90 file.
      !  - d is a free parameter that is encoded as vibcorrd in the .param files
      !
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! Original reference: M. Girod and B. Grammaticos, NPA 330 (1979).
      ! More practical:     T. Niksic et al. PRC 79, 034303 (2009)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! ATTENTION: hardcoded assumption that the (1:2,1:2) submatrices of
      !            M_3 and M2 as calculated in fission_MOI represent Q20 and Q22!
      M3_temp = M3(1:2,1:2,3)
      M3_inv  = inv_mat(M3_temp)

      temp = matmul(M3_inv, M2(1:2,1:2,3))
      Evib = - 0.25d0 * vibcorrd *  (temp(1,1) + temp(2,2))
    end select
  end function calc_vibrational_correction


  pure function inv_mat(A) result(Ainv)
    !
    ! Analytical inverse of a 2x2 matrix
    !

    real(KIND=dp), intent(in) :: A(2,2)
    real(KIND=dp)             :: Ainv(2,2)
    real(KIND=dp)             :: det

    ! Calculate the inverse (!) of the determinant
    det = 1/(A(1,1)*A(2,2) - A(1,2)*A(2,1))

    ! Calculate the inverse of the matrix
    Ainv(1,1) = +det * A(2,2)
    Ainv(2,1) = -det * A(2,1)
    Ainv(1,2) = -det * A(1,2)
    Ainv(2,2) = +det * A(1,1)

  end function inv_mat


  function classical_MOI(N,Z) result(MOI)
      !------------------------------------------------------------------------
      ! Calculate the classical rotational moment of inertia.
      !
      !  I_c = 2/5 * m_n * A * R**2/(hbar c**2)
      !   with R = (1.2 * A)**3
      !------------------------------------------------------------------------
      real(KIND=dp), intent(in) :: N,Z
      real(KIND=dp)             :: MOI
      real(KIND=dp)             :: A, R

      A   = N + Z
      R   = 1.2d0 * (A)**(1.d0/3.d0)
      MOI = 2.0d0/5.0d0 * sum(nucleonmass)/2 * (A)*R**2
      MOI = MOI/(hbarclum**2)
  end function classical_MOI


  function CompStabilisingFactor(PairE) result(stab)
    !---------------------------------------------------------------------------
    ! Calculate StabilisingGapFactor =  E_cut^2/E_pair^2 when needed, i.e.
    ! when PairingStabCut != 0. Otherwise set StabilisingGapFactor to zero.
    ! Taken (with minor modifications from MOCCav1, routine by M. Bender)
    !---------------------------------------------------------------------------
    integer                   :: it
    real(KIND=dp)             :: stab(2), cut(2)
    ! PE is the user's choice of pairing energy
    real(KIND=dp), intent(in) :: PairE(2)

    cut(1) = Estabn
    cut(2) = Estabp
    stab   = 0.0
    if (pairingtype.eq.0) return

    do it=1,2
      !-------------------------------------------------------------------------
      ! No stabilisation for this isospin. Set factor to zero.
      !-------------------------------------------------------------------------
      if ( abs(cut(it)) .lt. 1.d-10 ) cycle

      !-------------------------------------------------------------------------
      ! If pairing energy is non-zero, so just calculate the factor.
      ! If pairing energy is zero (meaning this is either the initial call or
      ! a failure), fall back on predefined value (0.1 MeV).
      !-------------------------------------------------------------------------
      if ( abs(PairE(it)) .gt. 1.d-8 ) then
        stab(it) = cut(it)**2/(PairE(it)**2)
      else
        if ( stab(it) .eq. 0.0_dp ) then
          stab(it) = 0.1_dp
           print '(" StabilisingFactor initialised to ",f12.6,  &
            &   " for it = ",i1, es15.5)', stab(it),it, PairE(it)
        endif
      endif
      if ( stab(it) .gt. 10.0 ) then
        print '(" WARNING: StabilisingFactor: for it = ",i1, &
        & " StabilisingGapFactor = ",1d16.8," for an energy of ",1d16.8)', &
        & it,stab(it),PairE(it)
      endif
    enddo

  end function CompStabilisingFactor


  subroutine verify_COM_motion()
    !---------------------------------------------------------------------------
    !  This routine performs calculations of the collective moments of inertia
    !  for the motion of the z-coordinate of the center of mass.
    !
    !  The collective coordinate for species q is thus
    !    Q_q = z_q / A
    !  to which corresponds a collective momentum  (in our convention)     
    !    P_q = - i \nabla_z (*)
    !  
    !  One can show that, analytically, the collective inertia associated with 
    !  movement of the centre-of-mass should be the TOTAL mass of the nucleus
    ! 
    !   M'_{0} = A m
    !
    !  This is what is always presented in the literature. Note however the 
    !  little accent M', indicating that this IS NOT the collective inertia
    !  associated with P_q. Rather, it is the collective mass associated with 
    !     P'_q = hbar P_q
    !  i.e. the 'physical' convention for the momentum. 
    !
    ! We calculate the collective inertia related to COM motion here in multiple
    ! ways:
    !
    !   (a) Analytically, printed as 'A m'
    !   (b) By using the Belyaev formula for the momentum in our convention (*)
    !       (and multiplying by hbar^2 afterward)
    !   (c) By using the perturbative cranking formula starting from Q_q
    !       (and converting convention again afterward)
    !
    ! The results of (b) and (c) are not necessarily close to (a) however: 
    ! for typical Skyrme interactions the effective mass m^*/m is not equal
    ! to one, spoiling the correspondence of the perturbative formulation. 
    ! The origin lies in the absence of Galileian invariance of the interaction 
    ! in the perturbative treatment, see 
    !
    !    K. Wen, and T. Nakatsukasa,  http://arxiv.org/abs/2112.13317
    ! 
    ! for a discussion. Ideally, we would calculate an "average" effective mass
    ! as these authors do and use it to correct our results.
    !
    !---------------------------------------------------------------------------
    
    use fission_MOI
    use densities
    
    1 format (' Pushing model                 M_0 (MeV/c^2)')
    2 format (25x, 'neutrons        protons          total', / &
    &         1x, 80('-'))
    3 format ('   Belyaev M_0      | ', 3f16.5)
    4 format ('   Pert. cranking   | ', 32x,f16.5)
    7 format (1x, 80('-'),/, &
    &         '   Am               | ', 3f16.5, /, 80('-'))
    
    real(KIND=dp), allocatable :: NablaMElements(:,:,:,:)
    real(KIND=dp) :: mat(2,2), Pmat(2,2), totalmass
      
    real(KIND=dp) :: Psp(nwt,nwt), Qsp(nwt,nwt), hbar
    real(KIND=dp) :: P20(nwt,nwt), Q20(nwt,nwt)
    
    ! Calculate hbar to make its use consistent
    hbar =  sqrt(hbm(1) * 2  * 0.5 * sum(nucleonmass))
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Calculate single-particle matrix elements of nabla_z 
    NablaMElements = compNablaMElements()
    Psp            =  NablaMElements(3,1,:,:) 
    if(pairingtype.eq.2) then
      ! Attention, in a HFB calculation these are the matrix elements in the 
      ! canonical basis, while Bogoliubov refers to the HF basis. So, we 
      ! transfer back to the HFbasis
      Psp = matmul(matmul(cantransfo, Psp), transpose(cantransfo))
    endif
    
    Psp(1:nwn,1:nwn)         =  Psp(1:nwn,1:nwn)         * hbar
    Psp(nwn+1:nwt,nwn+1:nwt) =  Psp(nwn+1:nwt,nwn+1:nwt) * hbar
        
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Calculate single-particle matrix elements of z-c.o.m. coordinate
    Qsp            = sqrt(4*pi/3) * Qlm_spme(1,0,.false.) * 10 
                    ! Q10 = sqrt(3/4pi) * z 
                    ! and the routine Qlm_spme uses units of b^1/2 => factor 10
    Qsp(1:nwn,1:nwn)         =      Qsp(1:nwn,1:nwn)        /(neutrons+protons)
    Qsp(nwn+1:nwt,nwn+1:nwt) =      Qsp(nwn+1:nwt,nwn+1:nwt)/(protons +neutrons)
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! TODO: Calculate the average effective mass for the pushing model.
    !       This used to work, but now fails thanks to the changes of the 
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Calculate average effective mass 
    !
    !   m^*/m_q  = 2m/N_q * int d^3r   rho_q(r) hbar^2/2m^*_q(r)
    ! 
    ! with hbar^2/2m^*_q =  hbar^2/2m_q(r) + F_Nm_Nm(r)
    ! 
    ! which is the logical generalization from Eq. (30) in 
    !    K. Wen, and T. Nakatsukasa,  http://arxiv.org/abs/2112.13317
    ! and an explicit factor of hbar^2.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !     avg_effmass = sum(F%F_Nm_Nm(:,:) * R%D_I_I(:,:) , 1) * dv
    !     avg_effmass(1) = avg_effmass(1) / (hbm(1) * neutrons)
    !     avg_effmass(2) = avg_effmass(2) / (hbm(2) * protons)
    !     avg_effmass    = avg_effmass + 1
    !     avg_effmass    = 1/avg_effmass
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    print 1
    print 2
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! perform the summations
    select case(pairingtype)
    case(0) ! HF
      Pmat                     =  Ksum_Mij_HF(Psp, Psp, 1, 1,(/1,3/))
      mat                      =  Ksum_Mij_HF(Qsp, Qsp, 1, 1,(/1,3/))
    case(1) ! BCS
      Pmat                     =  Ksum_Mij_BCS(Psp, Psp, 1, 1,(/1,3/))
      mat                      =  Ksum_Mij_BCS(Qsp, Qsp, 1, 1,(/1,3/))
    case(2) !HFB
      P20 = calc_Q20(Psp, bogoliubov, 1)        
      Pmat                     =  Ksum_Mij(P20, P20, 1, 1,(/1,3/))
      Q20 = calc_Q20(Qsp, bogoliubov, 1)        
      mat                      =  Ksum_Mij(Q20, Q20, 1, 1,(/1,3/))
    end select
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Calculating the inertia parameters
    totalmass   = 1.0d0/(sum(mat(1,:))) * sum(mat(2,:)) * 1.0d0/(sum(mat(1,:)))
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Printing 
    print 3,  Pmat(1,:), sum(Pmat(1,:))
    print 4,  totalmass * hbar**2
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Analytical result
    print 7,  neutrons * nucleonmass(1), protons*nucleonmass(2), &
    &         neutrons * nucleonmass(1)+ protons*nucleonmass(2) 

    print *

  end subroutine verify_COM_motion

  subroutine calcElectronEnergy()
      !NS: calculate kinetic energy of relativistic electron gas including exchange
      !(but latter in ultrarelativistic limit)
      
      real(KIND=dp) :: lamce, pfermi, xx, xx2, hi_x, E_rel, E_ultrarel, ne
      real(KIND=dp),parameter :: cc=2.99792458d23     !codata speed of light fm/s
      real(KIND=dp),parameter :: me=0.510998950d0 !codata electron mass in MeV
      real(KIND=dp),parameter :: hh=4.135667696d-21/(2.d0*pi) !codata h dirac MeV*s
      real(KIND=dp),parameter :: alphaem=7.2973525693d-3 !Codata fine structure

      if(protons .eq. 0.0d0) then
         ElectronEnergyKin   = 0.0d0
         ElectronChempotKin  = 0.0d0
         ElectronEnergyExch  = 0.0d0
         ElectronChempotExch = 0.0d0
      else
         ne=protons/(mv*dv)
         !Relativistic electrons
         lamce=hh*cc/me
         pfermi=(3.d0*(hh*2.d0*pi)**3.d0/(8.d0*pi)*ne)**(1.d0/3.d0)
         xx=pfermi*cc/me
         xx2=xx*xx
         hi_x=1.d0/(8.d0*pi*pi)*(xx*sqrt(1.d0+xx2)*(1.d0+2.d0*xx2)-log(xx+sqrt(1.d0+xx2)))
         E_rel= me/(lamce**3.d0)*hi_x
         
         !Ultrarelativistic electrons
         E_ultrarel=0.75d0*(3.d0*pi*pi)**(1.d0/3.d0)*hh*cc*ne**(4.d0/3.d0)
         
         ElectronEnergyKin=E_rel*mv*dv
         ElectronChempotKin=me*sqrt(1+xx2)

         !Electron exchange energy and chempot
         ElectronEnergyExch=E_ultrarel*alphaem/2.d0/pi*mv*dv
         ElectronChempotExch=4.d0/3.d0*E_ultrarel*alphaem/2.d0/pi/ne
      endif

  end subroutine calcElectronEnergy


 subroutine update_E_history()
    !---------------------------------------------------------------------------
    ! Update the history of the module with the values of various things
    ! currently in storage.
    !---------------------------------------------------------------------------
    integer :: i

    ! Move old values
    do i=4,1,-1
        Ehistory(i+1)    = Ehistory(i)
        Rhistory(i+1)    = Rhistory(i)
        Fhistory(i+1)    = Fhistory(i)
        Spwfhistory(i+1) = Spwfhistory(i)
    enddo
    Ehistory(1)    = TotalE
    Rhistory(1)    = Routhian
    Fhistory(1)    = FreeEner
    SpwfHistory(1) = SpwfEnergy
 end subroutine update_E_history
 

end submodule functional_energy
