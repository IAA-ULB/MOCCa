module functional
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
 ! Module containing the means to calculate (and print) the mean-field energy.
 ! Note that the actual coupling constants are contained in the constants.f90
 ! file. 
 !
 ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 ! Hephaestos keywords
 ! 
 ! Declaration     : [WAY TOO LONG TO INCLUDE HERE]
 ! Calccoef        : [WAY TOO LONG TO INCLUDE HERE]
 ! PrintCOEF_iso   : [WAY TOO LONG TO INCLUDE HERE]
 ! PrintCOEF_PN    : [WAY TOO LONG TO INCLUDE HERE]
 ! PrintCOEF_pair  : [WAY TOO LONG TO INCLUDE HERE]
 ! Calculation     : [WAY TOO LONG TO INCLUDE HERE]
 ! Total_EVEN      : [WAY TOO LONG TO INCLUDE HERE]
 ! Total_ODD       : [WAY TOO LONG TO INCLUDE HERE]
 ! TotalPAIR       : [WAY TOO LONG TO INCLUDE HERE]
 ! Print           : [WAY TOO LONG TO INCLUDE HERE]
 ! Calcfields      : [WAY TOO LONG TO INCLUDE HERE]
 ! SkyrmeAction    : [WAY TOO LONG TO INCLUDE HERE]
 ! PairingAction   : [WAY TOO LONG TO INCLUDE HERE]
 ! ERear           : [WAY TOO LONG TO INCLUDE HERE]
 ! CLEANING        : [WAY TOO LONG TO INCLUDE HERE]
 ! FIELDNUMBER     : [WAY TOO LONG TO INCLUDE HERE]
 ! WRITEPOTENTIALS : [WAY TOO LONG TO INCLUDE HERE]
 ! READPOTENTIALS  : [WAY TOO LONG TO INCLUDE HERE]
 ! NTR             : $NTR
 ! N2              : $N2
 ! N3              : $N3
 !==============================================================================
 
 use compilation
 use geninfo
 use densities
 use parameterization
 use pairing
 use timing
 use transform
 use Cranking

 implicit none
 
    !===========================================================================
    ! PARAMETERIZATION DEFINITION OPTIONS
    !===========================================================================
    !---------------------------------------------------------------------------
    ! Name of the parameterization
    character(len=20) :: name_param 
    ! Name of the functional file this code was compiled with
    character(len=20), parameter :: func_name = $FUNC_NAME
    !---------------------------------------------------------------------------
    ! Definition of global contributions to the energy
    real(KIND=dp) :: Kinetic(2), Skyrme(2), TotalE, SpwfEnergy, Ehistory(5)
    real(KIND=dp) :: tot_even(2), tot_odd(2)
    real(KIND=dp) :: COMCorrection(2,2), CoulombDirect, CoulombExchange
    ! Two definitions of the pairingenergy: one obtained by summing the gaps
    ! and one by integrating the particle-particle part of the functional
    real(KIND=dp) :: PairingEnergy(2), PairDenEnergy(2)
    ! Same thing, but with the added stabilization
    real(KIND=dp) :: PairE_stab(2), PairDenE_stab(2)
    !---------------------------------------------------------------------------
    ! Rotational correction
    real(KIND=dp) :: RotCorrection(3)
    !===========================================================================
    ! NUMERICAL OPTIONS
    !===========================================================================
    !---------------------------------------------------------------------------
    ! Small non-zero value that can be used in a .func file to safeguard against
    ! division by zero. 
    real(KIND=dp) :: eps =1d-20
    !---------------------------------------------------------------------------
    ! Factor in the preconditioning of the F_I_I potential
    real(KIND=dp) :: preconfactor = 1.0_dp
    !---------------------------------------------------------------------------
    ! Type of preconditioning to apply
    integer       :: potentialpreconditioning = 1
    !---------------------------------------------------------------------------
    ! Stabilisation factor for the pairing:
    !    f = E_cut^2 / E_pair^2
    real(KIND=dp) :: pairstabfactor(2) = 0.0
    !===========================================================================
    ! AUTOMATICALLY GENERATED DECLARATIONS
    !===========================================================================
    !---------------------------------------------------------------------------
    ! Declaration of
    ! a) Coupling constants
    ! b) Terms in the energy
    ! c) Fields
    ! d) Previous_values of the fields
    ! All automatically generated by Hephaestos
    !---------------------------------------------------------------------------
$DECLARATION
    !---------------------------------------------------------------------------
      
contains

 subroutine readfunctional(file_number)
    !---------------------------------------------------------------------------
    ! Initializes the functional
    ! a) read the details of the parameterization from file
    ! b) calculate the coupling constants
    !---------------------------------------------------------------------------
    
    integer(dp), intent(in), optional   :: file_number 

    namelist /func/ name_param
    
    if(present(file_number)) then
      read(unit=file_number, nml=func) 
    else
      read(unit=*, nml=func) 
    endif

    call readparameterization(name_param, func_name)
    call calcedfcoefs()

    !---------------------------------------------------------------------------
    ! Put the pairing routines pointers to the action of Delta
    delta_action_BCS => delta_action
    delta_action_HFB => delta_action
    
 end subroutine readfunctional
 
 subroutine calcedfcoefs()
     !--------------------------------------------------------------------------
     ! Calculate the coupling constants from the expressions passed into 
     ! Hephaestos.
     !--------------------------------------------------------------------------
$CALCCOEF
    
    nucleonmass = 0.5*(nucleonmass(1) + nucleonmass(2))
 
 end subroutine calcedfcoefs
 
 subroutine printfunctional
    !---------------------------------------------------------------------------
    ! Print all information on the functional.
    !---------------------------------------------------------------------------
    call printparameterization(name_param, func_name)
    call printedfcoefs
    
 end subroutine printfunctional
 
 subroutine printedfcoefs
    !---------------------------------------------------------------------------
    ! Print the values of the EFD coefs used.
    !---------------------------------------------------------------------------
    1 format (' - - - - - - - - - - -')
    2 format (5x, 58('_'))
    3 format (' Skyrme coupling constants ')
    4 format (38x, 'Isospin representation')
    5 format (36x, '  C_0            C_1    ')
    6 format (38x, 'BFH representation ')
    7 format (36x, '  C_0            C_q    ')
    8 format (38x, 'Pairing terms ')
    9 format (36x, '  neutron        proton ')
    
     print 1
     print 3
     print 2
     print 4
     print 5
     print 2
     ! This part is automatically generated by Hephaestos.
$PRINTCOEF_ISO
     print 2
     print 6
     print 7
$PRINTCOEF_PN
     print 2
     print 8
     print 9
$PRINTCOEF_PAIR
     print 2
 end subroutine printedfcoefs

 subroutine PrintEnergy()
    !---------------------------------------------------------------------------
    ! Print all of the information on the energy.
    !---------------------------------------------------------------------------
    use Coulombmod
    
    1 format (80('-'))
    5 format (30x, '       neutron        proton         total')
    6 format (15x, ' Kinetic Energy:', 3f15.6)
   61 format (15x, '     COM 1-body:', 3f15.6)
   62 format (15x, '     COM 2-body:', 3f15.6)
   63 format (15x, '   Rotational ', a1, ':', 30x, f15.6)
    7 format (15x, ' Coulomb Direct:', 3f15.6)
   71 format (15x, '   Dir. (point):', 3f15.6)
    8 format (15x, '       Exchange:', 3f15.6)
   81 format (15x, '   Exc. (point):', 3f15.6)  

    9 format (15x, 'Pairing (delta):', 3f15.6)
   91 format (15x, 'Pairing (densi):', 30x, f15.6)
   92 format ( 7x, 'Pair. (delta, no stab.):', 3f15.6)
   93 format ( 7x, 'Pair. (densi, no stab.):', 30x, f15.6)
   94 format ( 7x, 'Pair. (delta,    stab.):', 3f15.6)
   95 format ( 7x, 'Pair. (densi,    stab.):', 30x, f15.6)

   99 format (15x, '   Total energy:', 30x, f15.6)
  991 format (15x, '    (no corr.) :', 30x, f15.6)
  100 format (15x, '     from spwfs:', 30x, f15.6)
  101 format (15x, '    Free Energy:', 30x, f15.6)
  102 format (15x, '        Entropy:', 3f15.6)
  103 format (15x, '    E_fu - E_sp:', 30x, e15.6)
    real(KIND=dp) :: temp

    call printSkyrme

    print 1
    print 5
    print 6, Kinetic, sum(Kinetic)
    print 61, COMcorrection(1,:), sum(COMcorrection(1,:))
    if(any(COMcorrection(2,:).ne.0)) then
     print 62, COMcorrection(2,:), sum(COMcorrection(2,:))
    endif

    if(rotcorr .ne.  0) then
      print 63, 'X',  Rotcorrection(1)
      print 63, 'Y',  Rotcorrection(2)
      print 63, 'Z',  Rotcorrection(3)
    endif
  
    print *
    print 7, 0.0, CoulombDirect, CoulombDirect
    if(protonsize(1).ne.0 .and. (.not. nucleonsize_selfconsistent)) then
      temp = CoulombEnergy_Direct(D_I_I(:,2))
      print 71, 0.0, temp, temp
    endif
    print 8, 0.0, CoulombExchange, CoulombExchange
    if(protonsize(1).ne.0 .and. (.not. nucleonsize_selfconsistent)) then
      temp = CoulombEnergy_Exchange(D_I_I(:,2))
      print 81, 0.0, temp, temp
    endif
    print *
    if( abs(Estabp).lt.1d-10 .and. abs(Estabn).lt.1d-10) then
      print 9 , PairingEnergy, sum(PairingEnergy)
      print 91, sum(PairDenEnergy)
    else
      print 92, PairingEnergy, sum(PairingEnergy)
      print 93, sum(PairDenEnergy)
      print 94, PairE_stab, sum(PairE_stab)
      print 95, sum(PairdenE_stab)
    endif    
    print 1
    print  99, TotalE
    if(rotcorr.ne.0) then
        print 991, totalE - sum(rotcorrection)
    endif
    print 100, spwfenergy
    print 103, TotalE - spwfenergy
      
    if(inversetemp .ne. -1) then
        ! F = E - T * S
        print 101, TotalE - sum(entropy)/inversetemp
        print 102, entropy, sum(entropy)
    endif

    print 1
 end subroutine PrintEnergy
 
 subroutine CalcEnergy(iprint)
    !---------------------------------------------------------------------------
    ! Calculate all of the relevant energies.
    !---------------------------------------------------------------------------
    use momentsofinertia
    use Coulombmod

    integer :: i
    integer, intent(in) :: iprint

    call start_timer(T_energy)
    
    ! Kinetic energy
    Kinetic = CompKinetic()
    ! COM correction
    call CompCOMCorrection()
    ! Skyrme functional
    call compSkyrme()

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
      Skyrme = Skyrme - PairDenEnergy + PairdenE_stab
    endif


    if( all(protonsize.eq.0.0) .and. all(neutronsize.eq.0.0) ) then
      ! Direct contribution of the Coulomb potential
      CoulombDirect   = CoulombEnergy_Direct(D_I_I(:,2))
      ! Exchange contribution
      CoulombExchange = CoulombEnergy_Exchange(D_I_I(:,2)) 
    else
      ! Direct contribution of the Coulomb potential
      CoulombDirect   = CoulombEnergy_Direct(ChargeDensity)
      ! Exchange contribution
      CoulombExchange = CoulombEnergy_Exchange(ChargeDensity) 
    endif

    call calcrigid()
    if(iprint.eq. 1 .or. rotcorr .eq. 1) then
      ! Only calculate these things if we are going to print observables
      ! or we need a rotational correction.
      call start_timer(T_MOI)  
      call calcJ2andBelyaev()
      call stop_timer(T_MOI)  
      call calcRotationalCorrection()
    endif

    ! Saving history
    do i=4,1,-1
        Ehistory(i+1) = Ehistory(i)
    enddo
    Ehistory(1) = TotalE    

    ! Total energy
    TotalE = sum(Skyrme + Kinetic) + sum(COMCorrection)
    TotalE = TotalE + CoulombDirect + CoulombExchange + sum(Rotcorrection)

    ! Total energy from single-particle energies
    SpwfEnergy = calcspwfenergy()

    ! Entropy calculation when temperature is finite
    call calcentropy()

    call stop_timer(T_energy)

 end subroutine CalcEnergy
 
 subroutine CompSkyrme()
    !---------------------------------------------------------------------------
    ! Calculate the Skyrme part to the functional.
    !---------------------------------------------------------------------------
    real(KIND=dp) :: Edensity(mv,3)
    integer       :: m
    
$CALCULATION    

    tot_even = &
$TOTAL_EVEN

    tot_odd  = &
$TOTAL_ODD

    Skyrme = tot_even + tot_odd
    PairDenEnergy = &
$TOTALPAIR

 end subroutine CompSkyrme
 
 subroutine PrintSkyrme()
    !---------------------------------------------------------------------------
    ! Print all contributions to the Skyrme energy, automatically generated by
    ! Hephaestos. 
    !---------------------------------------------------------------------------
    
    1 format (80('-'))
    2 format (' Skyrme Energy ')
    3 format (35x, ' Isoscalar      Isovector   |   Total')
    !4 format (17x, 'Total Skyrme:', 3f15.6)
    5 format (17x, 'Total Skyrme:', 30x, f15.6)    
   51 format (17x, '   time-even:', 30x, f15.6)
   52 format (17x, '   time-odd :', 30x, f15.6)

     print 1
     print 2
     print 3
     print 1
$PRINT
     print 1
     !print 4, Skyrme, sum(Skyrme)
     print 5, sum(Skyrme)
     print 51, sum(tot_even)
     print 52, sum(tot_odd)
     print 1
 end subroutine PrintSkyrme

 function CompKinetic() result(kinetic)
    !---------------------------------------------------------------------------
    ! This subroutine computes the total kinetic energy,
    ! according to the following formula:
    !    E_k = -\hbar/2m \int d^3x \sum_{k} v_{k} \Psi_k^* \Delta \Psi_k
    !---------------------------------------------------------------------------
    ! Note that the 1-body c.o.m. correction is not taken into account here!
    !---------------------------------------------------------------------------
    use Constants

    integer          :: wave, it,k,i
    real(KIND=dp)    :: Inproduct
    real(KIND=dp)    :: Kinetic(2)
    
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Correctly set the pointers to the spwfs
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    select case(PairingType)
    case(0,1)
      ! HF or BCS Calculation
      DenPsi   => HFPsi    ; DenDPsi   => HFDPsi 
      DenddPsi => HFddPsi  ; DendddPsi => HFdddpsi
    case(2)
      ! HFB calculation
      DenPsi    => CanPsi   ; DenDPsi   => CanDPsi 
      DenddPsi  => CanddPsi ; DendddPsi => Candddpsi
    end select

    ! Kinetic Energy
    Kinetic = 0.0_dp
    do wave=1,nwt
        ! Isospin is neutron in the first half of blocks, proton in the rest
        it = 2
        if(wave.le.sum(HFBlocks(1:Blocks/2))) it = 1

        Inproduct = 0.0_dp
        do k=1,4          
                do i=1,mv
                       Inproduct = Inproduct + DenPsi(i,k,wave) *  & 
                       &  ( DenddPsi(i,1,k,wave) + &
                       &    DenddPsi(i,4,k,wave) + &
                       &    DenddPsi(i,6,k,wave))
                enddo
        enddo
        Kinetic(it)= Kinetic(it) + rho_can(wave)*Inproduct
    enddo
    Kinetic=-Kinetic * hbm * dv
    return
  end function CompKinetic
  
  subroutine CompCOMCorrection()
    !---------------------------------------------------------------------------
    ! M. Bender et al., Eur. Phys. J. A 7, 467-478 (2000)
    !
    ! For EV8-like symmetries, we have for the two-body part
    !
    !
    ! E_ph = + sum_km v^2_k v^2_m [ Re{nabla_x},k,m . Re{nabla_x},k,m          
    !                              +Im{nabla_y},k,m . Im{nabla_y},k,m     
    !                              +Re{nabla_z},k,m . Re{nabla_z},k,m ]   
    !
    !
    ! E_pp = -2 sum_k,m>0 v_k u_k v_m u_m {-Re(nabla_x)k,m . Re(nabla_x)-k,-m   
    !                                      +Im(nabla_y)k,m . Im(nabla_y)-k,-m   
    !                                      +Re(nabla_z)k,m . Re(nabla_z)-k,-m } 
    !
    ! Note: this routine will need quite some work to be generalized to 
    !       different symmetry combinations. 
    !
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
    integer       :: it, i,j
$NTR integer       :: B, ibar, jbar ii, jj, N, N2, N3, N4, si
    real(KIND=dp) :: NablaMElements(3,2,nwt,nwt),temp(3,2), fac
    real(KIND=dp) :: Butler_t, Butler_f
    
    COMCorrection = 0.0_dp
    select case(COM1Body)
    case(0)
      ! No contribution
    case(1,2)
      ! Deduce 1-body COM correction from the Kinetic Energy
      COMCorrection(1,:) = - Kinetic(:) * nucleonmass/ &
      &                 (neutrons * nucleonmass(1) + protons * nucleonmass(2))
    case(3)
      ! Deduce 1-body COM correction from the Kinetic Energy with Butlers 
      ! formula.

      Butler_t = (1.5 * (neutrons + protons))**(1./3.)
      Butler_f = 2./(Butler_t + 1./(3*Butler_t))

      COMCorrection(1,:) = - Kinetic(:) * nucleonmass * Butler_f/ &
      &                 (neutrons * nucleonmass(1) + protons * nucleonmass(2))
    end select    

    if(COM2body .eq. 1) then
      ! Calculate 2-body COMcorrection
      NablaMElements = compNablaMelements()
      
      COMCorrection(2,:) = 0.0
      !-------------------------------------------------------------------------
      ! particle-hole, and particle-particle part.
      !-------------------------------------------------------------------------
      temp = 0

      do i=1,nwt
         ! We sum over all possible (i,j) pairs, the matrix elements are 
         ! correctly calculated either way.
         it = 1
         if(i.gt.nwn) it = 2
         do j=1,nwt  
            ! v^2 v^2 part
            fac = rho_can(i)*rho_can(j) 
$TR         fac = fac / 4.0 ! rho_can is twice too large if T is conserved
            temp(1,it) = temp(1,it) + fac*NablaMElements(1,1,i,j)**2
            temp(2,it) = temp(2,it) + fac*NablaMElements(2,2,i,j)**2
            temp(3,it) = temp(3,it) + fac*NablaMElements(3,1,i,j)**2
$TR         ! uv uv part
$TR         fac = kappa_can(i)*kappa_can(j)
$TR         temp(1,it) = temp(1,it) + fac*NablaMElements(1,1,i,j)**2
$TR         temp(2,it) = temp(2,it) + fac*NablaMElements(2,2,i,j)**2
$TR         temp(3,it) = temp(3,it) + fac*NablaMElements(3,1,i,j)**2
         enddo
      enddo

$NTR      si = 0
$NTR      do B=1,8,4 ! This is essentially an isospin loop now
$NTR        N = HFblocks(B) ; if(N.eq.0) cycle
$NTR        N2 = HFblocks(B+1)
$NTR        N3 = HFblocks(B+2)
$NTR        N4 = HFblocks(B+3)
$NTR        it = 1  ; if(B .eq.5) it = 2
$NTR        ! We simply loop over all possible combinations of spwfs
$NTR        ! as the NablaMElements array is zero in the right places
$NTR        do i=1, N+N2+N3+N4
$NTR          ii   = si + i
$NTR          ibar = conjugp(ii) ; if(ibar .eq.0) cycle
$NTR          do j=1,  N+N2+N3+N4
$NTR            jj   = si +  j
$NTR            jbar = conjugp(jj) ; if(jbar .eq.0) cycle
$NTR            fac = -  kappa_can(ii)*kappa_can(jbar)
$NTR            temp(3,it) = temp(3,it) + fac*NablaMElements(3,1,ii,jj)        &
$NTR                                  &      *NablaMElements(3,1,ibar,jbar)
$NTR            temp(1,it) = temp(1,it) + fac*NablaMElements(1,1,ii,jj) &
$NTR                                  &      *NablaMElements(1,1,ibar,jbar)
$NTR            temp(2,it) = temp(2,it) - fac*NablaMElements(2,2,ii  ,jj) &
$NTR                                  &      *NablaMElements(2,2,ibar,jbar)
$NTR          enddo
$NTR        enddo
$NTR        si = si + N + N2 + N3 + N4
$NTR      enddo

      do it=1,2
        COMCorrection(2,it) = sum(temp(:,it))
      enddo
      
      ! In the case of Time-reversal conservation, we summed over only half 
      ! the states
$TR      COMCorrection(2,:) = 2*COMCorrection(2,:)

      ! Some constants
      COMCorrection(2,:) = COMCorrection(2,:) * hbm * nucleonmass/             & 
      &                 (neutrons * nucleonmass(1) + protons * nucleonmass(2))
     endif      

  end subroutine CompCOMCorrection

  subroutine calcRotationalCorrection()
    !---------------------------------------------------------------------------
    ! Calculate the rotational correction to the energy as
    !     E_crank = - \sum_{\mu} <J_mu^2>/(2 * I_{\mu})   
    !
    ! This is ill-defined for spherical nuclei, so we use the following
    ! prescription from     
    !    D. Pena-Arteaga, EPJA 52, 320 (2016).
    ! which is
    !    E_rot = E_crank * b * tanh(c|beta_2|)
    !---------------------------------------------------------------------------
    use momentsofinertia
    use moments  

    integer       :: i
    real(KIND=dp) :: B2, A, Q2(3), damp(3), compare(3), R
    type(moment), pointer :: rms

    Rotcorrection = 0.0
    if(Rotcorr .eq. 0) return

    A = neutrons+protons    
    Q2= calculatetotalql(2) 
    B2= abs(4*pi/3. /((1.2*A**(1./3.))**2 * A) * Q2(3))

    rms => FindMoment(-2,0, .false.)
    ! We calculate the classical moment of inertia along the axis
    do i=1, 3
!      compare(i) = 1./3. * 2./5. * sum(nucleonmass * rms%value(1:2)) 
      R = 1.2 * (neutrons+protons)**(1./3.)
      compare(i) = 1./3. * 2./5. * sum(nucleonmass)/2 * (neutrons+protons)*R**2
    enddo
    ! Putting it in correct units
    compare = compare/(hbarclum**2)

    ! Another possibility is to compare the moment of inertia to that one of the
    ! rigid rotor as calculated for the density in memory.
    !compare = compare/rigid
    select case(pairingtype)
    case(0,1)
      ! HF or BCS
      damp             = rotcorrb * tanh(rotcorrc * Belyaev(:,3)/compare)
      RotCorrection    =-J2(:,3)/(2*Belyaev(:,3))*damp
    
      ! Sanity check: no collective sense of rotational correction implemented
      !               yet for HF/BCStype calculations
      if(blocktype.ne.0) then
          print *, 'Rotational correction for odd nuclei not incorporated into BCS.'      
          stop
      endif
    case (2)
      if(inversetemp.lt.0) then
        damp          = rotcorrb * tanh(rotcorrc * Bely_coll(:,3)/compare)
        RotCorrection = - J2_coll(:,3)/(2*Bely_coll(:,3))*damp
      else
        damp          = rotcorrb * tanh(rotcorrc * Belyaev(:,3)/compare)
        RotCorrection = - J2(:,3)/(2*Belyaev(:,3))*damp
      endif
    end select

  end subroutine calcRotationalCorrection

  subroutine calcFields(calcall)
    !---------------------------------------------------------------------------
    ! Calculate all of the Skyrme potentials.
    !
    ! Calcall input decides whether or not to calculate ALL fields. 
    ! If Calcall is true, all of the potentials get recalculated.
    ! If Calcall is false, only potentials that are equal to zero get calculated.
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Includes preconditioning of F_I_I at the moment only.
    !---------------------------------------------------------------------------
    use Coulombmod , only : SolveCoulomb, CoulombPotential, Exchangepotential
    use Coulombmod , only : Foldedcoul,  FoldedExchange
    use moments
    
    integer                    :: it,i,j,k, maxit
    real(KIND=dp), allocatable :: update(:,:)
    logical, intent(in)        :: calcall
    logical                    :: rhoread

    call start_timer(T_fields)
  
    ! We need to determine if F_I_I was read from file or not. 
    ! If it was, it already includes Coulomb and constraining fields and we
    ! should not add them again.
    if(.not.calcall) then
        if(allocated(F_I_I))then
            rhoread = .true.
        else
            rhoread = .false.        
        endif
    else
        rhoread = .false.
    endif

$CALCFIELDS

    !-----------------------------------------------------------------------
    ! Solve for the Coulomb Potential
    call SolveCoulomb(D_I_I(:,2))

    if(.not. rhoread) then    
        !-----------------------------------------------------------------------
        ! Add the Coulomb contribution to the field corresponding to rho.
        ! The index juggling is ugly, but necessary. The Coulomb 
        ! potential is defined on a slightly larger box using boundary 
        ! conditions. A simple abstract statement might mess this up.
        if((all(protonsize.eq.0.0) .and. all(neutronsize.eq.0.0)) .or.         &
          &                             (.not. nucleonsize_selfconsistent)) then
          ! We simply put the coulomb potential. Note that this breaks 
          ! self-consistency if protons and neutrons are not treated as 
          ! point particles.
          do k=1,nz
            do j=1,ny
              do i=1,nx
                F_I_I(i+(j-1)*nx+(k-1)*ny*nx,2)=F_I_I(i+(j-1)*nx+(k-1)*ny*nx,2)&
                &                              + CoulombPotential(i,j,k)       &
                &                              + ExchangePotential(i,j,k)
              enddo
            enddo
          enddo

        else
          ! Use the folded coulombpotential, for full self-consistency.
          ! Note that both protons and neutrons feel a Coulomb force if their
          ! charge form factor is taken into account.
          if(.not. allocated(foldedcoul)) then
            print *, 'Nucleonsize_selfconsistent cannot be .false. if the protons are not point particles.'      
            stop
          endif 
          do it=1, 2
            do k=1,nz
              do j=1,ny
                do i=1,nx
                  F_I_I(i+(j-1)*nx+(k-1)*ny*nx,it)=  &
                  &                 F_I_I(i+(j-1)*nx+(k-1)*ny*nx,it)           &
                  &                              + FoldedCoul(i,j,k,it)        &
                  &                              + FoldedExchange(i,j,k,it)
                enddo
              enddo
            enddo
          enddo
        endif 
        !-----------------------------------------------------------------------
        ! Add the contribution from the constraints on the electric multipole 
        ! moments. 
        F_I_I =  F_I_I + Constraint_I_I

        !-----------------------------------------------------------------------
        ! Add the contribution of a cranking constraint to the 
        !    F_I_S and G_I_N  fields
$NTR    F_I_S = F_I_S + crank_spin_potential()     
$NTR    G_I_N = G_I_N + crank_current_potential() 
    endif

    !---------------------------------------------------------------------------
    ! Precondition the field corresponding to rho, F_I_I.
    if(.not.all(F_I_I_hist.eq.0.0_dp) .and. potentialpreconditioning.eq.1) then
      if(.not.allocated(update)) allocate(update(nx*ny*nz,2))
    
      update=  F_I_I - F_I_I_hist
      update=  PreconditionPotential(update,-preconfactor,1.0_dp,+1,+1,+1)
      F_I_I =  F_I_I_hist + update
    endif
    !---------------------------------------------------------------------------
    ! Precondition the field corresponding to s, F_I_S.
$NTR    if(.not.all(F_I_S_hist.eq.0.0_dp) .and. potentialpreconditioning.eq.1) then
$NTR      if(.not.allocated(update)) allocate(update(nx*ny*nz,2))
    
$NTR      ! X component
$NTR      update=  F_I_S(:,1,:) - F_I_S_hist(:,1,:)
$NTR      update=  PreconditionPotential(update,-preconfactor,1.0_dp,-1,+1,-1)
$NTR      F_I_S(:,1,:) =  F_I_S_hist(:,1,:) + update

$NTR      ! Y component
$NTR      update=  F_I_S(:,2,:) - F_I_S_hist(:,2,:)
$NTR      update=  PreconditionPotential(update,-preconfactor,1.0_dp,+1,-1,-1)
$NTR      F_I_S(:,2,:) =  F_I_S_hist(:,2,:) + update

$NTR      ! Z component
$NTR      update=  F_I_S(:,3,:) - F_I_S_hist(:,3,:)
$NTR      update=  PreconditionPotential(update,-preconfactor,1.0_dp,+1,+1,+1)
$NTR      F_I_S(:,3,:) =  F_I_S_hist(:,3,:) + update

$NTR    endif
    call stop_timer(T_fields)
 
  end subroutine calcFields 
  
  function sphamil(psi, dpsi, ddpsi, dddpsi, sx,sy,sz,iso, onthefly) &
                                                                  & result(hpsi)
    !---------------------------------------------------------------------------
    ! Apply the action of the single-particle hamiltonian to the 
    ! single-particle wave-functions.
    !---------------------------------------------------------------------------
    
    use derivatives
    
    ! Logical indicating if the derivatives need to be calculated before
    ! applying h.
    ! If false, the derivatives are passed in. If True, the derivatives are not
    ! passed in and need to be calculated.
    logical, intent(in)       :: onthefly 
    
    real(KIND=dp), intent(in)    :: psi(mv,4)  
    real(KIND=dp), intent(inout) :: dpsi(mv,3,4),ddpsi(mv,6,4), dddpsi(mv,10,4)
    integer, intent(in)       :: sx(4),sy(4),sz(4),   iso
    real(KIND=dp)             :: hpsi(mv,4)
    real(KIND=dp)             :: temp(mv,4)
    real(KIND=dp)             ::   dtemp(mv,3,4)
    real(KIND=dp)             ::  ddtemp(mv,3,3,4)
    real(KIND=dp)             :: dddtemp(mv,3,3,3,4)
    real(KIND=dp)             :: laptemp(mv,4)
    
    real(KIND=dp)             :: ReducedMass, Butler_t, Butler_f
    
    integer :: it, i,k
    
    call start_timer(T_sphamil)
    !---------------------------------------------------------------------------
    ! Determine the isospin index
    it = (iso + 3)/2
    !---------------------------------------------------------------------------
    ! Reduced mass in case of self-consistent 1-body COM correction
    Reducedmass = 1.0_dp
    select case(COM1Body)
    case(0,1)      
      Reducedmass = 1.0_dp
    case(2)
      Reducedmass = (1.0_dp-nucleonmass(it)/                                   &
      &                      (neutrons*nucleonmass(1)+protons*nucleonmass(2)))
    case(3)
      
      Butler_t = (1.5 * (neutrons + protons))**(1./3.)
      Butler_f = 2./(Butler_t + 1./(3*Butler_t))
      Reducedmass = (1.0_dp-nucleonmass(it) * Butler_f/                        &
      &                      (neutrons*nucleonmass(1)+protons*nucleonmass(2)))
    end select
    
    if(OnTheFly) then
      ! Calculate the derivatives
      do k=1,4
!-------------------------------------------------------------------------------
$N2        call Derive_tot(psi(:,k),sx(k),sy(k),sz(k),dpsi(:,:,k),ddpsi(:,:,k))
$N3        call Derive_tot(psi(:,k),sx(k),sy(k),sz(k),dpsi(:,:,k),ddpsi(:,:,k),&
$N3        &                                     dddpsi(:,:,k))
!-------------------------------------------------------------------------------
        enddo
    endif
    !---------------------------------------------------------------------------
    ! Action of the kinetic energy
    do k=1,4
        do i=1,mv
            hpsi(i,k) = - hbm(it)* reducedmass *       (ddpsi(i,1,k) &
            &                                         + ddpsi(i,4,k) &
            &                                         + ddpsi(i,6,k)) 
        enddo
    enddo
    !---------------------------------------------------------------------------
    ! Action of the Skyrme fields
    !
    ! Note that 
    ! a) Coulomb is included in the F_I_I field
    ! b) Every density contains the contributions from constraints on that 
    !    density. 
    ! c) The kinetic energy is NOT included in the F_N_N field, because 
    !    a constant is not in the Lagrange basis; so the current way of
    !    deriving stuff is not correct for a term of the form
    !          hbar^2_2m
    !---------------------------------------------------------------------------
$SKYRMEACTION

    call stop_timer(T_sphamil)

  end function sphamil
  
  function delta_action(psi, dpsi, ddpsi, dddpsi, sx,sy,sz,iso, onthefly)      &
  &                                                             result(deltapsi)
    !---------------------------------------------------------------------------
    !
    ! onthefly:
    !   Logical indicating if the derivatives need to be calculated before
    !   applying delta. If false, the derivatives are passed in. If True, the 
    !   derivatives are not passed in and need to be calculated.
    !---------------------------------------------------------------------------
    logical, intent(in)       :: onthefly 
    
    real(KIND=dp), intent(in)    :: psi(:,:)  
    real(KIND=dp), intent(inout) :: dpsi(:,:,:),ddpsi(:,:,:), dddpsi(:,:,:)
    integer, intent(in)        :: sx(:),sy(:),sz(:),   iso
    real(KIND=dp), allocatable :: deltapsi(:,:)
    real(KIND=dp)              ::    temp(mv,4)
    real(KIND=dp)              ::   dtemp(mv,3,4)
    real(KIND=dp)              ::  ddtemp(mv,3,3,4)
    real(KIND=dp)              :: dddtemp(mv,3,3,3,4)
    real(KIND=dp)              :: laptemp(mv,4)
    integer                    :: it,i
    
    !---------------------------------------------------------------------------
    ! Determine the isospin index
    it = (iso + 3)/2
    
    !---------------------------------------------------------------------------
    ! Zero the action of Delta. 
    ! This is the place to include contributions to the pairing that should 
    ! be coded manually
    allocate(deltapsi(mv, 4))
    deltapsi = 0.0
   
$PAIRINGACTION
   
  end function delta_action
  
  function calcspwfenergy() result(spwfenergy)
    !---------------------------------------------------------------------------
    ! Calculates the total energy from the single-particle energies. 
    !
    !---------------------------------------------------------------------------
    
    use wavefunctions
    use moments
    
    integer       :: wave
    real(KIND=dp) :: spwfenergy, e_rear
    
    ! Start by summing the single-particle energies
    spwfenergy = 0 ; e_rear = 0
    do wave=1,nwt
        if(pairingtype.lt.2) then
          spwfenergy = spwfenergy + rho_can(wave) * spenergies(wave)
        else
          spwfenergy = spwfenergy + rho_can(wave) * canenergies(wave)
        endif
    enddo
    
    ! Calculation of rearrangement energy (without Coulomb Exchange)
$EREAR   

    ! Add the rearrangement energy 
    spwfenergy = spwfenergy - e_rear 
    
    ! Add everything and don't forget about Coulomb exchange
    spwfenergy = 0.5 * spwfenergy + 0.5 * sum(kinetic) + CoulombExchange/3.d0
    
    ! Always add the 1-body COMcorrection. In case it is used iteratively, it
    ! is double counted along with the kinetic energy!
    if(COM1body.gt.0) then
        SpwfEnergy = SpwfEnergy  + sum(COMCorrection(1,:))/2.0_dp
    endif

    if(COM2body.gt.0) then
        SpwfEnergy = SpwfEnergy  + sum(COMCorrection(2,:))  
    endif
    
    ! Subtract contribution by multipole constraints
    SpwfEnergy = SpwfEnergy - sum(Constraint_I_I * D_I_I)*dv/2.0_dp

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
    
  end function calcspwfenergy
  
  subroutine output_Edensity(Edensity, N)
    !---------------------------------------------------------------------------
    ! Write the energydensity to a file with name N.
    !
    !
    !---------------------------------------------------------------------------
    character(len=*), intent(in)   :: N
    real(KIND=dp), intent(in), target ::  Edensity(nx*ny*nz,3)
    real(KIND=dp), pointer :: w(:,:,:,:)
    
    real(KIND=dp) :: r,x
    integer       :: i
    
    w(1:nx,1:ny,1:nz,1:3) => Edensity
    
    open(12, File=N)
    
    do i=1,nx
      x = dx/2 + (i-1)*dx
      r = sqrt(3*x**2)
      write(12, '(5f10.5)') r, w(i,i,i,1),  w(i,i,i,2),  w(i,i,i,3) 
    enddo
    close(12)
  end subroutine output_Edensity

  subroutine clean_potentials()
    !---------------------------------------------------------------------------
    ! Clean up the allocated potentials for multiple runs.
    !---------------------------------------------------------------------------
    $CLEANING

  end subroutine clean_potentials

  subroutine WritePotentials(chan)
    !---------------------------------------------------------------------------
    !  Subroutine writing the different potentials to file.
    !---------------------------------------------------------------------------
    integer, intent(in) :: chan
    integer             :: io

    ! Signalling how many fields have been stored.
    write(chan, iostat=io) $FIELDNUMBER

    ! Then, for every potential write the 
    ! * Name 
    ! * Value
    ! Note that the name is written as a length-30 string, padded with spaces.
    ! If not, the unformatted in/out cannot correctly determine the end of a
    ! string and comparisons can not be made.
$WRITEPOTENTIALS
  end subroutine WritePotentials

  subroutine ReadPotentials(chan, filenx, fileny, filenz)
    !---------------------------------------------------------------------------
    !  Subroutine writing the different potentials to file.
    !---------------------------------------------------------------------------
    integer, intent(in) :: chan, filenx, fileny, filenz
    integer             :: io, fieldnumber, fieldcount, it, filemv
    character(len=30)   :: fieldname

    filemv = filenx * fileny * filenz

    ! Checking how many fields have been stored
    read(chan, iostat=io) fieldnumber

    do fieldcount = 1,fieldnumber
        ! Read the fieldname
        read(chan, iostat=io) fieldname   
        ! Manually check for 

        select case(trim(fieldname))
$READPOTENTIALS
        CASE DEFAULT
          ! The potential is not in this program, forget about it
          read(chan, iostat=io)
        end select
    enddo

  end subroutine ReadPotentials

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

end module functional
