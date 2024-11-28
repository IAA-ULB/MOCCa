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
 ! Module containing the means to calculate (and print) the mean-field energy, 
 ! as well as all the mean-field potentials.
 !
 ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 ! Hephaestos keywords
 ! 
 ! NTERMS           : $NTERMS
 ! QUADRI           : $QUADRI
 ! Declaration      : [WAY TOO LONG TO INCLUDE HERE]
 ! Calccoef         : [WAY TOO LONG TO INCLUDE HERE]
 ! Calculation      : [WAY TOO LONG TO INCLUDE HERE]
 ! Total_EVEN       : [WAY TOO LONG TO INCLUDE HERE]
 ! Total_ODD        : [WAY TOO LONG TO INCLUDE HERE]
 ! Total_BI         : [WAY TOO LONG TO INCLUDE HERE]
 ! Total_TRI        : [WAY TOO LONG TO INCLUDE HERE]
 ! Total_QUAD       : [WAY TOO LONG TO INCLUDE HERE]
 ! Total_DD         : [WAY TOO LONG TO INCLUDE HERE]
 ! PAIRTOTAL_PROTON : [WAY TOO LONG TO INCLUDE HERE]
 ! PAIRTOTAL_NEUTRON: [WAY TOO LONG TO INCLUDE HERE]
 !
 ! PrintCOEF_iso    : [WAY TOO LONG TO INCLUDE HERE]
 ! PrintCOEF_pair   : [WAY TOO LONG TO INCLUDE HERE]
 ! Print            : [WAY TOO LONG TO INCLUDE HERE]
 ! Calcpotentials   : [WAY TOO LONG TO INCLUDE HERE]
 ! SkyrmeAction     : [WAY TOO LONG TO INCLUDE HERE]
 ! PairingAction    : [WAY TOO LONG TO INCLUDE HERE]
 ! ERear            : [WAY TOO LONG TO INCLUDE HERE]
 !
 ! POTENTIALPRECON  : [WAY TOO LONG TO INCLUDE HERE]
 ! CLEANING         : [WAY TOO LONG TO INCLUDE HERE]
 ! POTENTIALNUMBER  : [WAY TOO LONG TO INCLUDE HERE]
 ! WRITEPOTENTIALS  : [WAY TOO LONG TO INCLUDE HERE]
 ! READPOTENTIALS   : [WAY TOO LONG TO INCLUDE HERE]
 !
 ! TR               : $TR
 ! NTR              : $NTR
 ! N2               : $N2
 ! N3               : $N3
 ! N1DELTA          : $N1DELTA
 ! N2DELTA          : $N2DELTA
 ! N3DELTA          : $N3DELTA
 ! SYMDELTA         : $SYMDELTA
 !
 ! D2TEMPSPH        : $D2TEMPSPH
 ! D3TEMPSPH        : $D3TEMPSPH
 ! LAPTEMPSPH       : $LAPTEMPSPH
 !
 ! D1TEMPDELTA      : $D1TEMPDELTA
 ! D2TEMPDELTA      : $D2TEMPDELTA
 ! D3TEMPDELTA      : $D3TEMPDELTA
 ! LAPTEMPDELTA     : $LAPTEMPDELTA
 !
 ! TAUSCALAR        : $TAUSCALAR
 ! TAUTENSOR        : $TAUTENSOR
 !
 ! K2POT            : [WAY TOO LONG TO INCLUDE HERE]
 ! K4POT            : [WAY TOO LONG TO INCLUDE HERE]
 !------------------------------------------------------------------------------
 ! A density F_L_R is stored as
 !
 !      F_L_R (mv, [cartesian indices], [isospin indices])
 !             |     |                            |
 !             > spatial indices                  |
 !                   |                            |
 !                   > all cartesian indices      |
 !                                                > isospin indices
 !                                                  4 for normal densities
 !                                                     (n, p, 0, 1)
 !                                                  2 for pairing densities
 !                                                     (n,p)
 !
 !==============================================================================
 
 use compilation
 use geninfo
 use densities
 use parameterization
 use pairing
 use timing
 use transform
 use vectors
 use Cranking
 use pairing_strengths
 use HDF5

 implicit none
 
    !===========================================================================
    ! PARAMETERIZATION DEFINITION OPTIONS
    !===========================================================================
    !---------------------------------------------------------------------------
    ! Name of the parameterization
    character(len=30) :: name_param 
    ! Name of the parameterization used to generate the .wf file
    character(len=30) :: ini_name_param=''
    ! Name of the functional file this code was compiled with
    character(len=20), parameter :: func_name = $FUNC_NAME
    !---------------------------------------------------------------------------
    ! Definition of global contributions to the energy
    real(KIND=dp) :: Kinetic(2), Skyrme, TotalE, Ehistory(5)
    real(KIND=dp) :: ElectronEnergyKin, ElectronEnergyExch
    real(KIND=dp) :: ElectronChempotKin, ElectronChempotExch
    real(KIND=dp), parameter :: Qnp=1.29335236 !Mn-Mp
    real(KIND=dp) :: tot_even  , tot_odd
    real(KIND=dp) :: bilinear, trilinear, quadrilinear, densitydependent
    real(KIND=dp) :: COMCorrection(2,2), CoulombDirect, CoulombExchange
    ! Separation of 2-body Centre-of-mass correction into particle-hole 
    ! and pairing parts for diagnostic printing
    real(KIND=dp) :: COM2pp(2), COM2ph(2)
    ! Debugging quantities for the COM2body correction: a separation of all
    ! these energies according to Cartesian direction
    real(KIND=dp) :: COM2_pp_debug(3,2), COM2_ph_debug(3,2)
    ! Two definitions of the pairingenergy: one obtained by summing the gaps
    ! and one by integrating the particle-particle part of the functional
    real(KIND=dp) :: PairingEnergy(2), PairDenEnergy(2)
    ! Same thing, but with the added stabilization
    real(KIND=dp) :: PairE_stab(2), PairDenE_stab(2)
    !---------------------------------------------------------------------------
    ! Rotational correction
    real(KIND=dp) :: RotCorrection(3)
    ! and vibrational correction
    real(KIND=dp) :: VibCorrection(3)
    !---------------------------------------------------------------------------
    ! Value of the Routhian 
    real(KIND=dp) :: Routhian, RHistory(5) 
    ! Value of the free energy F = E - TS when finite temperature is active
    real(KIND=dp) :: FreeEner, FHistory(5)
    ! Value of the energy as calculated from the spwfs
    real(KIND=dp) :: SpwfEnergy, SpwfHistory(5)
    !===========================================================================
    ! NUMERICAL OPTIONS
    !===========================================================================
     ! level of compression in hdf5, 6 seems to be the best
    integer, parameter  :: comprlvl = 6
    !---------------------------------------------------------------------------
    ! Numerical parameter of the preconditioning of the potentials
    real(KIND=dp) :: preconfactor = 4.0_dp
    !---------------------------------------------------------------------------
    ! Stabilisation factor for the pairing:
    !    f = E_cut^2 / E_pair^2
    real(KIND=dp) :: pairstabfactor(2) = 0.0
    !===========================================================================
    ! Potential vectors that are predefined
    !===========================================================================
    type(PotentialVector) :: potentials,potentials_out
    type(PotentialVector) :: potentials_read
    type(PotentialVector), allocatable :: potential_iterates(:)
    type(PotentialVector), allocatable :: potential_updates(:)
    !===========================================================================
    ! SINGLE-PARTICLE POTENTIALS associated with other parts of the code
    ! W.R. 28/01/2023: moved here for consistency reasons.
    !===========================================================================
    !-----------------------------------------------------------------------------
    ! Contribution to the single-particle Hamiltonian by the constraints
    !  a) Electric multipole => Constraint_I_I => F_I_I
    !-----------------------------------------------------------------------------
    real(kind=dp), allocatable, target :: Constraint_I_I(:,:)
    !===========================================================================
    ! AUTOMATICALLY GENERATED DECLARATIONS
    !===========================================================================
    !---------------------------------------------------------------------------
    ! Declaration of
    ! a) Coupling constants
    ! b) Terms in the energy
    ! All automatically generated by Hephaestos
    !---------------------------------------------------------------------------
    real(KIND=dp) :: coupl_constant($NTERMS) = 0
$DECLARATION
    !---------------------------------------------------------------------------
    
   interface operator (+)
      !Overloading "+" to be used to add potential vectors.
      module procedure Add_potentialvector
   end interface

   interface operator (*)
      !Overloading "*" to be used to multiply potential vectors with scalars
      module procedure multiply_potentialvector
   end interface
  
contains

 subroutine readfunctional(file_number)
    !---------------------------------------------------------------------------
    ! Initializes the functional
    ! a) read the details of the parameterization from file
    ! b) calculate the coupling constants
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !   file_number : optional integer. If present, read from (open) channel
    !                 with this number. If absent, read from STDIN.
    !---------------------------------------------------------------------------
    integer(dp), intent(in), optional   :: file_number 
#if(USE_MPI>0)
    integer                             :: mpi_err
#endif

    namelist /func/ name_param
    
    if(MPI_RANK.eq.0) then
      ! Only the very first MPI rank reads stuff
      if(present(file_number)) then
        read(unit=file_number, nml=func) 
      else
        read(unit=*, nml=func) 
      endif
    endif
#if(USE_MPI > 0)
    ! broadcasting the name of the parameterization for consistency
    call MPI_Bcast(name_param, len(name_param), MPI_CHARACTER, 0, &
    &                                                   MPI_COMM_WORLD, mpi_err)
#endif

    ! Only the very first MPI rank goes on to read the .param file, but this
    ! is handled inside the readparameterization subroutine
    call readparameterization(name_param, func_name)

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Bookkeeping operations, including the calculation of the coupling 
    ! coefficients, that are to be executed by all MPI ranks 
    call calcedfcoefs()
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
    2 format (2x, 78('_'))
    3 format (' EDF coupling constants ')
    4 format (40x, 'Particle-hole terms')
   41 format (2x, 'Term', 35x, 'Isospin', 5x, '# #G', 5x,' Value ')
    5 format (40x, 'Pairing terms')

   97 format (2x, a38,'|', 2a2, 5x ,'|', 2i3,'|', 1x, f21.12)
   98 format (2x, a38,'|', 3a2, 3x ,'|', 2i3,'|', 1x, f21.12)
$QUADRI   99 format (2x, a38,'|', 4a2, 1x ,'|', 2i3,'|', 1x, f21.12)
    
     print 1
     print 3
     print 2
     print 4
     print 41
     print 2
     ! This part is automatically generated by Hephaestos.
$PRINTCOEF_PH
     print 2
     print 5
     print 2
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
  621 format (15x, '             ph:', 3f15.6)
  622 format (15x, '             pp:', 3f15.6)

   63 format (15x, '  Rotational  '  , a1, ':', 30x, f15.6)
  631 format (15x, '  Rotational  T:',          30x, f15.6)
  
   64 format (15x, '  Vibrational '   , a1, ':', 30x, f15.6)
  641 format (15x, '  Vibrational T:',          30x, f15.6)

   65 format (15x, '  Collective T: ',          30x, f15.6)
  
    7 format (15x, ' Coulomb Direct:', 3f15.6)
   !71 format (15x, '   Dir. (point):', 3f15.6)
    8 format (15x, '       Exchange:', 3f15.6)
   !81 format (15x, '   Exc. (point):', 3f15.6)  

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
  104 format (15x, '          dE   :', 30x, e15.6)
  105 format (15x, '       Routhian:', 30x, f15.6)
  106 format (15x, '          dR   :', 30x, e15.6)           

#if(PASTA > 0)
  107 format (30x, '         FOR PASTA CALCULATIONS    ')
  108 format (15x, '        e_pasta=(Total energy + electrons - Z[Mn-Mp])/A - Mn')
  109 format (15x, '        e_pasta:', 30x, f15.6)
  110 format (15x, '   Electron kin:', 30x, f15.6) 
  111 format (15x, '  Electron exch:', 30x, f15.6) 
  112 format (15x, 'Chempot_e total:', 30x, f15.6)
  113 format (15x, '      Chempot_n:', 30x, f15.6)
  114 format (15x, '      Chempot_p:', 30x, f15.6)
  115 format (15x, ' Chempot_p β-eq:', 30x, f15.6)
  116 format (15x, '       Pressure:', 30x, f15.6)
#endif

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Only the very first MPI rank needs to print to STDOUT
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    call printSkyrme

    print 1
    print 5
    print 6, Kinetic, sum(Kinetic)
    print 61, COMcorrection(1,:), sum(COMcorrection(1,:))
    if(any(COMcorrection(2,:).ne.0)) then
     print 62 , COMcorrection(2,:), sum(COMcorrection(2,:))
     print 621, COM2ph(:), sum(COM2ph(:))
     print 622, COM2pp(:), sum(COM2pp(:))

     !print *, 'DEBUG'
     !print ('(a3, 3f15.6)'), 'ph X',COM2_ph_debug(1,:)
     !print ('(a3, 3f15.6)'), 'ph Y',COM2_ph_debug(2,:)
     !print ('(a3, 3f15.6)'), 'ph Z',COM2_ph_debug(3,:)
     !print ('(a3, 3f15.6)'), 'pp X',COM2_pp_debug(1,:)
     !print ('(a3, 3f15.6)'), 'pp Y',COM2_pp_debug(2,:)
     !print ('(a3, 3f15.6)'), 'pp Z',COM2_pp_debug(3,:)
     !print *
    endif

    if(rotcorr .ne.  0) then  
      print 65, sum(Rotcorrection) + sum(Vibcorrection)
      print *
      print 631, sum(Rotcorrection)
      print 63, 'X',  Rotcorrection(1)
      print 63, 'Y',  Rotcorrection(2)
      print 63, 'Z',  Rotcorrection(3)
      print *
      print 641, sum(Vibcorrection)
      print 64, 'X',  Vibcorrection(1)
      print 64, 'Y',  Vibcorrection(2)
      print 64, 'Z',  Vibcorrection(3)
      print *
    endif
  
    print *
    print 7, 0.0, CoulombDirect, CoulombDirect
    !if(protonsize(1).ne.0 .and. (.not. nucleonsize_selfconsistent)) then
    !  temp = CoulombEnergy_Direct(Density, Potentials)
    !  print 71, 0.0, temp, temp
    !endif
    print 8, 0.0, CoulombExchange, CoulombExchange
    !if(protonsize(1).ne.0 .and. (.not. nucleonsize_selfconsistent)) then
    !  temp = CoulombEnergy_Exchange(Density)
    !  print 81, 0.0, temp, temp
    !endif

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
    print 100, spwfenergy
    print 103, TotalE - spwfenergy

    if(rotcorr.ne.0) then
        print 991, totalE - sum(rotcorrection)      &
        &                 - sum(COMcorrection(2,:)) & 
        &                 - sum(vibcorrection)
    endif
      
    if(inversetemp .ne. -1) then
        ! F = E - T * S
        print 101, FreeEner
        print 102, entropy, sum(entropy)
    endif
    print 105, Routhian

    print 104,  TotalE   - Ehistory(1)
    print 106,  Routhian - Rhistory(1)
#if(PASTA==1)
    print 1
    print 107
    print 1
    print 108
    print 109, (TotalE+ElectronEnergyKin+ElectronEnergyExch  &
    &                   -protons*Qnp)/dble(protons+neutrons)
    print 110, ElectronEnergyKin
    print 111, ElectronEnergyExch
    print 112, ElectronChempotKin+ElectronChempotExch
    print 113, FermiEnergy(1)
    print 114, FermiEnergy(2)
    print 115, FermiEnergy(1)-ElectronChempotKin-            &
    &                     ElectronChempotExch+Qnp
    print 116, (-TotalE-ElectronEnergyKin-ElectronChempotExch&
    &                          +dble(neutrons)*FermiEnergy(1)&
    &                          +dble(protons)*(FermiEnergy(2)&
    &       +ElectronChempotKin+ElectronChempotExch))/(mv*dv)
#endif

    print 1
end subroutine PrintEnergy

subroutine save_potential_history(F_in, F_out)
    !---------------------------------------------------------------------------
    ! Update the history of the mean-field densities, before calculating an
    ! updated version.
    !---------------------------------------------------------------------------
    integer :: i
    type(PotentialVector), intent(in) :: F_in, F_out
  
    if(.not.allocated(Potential_updates)) allocate(Potential_updates(memory))
    if(.not.allocated(Potential_iterates)) allocate(Potential_iterates(memory))

    do i=1,memory-1
      Potential_updates(memory-i+1)  = Potential_updates(memory-i)
      Potential_iterates(memory-i+1) = Potential_iterates(memory-i)
    enddo
    Potential_iterates(1)  = F_in
    Potential_updates(1)   = F_out + (-1.0d0) * F_in
!  print *, 'SAVING', maxval(abs(Potential_updates(1)%F_I_I))
!  print *,  maxval(abs(F_in%F_I_I))
!  print *,  maxval(abs(F_out%F_I_I))
end subroutine save_potential_history

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

 function calcRouth_onthefly(Rin, Fin) result(routh)
    !---------------------------------------------------------------------------
    ! This routine calculates the total Routhian for a given set of densities.
    ! 
    ! Note: 
    !   (*) It only includes contributions that are included in the calculation
    !       in a variational way, i.e. any rotational correction or two-body
    !       COM-correction is not included. 
    ! 
    !   (*) This routine has NO side-effects, i.e. it calculates the total 
    !       routhian and ONLY the total Routhian. No module-level variables
    !       are set.
    !
    ! ATTENTION: testing routine, results not guaranteed!
    !
    !
    !---------------------------------------------------------------------------
    use Coulombmod
    use Densities

    real(KIND=dp)                     :: routh
    type(DensityVector), intent(in)   :: Rin
    type(PotentialVector), intent(in) :: Fin
    
    real(KIND=dp) :: K(2), TE, TO, PE(2), S, C, CE, COM(2,2)

    routh = 0.0d0

    ! Kinetic energy
    K = CompKinetic_density(Rin)
    ! Note that this estimate of the energy does NOT include the 2-body COM
    ! To save CPU cycles, we override its calculation
#if(PASTA == 0)
    call compComCorrection(K, .true., COM)
#else
    ! A waste of CPU time for pasta calculations
    COM = 0.0d0
#endif    
    ! Skyrme energy
    call compSkyrme(Rin, Fin, S, TE, TO, PE)
    ! Coulomb energies
    C  = CoulombEnergy_Direct(Rin,Fin)
    CE = CoulombEnergy_Exchange(Rin) 

    ! The total energy is comprised of 
    !      Kinetic part + Skyrme part + corrections + Coulomb energy
    routh =sum(S + K) + sum(COM) + C + CE
    ! Note again: no 2-body COM or rotational correction in here

    ! Addition of the various constraints
    ! TODO: correct this!
    !routh = routh  + sum(Constraint_I_I * Rin%D_I_I)*dv/2.0_dp
    ! Missing = cranking contribution and correct multipole contribution


 end function calcRouth_onthefly
 
 subroutine CalcEnergy(Rin, Fin, calc_expensive)
    !---------------------------------------------------------------------------
    ! Calculate (i)   the energy
    !           (ii)  the Routhian
    !           (iii) the free energy (when T!= 0)
    ! for a given set of mean-field densities.
    !
    ! Input: 
    !    Rin           : density-vector 
    !    calc_expensive: controls the calculation of the numerically expensive 
    !                    parts of the total energy. 
    !                    Right now these are:
    !                      (i) the 2-body centre-of-mass correction 
    !                     (ii) the rotational correction 
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
    use momentsofinertia
    use Coulombmod
    use densities

    logical, intent(in) :: calc_expensive
    type(DensityVector), intent(in)   :: Rin
    type(PotentialVector), intent(in) :: Fin

    call start_timer(T_energy)
    
    !---------------------------------------------------------------------------
    ! First we calculate all the individual terms/parts

    ! Kinetic energy
    Kinetic = CompKinetic_density(Rin)
    ! COM correction (attention to input!)
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
    CoulombDirect   = CoulombEnergy_Direct(Rin, Fin)
    
    ! Exchange contribution
    CoulombExchange = CoulombEnergy_Exchange(Rin) 

    call calcrigid(Rin)
    if(calc_expensive) then
#if(PASTA == 0)
      ! Collective correction
      call start_timer(T_MOI)  
      call calcJ2andBelyaev()
      call stop_timer(T_MOI)  
      call calcRotationalCorrection()
#else
      Vibcorrection = 0.0d0
      Rotcorrection = 0.0d0
#endif
    endif
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
    !      Kinetic part + Skyrme part + corrections + Coulomb energy
    TotalE = Skyrme + sum(Kinetic) + sum(COMCorrection)
    TotalE = TotalE + CoulombDirect + CoulombExchange 
    ! Plus schematic corrections for the collective energy
    TotalE = TotalE + sum(Rotcorrection) + sum(Vibcorrection)

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

    S = tot_even + tot_odd

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
    
    1 format (80('-'))
    3 format (' Skyrme Energy',16x, 'Isospin  1 2 3 4', 19x, 'Energy [MeV]')
    !4 format (17x, 'Total Skyrme:', 3f15.6)
    5 format (17x, 'Total Skyrme:', 31x, f15.6)    

   51 format (17x, '   time-even:', 31x, f15.6)
   52 format (17x, '   time-odd :', 31x, f15.6)

   53 format (17x,     '   bilinear :', 31x, f15.6)
   54 format (17x,     '   trilinear:', 31x, f15.6)
   55 format (17x,     'quadrilinear:', 31x, f15.6)
   56 format (12x, 'density-dependent:', 31x, f15.6)
   
   97 format ( a38, 2a2, 19x, f15.6)
   98 format ( a38, 3a2, 17x, f15.6)
$QUADRI   99 format ( a38, 4a2, 15x, f15.6)

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

    ! Kinetic Energy
    Kinetic = 0.0_dp
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
#if(USE_MPI > 0)
    ! Sum the contributions across all MPI ranks
    call MPI_ALLREDUCE(MPI_IN_PLACE, Kinetic, 2, MPI_REAL8, MPI_SUM,           &
    &                                                   MPI_COMM_WORLD, mpi_err)
#endif

    Kinetic=-Kinetic * hbm * dv
    return
  end function CompKinetic_spwfs

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

  function effmass_pot(F_in) result(em_pot)
    !---------------------------------------------------------------------------
    ! From a given set of mean-field potentials F_in, determine the F^{(N,N)}
    ! potential in infinite nuclear matter with the same density at each 
    ! mesh point. This is abstracted in a routine since 
    !
    ! (a) F_Nm_Nm is not always represented in a calculation
    ! (b) This is "cooking" in the sense that there is probably no unique 
    !     recipe; it might be changed in the future.
    !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - --
    ! Input:
    !  F_in   : potentialvector
    ! Output:
    !  em_pot : the potential associated with F^{N,N} in infinite nuclear matter
    !---------------------------------------------------------------------------
    type(PotentialVector), intent(in) :: F_in
    real(KIND=dp) :: em_pot(mv,4)
    
    !----------------------------------------------------------------------------
    ! F_Nm_Nm is represented explicitly, just use that
    !----------------------------------------------------------------------------
$TAUSCALAR em_pot = F_in%F_Nm_Nm
    !----------------------------------------------------------------------------
    ! F_N_N is represented; we take the scalar part of this tensor
    !----------------------------------------------------------------------------
$TAUTENSOR em_pot = (F_in%F_N_N(:,1,1,:)+F_in%F_N_N(:,2,2,:)+F_in%F_N_N(:,3,3,:))/3
    
  end function effmass_pot

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

      Butler_t = (1.5 * (neutrons + protons))**(1./3.)
      Butler_f = 2./(Butler_t + 1./(3*Butler_t))

      COMCorr(1,:) = - Kin(:) * nucleonmass * Butler_f/ &
      &                 (neutrons * nucleonmass(1) + protons * nucleonmass(2))
    end select    

    call stop_timer(T_com1)
    
    ! The calculations is not yet implemented for Hartree-Fock calculations
$NTR    if(COM2body .ne. 0 .and. pairingtype .eq. 0) then
$NTR      call stp('Two-body COM not implemented yet for Hartree-Fock &
$NTR             & calculations with time-reversal breaking.')
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
$TR         fac = fac / 4.0 ! rho_can is twice too large if T is conserved
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

  subroutine calcRotationalCorrection()
    !---------------------------------------------------------------------------
    ! Calculate a phenomenological collective correction:
    !
    !    E_corr = - \sum_{\mu} (f^rot_mu  +  f^vib_mu ) <J_mu^2>/(2 * I_{\mu})   
    !
    ! where
    !      * the sum is over all three Cartesian directions
    !      *  < J^2_{mu} > is the expectation value of the angular momentum
    !                      squared in a given direction
    !      *  I_mu is the Belyaev moment of inertia along a given direction
    !     
    ! This incorporates more than 'just' the rotational correction: the factors 
    ! f have different interpretation:
    !
    !      * f^rot_mu is a cutoff function for the rotational correction
    !      * f^vib_mu is a modification of the rotational correction, in order
    !        to mimic a vibrational correction.
    ! 
    ! We take for both f-values
    ! 
    !      f^rot_mu = b tanh( c B_mu )  
    !      f^vib_mu = d B_mu exp ( -l  (B - b_vib)**2  ) 
    !
    ! where 
    !
    !      B_mu = I_mu / I_c 
    !
    ! is the ratio between the calculated Belyaev moment of inertia and (one
    ! third of) the classical moment of inertia: 
    !
    !      I_c = 2/15 * m_n * A * (1.2 * A)**2/(hbar c**2)
    !
    ! All of this is determined by five parameters: b, c, d, l and B_vib.
    !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
    ! 
    ! References:
    !   G. Scamps, S. Goriely, E. Olsen, M. Bender and W. Ryssens, PRC XX (2021) 
    ! & S. Goriely, M. Samyn and J. M. Pearson, PRC 75, 065312 (2007).
    !---------------------------------------------------------------------------
    use momentsofinertia
    use moments  

    integer       :: i
    real(KIND=dp) :: B(3), A,  compare(3), R, f_rot(3), f_vib(3)
    real(KIND=dp) :: bely(3), J2_temp(3)

    Rotcorrection = 0.0
    Vibcorrection = 0.0
    if(Rotcorr .eq. 0) return
  
    !---------------------------------------------------------------------------
    ! Selecting the right quantities to use for the moment of inertia and <J^2>
    select case(pairingtype)
    case(0,1)
      ! HF or BCS
      Bely    = Belyaev(:,3)
      J2_temp = J2(:,3)    
      ! Sanity check: no collective sense of rotational correction implemented
      !               yet for HF/BCStype calculations
      if(blocktype.ne.0) then
        call stp('Rotational correction for odd nuclei not incorporated into HF/BCS.')
      endif
    case (2)
      ! HFB
      if(inversetemp.lt.0) then
        Bely    = Bely_coll(:,3)
        J2_temp = J2_coll(:,3)
      else
        Bely    = Belyaev(:,3)
        J2_temp = J2(:,3)
      endif
    end select

    !---------------------------------------------------------------------------
    ! We calculate the classical moment of inertia along every Cartesian axis
    A = neutrons+protons    

    do i=1, 3
      R = 1.2 * (neutrons+protons)**(1./3.)
      compare(i) = 1./3. * 2./5. * sum(nucleonmass)/2 * (neutrons+protons)*R**2
    enddo
    ! Putting it in correct units
    compare = compare/(hbarclum**2)
    
    !---------------------------------------------------------------------------
    ! Actual calculation
    B = Bely/compare
    f_rot         = rotcorrb * tanh(rotcorrc * B)
    f_vib         = vibcorrd * B * exp( - vibcorrl * (B - vibcorrb)**2)
    RotCorrection = - f_rot * J2_temp/(2*Bely)
    
    VibCorrection = - f_vib * J2_temp/(2*Bely)

  end subroutine calcRotationalCorrection

  function calcPotentials(R, Fread, Coulomb_guess) result(F)
    !---------------------------------------------------------------------------
    ! Calculate all of the Skyrme potentials and the Coulomb potential
    ! corresponding to the input densityvector R.
    !
    ! Input: 
    !    R        : density-vector for the calculation of the fields
    !    Fread    : optional Potentialvector.
    !               If present, this potentialvector is used to initialize 
    !               the calculated potentialvector; only the fields that 
    !               are detected to be completely zero are then recalculated
    !    Coulomb_guess : optional initial guess for the Coulomb potential
    ! Output: 
    !    F        : a fresh potentialvector 
    !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! 
    ! Notes
    ! (1)  this function does not do any preconditioning, mixing or anything 
    !      else to the calculated fields. They are calculated "as is" from
    !      the density vector passed in.
    !
    ! (2)  The Coulomb potential is part of the field-vector, and hence 
    !      gets calculated here as well. 
    !---------------------------------------------------------------------------
    use Coulombmod , only : SolveCoulomb
    use Coulombmod , only : Coulomb_read_from_file
    use Coulombmod , only : coul_offset_x, coul_offset_y, coul_offset_z
    use pairing_strengths, only : vmicro, vmicro_stored
    use moments

    type(DensityVector), intent(in)             :: R
    type(PotentialVector), intent(in), optional :: Fread
    type(PotentialVector)                       :: F
    real(KIND=dp), intent(in), optional         :: Coulomb_guess(:,:,:)    
    integer                                     :: it,i,j,k, ox, oy, oz

    call start_timer(T_potentials)

    ! Copy all the relevant information if this information was presented
    if(present(Fread)) F = Fread

    ! Signal that microscopic pairing strengths have to be recalculated
    vmicro_stored = .false.

$CALCPOTENTIALS
    
    !---------------------------------------------------------------------------
    ! Additions to the potential F_I_I associated with the density
    ! (1) Coulomb potential, direct and exchange
    ! (2) Constraints
    !
    ! and to F_I_S and G_I_N: 
    ! (1) cranking potential
    !---------------------------------------------------------------------------
    if((.not. present(Fread)) .or. (.not. Coulomb_read_from_file)) then
      if(present(Coulomb_guess)) then
        call SolveCoulomb(R,F,Coulomb_guess)
      else
        ! Start solving from a zero'd initial Coulomb potentials
        call SolveCoulomb(R,F)
      endif
    endif
    if(.not. present(Fread)) then    
        !-----------------------------------------------------------------------
        ! Add the Coulomb contribution to the potential F_I_I. The index 
        ! juggling is ugly, but necessary, because the Coulomb potential is
        ! defined on a larger mesh.
        !-----------------------------------------------------------------------
        if((all(protonsize.eq.0.0) .and. all(neutronsize.eq.0.0)) .or.         &
          &                             (.not. nucleonsize_selfconsistent)) then
          ox = coul_offset_x ; oy = coul_offset_y ; oz = coul_offset_z
          do k=1,nz
            do j=1,ny
              do i=1,nx
                F%F_I_I(meshindex(i,j,k),2)=F%F_I_I(meshindex(i,j,k),2)        &
                &                       + F%CoulombPotential(i+ox,j+oy,k+oz)   &
                &                       + F%ExchangePotential(i,j,k)
              enddo
            enddo
          enddo

        else
        !-----------------------------------------------------------------------
        ! Use the folded coulombpotential, for full self-consistency.
        ! Note that both protons and neutrons feel a Coulomb force if their
        ! charge form factor is taken into account.
        !-----------------------------------------------------------------------
          if(.not. allocated(F%foldedcoul)) then
            call stp('Nucleonsize_selfconsistent cannot be .false. if the protons are not point particles.')
          endif 
          do it=1, 2
            do k=1,nz
              do j=1,ny
                do i=1,nx
                  F%F_I_I(meshindex(i,j,k),it)=F%F_I_I(meshindex(i,j,k),it)    &
                  &                              + F%FoldedCoul(i,j,k,it)      &
                  &                              + F%FoldedExchange(i,j,k,it)
                enddo
              enddo
            enddo
          enddo
        endif 
        !-----------------------------------------------------------------------
        ! Add the contribution from the constraints on the electric multipole 
        ! moments. 
        Constraint_I_I = constraints_sph_elmult(.false.)
        F%F_I_I(:,1:2) = F%F_I_I(:,1:2) + Constraint_I_I(:,1:2)
        !-----------------------------------------------------------------------
        ! We added stuff to the proton and neutron potentials, we should be 
        ! consistent with the isospin 0 and 1 potentials
        F%F_I_I(:,3) = F%F_I_I(:,1) + F%F_I_I(:,2)
        F%F_I_I(:,4) = F%F_I_I(:,1) - F%F_I_I(:,2)
        !-----------------------------------------------------------------------
        ! Add the contribution of a cranking constraint to the 
        !    F_I_S and G_I_N potentials
$NTR    F_I_S = F_I_S + crank_spin_potential()
$NTR    G_I_N = G_I_N + crank_current_potential()
    endif
    call stop_timer(T_potentials)

  end function calcPotentials
  
  function precondition_potentials(F_in, F_out) result(F)
    !---------------------------------------------------------------------------
    ! Precondition the change in potentials from one iteration to the next. 
    ! All code generated by Hephaestos. 
    !
    ! Input:
    !   F_in : input potentials for an SCF iteration
    !   F_out: output potentials for an SCF iteration
    !
    ! Output:
    !   F    : set of preconditioned potentials. For every 
    !            F_a = F_in,a + P_a [F_out,a - F_in,a]
    !          where P_a is the preconditioning operator for µ
    !          mean-field potential a.
    !------------------------------------------ ---------------------------------
    type(PotentialVector), intent(in) :: F_in, F_out
    type(PotentialVector)             :: F
    real(KIND=dp), allocatable :: update(:,:)

    call start_timer(T_potentials)
    call start_timer(T_pot_precon)

    ! Everything which is not specifically preconditioned below just gets 
    ! explicitly copied from the new values.
    F = F_out

$POTENTIALPRECON
    call stop_timer(T_pot_precon)
    call stop_timer(T_potentials)
 
  end function precondition_potentials
  
  pure function pow( f, alpha) result(pf)
    !---------------------------------------------------------------------------
    ! Safely take powers of a density f, avoiding negative powers of numbers
    ! that might be accidentally 0 or negative below machine precision. This is
    ! achieved by adding a small (positive value) to the density. 
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in) :: f(mv), alpha
    real(KIND=dp)             :: pf(mv)
    
    if(alpha .lt. 0) then
      pf = (f + eps)**(alpha)
    else
      pf = (f)**(alpha)
    endif
  end function pow

  function INM_k2_pot(rho) result(pot)
    !-----------------------------------------------------------------
    ! Imagine homogeneous and unpolarised infinite nuclear matter:
    ! if one restricts itselfs to fourth order in gradients, the
    ! single-particle energies are:
    !
    !       e(k) = U_0 + U_2 k^2 + U_4 k^4
    !
    ! where the U_i are k-indepedent but possibly rho dependent.
    ! This routine calculates the potential U_2 in this expression
    ! as a function of the density on the mesh, starting only from
    ! the density D_I_I and using the local density approximation
    ! to get higher order densities.
    !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !   rho : density at every point on the mesh
    ! Output:
    !   pot : U_2 at every point on the mesh, assuming
    !         homogeneous unpolarised INM at that specific density.
    !-----------------------------------------------------------------
    real(KIND=dp), intent(in) :: rho(mv,4)
    real(KIND=dp)             :: pot(mv,4)
    real(KIND=dp)             :: kfn(mv), kfp(mv)

    kfn=(3.d0*pi**2*rho(:,1))**(1.0d0/3.0d0) ! Neutron density
    kfp=(3.d0*pi**2*rho(:,2))**(1.0d0/3.0d0) ! Proton  density

    !- - - - - - - - - - - - - - - - - - - - - -
    ! Initialize to zero
    pot = 0.0d0

    !- - - - - - - - - - - - - - - - - - - - - -
    ! Calculation of isospin 0 and 1
$K2POT

    !- - - - - - - - - - - - - - - - - - - - - -
    ! Recombine to proton and neutron potentials
    pot(:,1) = pot(:,3) + pot(:,4)
    pot(:,2) = pot(:,3) - pot(:,4)
  end function INM_k2_pot

  function INM_k4_pot(rho) result(pot)
    !-----------------------------------------------------------------
    ! Imagine homogeneous and unpolarised infinite nuclear matter:
    ! if one restricts itselfs to fourth order in gradients, the
    ! single-particle energies are:
    !
    !       e(k) = U_0 + U_2 k^2 + U_4 k^4
    !
    ! where the U_i are k-indepedent but possibly rho dependent.
    ! This routine calculates the potential U_4 in this expression
    ! as a function of the density on the mesh, starting only from
    ! the density D_I_I and using the local density approximation
    ! to get higher order densities.
    !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !   rho : density at every point on the mesh
    ! Output:
    !   pot : U_4 at every point on the mesh, assuming
    !         homogeneous unpolarised INM at that specific density.
    !-----------------------------------------------------------------
    real(KIND=dp), intent(in) :: rho(mv,4)
    real(KIND=dp)             :: pot(mv,4)
    real(KIND=dp)             :: kfn(mv), kfp(mv)

    kfn=(3.d0*pi**2*rho(:,1))**(1.0d0/3.0d0) ! Neutron density
    kfp=(3.d0*pi**2*rho(:,2))**(1.0d0/3.0d0) ! Proton  density

    !- - - - - - - - - - - - - - - - - - - - - -
    ! Initialize to zero
    pot = 0.0d0

    !- - - - - - - - - - - - - - - - - - - - - -
    ! Calculation of isospin 0 and 1
$K4POT

    !- - - - - - - - - - - - - - - - - - - - - -
    ! Recombine to proton and neutron potentials
    pot(:,1) = pot(:,3) + pot(:,4)
    pot(:,2) = pot(:,3) - pot(:,4)
end function INM_k4_pot

  function apply_sphamil(psi, dpsi, ddpsi, &
$N3                                 dddpsi, &
&                                        sx,sy,sz,iso, onthefly, F) result(hpsi)
    !---------------------------------------------------------------------------
    ! Apply the single-particle hamiltonian to a single-particle wavefunction.
    ! - - - - - - - - - - - -- - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Input:
    !          psi : spwf to act on with h
    !  d/dd/dddpsi : arrays containing the first, second and third derivatives 
    !                of the spwf. ddpsi does not need to be a "full" matrix 
    !                when dealing with standard NLO EDFs. dddpsi is only used 
    !                when dealing with N3LO EDFs.
    ! sx/sy/sz     : signs under reflection symmetry for this particular spwf
    !                not referenced when onthefly = .false.
    ! onthefly     : if .true., recalculate the derivatives of psi and store
    !                them in the array psi.
    ! F            : a set of mean-field potentials determining the 
    !                singleparticle hamiltonian
    ! - - - - - - - - - - - -- - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Output:
    !    hpsi      : h | psi >
    !  d/dd/dddpsi : arrays containing the derivatives of psi
    !                if onthefly = .true.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Important notes:
    ! - the input values of the derivative matrices dpsi/ddpsi/dddpsi are not
    !   relevant if onthefly=.true. These arrays will be overwritten on output
    !   in case the derivatives can serve afterwards.
    ! - the array of third derivatives is only necessary for N3LO calculation, a
    !   and Hephaestos comments them out of the interface of this function when
    !   possible. This is why the function declaration above is spread across
    !   a few lines.
    !---------------------------------------------------------------------------

    use derivatives

    ! Logical indicating if the derivatives need to be calculated before
    ! applying h.
    ! If false, the derivatives are passed in. If True, the derivatives are not
    ! passed in and need to be calculated.
    logical, intent(in)               :: onthefly 
    real(KIND=dp), intent(in)         :: psi(mv,4)
    real(KIND=dp), intent(inout)      :: dpsi(mv,3,4),ddpsi(mv,6,4)
$N3 real(KIND=dp), intent(inout)      :: dddpsi(mv,10,4)
    integer, intent(in)               :: sx(4),sy(4),sz(4), iso
    type(PotentialVector), intent(in) :: F

    integer                   :: sym(4)
    real(KIND=dp)             :: hpsi(mv,4)
    real(KIND=dp)             :: temp(mv,4)
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Declaration of temporary spinors 
    !
    ! Technical note: these arrays are declared with the spinor indices (1-4)
    !                 BEFORE the derivative indices (3). This is to aid the 
    !                 memory locality of operations in this particular function
    !                 and is OPPOSITE the conventions of the rest of the code.
    real(KIND=dp)             ::   dtemp(mv,4,3)
$D2TEMPSPH    real(KIND=dp)   ::  ddtemp(mv,4,3,3)
$D3TEMPSPH    real(KIND=dp)   :: dddtemp(mv,4,3,3,3)
$LAPTEMPSPH   real(KIND=dp)   :: laptemp(mv,4)
    
    real(KIND=dp)             :: ReducedMass, Butler_t, Butler_f
    
    integer :: it, i,k
    
    call start_timer(T_sphamil)
    !---------------------------------------------------------------------------
    ! Determine the isospin index
    it = (iso + 3)/2
    !---------------------------------------------------------------------------
    ! Reduced mass in case of self-consistent 1-body COM correction
    ! If doing pasta calculations, just skip.
    Reducedmass = 1.0_dp
#if(PASTA == 0)
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
#endif
    !---------------------------------------------------------------------------

    if(OnTheFly) then
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      ! Calculate the derivatives of this spwf on the fly
#if(USE_Periodic==0)
      ! Original Lagrange mesh boundary conditions
      do k=1,4
$N2        call Derive_tot(psi(:,k),sx(k),sy(k),sz(k),dpsi(:,:,k),ddpsi(:,:,k))
$N3        call Derive_tot(psi(:,k),sx(k),sy(k),sz(k),dpsi(:,:,k),ddpsi(:,:,k),&
$N3        &                                     dddpsi(:,:,k))
      enddo
#else  
      !NS: for periodic boundary conditions
      do k=1,2
$N2        call Derive_tot_periodic(psi(:,   (2*k-1):2*k), &
$N2             &                         sx((2*k-1):2*k), &
$N2             &                         sy((2*k-1):2*k), & 
$N2             &                         sz((2*k-1):2*k), & 
$N2             &                   dpsi(:,:,(2*k-1):2*k), &
$N2             &                  ddpsi(:,:,(2*k-1):2*k))
      enddo
#endif
    endif
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
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
    ! Action of the Skyrme potentials
    !
    ! Note that 
    ! a) Coulomb is included in the F_I_I potential
    ! b) Every density contains the contributions from constraints on that 
    !    density. 
    ! c) The kinetic energy is NOT included in the F_N_N potential, because 
    !    a constant is not in the Lagrange basis; so the current way of
    !    deriving stuff is not correct for a term of the form
    !          hbar^2_2m
    !---------------------------------------------------------------------------
$SKYRMEACTION

    call stop_timer(T_sphamil)

  end function apply_sphamil
  
  function delta_action(        psi,   &
$N1DELTA                   &   dpsi,   &
$N2DELTA                   &  ddpsi,   &
$N3DELTA                   & dddpsi,   &
$SYMDELTA                  & sx,sy,sz, &
&                                         iso, onthefly, F) result(deltapsi)
    !---------------------------------------------------------------------------
    !
    ! onthefly:
    !   Logical indicating if the derivatives need to be calculated before
    !   applying delta. If false, the derivatives are passed in. If True, the 
    !   derivatives are not passed in and need to be calculated.
    !---------------------------------------------------------------------------
    logical, intent(in)               :: onthefly 
    integer, intent(in)               :: iso
    type(PotentialVector), intent(in) :: F
    real(KIND=dp), intent(in)         :: psi(:,:)  
$N1DELTA    real(KIND=dp), intent(inout) :: dpsi(:,:,:)
$N2DELTA    real(KIND=dp), intent(inout) :: ddpsi(:,:,:)
$N3DELTA    real(KIND=dp), intent(inout) :: dddpsi(:,:,:)
$SYMDELTA    integer, intent(in)          :: sx(:),sy(:),sz(:)
    real(KIND=dp), allocatable   :: deltapsi(:,:)
    real(KIND=dp)                ::    temp(mv,4)
$D1TEMPDELTA    real(KIND=dp)    ::   dtemp(mv,3,4)
$D2TEMPDELTA    real(KIND=dp)    ::  ddtemp(mv,3,3,4)
$D3TEMPDELTA    real(KIND=dp)    :: dddtemp(mv,3,3,3,4)
$LAPTEMPDELTA   real(KIND=dp)    :: laptemp(mv,4)
    integer                      :: it,i
    
    !---------------------------------------------------------------------------
    ! Determine the isospin index
    it = (iso + 3)/2
    
    if(onthefly) then
      call stp('On the fly calculation of derivatives in delta_action not implemented.')
    endif
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
    spwfenergy = 0 
    do wave=1,nwt
        if(pairingtype.lt.2) then
          spwfenergy = spwfenergy + rho_can(wave) * spenergies(wave)
        else
          spwfenergy = spwfenergy + rho_can(wave) * canenergies(wave)
        endif
    enddo
    !
    spwfenergy = 0.5 * spwfenergy
    
    ! Calculation of rearrangement energy (without Coulomb Exchange)
    e_rear = 0                    
    e_rear = e_rear - 0.5d0*  trilinear
    e_rear = e_rear -         quadrilinear
$EREAR   
   
   
    spwfenergy = spwfenergy + e_rear
    ! Add kinetic and CoulombExchange contributions
    spwfenergy = spwfenergy + 0.5 * sum(kinetic) + CoulombExchange/3.d0
    
    ! Always add the 1-body COMcorrection. In case it is used iteratively, it
    ! is double counted along with the kinetic energy!
    if(COM1body.gt.0) then
        SpwfEnergy = SpwfEnergy  + sum(COMCorrection(1,:))/2.0_dp
    endif

    if(COM2body.gt.0) then
        SpwfEnergy = SpwfEnergy  + sum(COMCorrection(2,:))  
    endif
    
    ! Subtract contribution by multipole constraints
    Constraint_I_I = constraints_sph_elmult(.false.)
    SpwfEnergy = SpwfEnergy - &
    &                sum(Constraint_I_I(:,1:2) * Density%D_I_I(:,1:2))*dv/2.0_dp

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
    Spwfenergy = Spwfenergy + sum(vibcorrection)
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

  end subroutine clean_potentials
  
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

  subroutine WritePotentials(chan, F)
    !---------------------------------------------------------------------------
    !  Subroutine writing the different potentials to file.
    !---------------------------------------------------------------------------
    type(PotentialVector), intent(in) :: F
    integer, intent(in) :: chan
    integer             :: io

    ! Signalling how many potentials have been stored.
    write(chan, iostat=io) $POTENTIALNUMBER

    ! Then, for every potential write the 
    ! * Name 
    ! * Value
    ! Note that the name is written as a length-30 string, padded with spaces.
    ! If not, the unformatted in/out cannot correctly determine the end of a
    ! string and comparisons can not be made.
$WRITEPOTENTIALS
  end subroutine WritePotentials
 
 
 subroutine WritePotentials_hdf5(file_id, F)
   !---------------------------------------------------------------------------
   !  Subroutine writing the different potentials to hdf5 file.
   !---------------------------------------------------------------------------
   type(PotentialVector), intent(in) :: F
   integer(hid_t), intent(in) :: file_id

   ! Then, for every potential write the 
   ! * Name 
   ! * Value
   ! Note that the name is written as a length-30 string, padded with spaces.
   ! If not, the unformatted in/out cannot correctly determine the end of a
   ! string and comparisons can not be made.
$WRITEPOTENTIALS_HDF5
 end subroutine WritePotentials_hdf5

  subroutine hdf5_writepot(id, name, dset,n)
    ! writes double precision dataset array with some name in the hdf5 file
    use HDF5
    character(len=30), intent(in) :: name
    integer(hid_t), intent(in) :: id
    integer       , intent(in) :: n
    real(kind=dp),  intent(in) :: dset(n) 
    integer(hid_t)             :: space_id, dset_id, plist_id  
    integer                    :: error
    integer(hsize_t), dimension(1) :: dims,data_dims!, chdims

    dims=(/n/)
    data_dims(1)=n
    ! Create dataspace for data_set 
    call h5screate_simple_f(1, dims, space_id, error)
    ! create property list
    call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
    ! create chunks with property list for compression, as of now size of chunk  
    ! is just equal to the size of array (for some reason work better). 
    ! Modify for MPI reading?
    call h5pset_chunk_f(plist_id, 1, dims, error)
    ! shuffling for better compression?
    call h5pset_shuffle_f(plist_id, error)
    ! zlib compression with deflate
    call h5pset_deflate_f(plist_id, comprlvl, error)
    ! Create dataset with default properties "dset_id" is returned
    call h5dcreate_f(id,'potentials/'//trim(name),H5T_NATIVE_DOUBLE, space_id, &
                      dset_id, error, plist_id)
    ! Write dataset 
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, dset, data_dims, error)
    ! Close access to dataset 
    call h5dclose_f(dset_id, error)
    ! Close access to data space 
    call h5sclose_f(space_id, error)
    ! close access to plist
    call h5pclose_f(plist_id, error)

    if (error.ne.0) then
      call stp('ERROR: writting potentials in hdf5 format')
    endif

  end subroutine hdf5_writepot

  function ReadPotentials(chan, filenx, fileny, filenz, symtransfo_needed) &
  & result(F)
    !---------------------------------------------------------------------------
    ! Subroutine that reads the different mean-field potentials from a 
    ! wavefunction file created by a previous run of the code.
    ! Notes:
    !  1. this function does not rely on MPI I/O and simply reads everything 
    !     with rank 0 and then does a bunch of MPI_BCASTS.
    !  2. this function should not be confused with the readpotentials_separate
    !     routine below
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !   chan                  : integer, channel number for input
    !   filenx, fileny,filenz : integers, number of mesh points in every 
    !                           direction for the quantities on file
    !   symtransfo_needed     : logical, if a symmetry transformation is 
    !                           needed (.true.) or not (.false.)
    !                  .false.: use the potentials as read from 
    !                           file, transforming only the number of mesh 
    !                           points if needed. 
    !                  .true. : use the potentials from file for further 
    !                           calculations. This means just reading them here
    !                           and trusting the rest of the program to do the
    !                           the rest.
    ! Output:
    !   F                      : a potential-vector, read from file
    !---------------------------------------------------------------------------
    integer, intent(in)   :: chan, filenx, fileny, filenz
    logical, intent(in)   :: symtransfo_needed
    type(PotentialVector) :: F, F_temp
    integer               :: io, potnumber, potcount, it, filemv
    character(len=30)     :: potname

#if(USE_MPI > 0)
    integer             :: mpi_err
#endif

    filemv = filenx * fileny * filenz

    ! Checking how many potentials have been stored on file
    if(MPI_RANK .eq. 0) read(chan, iostat=io) potnumber
#if(USE_MPI > 0)
    call MPI_BCAST(potnumber, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)
#endif

    do potcount = 1,potnumber
        ! Read the name of the potential stored
        if(MPI_RANK .eq. 0) read(chan, iostat=io) potname
#if(USE_MPI > 0)
        call MPI_BCAST(potname,30, MPI_CHARACTER, 0, MPI_COMM_WORLD, mpi_err)
#endif
        ! Then select which potential we are going to be reading
        select case(trim(potname))
$READPOTENTIALS
        CASE DEFAULT
          ! The potential is not in this program, forget about it
          if(MPI_RANK .eq. 0 ) read(chan, iostat=io)
        end select
    enddo

  end function ReadPotentials

  function ReadPotentials_hdf5(file_id, filenx, fileny, filenz, symtransfo_needed) &
   & result(F)
    !---------------------------------------------------------------------------
    ! Subroutine that reads the different mean-field potentials from a 
    ! wavefunction file created by a previous run of the code.
    ! Notes:
    !  1. now only sequental reading
    !  2. this function should not be confused with the readpotentials_separate
    !     routine below
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !   chan                  : id of the group potentials
    !   filenx, fileny,filenz : integers, number of mesh points in every 
    !                           direction for the quantities on file
    !   symtransfo_needed     : logical, if a symmetry transformation is 
    !                           needed (.true.) or not (.false.)
    !                  .false.: use the potentials as read from 
    !                           file, transforming only the number of mesh 
    !                           points if needed. 
    !                  .true. : use the potentials from file for further 
    !                           calculations. This means just reading them here
    !                           and trusting the rest of the program to do the
    !                           the rest.
    ! Output:
    !   F                      : a potential-vector, read from file
    !---------------------------------------------------------------------------
    integer(hid_t), intent(in) :: file_id
    integer(hid_t)             :: group_id
    character(len=11)          :: groupname
    integer, intent(in)        :: filenx, fileny, filenz
    logical, intent(in)        :: symtransfo_needed
    type(PotentialVector)      :: F, F_temp
    real(kind=dp), allocatable :: Ftmp(:) 
    integer                    :: filemv, it, h5ferr
#if(USE_MPI > 0)
    integer                    :: mpi_err
#endif

    filemv = filenx * fileny * filenz
    
    if(MPI_RANK .eq. 0) then
      groupname='/potentials'
      call h5gopen_f(file_id,groupname,group_id,h5ferr)
    endif

$READPOTENTIALS_HDF5

    if(MPI_RANK.eq.0) then
      !close the potentials group
      call h5gclose_f(group_id, h5ferr)
    endif

  end function ReadPotentials_hdf5

  subroutine hdf5_readpot(id, name, dset, n)
    ! reads double precision potential of length n with some name in the hdf5 file
    use HDF5
    character(len=*), intent(in) :: name
    integer(hid_t), intent(in)   :: id
    integer, intent(in)          :: n
    real(kind=dp), intent(inout) :: dset(n)
    integer(hid_t)               :: dset_id !identifiers
    integer                      :: error
    integer(hsize_t), dimension(1) :: dims, data_dims

    dims=(/n/)
    data_dims(1)=n
    ! open dataset, "dset_id" is returned
    call h5dopen_f(id, name, dset_id, error)
    ! read dataset 
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, dset, data_dims, error)
    ! Close access to dataset 
    call h5dclose_f(dset_id, error)

    if (error.ne.0) then
      call stp('ERROR: reading potentials in hdf5 format')
    endif

  end subroutine hdf5_readpot

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

  function readpotentials_separate(chan, ifn) result(F)
    !---------------------------------------------------------------------------
    ! Read mean-field potentials from a separate user-provided file that is 
    ! NOT a wavefunction file. This routine should not be confused with the 
    ! readpotentials subroutine from the functional module.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Input: 
    !  * chan : integer, channel number to read the file
    !  * ifn  : input filename (will be checked for existence)
    !
    ! Output:
    !      F  : a set of mean-field potentials
    !
    ! Caution: this routine is currently foreseen for a specific application, 
    !          limited to maximally symmetric calculations and .func files
    !          for which F_Nm_Nm and G_I_NS potentials are defined.
    !
    ! 
    !---------------------------------------------------------------------------
    use Coulombmod ! module explicitly 'used' in order to be able to place the 
                   ! values of the direct and exchange Coulomb potentials 
                   ! correctly on the mesh
  
    type(PotentialVector)        :: F
    integer, intent(in)          :: chan
    character(len=*), intent(in) :: ifn

    logical :: exists
    integer :: i,j,k,io, it, mu, nu, ox, oy, oz, headercount, mpi_err
    real(KIND=dp), allocatable :: Vc(:), Ec(:)
    real(KIND=dp)              :: x,y,z
    character(len=200)         :: temp
     
    inquire(file=ifn, exist=exists)
    if(.not.exists) then
      print *, 'Input file specified does not exist!'
      stop
    endif
    
    open (chan,file=ifn)
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! We need to skip any header lines (indicated by #).
    ! For a Tantalus-created file, there are 14 of them by default but other 
    ! people might write a different amount
    io = 0; headercount = -1
    do while(io.eq.0) 
      headercount = headercount + 1
      read(chan, iostat=io, fmt='(a200)') temp
      if(temp(1:1) .ne. '#') io = 1
    enddo  
    ! We've found an error; we have counted the number of header lines!
    rewind(chan)
    ! ... and now we skip this number of lines    
    do i=1,headercount
        read(chan, fmt=('()'))
    enddo
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Allocate the relevant potentials    
    allocate(F%F_I_I   (nx*ny*nz,4))     ; F%F_I_I    = 0.0d0
    allocate(F%FP_I_I  (nx*ny*nz,2))     ; F%FP_I_I   = 0.0d0
$TAUSCALAR    allocate(F%F_Nm_Nm(nx*ny*nz,4))   ; F%F_Nm_Nm = 0.0d0
$TAUTENSOR    allocate(F%F_N_N(nx*ny*nz,3,3,4)) ; F%F_N_N   = 0.0d0
    allocate(F%G_I_NS (nx*ny*nz,3,3,4)) ; F%G_I_NS  = 0.0d0
    allocate(Vc(nx*ny*nz))              ; VC        = 0.0d0
    allocate(Ec(nx*ny*nz))              ; EC        = 0.0d0

    ! We assume the points on the file are correctly ordered in 
    ! FORTRAN fashion, such that we do not have to worry about looping 
    ! separately over x/y/z and can just loop once over all mesh points.
    ! This also means the coordinate information is not used.
    do i=1,nx*ny*nz
      read(chan, fmt='(7es25.12)', iostat=io, advance='no')        & 
      &                            x,y,z,                          & !unused
      &                            F%F_I_I(i,1),F%F_I_I(i,2),      & ! U(r)
      &                            Vc(i),  Ec(i)                     ! Coulomb
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! It is easy to read the kinetic potential from ETFSI calculations if 
      ! D_Nm_Nm is defined as a contracted density
$TAUSCALAR      read(chan, fmt='(2es25.12)', iostat=io, advance='no')    & 
$TAUSCALAR      &                            F%F_Nm_Nm(i,1), F%F_Nm_Nm(i,2)    ! kinetic
!      ! If not, then we have to do some reorganisation
$TAUTENSOR      read(chan, fmt='(2es25.12)', iostat=io, advance='no')    & 
$TAUTENSOR      &                            F%F_N_N(i,1,1,1), F%F_N_N(i,1,1,2)! kinetic
$TAUTENSOR      F%F_N_N(i,2,2,:) = F%F_N_N(i,1,1,:)/3
$TAUTENSOR      F%F_N_N(i,3,3,:) = F%F_N_N(i,1,1,:)/3
$TAUTENSOR      F%F_N_N(i,1,1,:) = F%F_N_N(i,1,1,:)/3
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! Pairing field
      read(chan, fmt='(7es25.12)', iostat=io, advance='no')                    &
      &                            F%FP_I_I(i,1),F%FP_I_I(i,2)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! All components of the spin-orbit field
      do mu=1,3
        do nu=1,3
          read(chan, fmt='(2es25.12)', advance='no', iostat=io) &
          &               F%G_I_NS(i,mu,nu,1), F%G_I_NS(i,mu,nu,2)
        enddo
      enddo
      
      read(chan, *) ! Advance to new line
      
      if(io.ne.0) then
        print *, 'Problem encountered reading potential file ', ifn
        print *, 'IOSTAT = ', io
        stop
      endif
    enddo

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Make sure the Coulomb module is configured with the right array dimensions
    call setupCoulomb(F)
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! The index juggling is ugly, but necessary, because the Coulomb 
    ! potential has a different size than the Lagrange mesh.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ox = coul_offset_x ; oy = coul_offset_y ; oz = coul_offset_z
      
    do k=1,nz
      do j=1,ny
        do i=1,nx
          F%CoulombPotential(i+ox,j+oy,k+oz)  = Vc(meshindex(i,j,k))
          F%ExchangePotential(i+ox,j+oy,k+oz) = Ec(meshindex(i,j,k))
        enddo
      enddo
    enddo
    
    if((all(protonsize.eq.0.0) .and. all(neutronsize.eq.0.0)) .or.         &
    &                             (.not. nucleonsize_selfconsistent)) then
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      ! No finite size effects; correction is simple
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -      
      F%F_I_I(:,2) = F%F_I_I(:,2) + Vc(:) + Ec(:)
    else
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! Finite size effects taken into account
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! Calculate folded potentials from the read-in potentials
      call obtain_folded_potentials(F)
      ! ... and correct F_I_I for them with ugly index juggling
      do it=1, 2
        do k=1,nz
          do j=1,ny
            do i=1,nx

              F%F_I_I(meshindex(i,j,k),it)= F%F_I_I(meshindex(i,j,k),it)       &
              &                           + F%FoldedCoul(i,j,k,it)             &
              &                           + F%FoldedExchange(i,j,k,it)
            enddo
          enddo
        enddo
      enddo
    endif

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Make sure isospin combinations are made correctly for all potentials
    F%F_I_I(:,3) = F%F_I_I(:,1) + F%F_I_I(:,2)
    F%F_I_I(:,4) = F%F_I_I(:,1) - F%F_I_I(:,2)

$TAUSCALAR    F%F_Nm_Nm(:,3)   = F%F_Nm_Nm(:,1)   + F%F_Nm_Nm(:,2)
$TAUSCALAR    F%F_Nm_Nm(:,4)   = F%F_Nm_Nm(:,1)   - F%F_Nm_Nm(:,2)
$TAUTENSOR    F%F_N_N(:,:,:,3) = F%F_N_N(:,:,:,1) + F%F_N_N(:,:,:,2)
$TAUTENSOR    F%F_N_N(:,:,:,4) = F%F_N_N(:,:,:,1) - F%F_N_N(:,:,:,2)

    do mu=1,3
      do nu=1,3
        F%G_I_NS(:,mu,nu,3) = F%G_I_NS(:,mu,nu,1) + F%G_I_NS(:,mu,nu,2)
        F%G_I_NS(:,mu,nu,4) = F%G_I_NS(:,mu,nu,1) - F%G_I_NS(:,mu,nu,2)
      enddo
    enddo

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Close channel after succesfull IO operations.
    close(chan)
  end function readpotentials_separate

end module functional
