module functional
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
 ! Calcpotentials_PERTURBED   : [WAY TOO LONG TO INCLUDE HERE]
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
 ! WRITEPOTENTIALS_HDF5 : [WAY TOO LONG TO INCLUDE HERE]
 ! READPOTENTIALS_HDF5 : [WAY TOO LONG TO INCLUDE HERE]
 ! 
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
 ! DTEMPSPH         : $DTEMPSPH
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
 ! TAUPRESENT       : $TAUPRESENT
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
 use compilation, only: dp
 use geninfo, only: dv, dx, inversetemp, meshgrid, meshindex,          &
                    mpi_rank, mv, neutrons, nx, ny, nz, pi, protons,   &
                    reduZ, rps, simulate_spherical_bc,                 &
                    store_derivatives, stp
 use densities, only: CompNablaMelements, denddpsi, denpsi,              &
                     Density, print_maxval,                            &
                     sx_rho, sx_rho_antisym, sy_rho, sy_rho_antisym,    &
                     sz_rho, sz_rho_antisym
 ! The raw EDF couplings (t0, x3, sigma_b, t1n2, ...) are functional-
 ! dependent: which names exist and are referenced by the generated code
 ! differs per .func file.  Hephaestos therefore fills in the only:-list
 ! of the coupling names here via a dedicated marker, and, in the
 ! generated parameterization, declares the full union of couplings over
 ! all .func files so that these names exist in every build.
 use parameterization, only: Cc, COM1body, COM2body,              &
                            Estabn, Estabp, hbarclum, hbm,          &
                            neutronsize, nucleonmass,               &
                            nucleonsize_selfconsistent,              &
                            printparameterization, protonsize,     &
                            readparameterization, RotCorr, RotcorrB, &
                            rotcorrC, Vibcorr, vibcorrB, vibcorrD,   &
                            vibcorrL, pow, lr_pow,                   &
                            $COUPLINGS
 use pairing, only: Bogoliubov, CalcEntropy, calcpairingenergy,           &
                   CanTransfo, entropy, FermiEnergy, kappa_can,         &
                   PairingType, rho_can, rho_pairing
 use BCS, only: delta_action_BCS
 use HFB, only: conjugp, delta_action_HFB
 use timing, only: start_timer, stop_timer, T_COM, T_COM1, T_COM2,      &
                   T_COM2_summation, T_energy, T_MOI, T_pot_precon,     &
                   T_potentials, T_sphamil
 use transform, only: changeboxsize_function, sym_transfo_needed
 use vectors, only: DensityVector, PotentialVector, memory
 use preconditioning, only: PreconditionPotential,                   &
                          KerkerPreconditionPotential
 use wavefunctions, only: hfblocks, hfblocks_global, hfddpsi, hfpsi,    &
                         nwn, nwt, nwt_local, spwf_map
 use Cranking, only: crank_current_potential, crank_spin_potential,     &
                     crankenergy_cut
 use pairing_strengths, only: vmicro, vmicro_stored
#if(USE_HDF5>0)
 use HDF5, only: hid_t
#endif

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
    real(KIND=dp), parameter :: Qnp=1.29335236d0 !Mn-Mp
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
    real(KIND=dp) :: VibCorrection
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
    !---------------------------------------------------------------------------
    ! Numerical parameter of the preconditioning of the Skyrme potentials
#if(PASTA > 0)
    ! The 'optimal' default value for potential preconditioning in pasta 
    ! calculations is higher than for finite nuclei; this is not currently
    ! understood. 
    real(KIND=dp) :: preconfactor = 4.0_dp
#else
    real(KIND=dp) :: preconfactor = 1.0_dp
#endif
    ! Kerker parameter for the preconditioning of the Coulomb potential
    ! Nuclei => don't do kerker by default
    ! Pasta  => do Kerker by default
    real(KIND=dp) :: kerker_k0    = 0.0d0 ! typical screening length ~ 100 fm
    real(KIND=dp) :: mixstepsize  = 1.0d0
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

    interface INM_k2_pot
      module procedure INM_k2_pot_real
      module procedure INM_k2_pot_complex
    end interface

    interface INM_k4_pot
      module procedure INM_k4_pot_real
      module procedure INM_k4_pot_complex
    end interface 
  interface
    !-------------------------------------------------------------------
    ! The following procedures are implemented in submodules (files
    ! src_orig/functional_{energy,potentials,fields,io}.f90).
    !-------------------------------------------------------------------
    ! ---- functional_energy
    module subroutine CalcEnergy(Rin, Fin, calc_expensive)
    logical, intent(in) :: calc_expensive
    type(DensityVector), intent(in)   :: Rin
    type(PotentialVector), intent(in) :: Fin
  end subroutine CalcEnergy
    module subroutine CompSkyrme(R, F, S, T_even, T_odd, PE)
    type(DensityVector), intent(in)   :: R
    type(PotentialVector), intent(in) :: F
    real(KIND=dp), intent(out)        :: S, T_even, T_odd, PE(2)
  end subroutine CompSkyrme
    module subroutine PrintSkyrme()
  end subroutine PrintSkyrme
    module function CompKinetic_spwfs() result(kinetic)
    real(KIND=dp)    :: kinetic(2)
  end function CompKinetic_spwfs
    module function CompKinetic_density(Rin) result(kinetic)
    real(KIND=dp)                   :: kinetic(2)
    type(DensityVector), intent(in) :: Rin
  end function CompKinetic_density
    module function calculate_epasta(nucE) result(Epasta)
    real(KIND=dp), intent(in) :: nucE
    real(KIND=dp)             :: Epasta
  end function calculate_epasta
    module function calcspwfenergy() result(spwfenergy)
    real(KIND=dp) :: spwfenergy
  end function calcspwfenergy
    module subroutine CompCOMCorrection(kin, override_2body, Comcorr)
    real(KIND=dp), intent(out)    :: Comcorr(2,2)
    logical, intent(in)           :: override_2body
    real(KIND=dp), intent(in)     :: kin(2)
  end subroutine CompCOMCorrection
    module function calc_rotational_correction(rot_option) result(Erot)
    integer, intent(in)        :: rot_option
    real(KIND=dp)              :: Erot(3)
  end function calc_rotational_correction
    module function calc_vibrational_correction(vib_option) result(Evib)
    integer, intent(in)        :: vib_option
    real(KIND=dp)              :: Evib
  end function calc_vibrational_correction
    module pure function inv_mat(A) result(Ainv)
    real(KIND=dp), intent(in) :: A(2,2)
    real(KIND=dp)             :: Ainv(2,2)
  end function inv_mat
    module function classical_MOI(N,Z) result(MOI)
    real(KIND=dp), intent(in) :: N,Z
    real(KIND=dp)             :: MOI
  end function classical_MOI
    module function CompStabilisingFactor(PairE) result(stab)
    real(KIND=dp)             :: stab(2)
    real(KIND=dp), intent(in) :: PairE(2)
  end function CompStabilisingFactor
    module subroutine verify_COM_motion()
  end subroutine verify_COM_motion
    module subroutine calcElectronEnergy()
  end subroutine calcElectronEnergy
    module subroutine update_E_history()
  end subroutine update_E_history
    
    ! ---- functional_potentials
    module function calcPotentials(R, Fread) result(F)
    type(DensityVector), intent(in)             :: R
    type(PotentialVector), intent(in), optional :: Fread
    type(PotentialVector)                       :: F
  end function calcPotentials
    module subroutine calc_perturbed_potentials( R, dRs, dRa, dR_pp_plus,        &
    &                                     dR_pp_minus, dFs, dFa,         &
    &                                     dF_pp_plus, dF_pp_minus)
    type (DensityVector), intent(in) :: R, dRs, dRa, dR_pp_plus, dR_pp_minus
    type (PotentialVector)           :: dFs, dFa, dF_pp_plus, dF_pp_minus
  end subroutine calc_perturbed_potentials
    module function calc_perturbed_potentials_oneoff(R,R_pert) result (F)
    type (DensityVector), intent(in) :: R, R_pert
    type (PotentialVector)           :: F
  end function calc_perturbed_potentials_oneoff
    module subroutine combine_potentials(F)
    type(PotentialVector), intent(inout) :: F
  end subroutine combine_potentials
    module subroutine add_constraint_potential(F)
    type(PotentialVector), intent(inout) :: F
  end subroutine add_constraint_potential
    module subroutine add_coulomb_potential(F)
    type(PotentialVector), intent(inout) :: F
  end subroutine add_coulomb_potential
    module subroutine add_cranking_potentials(F)
    type(PotentialVector), intent(inout) :: F
  end subroutine add_cranking_potentials
    module function precondition_potentials(F_in, F_out) result(F)
    type(PotentialVector), intent(in) :: F_in, F_out
    type(PotentialVector)             :: F
  end function precondition_potentials
    module subroutine clean_potentials()
  end subroutine clean_potentials
    module subroutine save_potential_history(F_in, F_out)
    type(PotentialVector), intent(in) :: F_in, F_out
  end subroutine save_potential_history
    
    ! ---- functional_fields
    module function Add_potentialvector(F1, F2) result(F)
    type(PotentialVector), intent(in) :: F1, F2
    type(PotentialVector)             :: F
  end function Add_potentialvector
    module function multiply_potentialvector(a, F1) result(F)
    type(PotentialVector), intent(in) :: F1
    real(KIND=dp), intent(in)         :: a
    type(PotentialVector)             :: F
  end function multiply_potentialvector
    module function PVectorInproduct(F1, F2) result(x)
    type(PotentialVector), intent(in) :: F1,F2
    real(KIND=dp) :: x
  end function PVectorInproduct
    
    ! ---- functional_io
    module subroutine WritePotentials(chan, F)
    type(PotentialVector), intent(in) :: F
    integer, intent(in) :: chan
  end subroutine WritePotentials
#if(USE_HDF5>0 && $FAM == 0)
    module subroutine write_hdf5_potentials(file_id, F)
    type(PotentialVector), intent(in) :: F
    integer(hid_t), intent(in)        :: file_id
  end subroutine write_hdf5_potentials
#endif
    module function ReadPotentials(chan, filenx, fileny, filenz, sym_transfo_needed) result(F)
    integer, intent(in)   :: chan, filenx, fileny, filenz
    logical, intent(in)   :: sym_transfo_needed
    type(PotentialVector) :: F
  end function ReadPotentials
#if(USE_HDF5 > 0)
    module function read_potentials_hdf5(file_id, filenx, fileny, filenz, sym_transfo_needed) result(F)
    integer(hid_t), intent(in)   :: file_id
    integer, intent(in)          :: filenx, fileny, filenz
    logical, intent(in)          :: sym_transfo_needed
    type(PotentialVector)      :: F
  end function read_potentials_hdf5
#endif
    module subroutine hdf5_readpot(id, name, dset, n)
    character(len=*), intent(in) :: name
    integer(hid_t), intent(in)   :: id
    integer, intent(in)          :: n
    real(kind=dp), intent(inout) :: dset(n)
  end subroutine hdf5_readpot
    module function readpotentials_separate(chan, ifn) result(F)
    type(PotentialVector)        :: F
    integer, intent(in)          :: chan
    character(len=*), intent(in) :: ifn
  end function readpotentials_separate

    ! ---- functional_coefs
    module subroutine calcedfcoefs()
  end subroutine calcedfcoefs

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

    1 format (91('-'))
    5 format (40x, '    neutron             proton              total')
    6 format (15x, ' Kinetic Energy:', 3f20.6)
   61 format (15x, '     COM 1-body:', 3f20.6)
   62 format (15x, '     COM 2-body:', 3f20.6)
  621 format (15x, '             ph:', 3f20.6)
  622 format (15x, '             pp:', 3f20.6)

   63 format (15x, '  Rotational  '  , a1, ':', 40x, f20.6)
  631 format (15x, '  Rotational  T:',          40x, f20.6)
  
   65 format (15x, '  Collective T: ',          40x, f20.6)
  
    7 format (15x, ' Coulomb Direct:', 3f20.6)
   !71 format (15x, '   Dir. (point):', 3f20.6)
    8 format (15x, '       Exchange:', 3f20.6)
   !81 format (15x, '   Exc. (point):', 3f20.6)  

   !    64 format (15x, '  Vibrational '   , a1, ':', 30x, f15.6)
  641 format (15x, '  Vibrational T:',          30x, f15.6)

    9 format (15x, 'Pairing (delta):', 3f20.6)
   91 format (15x, 'Pairing (densi):', 40x, f20.6)
   92 format ( 7x, 'Pair. (delta, no stab.):', 3f20.6)
   93 format ( 7x, 'Pair. (densi, no stab.):', 40x, f20.6)
   94 format ( 7x, 'Pair. (delta,    stab.):', 3f20.6)
   95 format ( 7x, 'Pair. (densi,    stab.):', 40x, f20.6)

   99 format (15x, '   Total energy:', 40x, f20.6)
  991 format (15x, '    (no corr.) :', 40x, f20.6)
  100 format (15x, '     from spwfs:', 40x, f20.6)
  101 format (15x, '    Free Energy:', 40x, f20.6)
  102 format (15x, '        Entropy:', 3f20.6)
  103 format (15x, '    E_fu - E_sp:', 40x, e20.6)
  104 format (15x, '          dE   :', 40x, e20.6)
  105 format (15x, '       Routhian:', 40x, f20.6)
  106 format (15x, '          dR   :', 40x, e20.6)           

#if(PASTA > 0)
  107 format (40x, '         FOR PASTA CALCULATIONS    ')
  108 format (15x, '        e_pasta=(Total energy + electrons - Z[Mn-Mp])/A - Mn')
  109 format (15x, '        e_pasta:', 40x, f20.6)
  110 format (15x, '   Electron kin:', 40x, f20.6) 
  111 format (15x, '  Electron exch:', 40x, f20.6) 
  112 format (15x, 'Chempot_e total:', 40x, f20.6)
  113 format (15x, '      Chempot_n:', 40x, f20.6)
  114 format (15x, '      Chempot_p:', 40x, f20.6)
  115 format (15x, ' Chempot_p β-eq:', 40x, f20.6)
  116 format (15x, '       Pressure:', 40x, f20.6)
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

    if(rotcorr .ne.  0 .or. vibcorr .ne. 0) then
      print 65, sum(Rotcorrection) + Vibcorrection
      print *
      print 631, sum(Rotcorrection)
      print 63, 'X',  Rotcorrection(1)
      print 63, 'Y',  Rotcorrection(2)
      print 63, 'Z',  Rotcorrection(3)
      print *
      print 641, Vibcorrection
!       print 64, 'X',  Vibcorrection(1)
!       print 64, 'Y',  Vibcorrection(2)
!       print 64, 'Z',  Vibcorrection(3)
      print *
    endif

    print *
    print 7, 0.0d0, CoulombDirect, CoulombDirect
    !if(protonsize(1).ne.0 .and. (.not. nucleonsize_selfconsistent)) then
    !  temp = CoulombEnergy_Direct(Density, Potentials)
    !  print 71, 0.0d0, temp, temp
    !endif
    print 8, 0.0d0, CoulombExchange, CoulombExchange
    !if(protonsize(1).ne.0 .and. (.not. nucleonsize_selfconsistent)) then
    !  temp = CoulombEnergy_Exchange(Density)
    !  print 81, 0.0d0, temp, temp
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
        &                 - vibcorrection
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
    print 109, calculate_epasta(totalE)
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

#if(PASTA==1)

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
#endif


  function transfer_coulomb_mesh(F, exchange) result(pot)
    !---------------------------------------------------------------------------
    ! Restrict a Coulomb potential defined on the large grid defined as a 3D 
    ! array to the grid of the rest of the code. This is abstracted because of 
    ! the "if" statements  below.
    ! 
    ! Input:
    !   F  : potentialvector containing the relevant fields
    !   exchange : logical, determines whether to ask for the exchange or 
    !              direct coulomb potential
    ! Output:
    !   pot: the requested coulomb field on the grid of the Skyrme potentials.
    !---------------------------------------------------------------------------
    use coulombmod

    Type(PotentialVector), intent(in) :: F
    logical, intent(in)               :: exchange
    integer                           :: i,j,k, ox, oy, oz, it
#if($FAM == 0)
    ! Coulomb potentials are real for mean-field calculations
    real(KIND=dp)                     :: pot(mv,4) 
#else
    ! Coulomb potentials are complex for FAM calculations
    complex(KIND=dp)                  :: pot(mv,4)
#endif
    pot = 0.0d0    
    if((all(protonsize.eq.0.0d0) .and. all(neutronsize.eq.0.0d0)) .or.     &
      &                             (.not. nucleonsize_selfconsistent)) then

      !------------------------------------------------------------------------
      ! Protons and neutrons are treated as point particles
      !------------------------------------------------------------------------
      ox = coul_offset_x ; oy = coul_offset_y ; oz = coul_offset_z
      do k=1,nz
        do j=1,ny
          do i=1,nx
            if(.not. exchange) then
              pot(meshindex(i,j,k),2)=  F%CoulombPotential(i+ox,j+oy,k+oz)
            else
              pot(meshindex(i,j,k),2)=  F%ExchangePotential(i,j,k)
            endif
          enddo
        enddo
      enddo
    else
      !-----------------------------------------------------------------------
      ! Use the folded coulombpotential, for full self-consistency.
      ! Note that both protons and neutrons feel a Coulomb force if their
      !  charge form factor is taken into account.
      !-----------------------------------------------------------------------
      if(.not. allocated(F%foldedcoul)) then
        call stp('Nucleonsize_selfconsistent cannot be .false. if the protons are not point particles.')
      endif
      do it=1, 2
        do k=1,nz
          do j=1,ny
            do i=1,nx
              if(.not. exchange) then
                pot(meshindex(i,j,k),it)= F%FoldedCoul(i,j,k,it)
              else
                pot(meshindex(i,j,k),it)= F%FoldedExchange(i,j,k,it)
              endif
            enddo
          enddo
        enddo
      enddo
    endif
    
    ! Make sure the potentials are consistent among neutron/proton and isospin 0/1 channels.
    pot(:,3) = pot(:,1) + pot(:,2)
    pot(:,4) = pot(:,1) - pot(:,2)
    
  end function transfer_coulomb_mesh
  

  subroutine set_coul(coul, F)
    !---------------------------------------------------------------------------
    ! Set the Coulombpotential of potentialvector F to the input function coul, 
    ! at least within the zone that corresponds to the mesh of the other
    ! potentials. This subroutine achieves more or less the inverse of 
    ! transfer_coulomb_mesh.
    !
    ! Input:
    !   coul :  Coulomb potential, defined on a (nx*ny*nz,4) grid
    !      F : potentialvector
    ! Output:
    !      F : potentialvector with its direct Coulomb potential equal 
    !          coul (with the exception of the points on the boundary.)
    !---------------------------------------------------------------------------
    use coulombmod
    
    real(KIND=dp), intent(in)            :: coul(nx*ny*nz,4)
    Type(PotentialVector), intent(inout) :: F
    integer                              :: i,j,k, ox, oy, oz, it

    if((all(protonsize.eq.0.0d0) .and. all(neutronsize.eq.0.0d0)) .or.     &
      &                             (.not. nucleonsize_selfconsistent)) then
      !------------------------------------------------------------------------
      ! Protons and neutrons are treated as point particles
      !------------------------------------------------------------------------
      ox = coul_offset_x ; oy = coul_offset_y ; oz = coul_offset_z
      do k=1,nz
        do j=1,ny
          do i=1,nx
            F%CoulombPotential(i+ox,j+oy,k+oz) = coul(meshindex(i,j,k),2) 
          enddo
        enddo
      enddo
    else
      !-----------------------------------------------------------------------
      ! Use the folded coulombpotential, for full self-consistency.
      ! Note that both protons and neutrons feel a Coulomb force if their
      !  charge form factor is taken into account.
      !-----------------------------------------------------------------------
      if(.not. allocated(F%foldedcoul)) then
        call stp('Nucleonsize_selfconsistent cannot be .false. if the protons are not point particles.')
      endif
      do it=1,2
        do k=1,nz
          do j=1,ny
            do i=1,nx
              F%FoldedCoul(i,j,k,it) =coul(meshindex(i,j,k),it) 
            enddo
          enddo
        enddo
      enddo
    endif
  
  end subroutine set_coul


  function INM_k2_pot_real(rho) result(pot)
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
  end function INM_k2_pot_real


  function INM_k2_pot_complex(rho) result (pot)
    ! Same as above, but for complex densities.
    complex(KIND=dp), intent(in) :: rho(mv,4)
    real(KIND=dp)                :: pot(mv,4)
    real(KIND=dp)                :: kfn(mv), kfp(mv)

    kfn=(3.d0*pi**2*DBLE(rho(:,1)))**(1.0d0/3.0d0) ! Neutron density
    kfp=(3.d0*pi**2*DBLE(rho(:,2)))**(1.0d0/3.0d0) ! Proton  density

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
  end function INM_k2_pot_complex


  function INM_k4_pot_real(rho) result(pot)
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
end function INM_k4_pot_real


function INM_k4_pot_complex(rho) result(pot)
  ! 
  complex(KIND=dp), intent(in) :: rho(mv,4)
  real(KIND=dp)                :: pot(mv,4)
  real(KIND=dp)             :: kfn(mv), kfp(mv)

  kfn=(3.d0*pi**2*DBLE(rho(:,1)))**(1.0d0/3.0d0) ! Neutron density
  kfp=(3.d0*pi**2*DBLE(rho(:,2)))**(1.0d0/3.0d0) ! Proton  density

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
end function INM_k4_pot_complex


  function apply_sphamil(psi, dpsi, ddpsi, &
$N3                                 dddpsi, &
&                                        sx,sy,sz,iso, onthefly, Fin) result(hpsi)
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
    ! Fin          : a set of mean-field potentials determining the
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
    type(PotentialVector), intent(in) :: Fin
    type(PotentialVector)             :: F

    real(KIND=dp)             :: hpsi(mv,4)
    real(KIND=dp)             :: temp(mv,4)
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Declaration of temporary spinors
    !
    ! Technical note: these arrays are declared with the spinor indices (1-4)
    !                 BEFORE the derivative indices (3). This is to aid the
    !                 memory locality of operations in this particular function
    !                 and is OPPOSITE the conventions of the rest of the code.
$DTEMPSPH     integer         :: sym(4)
$DTEMPSPH     real(KIND=dp)   ::   dtemp(mv,4,3)
$D2TEMPSPH    real(KIND=dp)   ::  ddtemp(mv,4,3,3)
$D3TEMPSPH    real(KIND=dp)   :: dddtemp(mv,4,3,3,3)
$LAPTEMPSPH   real(KIND=dp)   :: laptemp(mv,4)

    real(KIND=dp)             :: ReducedMass
#if(PASTA == 0)
    real(KIND=dp)             :: Butler_t, Butler_f
#endif
    integer :: it, i,k

    call start_timer(T_sphamil)
    !---------------------------------------------------------------------------
    ! Dirty trick: the potential vector Fin contains all relevant information on
    !              the potentials, namely F_I_I and the coulombpotential.
    !              Below however, we only act with F_I_I; for the purpose of this
    !              routine we add these potentials together.
    !---------------------------------------------------------------------------
    ! TODO: add in constraint_I_I to the definition of a potentialvector
    F = Fin
    call combine_potentials(F)
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
      Butler_t = (1.5d0 * (neutrons + protons))**(1.d0/3.d0)
      Butler_f = 2.d0/(Butler_t + 1.d0/(3*Butler_t))
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
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    ! ----------
    !         psi :
    !           spwf to act on with delta
    !  d/dd/dddpsi :
    !           arrays containing the first, second and third derivatives
    !           of the spwf. 
    ! sx/sy/sz :
    !           signs under reflection symmetry for this particular spwf
    !           not referenced when onthefly = .false.
    ! iso :
    !   isospin of the spwf: +1 for protons, -1 for neutrons
    !
    ! onthefly:
    !   Logical indicating if the derivatives need to be calculated before
    !   applying delta. If false, the derivatives are passed in. If True, the
    !   derivatives are not passed in and need to be calculated.
    !   NOT IMPLEMENTED YET 
    !
    ! F :
    !   potentialvector containing the pairing fields
    !
    ! Output:
    ! --------
    !   deltapsi :
    !     delta | psi >
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
    ! Start by zeroing the array
    allocate(deltapsi(mv,4))
    deltapsi = 0.0

$PAIRINGACTION

  end function delta_action


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


end module functional
