submodule (densities) densities_offdiag
!==============================================================================
!
! Implementation submodule of module densities: densities_offdiag
! The interfaces live in module densities (src_orig/densities.f90).
!
!==============================================================================


 implicit none

contains

module subroutine densit_offdiag(rho, kappa_plus, kappa_minus, Rs, Ra, R_pp_plus, R_pp_minus)
    !----------------------------------------------------------------------------
    ! Calculate normal and anomalous densities through a double sum across spwfs
    ! by summing symmetric and antisymmetric parts.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input :
    !   rho        :  complex matrix
    !   kappa_plus : (antisymmetric) complex matrix
    !   kappa_minus: (antisymmetric) complex matrix 
    !
    ! Output:
    !   Rs   : density vector containing the symmetric part of the densities
    !   Ra   : density vector containing the antisymmetric part of the densities
    !   R_pp_plus  : density vector containing the pairing densities for
    !                 kappa_plus
    !   R_pp_minus : density vector containing the pairing densities for
    !                 kappa_minus
    !----------------------------------------------------------------------------
    complex(KIND=dp), intent(in)     :: rho(:,:), kappa_plus(:,:), kappa_minus(:,:)
    type(DensityVector), intent(out) :: Rs, Ra, R_pp_plus, R_pp_minus

    call start_timer(T_den_perturbed)

    ! building the ph densities
    Rs = densit_offdiag_ph_symmetric(rho)      ! symmetric ph densities
    Ra = densit_offdiag_ph_antisymmetric(rho)  ! antisymmetric ph densities
    R_pp_plus = densit_offdiag_pp(kappa_plus)  ! pp densities for kappa_plus 
    R_pp_minus= densit_offdiag_pp(kappa_minus) ! pp densities for kappa_minus
    call stop_timer(T_den_perturbed)

    !call print_maxval('D_I_I'  , Rs%D_I_I  , Ra%D_I_I)
    !call print_maxval('D_Nm_Nm', Rs%D_Nm_Nm, Ra%D_Nm_Nm)
    !call print_maxval('D_I_Sx', Rs%D_I_S(:,1,:), Ra%D_I_S(:,1,:))
    !call print_maxval('D_I_Sy', Rs%D_I_S(:,2,:), Ra%D_I_S(:,2,:))
    !call print_maxval('D_I_Sz', Rs%D_I_S(:,3,:), Ra%D_I_S(:,3,:))
    !call print_maxval('C_I_Nx', Rs%C_I_N(:,1,:), Ra%C_I_N(:,1,:))
    !call print_maxval('C_I_Ny', Rs%C_I_N(:,2,:), Ra%C_I_N(:,2,:))
    !call print_maxval('C_I_Nz', Rs%C_I_N(:,3,:), Ra%C_I_N(:,3,:))
    !call print_maxval('C_I_NSxy', Rs%C_I_NS(:,1,2,:), Ra%C_I_NS(:,1,2,:))

end subroutine densit_offdiag

module function densit_offdiag_restricted(rho, kappa) result(R)
  !----------------------------------------------------------------------------
  ! Calculate the mean-field densities for the restricted set of density
  ! matrices, i.e., diagonal rho and kappa in the canonical basis.
  !
  ! TODO: document
  !
  ! Input :
  !   rho      : real matrix
  !   kappa    : real matrix
  !
  ! Output:
  !   R        : densityvector, values for the ph and pp densities
  !----------------------------------------------------------------------------
    
  real(KIND=dp), intent(in)     :: rho(:,:), kappa(:,:)
  complex(KIND=dp), allocatable :: rho_temp(:,:), kappa_plus_temp(:,:), kappa_minus_temp(:,:)
  type(DensityVector)           :: R
  type(DensityVector)           :: Rs, Ra, R_pp_plus, R_pp_minus

  allocate(rho_temp(nwt,nwt), kappa_plus_temp(nwt,nwt), kappa_minus_temp(nwt,nwt))
  rho_temp         = rho 

  kappa_plus_temp  = kappa
  kappa_minus_temp = 0.0d0

  call densit_offdiag(rho_temp, kappa_plus_temp, kappa_minus_temp, Rs, Ra, R_pp_plus, R_pp_minus)
  ! Combine the correct densities
  R = Rs + R_pp_plus

end function densit_offdiag_restricted

module function densit_offdiag_pp(kappa) result(R_total)
    !----------------------------------------------------------------------------
    ! Calculate the pairing mean-field densities, based on arbitrary
    ! anomalous density matrix kappa.
    !
    ! Input :
    !   kappa    : complex matrix
    !   R        : densityvector, values for the ph densities will be maintained
    !
    ! Output:
    !   R        : densityvector, this time with values for the pp densities
    !              added. 
    !----------------------------------------------------------------------------
    COMPLEX(KIND=dp), intent(in)       :: kappa(:,:)
    type(DensityVector)                :: R, R_total

    integer                            :: wave, wave2, i, B, it, N, si, N2, T
    INTEGER                            :: wave_global, wave2_global, der_index
    complex(KIND=dp)                   :: weight
    complex(KIND=dp), allocatable      :: kappa_cut(:,:)
    ! TODO: this for sure declares too much
$SPWF_DECLARATION
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Allocation and initialization
$INITIALIZATION

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Zero the current density
$ZEROING
    ! Ensure allocation of charge density to avoid trouble when combining with
    !  other density vectors
    allocate(R%chargedensity(nx,ny,nz)) ; R%chargedensity = 0.0d0

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Make a copy to ensure everything is allocated and zero'd
    R_total = R

    call start_timer(T_den_perturbed_pp)
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! PAIRING DENSITIES
    select case (PairingType) 
    case(0)
      !----------------------------------------
      ! HF calculation, no need to sum pairing densities
      !----------------------------------------
    case(1)
      call stp('Linear response calculations for HF+BCS not implemented yet.')
    case(2)
      !-------------------------------------------------------------------------
      ! HFB calculations: full summations.
      !-------------------------------------------------------------------------
      ! Sum the pairing densities in the HFbasis. This could be done in the 
      ! canonical basis, but this would surely be less straightforward because
      ! of the presence of pairing cutoffs.
      DenPsi    => HFPsi   ; DenDPsi   => HFdPsi 
      DenddPsi  => HFddPsi ; DendddPsi => HFdddpsi

      ! a) start out by "just" copying kappa
      allocate(kappa_cut(nwt,nwt))
      kappa_cut = kappa

      ! ... and then generating the entirety of kappa_cut
      si = 0
      do B=1,8,2  ! <------- this loop ranges over the global set of spwfs
        N = HFBlocks_global(B) ;  if (N.eq.0) cycle
        N2= HFBlocks_global(B+1)
        T = N+N2
        it = 2
        if( B.le. 4) it = 1

        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ! Calculation of the pairing cutoffs * kappa
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        do wave=1,T
          do wave2=1,T
            kappa_cut(si+wave,si+wave2) = kappa_cut(si+wave,si+wave2) &
            &                     *Pcutoffs(si+wave)*Pcutoffs(si+wave2)
          enddo
        enddo
        si = si + N + N2
      enddo

!$OMP PARALLEL PRIVATE(si, B, N, N2, T, it,           &
!$OMP                  wave      , wave2          ,   &
!$OMP                  wave_global, wave2_global  ,   &
!$OMP                  i, weight,                     &
$OMP_DENSITY_VARS
!$OMP                 )                               &
!$OMP          FIRSTPRIVATE (R)                       &
!$OMP          SHARED ( kappa_cut, HFBlocks, R_total, &
!$OMP                   spwf_map, denpsi, dendpsi,    &
!$OMP                   mv, denddpsi) DEFAULT(PRIVATE)

      ! with kappa_cut in hand, we can turn to the summation of the densities.
      do B=1,8,2  ! <------- this loop ranges over the LOCAL set of spwfs
        N = HFBlocks(B) ;  if (N.eq.0) cycle
        N2= HFBlocks(B+1)
        T = N+N2

        si = sum(HFBlocks(1:B-1)) ! OMP requires the direct calculation of offset
        it = 2
        if( B.le. 4) it = 1
!$OMP DO SCHEDULE(STATIC)
!       Note: these loops cannot be OMP COLLAPSEd, as they form a
!       non-rectangular iteration space. Some compilers would accept this,
!       but not all.
        do wave=1,N                         ! local index of the spwf
$TR          do wave2=wave,N
$NTR          do wave2=N+1,N+N2      
                wave_global  = spwf_map(si+wave)   ! inside the wave2_loop to make the COLLAPSE statement work
                wave2_global = spwf_map(si+wave2)
                ! This factor two is the antisymmetry of \kappa_cut
                weight=2*kappa_cut(wave_global,wave2_global)
$TR             if(wave.ne.wave2) weight = 2 * weight
            do i=1,mv
$HFBEXPRESSION  
            enddo
          enddo
        enddo
!$OMP END DO
      enddo
!$OMP CRITICAL
    R_total = R_total + R
!$OMP END CRITICAL
!$OMP END PARALLEL
    end select

    call stop_timer(T_den_perturbed_pp)
end function densit_offdiag_pp

module function densit_offdiag_ph_antisymmetric(rho) result(R_total)
    !------------------------------ ---------------------------------------------
    ! Calculate the antisymmetric part(*) of the particle-hole mean-field densities, 
    ! based on arbitrary matrix rho
    !
    ! TODO: explain "antisymmetric part"
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input :
    !   rho      real/complex matrix
    !
    ! Output:
    !   R_total  densityvector   values of the mean-field densities.
    !----------------------------------------------------------------------------
    complex(KIND=dp), intent(in) :: rho(:,:)
    type(DensityVector)          :: R, R_total

    complex(KIND=dp)          :: weight_asym
    integer                   :: wave_i       , wave_j
    integer                   :: wave_global_i, wave_global_j
    integer                   :: it_i, it_j, it, der_index_i,der_index_j, si, B, N

    integer :: i

$SPWF_DECLARATION
    call start_timer(T_den_perturbed_asym)

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Allocation and initialization
$INITIALIZATION
    ! Explicit allocation of this array; below we will be incrementing densities!
    if(.not. allocated(R%chargedensity)) allocate(R%chargedensity(nx,ny,nz))

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Zero the current density
$ZEROING
    R%chargedensity = 0.0d0
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Correctly set the pointers to the spwfs
    ! This should always be the HF basis!
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    DenPsi   => HFPsi    ; DenDPsi   => HFDPsi
    DenddPsi => HFddPsi  ; DendddPsi => HFdddpsi
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Make a copy to ensure everything is allocated
    R_total = R
    call start_timer(T_den_ph)
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! PARTICLE-HOLE DENSITIES
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!$OMP PARALLEL PRIVATE(si, B, N, it,                &
!$OMP                  wave_i       ,wave_j       , &
!$OMP                  wave_global_i,wave_global_j, &
!$OMP                  der_index_i  ,der_index_j,i, &
!$OMP                  weight_asym,                 &
$OMP_DENSITY_VARS
!$OMP                 )                             &
!$OMP          FIRSTPRIVATE (R)                     &
!$OMP          SHARED ( rho, HFBlocks, R_total,     &
!$OMP                   spwf_map, denpsi, dendpsi,  &
!$OMP                   mv, denddpsi) DEFAULT(PRIVATE)
    ! Note: R is FIRSTPRIVATE, to ensure it is correctly ZERO'd in EVERY
    !       thread.

    si = 0
    do B=1,8
      N = HFBLocks(B) ; if(N.eq.0) cycle
      si = sum(HFBlocks(1:B-1)) ! offset has to be explicitly calculated for OpenMP
      ! Isospin is neutron in the first half of blocks, proton in the rest
      it = 1
      if(B .ge. 5) it = 2
!$OMP DO COLLAPSE(2) SCHEDULE(STATIC)
      do wave_i=si+1,si+N                ! Loop over the local spwf index
        do wave_j=si+1,si+N
          ! Note: the i-based indexing could be moved to the i-loop, but not all
          !       compiler versions then accept the COLLAPSE statement in the
          !       OpenMP clause.
          wave_global_i = spwf_map(wave_i) ! Global spwf index
          der_index_i = wave_i
          wave_global_j = spwf_map(wave_j) ! Global spwf index
          der_index_j = wave_j
          !----------------------------------------------------------------------------
          ! The summation weights for particle-hole densities
          weight_asym = 0.5d0*( &
          &      rho(wave_global_i, wave_global_j) - rho(wave_global_j, wave_global_i))
$TR       weight_asym = 2 * weight_asym ! <------ factor two for time-reversal symmetry 

          do i=1,mv
$EXPRESSION_OFFDIAG_ANTISYMMETRIC
          enddo
        enddo
      enddo
!$OMP END DO
    enddo

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Calculation of the 'derived' densities, densities obtainable by
    ! deriving other ones.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    do it=1,2
$DERIVATION_OFFDIAG_ANTISYMMETRIC
    enddo

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Calculate the densities in isospin representation
$ISOSPINCOUPL_ANTISYMMETRIC

!$OMP CRITICAL
    R_total = R_total + R
!$OMP END CRITICAL
!$OMP END PARALLEL
    call stop_timer(T_den_ph)
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Construct the charge density
    call construct_charge_density(R_total, sx_rho_antisym, sy_rho_antisym, sz_rho_antisym)

    call stop_timer(T_den_perturbed_asym)

end function densit_offdiag_ph_antisymmetric

end submodule densities_offdiag
