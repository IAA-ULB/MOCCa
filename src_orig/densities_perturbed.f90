submodule (densities) densities_perturbed
!==============================================================================
!
! Implementation submodule of module densities: densities_perturbed
! The interfaces live in module densities (src_orig/densities.f90).
!
!==============================================================================


 implicit none

contains

module function densit_offdiag_ph_symmetric(rho) result(R_total)
    !------------------------------ ---------------------------------------------
    ! Calculate the symmetric part(*) of the particle-hole mean-field densities, 
    ! based on arbitrary matrix rho.
    !
    ! TODO: explain "symmetric part"
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input :
    !   rho      real/complex matrix
    !
    ! Output:
    !   R_total     densityvector   values of the mean-field densities.
    !----------------------------------------------------------------------------

    complex(KIND=dp), intent(in) :: rho(:,:)
    type(DensityVector)          :: R, R_total

    complex(KIND=dp)          :: weight_sym
    integer                   :: wave_i       , wave_j
    integer                   :: wave_global_i, wave_global_j, B, si, N
    integer                   :: it_i, it_j, it, der_index_i,der_index_j

    integer :: i

$SPWF_DECLARATION
    call start_timer(T_den_perturbed_sym)

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
    ! Make a copy to ensure everything is allocated and zero'd
    R_total = R
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Correctly set the pointers to the spwfs
    ! This should always be the HF basis!
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    DenPsi   => HFPsi    ; DenDPsi   => HFDPsi
    DenddPsi => HFddPsi  ; DendddPsi => HFdddpsi
    call start_timer(T_den_ph)
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! PARTICLE-HOLE DENSITIES
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!$OMP PARALLEL PRIVATE(si, B, N, it,                &
!$OMP                  wave_i       ,wave_j       , &
!$OMP                  wave_global_i,wave_global_j, &
!$OMP                  der_index_i  ,der_index_j,i, &
!$OMP                  weight_sym,                  &
$OMP_DENSITY_VARS
!$OMP                 )                             &
!$OMP          FIRSTPRIVATE (R)                     &
!$OMP          SHARED ( rho, HFBlocks, R_total,     &
!$OMP                   spwf_map, denpsi, dendpsi,  &
!$OMP                   mv, denddpsi) DEFAULT(PRIVATE)
    ! Note: R is FIRSTPRIVATE, to ensure it is correctly ZERO'd in EVERY
    !       thread.
    do B=1,8
      N = HFBlocks(B) ; if(N.eq.0) cycle
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
          wave_global_i = spwf_map(wave_i)  ! Global spwf index
          der_index_i = wave_i
          wave_global_j = spwf_map(wave_j)  ! Global spwf index
          der_index_j = wave_j
          !----------------------------------------------------------------------------
          ! The summation weights for particle-hole densities
          weight_sym = 0.5d0*( &
          &      rho(wave_global_i, wave_global_j) + rho(wave_global_j, wave_global_i))
$TR       weight_sym = 2 * weight_sym ! <------ factor two for time-reversal symmetry

          do i=1,mv
$EXPRESSION_OFFDIAG_SYMMETRIC
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
$DERIVATION_OFFDIAG_SYMMETRIC
    enddo
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Calculate the densities in isospin representation
$ISOSPINCOUPL_SYMMETRIC

!$OMP CRITICAL
    R_total = R_total + R
!$OMP END CRITICAL
!$OMP END PARALLEL
    call stop_timer(T_den_ph)

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Construct the charge density
    call construct_charge_density(R_total, sx_rho, sy_rho, sz_rho)

    call stop_timer(T_den_perturbed_sym)

end function densit_offdiag_ph_symmetric

end submodule densities_perturbed
