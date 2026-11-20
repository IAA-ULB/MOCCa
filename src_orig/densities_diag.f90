submodule (densities) densities_diag
!==============================================================================
!
! Implementation submodule of module densities: densities_diag
! The interfaces live in module densities (src_orig/densities.f90).
!
!==============================================================================


 implicit none

contains

module function densit(rho, kappa) result(R)
    !---------------------------------------------------------------------------
    ! Calculate all of the mean-field densities, both normal and pairing. 
    ! This includes the calculation of the charge density.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Input :
    !   rho      real vector     diagonal matrix elements of the density in
    !                            either the HF or canonical basis
    !   kappa    real matrix     matrix elements of the anomalous density
    !                            in the HF-basis.
    ! Output:
    !   R        densityvector   values of the mean-field densities.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Remarks: 
    ! * the full-matrix storage of kappa in the BCS case is wasteful, but is 
    !   necessary to not have to change the interface of this function.
    ! * this routine assumes the right set of SPWFs are in the right location,
    !   meaning that (if necessary) the canonical basis has already been 
    !   constructed.  
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in) :: rho(:), kappa(:,:)
    type(DensityVector)       :: R

    integer                    :: i, k, it, wave, wave2, B, N, si, N2, T
    integer                    :: wave_global, wave2_global, der_index
    real(KIND=dp)              :: weight
    real(KIND=dp), allocatable :: kappa_cut(:,:)

$SPWF_DECLARATION

#if(USE_MPI>0)
    integer      :: mpi_err
#endif
    call start_timer(T_densities)

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Allocation and initialization
$INITIALIZATION


    if(.not.allocated(R%divJ)) then
      allocate(R%divJ(nx*ny*nz,4))
    endif
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Zero the current density
$ZEROING
    R%divJ = 0.0d0
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
      if(.not. efficientHFB) then
        DenPsi    => CanPsi   ; DenDPsi   => CanDPsi 
        DenddPsi  => CanddPsi ; DendddPsi => Candddpsi
      else
        DenPsi    => HFPsi    ; DenDPsi   => HFDPsi 
        DenddPsi  => HFddPsi  ; DendddPsi => HFdddpsi      
      endif
    end select
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    call start_timer(T_den_ph)
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! PARTICLE-HOLE DENSITIES
    do wave=1,nwt_local              ! Loop over the local spwf index
        wave_global = spwf_map(wave) ! Global spwf index
    
        ! Isospin is neutron in the first half of blocks, proton in the rest
        it = 2
        if(wave_global.le.nwn) it = 1
        
        ! For ordinary densities
        weight  = rho(wave_global) 

        !---------------------------------------------------------------------------
        ! Perform some gymnastics to see where we are getting the information
        ! on derivatives: in memory (target_index = wave) or by calculation
        ! on the fly.
        if(store_derivatives) then
          ! We have centrally stored derivative information
          der_index = wave
        else
          der_index = 1
#if(USE_Periodic==0)
          do k=1,4
           call Derive_tot(denPsi(:,k,wave), sx(k,wave), sy(k,wave), sz(k,wave), &
           &                                           dendPsi(:,:,k,der_index), &
           &                                           denddPsi(:,:,k,der_index))
          enddo
#else
          do k=1,2
           call Derive_tot_periodic(denPsi   (:,   (2*k-1):2*k,wave),     &
           &                        sx       ((2*k-1):2*k,wave),          &
           &                        sy       ((2*k-1):2*k,wave),          &
           &                        sz       ((2*k-1):2*k,wave),          &
           &                        dendpsi  (:,:,(2*k-1):2*k,der_index), &
           &                        denddpsi (:,:,(2*k-1):2*k,der_index))
          enddo
#endif
        endif


        !---------------------------------------------------------------------------

        do i=1,mv
$EXPRESSION
$DERIVATION_SUM_SPWF_PH
        enddo

        ! Separate expression for the sum of divJ - see comments in that function
        R%divJ(:,it) = R%divJ(:,it) + weight * divJ_spwf(der_index)
    enddo
    call stop_timer(T_den_ph)

    call start_timer(T_den_pp)
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! PAIRING DENSITIES
    select case (PairingType) 
    case(0)
      !----------------------------------------
      ! HF calculation, no need to get the gaps
      !----------------------------------------
    case(1)
      !----------------------------------------------
      ! BCS calculation, the sums are over i == ibar.
      !----------------------------------------------
      do wave=1,nwt_local              ! local index of the spwf
          wave_global = spwf_map(wave) ! global index of the spwf
          ! Isospin is neutron in the first half of blocks, proton in the rest
          it = 2
          if(wave_global.le.nwn) it = 1
          weight  = 2 * kappa_can(wave_global) * Pcutoffs(wave_global)**2
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
          ! Note about the factor two in the weight:
          ! 
          ! \tilde{\rho}(r) = \sum_{ij} \kappa_{ij} \phi_j(r) \phi_i(r)
          ! 
          ! Now kappa_ij is only non-zero if both spwfs have opposite signature, 
          ! meaning that the summation can be split into two parts.
          !
          ! \tilde{\rho}(r) = \sum_{i>0,j<0} \kappa_{ij} \chi_{ji}(r)
          !                 + \sum_{i<0,j>0} \kappa_{ij} \chi_{ji}(r)
          !
          ! with chi_ji(r) the spatial contribution of spwfs i and j to whatever
          ! the pairing density is. These two parts are identical because of
          ! the skew symmetry of chi and kappa.
          !
          ! Historically in MOCCa and this code, this factor two was absorbed 
          ! in the definition of the pairing strengths, but no longer.
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

          do i=1,mv
$BCSEXPRESSION
$DERIVATION_SUM_SPWF_BCS
          enddo
      enddo
    case(2)
      !-------------------------------------------------------------------------
      ! HFB calculations: full summations.
      !-------------------------------------------------------------------------
      ! Sum the pairing densities in the HFbasis. This could be done in the 
      ! canonical basis, but this would surely be less straightforward because
      ! of the presence of pairing cutoffs.
      DenPsi    => HFPsi   ; DenDPsi   => HFdPsi 
      DenddPsi  => HFddPsi ; DendddPsi => HFdddpsi

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! Tricky loop-structure due to MPI implementation. 
      ! The previous implementation for kappa_cut did not generalize easily, 
      ! since it was constructed for each symmetry-block separately.
      ! The new implementation just constructs kappa_cut completely before 
      ! summing any of the pairing densities.
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      ! a) start out by "just" copying kappa
      allocate(kappa_cut(nwt,nwt))
      kappa_cut = kappa
      if((.not. diagsphamil)) then
        kappa_cut = transform_mat(kappa_cut, HFTransfo)
      endif

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
        ! W.R. Nasty bug 07/04/'22
        ! The multiplication by the cutoffs needs to happen for ALL matrix
        ! elements of kappa, not just the ones that are used to sum the 
        ! pairing densities. The reason is that the HFTransfo transformation 
        ! does see all of them, at least as it is coded at the moment.
        ! We could use the symmetries of kappa to reduce the workload, but 
        ! this is only O(nwt**2) effort and the implementation would depend on 
        ! whether or not Timereversal is broken.
        do wave=1,T
          do wave2=1,T
            kappa_cut(si+wave,si+wave2) = kappa_cut(si+wave,si+wave2) &
            &                     *Pcutoffs(si+wave)*Pcutoffs(si+wave2)
          enddo
        enddo
        si = si + N + N2
      enddo

      ! transform back to the HFbasis
      if((.not. diagsphamil)) then
        kappa_cut = transform_mat(kappa_cut, transpose(HFTransfo))
      endif

      ! with kappa_cut in hand, we can turn to the summation of the densities.
      si = 0
      do B=1,8,2  ! <------- this loop ranges over the LOCAL set of spwfs
        N = HFBlocks(B) ;  if (N.eq.0) cycle
        N2= HFBlocks(B+1)
        T = N+N2
        it = 2          
        if( B.le. 4) it = 1

        do wave=1,N                         ! local index of the spwf
          wave_global = spwf_map(si+wave)   ! global index of the spwf
$TR          do wave2=wave,N               
$NTR          do wave2=N+1,N+N2      
                wave2_global = spwf_map(si+wave2)

            weight=2*kappa_cut(wave_global,wave2_global)
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
            ! Note about the factor two in the weight:
            ! 
            ! \tilde{\rho}(r) = \sum_{ij} \kappa_{ij} \phi_j(r) \phi_i(r)
            !  
            ! Now kappa_ij is only non-zero if both spwfs have opposite signature, 
            ! meaning that the summation can be split into two parts.
            !
            ! \tilde{\rho}(r) = \sum_{i>0,j<0} \kappa_{ij} \chi_{ji}(r)
            !                 + \sum_{i<0,j>0} \kappa_{ij} \chi_{ji}(r)
            !
            ! with chi_ji(r) the spatial contribution of spwfs i and j to whatever
            ! the pairing density is. These two parts are identical because of
            ! the skew symmetry of chi and kappa.
            !
            ! Historically in MOCCa and this code, this factor two was absorbed 
            ! in the definition of the pairing strengths, but no longer.
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
            ! In the HFB case, there is an EXTRA degeneracy that we exploit.
            ! When signature is conserved we know that 
            !
            !         ( 0          kappa^{+-})
            ! kappa = (                      )
            !         ( kappa^{-+}  0        )
            !
            ! and using this structure is exactly the factor two used above.
            !
            ! If there is time-reversal symmetry, we know in addition that 
            ! kappa^{+-} itself is symmetric, hence the factor two below this 
            ! comment. Note that it doesn't get applied for the diagonal
            ! component.
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
$TR            if(wave.ne.wave2) weight = 2 * weight

            do i=1,mv
$HFBEXPRESSION
$DERIVATION_SUM_SPWF_HFB
            enddo
          enddo
        enddo
        si = si + N + N2
      enddo
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      ! We (possibly) changed these pointers, and now make sure they are still
      ! pointing the right way
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      select case(PairingType)
      case(0,1)
        ! HF or BCS Calculation
        DenPsi   => HFPsi    ; DenDPsi   => HFDPsi 
        DenddPsi => HFddPsi  ; DendddPsi => HFdddpsi
      case(2)
        ! HFB calculation
        if(.not. efficientHFB) then
          DenPsi    => CanPsi   ; DenDPsi   => CanDPsi 
          DenddPsi  => CanddPsi ; DendddPsi => Candddpsi
        else
          DenPsi    => HFPsi    ; DenDPsi   => HFDPsi 
          DenddPsi  => HFddPsi  ; DendddPsi => HFdddpsi      
        endif
      end select
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      ! clean up explicitly
      deallocate(kappa_cut)
    end select
    call stop_timer(T_den_pp)
    
#if(USE_MPI > 0)
   call start_timer(T_allreduce)
   ! Sum the density over all processes
$MPIDEN
   call stop_timer(T_allreduce)
#endif
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Calculation of the 'derived' densities, densities obtainable by 
    ! deriving other ones. 
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    call start_timer(T_den_der)
    do it=1,2
$DERIVATION
    enddo  
    call stop_timer(T_den_der)

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Calculate the densities in isospin representation 
$ISOSPINCOUPL

    R%divJ(:,3) = R%divJ(:,1) + R%divJ(:,2)
    R%divJ(:,4) = R%divJ(:,1) - R%divJ(:,2)
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Construct the charge density
    call construct_charge_density(R, sx_rho, sy_rho, sz_rho)

    call stop_timer(T_densities)
end function densit


end submodule densities_diag
