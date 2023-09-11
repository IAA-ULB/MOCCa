module densities
!===============================================================================
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
!=============================================================================== 
! Module that defines, calculates and generally deals with all densities. 
! 
!=============================================================================== 
! Hephaestos keywords
!
! DECLARATION     : [WAY too long to include here]
! INITIALIZATION  : [WAY too long to include here]
! ZEROING         : [WAY too long to include here]
! EXPRESSION      : [WAY too long to include here]
! BCSEXPRESSION   : [WAY too long to include here]
! HFBEXPRESSION   : [WAY too long to include here]
! DERIVATION      : [WAY too long to include here]
! ISOSPINCOUPL    : [WAY too long to include here]
! CLEANING        : [WAY too long to include here]
! MPIDEN          : [WAY too long to include here]
!
! TR              : $TR
! NTR             : $NTR 
! 
! SX_RHO          : $SX_RHO 
! SY_RHO          : $SY_RHO
! SZ_RHO          : $SZ_RHO
!
! SX_SX/SY/SZ     : $SX_SX, $SX_SY, $SX_SZ
! SY_SX/SY/SZ     : $SY_SX, $SY_SY, $SY_SZ
! SY_SX/SY/SZ     : $SZ_SX, $SZ_SY, $SZ_SZ
!
! PBROKEN         : $PBROKEN
!
!===============================================================================
!
! A density D_L_R is stored as
!
!      D_L_R (mv, [cartesian indices], [isospin indices])
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
!=============================================================================== 
! Some technical notes:
!
! a) The densities are represented as vectors on the mesh, instead of 3D arrays. 
!    Advantages: 
!      *) Two less indices, which is needed to naively store densities
!         and their derivatives for N2/3LO. Fortran 90 was limited to 
!         arrays of rank 7, while Fortran 2008 lifted this limit to 15
!         many compilers (i.e. gfortran) do not yet support it.
!      *) Less involved coding in the functional.f90 file, since there
!         are less ':' to be put in the sums, and all of the spatial
!         indices are combined into one.
!      *) It MIGHT result in faster summing of densities, but this has never 
!         been a problem on the mean-field level.
!
!    Disadvantages: 
!      *) This needed a wrapper routine for the derivative functions Derive_grad
!         and derive_lap. Currently everything is accomplished by defining 
!         3D pointers, which is hopefully more efficient than a call to RESHAPE.
!         If this ever takes up a significant fraction of computation time, one
!         can think of explicitly writing the 1D derivative code.
!
! b) Note that the storage scheme for derivatives is not yet implemented on the 
!    level of densities, only on the level of derivatives of densities.
!    Thus
!           D_N_N is fully stored with indices (nx*ny*nz,3,3,2)
!    But 
!           Der_Der_D_I_I is stored as (nx*ny*nz,7,2)
!
!===============================================================================
use compilation
use geninfo
use wavefunctions
use pairing
use derivatives 
use preconditioning 
use basis_transform
use timing

implicit none

    !---------------------------------------------------------------------------
    ! Type declaration of the various densities
$DECLARATION   
    !---------------------------------------------------------------------------
    ! Separate, manual, declaration of Div.J(r) as calculated from the spwfs
    ! for the more accurate calculation of its multipole moments
    real(KIND=dp), allocatable :: divJ(:,:)
    !---------------------------------------------------------------------------
    ! Density-mixing parameter default value.
    ! This can be set in the scfiteration namelist in the scfiteration model.
    real(KIND=dp) :: denmix = 0.75_dp
    !---------------------------------------------------------------------------
    ! Type of density mixing to perform. (Default = None)
    ! This can be set in the scfiteration namelist in the scfiteration model.
    integer       :: densitymixing = 0
    !---------------------------------------------------------------------------
    ! Previous value(s) of the density rho and s
    real(KIND=dp), allocatable :: D_I_I_hist(:,:,:)
$NTR    real(KIND=dp), allocatable :: D_I_S_hist(:,:,:,:)
    !---------------------------------------------------------------------------
    ! Charge density of the protons, possibly including the correction for the 
    ! finite size of the proton. It is stored here, as both the moments module 
    ! and the coulomb module need it, even though Coulomb depends on the moments
    ! module.
    real(KIND=dp), allocatable :: chargedensity(:,:,:)
    !---------------------------------------------------------------------------
    ! The amount of iterations to keep in memory for the density mixing and 
    ! estimation of the convergence rate
    integer           :: memory = 3
    !---------------------------------------------------------------------------
    ! As several other modules deal with the density D_I_I and its derivatives
    ! in various forms,  Hephaestos fills in here the appropriate symmetries.
    integer, parameter :: sx_rho = $SX_RHO
    integer, parameter :: sy_rho = $SY_RHO
    integer, parameter :: sz_rho = $SZ_RHO
    ! and similar for the vector spin density s, which is needed in the 
    ! preconditioning of the functionals
    integer, parameter :: sx_s(3) = (/$SX_SX,$SX_SY,$SX_SZ/)
    integer, parameter :: sy_s(3) = (/$SY_SX,$SY_SY,$SY_SZ/)
    integer, parameter :: sz_s(3) = (/$SZ_SX,$SZ_SY,$SZ_SZ/)
    
contains
 
 subroutine ConstructCanonicalBasis()
    !---------------------------------------------------------------------------
    ! Construction of the canonical basis, including all possible transformation
    ! matrices if we are doing HFB calculations in a single basis. 
    ! 
    ! a. Diagonalize the density matrix rho to obtain the canonical 
    !    transformation.
    ! b. Apply the canonical transformation, either in-place or allocating 
    !    a second set of spwfs.
    ! c. If the canonical basis was constructed in-place, we transform all
    !    relevant matrices as well. (Including the canonical transformation, 
    !    which becomes trivial.)
    !---------------------------------------------------------------------------
    integer :: ifail, wave

    call start_timer(T_den_can)

    if(.not.allocated(Canenergies)) allocate(Canenergies(nwt)) 

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! a. construct the transformation CanTransfo that brings us into the 
    !    canonical basis by diagonalizing rho
    call Canonical(rho_pairing, kappa_pairing, rho_can, kappa_can,             &
    &              cantransfo,cancuttransfo,ifail)
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! b. use this transformation to construct the physical wavefunctions
    if(efficientHFB) then
       ! Saving memory
       call transform_spwfs_inplace(hfpsi, cantransfo)
       if(allocated(momentum_updates)) then
         call transform_spwfs_inplace(momentum_updates, cantransfo)
       endif
    else
       ! Full-on transformation, we have memory to burn
       call transform_spwfs(hfpsi, canpsi, cantransfo)
    endif
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! c. Transform all relevant matrices into the new basis
    if(efficientHFB) then
      rho_pairing   = transform_mat(rho_pairing, cantransfo)
      ! Kappa transforms differently from rho, but in case of real 
      ! matrices this is largely irrelevant
      kappa_pairing = transform_mat(kappa_pairing, cantransfo)

      current_sph   = transform_mat(current_sph, cantransfo)
      do wave=1,nwt
        canenergies(wave) = current_sph(wave,wave)
      enddo 
      HFBgaps    = transform_mat(HFBgaps, cantransfo)

      HFtransfo   = transform_vec(HFtransfo  , cantransfo)
      Bogoliubov  = transform_bogo(Bogoliubov, cantransfo)

      ! The new basis IS the canonical basis, hence the canonical transformation
      ! is trivial.
      cantransfo = 0.0d0
      do wave=1,nwt
        cantransfo(wave,wave) = 1.0d0
      enddo 
    else
      canenergies = transform_diag(current_sph, cantransfo) 
    endif
    call stop_timer(T_den_can)

    ! If we are not doing HFB efficiently, we should rederive the spwfs
    if(.not. efficientHFB) call derivecan()

 end subroutine ConstructCanonicalBasis

 subroutine densit(SaveRho)
    !---------------------------------------------------------------------------
    ! Calculate all of the densities. 
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Input:
    !    SaveRho : logical. If .true., save the values of the densities on input
    !                       to their respective histories. If .false., do not
    !                       keep this information.
    !---------------------------------------------------------------------------
    integer      :: i, it, wave, wave2, B, N, si, N2, T
    integer      ::  wave_global, wave2_global
    real(KIND=dp):: weight
    logical      :: SaveRho
    real(KIND=dp), allocatable :: kappa_cut(:,:)
#if(USE_MPI>0)
    integer      :: mpi_err
#endif


    call start_timer(T_densities)

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Allocation and initialization
$INITIALIZATION

    if(.not.allocated(divJ)) then
      allocate(divJ(nx*ny*nz,4))    
    endif
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Save old density for next iteration and mixing.
    ! Note that this is only necessary at the moment for the ordinary rho
    ! density, it is the one that can make calculations unstable.
    if(.not. allocated(D_I_I_hist)) then
        allocate(D_I_I_hist(nx*ny*nz,2,memory)) ; D_I_I_hist = 0.0_dp
$NTR    allocate(D_I_S_hist(nx*ny*nz,3,2,memory)) ; D_I_S_hist = 0.0_dp
    endif   
    if(SaveRho) then
      do i=1,memory-1
          D_I_I_hist(:,:,memory-i+1) = D_I_I_hist(:,:,memory-i)
$NTR      D_I_S_hist(:,:,:,memory-i+1) = D_I_S_hist(:,:,:,memory-i)
      enddo
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
      ! Saving the input density for mixing    
      D_I_I_hist(:,:,1) = D_I_I
$NTR      D_I_S_hist(:,:,:,1) = D_I_S
    endif
    
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Zero the current density
$ZEROING

    divJ = 0.0d0
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
        weight  = rho_can(wave_global) 

        do i=1,mv
$EXPRESSION
        enddo
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
      kappa_cut = kappa_pairing
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
        ! pairing densities. The reason is that the HFTransfo multiplication 
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
    end select
    call stop_timer(T_den_pp)
    
#if(USE_MPI > 0)
   ! Sum the density over all processes
$MPIDEN
#endif
    
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! The asked for mixing+preconditioning scheme.
    call MassageDensity()
    
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

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Sum DivJ from the spwfs separately
    divJ = sum_divJ_spwf()

    call stop_timer(T_densities)

!#if(USE_MPI > 0)      
!      call stp('End of densit')
!#endif      

end subroutine densit

function sum_divJ_spwf() result(divJ)
    !---------------------------------------------------------------------------
    ! Calculate the 
    !            nabla cdot J
    ! where J is the vector component of the spin-current density J_munu. 
    ! The contribution from a single spwf is 
    !
    !  sum_{mu nu kappa} eps_munukappa
    !            rho_ii  Im [ \nabla_mu Psi_i^* \nabla_kappa \sigma_nu \Psi_ii]
    !
    ! where eps_munukappa is a Levi-Civita symbol.
    !
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! The reason this thing is calculated separately is because the sum of 
    ! spwfs is represented much more accurately on the mesh than the sum of
    ! derivatives of J on the mesh. We use this summation for the calculation
    ! of multipole moments of J, but it should in W.R.'s opinion NOT be used
    ! in the calculation of any energy for consistency reasons.
    !---------------------------------------------------------------------------
    real(KIND=dp) :: divJ(nx*ny*nz,4)
    real(KIND=dp) :: temp(nx*ny*nz,4)
    real(KIND=dp) :: weight
    integer       :: wave,it

    divJ = 0.0d0
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
    do wave=1,nwt
        ! Isospin is neutron in the first half of blocks, proton in the rest
        it = 2
        if(wave.le.nwn) it = 1
        
        ! For ordinary densities
        weight  = rho_can(wave) 

        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ! I code this with calls to the Pauli and ImagMultiplySpinor functions 
        ! to make no mistakes
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ! x y z
        temp = Pauli(DenDPsi(:,2,:,wave), 3)
        divJ(:,it) = divJ(:,it) &
        &            + weight * ImagMultiplySpinor(DendPsi(:,1,:,wave), temp) 
        ! y x z
        temp = Pauli(DenDPsi(:,1,:,wave), 3)
        divJ(:,it) = divJ(:,it) &
        &            - weight * ImagMultiplySpinor(DendPsi(:,2,:,wave), temp) 
        ! x z y 
        temp = Pauli(DenDPsi(:,3,:,wave), 2)
        divJ(:,it) = divJ(:,it) &
        &            - weight * ImagMultiplySpinor(DendPsi(:,1,:,wave), temp) 
        ! z x y 
        temp = Pauli(DenDPsi(:,1,:,wave), 2)
        divJ(:,it) = divJ(:,it) &
        &            + weight * ImagMultiplySpinor(DendPsi(:,3,:,wave), temp) 
        ! y z x 
        temp = Pauli(DenDPsi(:,3,:,wave), 1)
        divJ(:,it) = divJ(:,it) &
        &            + weight * ImagMultiplySpinor(DendPsi(:,2,:,wave), temp) 
        ! z y x 
        temp = Pauli(DenDPsi(:,2,:,wave), 1)
        divJ(:,it) = divJ(:,it) &
        &            - weight * ImagMultiplySpinor(DendPsi(:,3,:,wave), temp) 
    enddo
    ! Taking isospin combinations
    divJ(:,3) = divJ(:,1) + divJ(:,2)
    divJ(:,4) = divJ(:,1) - divJ(:,2)
  
end function sum_divJ_spwf

subroutine MassageDensity()
    !---------------------------------------------------------------------------
    ! Operate on the density before feeding it into the rest of the program.
    !---------------------------------------------------------------------------
    real(KIND=dp), target :: resid(nx*ny*nz,2)
$NTR real(KIND=dp), target :: sresid(nx*ny*nz,3,2)
    if(all(D_I_I_hist.eq.0.0)) return
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Compute the residual
    resid = D_I_I - D_I_I_hist(:,:,1)
$NTR    sresid = D_I_S - D_I_S_hist(:,:,:,1)
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Perform mixing
    select case(densitymixing) 
    case(0)
        !-----------------------------------------------------------------------
        ! Precondition the potentials instead of the densities. 
        ! So do nothing to the densities.
    case(1)
        !-----------------------------------------------------------------------
        ! Simple linear mixing at the moment.
        D_I_I = D_I_I_hist(:,:,1) + (1-denmix) * resid
$NTR    D_I_S = D_I_S_hist(:,:,:,1) + (1-denmix) * sresid
    end select
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Safeguard
    where(D_I_I.lt.1d-10) D_I_I = 0 
end subroutine MassageDensity

function couple_iso(density, iso) result(coupled)
    !---------------------------------------------------------------------------
    ! All densities are calculated and stored in proton-neutron format, but 
    ! terms in the EDF are generally calculated in isospin formalism. 
    !
    ! This function takes as argument a density (with ONLY an proton/neutron)
    ! index and produces the isoscalar or isovector combination for on-the-fly
    ! resummation and easy Hephaestos code generation.
    !
    !---------------------------------------------------------------------------
    integer, intent(in)       ::  iso
    real(KIND=dp), intent(in) ::  density(mv,2)
    real(KIND=dp)             ::  coupled(mv)
    
    select case(iso)
    case(0)
      ! Isoscalar = neutron + proton
      coupled = sum(density(:,2))
    case(1)
      ! Isovector = neutron - proton
      coupled = density(:,1) - density(:,2)
    case DEFAULT
      call stp('Invalid iso argument to couple_iso.')
    end select

end function couple_iso

subroutine Transfer_psi(psi,send_rank, calc_rank, wave)
    !---------------------------------------------------------------------------
    ! Transfer the wavefunction  
    !        denpsi(:,:,wave)   
    ! from the sending MPI_rank to a rank fit for calculations. 
    !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input : 
    !   send_rank : MPI rank storing the requested wavefunction
    !   calc_rank : MPI rank supposed to be doing calculations with the 
    !               requested wavefunction.
    !   wave      : LOCAL index of the requested spwf on the send_rank
    ! Output:
    !  psi        : the requested spwf, but only on CALC_RANK. For all 
    !               other MPI ranks, the result will be unallocated. 
    !---------------------------------------------------------------------------
    real(KIND=dp),  intent(out) :: psi(mv,4)
    integer, intent(in)         :: send_rank, calc_rank, wave
#if(USE_MPI>0)
    integer                     :: mpi_err

    if((MPI_RANK.eq. calc_rank) .AND. (send_rank.eq.calc_rank)) then
        ! nothing to send or receive
        psi   = DenPsi(:,:,wave)
    elseif(MPI_RANK.eq.calc_rank) then
        ! calc_rank receives
        call MPI_RECV(            psi, 4*mv, MPI_REAL8, send_rank, 2, &
        &                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_err)
    else if(MPI_RANK .eq. send_rank)  then
        ! send_rank sends the wavefunction
        call MPI_SEND(denpsi(:,:,wave), 4*mv, MPI_REAL8, calc_rank, 2,&
        &                                           MPI_COMM_WORLD, mpi_err)
    endif
#else 
    psi   = DenPsi(:,:,wave)
#endif

end subroutine Transfer_psi

subroutine Transfer_derpsi(derpsi,send_rank, calc_rank, wave, direction, TR)
    !---------------------------------------------------------------------------
    ! Transfer the (possibly time-reversed) derivative of a wavefunction  
    !        dendpsi(:,:,direction,wave)   
    ! from the sending MPI_rank to a rank fit for calculations. 
    !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input : 
    !   send_rank : MPI rank storing the requested wavefunction
    !   calc_rank : MPI rank supposed to be doing calculations with the 
    !               requested wavefunction.
    !   wave      : LOCAL index of the requested spwf on the send_rank
    !   direction : direction of the derivative (x/y/z = 1/2/3) of the 
    !               requested spwf
    !   TR        : logical, whether or not to apply a time-reversal operation
    !               before returning.
    ! Output:
    !  derpsi     : the requested derivative of an spwf, but only on CALC_RANK.
    !               For all other MPI ranks, the array is not changed.
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(out) :: derpsi(mv,4)
    integer, intent(in)        :: send_rank, calc_rank, wave
    integer, intent(in)        :: direction
    logical, intent(in)        :: TR
#if(USE_MPI>0)
    integer                    :: mpi_err

    if((MPI_RANK.eq. calc_rank) .AND. (send_rank.eq.calc_rank)) then
        ! nothing to send or receive
        derpsi   = DendPsi(:,:,direction,wave)
        if(TR)   derpsi = TimeReverse(derpsi)
    elseif(MPI_RANK.eq.calc_rank) then
        ! calc_rank receives
        call MPI_RECV(                     derpsi, 4*mv, MPI_REAL8,send_rank,2,&
        &                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_err)
        ! and time-reverses if needed
        if(TR)   derpsi = TimeReverse(derpsi)
    else if(MPI_RANK .eq. send_rank)  then
        ! ranki sends the wavefunction
        call MPI_SEND(DendPsi(:,:,direction,wave), 4*mv, MPI_REAL8,calc_rank,2,&
        &                                           MPI_COMM_WORLD, mpi_err)
    endif
#else 
    derpsi   = DendPsi(:,:,direction,wave)
    if(TR)   derpsi = TimeReverse(derpsi)
#endif

end subroutine Transfer_derpsi

function CompNablaMelements() result(NablaMelements)
    !---------------------------------------------------------------------------
    ! Computes the matrix elements of Nabla
    !
    !   < Psi_i | \nabla | \Psi_j >
    !
    ! In the basis from which the densities are constructed:
    !    (a) HF-basis for HF and BCS calculations
    !    (b) Canonical basis for HFB calculations
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! The current MPI parallelization strategy to calculate a given matrix 
    ! element is simple: transfer all relevant spwfs to an indicated rank.
    !    1) find out the MPI ranks that store psi_i and psi_j, rank_i and rank_j
    !    2) send both to the designated rank
    !    3) which integrates psi_i \nabla_mu psi_j 
    !    4) at this point we use the symmetry of the operator, I.e.
    !          <psi_i |\nabla_mu |psi_j> = <psi_j| \nabla_mu| psi_i>^*
    !
    ! This does not load balance very well, but at least we take the effort of
    ! load balancing the protons and neutrons. 
    !    Sequential : everything done by rank 0 
    !    Parallel   : everything done by rank 0 (neutrons)
    !                 and the last rank. These should in all realistic cases 
    !                 only store neutron/proton spwfs respectively and should
    !                 hence never communicate.
    !    5) MPI_BCAST creates the complete array for all ranks
    ! 
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Some notes:
    ! 1) These matrix elements are, in general, complex numbers!
    ! 2) Symmetries restrict them in weird ways:
    !    Same p             => <Nx> = <Ny> = <Nz> = 0
    !    Same signature     => <Nx> = <Ny> = 0
    !    Different signature=> <Nz> = 0
    !    Different isospin  => <Nx> = <Ny> = <Nz> = 0
    ! 3) We must not forget the time-reversed states when time reversal is
    !    conserved!
    ! 4) This correction is a nightmare to efficiently parallelize using MPI.
    !    The reason is that all matrix elements that are not zero involve 
    !    spwfs with different parity and often opposite signature. These spwfs
    !    are unlikely to be found on the same MPI rank => massive communication 
    !    costs. I have at the moment not made any effort to be clever. The 
    !    implementation below is naive but works and can serve as benchmark
    !    to future implementations that are more clever.
    ! 5) This routine is NOT ready for generalization beyond a CR4 (PT-broken)
    !    context.
    !
    ! Luckily, these matrix elements are not necessary in the most demanding 
    ! case with respect to scaling with the number of spwfs: nuclear pasta 
    ! calculations describe essentially infinite systems which are well-located
    ! in space.
    !---------------------------------------------------------------------------
    integer       :: i,         j, B, N, si, N2, N3, N4, wave, wave2
    integer       :: ranki, rankj, wave_global, wave2_global
    integer       :: designated_rank(2), calc_rank, T
    real(KIND=dp) :: NablaMElements(3,2,nwt,nwt), psi(mv,4)
    real(KIND=dp) :: derx(mv,4), dery(mv,4), derz(mv,4)
#if(USE_MPI>0)
    integer       :: mpi_err
#endif
    !---------------------------------------------------------------------------
    NablaMElements= 0.0_dp

    designated_rank(1) = 0
    designated_rank(2) = Ncores - 1

#if(USE_MPI>0)
    ! Verify that both designated ranks should have no communications; if they
    ! do this calculation will take very long.
    if((MPI_RANK.eq.designated_rank(1)) &
    &          .and.                    &
    &        (sum(HFBlocks(5:8)).gt. 0)) then
      call stp('This load balancing cannot be used in CompNablaMelements.')
    else if((MPI_RANK.eq.designated_rank(2)) &
    &          .and.                         &
    &            (sum(HFBlocks(1:4)).gt. 0)) then
      call stp('This load balancing cannot be used in CompNablaMelements.') 
    endif
#endif
    si = 0 
    do B=1,8,4 ! This loop is essentially over protons vs neutrons
      N = HFblocks_global(B)  ;  N2= HFBlocks_global(B+1)
      N3= HFBlocks_global(B+2);  N4= HFBlocks_global(B+3)

      calc_rank = designated_rank(mod(B,5) + 1)

      !-------------------------------------------------------------------------
      ! Re < p_z > : only non-zero matrix elements when equal signature
      !              and opposite parity (if it is conserved)
      !-------------------------------------------------------------------------
      ! Block 1 with block 3
      do i=1,N
        wave_global = si+i  ! this is the global index of the spwf
        ! .... and we have to find the MPI RANK and the index it is on
        ranki = rank_map(wave_global)
        wave  = spwf_inverse(wave_global)
        ! sending, receing, copying etc to the array psi for calc_rank only 
        call transfer_psi(psi, ranki, calc_rank, wave)

        do j=1,N3
          wave2_global = si+N+N2+j
          rankj        = rank_map(wave2_global)
          wave2        = spwf_inverse(wave2_global)

          ! sending, receing, copying etc to the array psi for calc_rank only 
          call transfer_derpsi(derz, rankj, calc_rank, wave2, 3,.false.)

          if(MPI_RANK .eq. calc_rank) then
            ! Perform the calculation with the rank designated to calculate
            NablaMElements(3,1,wave_global,wave2_global) = dv*                 &
            & sum(                  derz(:,1) * psi(:,1) + derz(:,2) * psi(:,2)&
            &                    +  derz(:,3) * psi(:,3) + derz(:,4) * psi(:,4))

            NablaMElements(3,1,wave2_global,wave_global) =&
            &                       NablaMElements(3,1,wave_global,wave2_global)
          endif
        enddo
      enddo
      
      ! Block 1 with block 1 (parity broken)
$PBROKEN   do i=1,N
$PBROKEN     wave_global = si+i  ! this is the global index of the spwf
$PBROKEN     ! .... and we have to find the MPI RANK and the index it is on
$PBROKEN     ranki = rank_map(wave_global)
$PBROKEN     wave  = spwf_inverse(wave_global)
$PBROKEN     call transfer_psi(psi, ranki, calc_rank, wave)
$PBROKEN
$PBROKEN     do j=1,N
$PBROKEN       wave2_global = si+j
$PBROKEN       rankj        = rank_map(wave2_global)
$PBROKEN       wave2        = spwf_inverse(wave2_global)
$PBROKEN
$PBROKEN       ! sending, receing, copying etc to the array psi for calc_rank only 
$PBROKEN       call transfer_derpsi(derz, rankj, calc_rank, wave2, 3,.false.)
$PBROKEN
$PBROKEN       if(MPI_RANK .eq. calc_rank) then
$PBROKEN         NablaMElements(3,1,wave_global,wave2_global) = dv*            &
$PBROKEN         & sum(             derz(:,1) * psi(:,1) + derz(:,2) * psi(:,2)&
$PBROKEN         &               +  derz(:,3) * psi(:,3) + derz(:,4) * psi(:,4))
$PBROKEN
$PBROKEN         NablaMElements(3,1,wave2_global,wave_global) = &
$PBROKEN         &                  NablaMElements(3,1,wave_global,wave2_global)
$PBROKEN        endif
$PBROKEN     enddo
$PBROKEN   enddo

      ! Block 2 with block 4
      do i=1,N2
        wave_global  = si+N+i
        ranki        = rank_map(wave_global)
        wave         = spwf_inverse(wave_global)
        ! sending, receing, copying etc to the array psi for calc_rank only 
        call transfer_psi(psi, ranki, calc_rank, wave)

        do j=1,N4
          wave2_global = si+N+N2+N3+j
          rankj        = rank_map(wave2_global)
          wave2        = spwf_inverse(wave2_global)

          ! sending, receing, copying etc to the array psi for calc_rank only 
          call transfer_derpsi(derz, rankj, calc_rank, wave2, 3, .false.)

          if(MPI_RANK.eq.calc_rank) then
          ! Re < p_z > 
            NablaMElements(3,1,wave_global,wave2_global) = dv*                 &
            & sum(                  derz(:,1) * psi(:,1) + derz(:,2) * psi(:,2)&
            &                    +  derz(:,3) * psi(:,3) + derz(:,4) * psi(:,4))

            NablaMElements(3,1,wave2_global,wave_global) = &
            &                       NablaMElements(3,1,wave_global,wave2_global)
          endif
        enddo
      enddo
      
      ! Block 2 with block 2 (parity broken)
$PBROKEN   do i=1,N2
$PBROKEN     wave_global = si+i  ! this is the global index of the spwf
$PBROKEN     ! .... and we have to find the MPI RANK and the index it is on
$PBROKEN     ranki = rank_map(wave_global)
$PBROKEN     wave  = spwf_inverse(wave_global)
$PBROKEN     call transfer_psi(psi, ranki, calc_rank, wave)
$PBROKEN 
$PBROKEN     do j=1,N2
$PBROKEN       wave2_global = si+N+j
$PBROKEN       rankj        = rank_map(wave2_global)
$PBROKEN       wave2        = spwf_inverse(wave2_global)
$PBROKEN
$PBROKEN       call transfer_derpsi(derz, rankj, calc_rank, wave2, 3, .false.)
$PBROKEN       if(MPI_RANK.eq.calc_rank) then
$PBROKEN         ! Re < p_z > 
$PBROKEN         NablaMElements(3,1,wave_global,wave2_global) = dv*            &
$PBROKEN         & sum(             derz(:,1) * psi(:,1) + derz(:,2) * psi(:,2)&
$PBROKEN         &               +  derz(:,3) * psi(:,3) + derz(:,4) * psi(:,4))
$PBROKEN
$PBROKEN         NablaMElements(3,1,wave2_global,wave_global) = &
$PBROKEN         &                  NablaMElements(3,1,wave_global,wave2_global)
$PBROKEN       endif
$PBROKEN     enddo
$PBROKEN   enddo

      !-------------------------------------------------------------------------
      ! Re < p_x >, Im <p_y> 

      do i=1,N
        wave_global  = si+i
        ranki        = rank_map(wave_global)
        wave         = spwf_inverse(wave_global)        
        ! sending, receing, copying etc to the array psi for calc_rank only 
        call transfer_psi(psi, ranki, calc_rank, wave)
           ! Block 1 with block 4  (T-broken)    
$NTR       do j=1, N4
$NTR          wave2_global = si+N+N2+N3+j
$NTR          Derx  =  DendPsi(:,1,:,wave2)
$NTR          Dery  =  DendPsi(:,2,:,wave2)
           ! Block 1 with block 3  (T-conserved)    
$TR        do j=1, N3
$TR           wave2_global = si+N+j
$TR           rankj        = rank_map(wave2_global)
$TR           wave2        = spwf_inverse(wave2_global)
$TR           ! There is a timereversal operation hidden behind the .true. in
$TR           ! the lines below!
$TR           call transfer_derpsi(derx, rankj, calc_rank, wave2, 1, .true.)
$TR           call transfer_derpsi(dery, rankj, calc_rank, wave2, 2, .true.)

          if(MPI_RANK.eq.calc_rank) then

            NablaMElements(1,1,wave_global,wave2_global) = dv*                 &
            & sum(                 derx(:,1) * psi(:,1) + derx(:,2) * psi(:,2) &
            &                   +  derx(:,3) * psi(:,3) + derx(:,4) * psi(:,4))
  
            NablaMElements(2,2,wave_global,wave2_global) = dv*                 &
            & sum(                  dery(:,2) * psi(:,1) - dery(:,1) * psi(:,2)&
            &                    -  dery(:,3) * psi(:,4) + dery(:,4) * psi(:,3))

            ! Attention to the signs here
            NablaMElements(1,1,wave2_global,wave_global) = &
            &                     - NablaMElements(1,1,wave_global,wave2_global)
            NablaMElements(2,2,wave2_global,wave_global) = &
            &                       NablaMElements(2,2,wave_global,wave2_global)
         endif
        enddo
      enddo
      
$PBROKEN   do i=1,N
$PBROKEN     wave_global  = si+i
$PBROKEN     ranki        = rank_map(wave_global)
$PBROKEN     wave         = spwf_inverse(wave_global)        
$PBROKEN     ! sending, receing, copying etc to the array psi for calc_rank only 
$PBROKEN     call transfer_psi(psi, ranki, calc_rank, wave)
                ! Block 1 with block 2  (T-broken, P-broken)    
$PBROKEN $NTR do j=1, N2
$PBROKEN $NTR   wave2_global = si+N+j
$PBROKEN $NTR   rankj        = rank_map(wave2_global)
$PBROKEN $NTR   wave2        = spwf_inverse(wave2_global)
$PBROKEN $NTR   call transfer_derpsi(derx, rankj, calc_rank, wave2, 1, .false.)
$PBROKEN $NTR   call transfer_derpsi(dery, rankj, calc_rank, wave2, 2, .false.)
                ! Block 1 with block 1  (T-conserved, P-broken)    
$PBROKEN $TR  do j=1, N
$PBROKEN $TR    wave2_global = si+j
$PBROKEN $TR    rankj        = rank_map(wave2_global)
$PBROKEN $TR    wave2        = spwf_inverse(wave2_global)
$PBROKEN $TR    call transfer_derpsi(derx, rankj, calc_rank, wave2, 1, .true.)
$PBROKEN $TR    call transfer_derpsi(dery, rankj, calc_rank, wave2, 2, .true.)

$PBROKEN        if(MPI_RANK .eq. calc_rank) then
$PBROKEN          NablaMElements(1,1,wave_global,wave2_global) = dv*           &
$PBROKEN          & sum(         derx(:,1) * psi(:,1) + derx(:,2) * psi(:,2)   &
$PBROKEN          &           +  derx(:,3) * psi(:,3) + derx(:,4) * psi(:,4))

$PBROKEN          NablaMElements(2,2,wave_global,wave2_global) = dv*           &
$PBROKEN          & sum(         dery(:,2) * psi(:,1) - dery(:,1) * psi(:,2)   &
$PBROKEN          &           -  dery(:,3) * psi(:,4) + dery(:,4) * psi(:,3))

$PBROKEN          NablaMElements(1,1,wave2_global,wave_global) =               &
$PBROKEN          &               - NablaMElements(1,1,wave_global,wave2_global)
$PBROKEN          NablaMElements(2,2,wave2_global,wave_global) =               &
$PBROKEN          &                 NablaMElements(2,2,wave_global,wave2_global)
$PBROKEN        endif
$PBROKEN     enddo
$PBROKEN   enddo

      do i=1,N2
        wave_global  = si+N+i
        ranki        = rank_map(wave_global)
        wave         = spwf_inverse(wave_global)
        ! sending, receing, copying etc to the array psi for calc_rank only 
        call transfer_psi(psi, ranki, calc_rank, wave)
           ! Block 2 with block 3  (T-broken)    
$NTR       do j=1, N3
$NTR          wave2 = si+N+N2+j
$NTR          Derx  =  DendPsi(:,1,:,wave2)
$NTR          Dery  =  DendPsi(:,2,:,wave2)
           ! Block 2 with block 4  (T-conserved)    
$TR        do j=1, N4
$TR           wave2_global = si+N+N2+N3+j
$TR           rankj        = rank_map(wave2_global)
$TR           wave2        = spwf_inverse(wave2_global)
$TR           ! There is a timereversal operation hidden behind the .true. in
$TR           ! the lines below!
$TR           call transfer_derpsi(derx, rankj, calc_rank, wave2, 1, .true.)
$TR           call transfer_derpsi(dery, rankj, calc_rank, wave2, 2, .true.)

          if(MPI_RANK.eq.calc_rank) then
            NablaMElements(1,1,wave_global,wave2_global) = dv*                 &
            & sum(                 derx(:,1) * psi(:,1) + derx(:,2) * psi(:,2) &
            &                   +  derx(:,3) * psi(:,3) + derx(:,4) * psi(:,4))

            NablaMElements(2,2,wave_global,wave2_global) = dv*                 &
            & sum(                  dery(:,2) * psi(:,1) - dery(:,1) * psi(:,2)&
            &                    -  dery(:,3) * psi(:,4) + dery(:,4) * psi(:,3))
            ! Attention to the signs below
            NablaMElements(1,1,wave2_global,wave_global) = &
            &                     - NablaMElements(1,1,wave_global,wave2_global)
            NablaMElements(2,2,wave2_global,wave_global) = &
            &                       NablaMElements(2,2,wave_global,wave2_global)
          endif
        enddo
      enddo

$PBROKEN   do i=1,N2
$PBROKEN      wave_global  = si+N+i
$PBROKEN      ranki        = rank_map(wave_global)
$PBROKEN      wave         = spwf_inverse(wave_global)
              ! Block 2 with block 1 (P-broken, T-broken)
              ! This one is likely superfluous, as we have this 
              ! combination already once above
$PBROKEN $NTR do j=1, N
$PBROKEN $NTR    wave2_global = si+j
$PBROKEN $NTR    rankj        = rank_map(wave2_global)
$PBROKEN $NTR    wave2        = spwf_inverse(wave2_global)
$PBROKEN $NTR    call transfer_derpsi(derx, rankj, calc_rank, wave2, 1, .false.)
$PBROKEN $NTR    call transfer_derpsi(dery, rankj, calc_rank, wave2, 2, .false.)

              ! Block 2 with block 2 (P-broken, T-conserved)
$PBROKEN $TR   do j=1, N2
$PBROKEN $TR    wave2_global = si+N+j
$PBROKEN $TR    rankj        = rank_map(wave2_global)
$PBROKEN $TR    wave2        = spwf_inverse(wave2_global)
$PBROKEN $TR    call transfer_derpsi(derx, rankj, calc_rank, wave2, 1, .false.)
$PBROKEN $TR    call transfer_derpsi(dery, rankj, calc_rank, wave2, 2, .false.)

$PBROKEN        if(MPI_RANK.eq.calc_rank) then
$PBROKEN          NablaMElements(1,1,wave_global,wave2_global) = dv*           &
$PBROKEN          & sum(         derx(:,1) * psi(:,1) + derx(:,2) * psi(:,2)   &
$PBROKEN          &           +  derx(:,3) * psi(:,3) + derx(:,4) * psi(:,4))

$PBROKEN          NablaMElements(2,2,wave_global,wave2_global) = dv*           &
$PBROKEN          & sum(          dery(:,2) * psi(:,1) - dery(:,1) * psi(:,2)  &
$PBROKEN          &            -  dery(:,3) * psi(:,4) + dery(:,4) * psi(:,3))

$PBROKEN          NablaMElements(1,1,wave2_global,wave_global) = &
$PBROKEN          &               - NablaMElements(1,1,wave_global,wave2_global)
$PBROKEN          NablaMElements(2,2,wave2_global,wave_global) = &
$PBROKEN          &                 NablaMElements(2,2,wave_global,wave2_global)
$PBROKEN        endif
$PBROKEN      enddo
$PBROKEN    enddo
      T = N + N2 + N3 + N4 
      si = si + T 
    enddo

#if(USE_MPI>0)
    ! Broadcast all calculated matrix elements to all ranks. 
    ! Note that this is an explicit, separate, loop such that the parallel 
    ! structure of the loop above does not get disturbed
    si = 0 
    do B=1,8,4 ! This loop is essentially over protons vs neutrons
      N = HFblocks_global(B)  ;  N2= HFBlocks_global(B+1)
      N3= HFBlocks_global(B+2);  N4= HFBlocks_global(B+3)
      T = N + N2 + N3 + N4 
      calc_rank = designated_rank(mod(B,5) + 1)

      ! Broadcast all matrix elements from the calculating rank t all the rest
      call MPI_BCAST(NablaMElements(:,:,si+1:si+T, si+1:si+T), 6*T**2, &
      &              MPI_REAL8, calc_rank, MPI_COMM_WORLD, mpi_err)
      si = si + T 
    enddo
#endif
    !---------------------------------------------------------------------------
end function CompNablaMelements

subroutine print_boxsize_check()
  !-----------------------------------------------------------------------------
  ! Print the maximum values of the densities D_I_I and DP_I_I at the edges of 
  ! the box, to check that we are not dealing with a too small box.
  !-----------------------------------------------------------------------------

  real(KIND=dp), pointer :: rho3D(:,:,:,:)
  real(KIND=dp), pointer :: rhoP_3D(:,:,:,:)

  1 format ('------------------------- Box Size Check ---------------------------')
 11 format (' Normal  density        rho  =  D_I_I')
 12 format (' Pairing density \tilde{rho} = DP_I_I')
  2 format (' Xmax = (nx+0.5)dx = ', f10.3, ' fm,  max(rho(X=Xmax)) = ', es12.3 )
! 21 format (' Xmin =-(nx+0.5)dx = ', f6.3, 'fm,  max(rho(X=Xmin)) = ', e12.3 )

  3 format (' Ymax = (ny+0.5)dx = ', f10.3, ' fm,  max(rho(Y=Ymax)) = ', es12.3 )
! 31 format (' Ymin =-(ny+0.5)dx = ', f5.3, 'fm,  max(rho(Y=Ymax)) = ', e12.3 )

  4 format (' Zmax = (nz+0.5)dx = ', f10.3, ' fm,  max(rho(Z=Zmax)) = ', es12.3 )
$PBROKEN 41 format (' Zmin =-(nz+0.5)dx = ', f10.3, ' fm,  max(rho(Z=Zmin)) = ', es12.3 )
  
  rho3D(1:nx,1:ny,1:nz,1:2)   => D_I_I
  
  print 1
  print 11
  print 2, meshX(nx) , maxval(sum(rho3D(nx,:,:,:),3))
  print 3 , meshY(ny), maxval(sum(rho3D(:,ny,:,:),3))
  print 4 , meshZ(nz), maxval(sum(rho3D(:,:,nz,:),3))
$PBROKEN  print 41, meshZ(1) , maxval(sum(rho3D(:,:,1,:),3))

  if(pairingtype .ne. 0) then
    rhoP_3D(1:nx,1:ny,1:nz,1:2) => DP_I_I

    print 12
    print 2, meshX(nx) , maxval(abs(sum(rhoP_3D(nx,:,:,:),3)))
    print 3 , meshY(ny), maxval(abs(sum(rhoP_3D(:,ny,:,:),3)))
    print 4 , meshZ(nz), maxval(abs(sum(rhoP_3D(:,:,nz,:),3)))
$PBROKEN  print 41, meshZ(1) , maxval(abs(sum(rhoP_3D(:,:,1,:),3)))
  endif
  
end subroutine print_boxsize_check

subroutine clean_densities
$CLEANING
    if(allocated(D_I_I_hist))    deallocate(D_I_I_hist)
    if(allocated(chargedensity)) deallocate(chargedensity)
    nullify(DenPsi)
    nullify(DendPsi)
    nullify(DenddPsi)
    nullify(DendddPsi)

end subroutine clean_densities

end module densities
