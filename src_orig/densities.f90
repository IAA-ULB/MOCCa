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
! DISABLED KEYWORD, still present in Hephaestos
! CLEANING        : [WAY too long to include here]
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
use vectors
use wavefunctions
use pairing
use derivatives 
use preconditioning 
use basis_transform
use timing

implicit none

    interface operator (+)
      !Overloading "+" to be used to add density vectors.
      module procedure Add_densityvector
    end interface
    
    interface operator (*)
      !Overloading "*" to be used to multiply density vectors with scalars
      module procedure multiply_densityvector
    end interface
    !---------------------------------------------------------------------------
    ! Different density-vectors that can be used throughout the code
    type(DensityVector), target :: Density
    type(DensityVector), target :: Density_out
    type(DensityVector), target, allocatable :: DensityHistory(:)
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
 
 subroutine save_density_history(R)
  !-----------------------------------------------------------------------------
  ! Update the history of the mean-field densities with R.
  !-----------------------------------------------------------------------------
  integer :: i
  type(DensityVector), intent(in) :: R
  
  if(.not.allocated(DensityHistory)) allocate(DensityHistory(memory))
    
  do i=1,memory-1
      DensityHistory(memory-i+1) = DensityHistory(memory-i)
  enddo
  DensityHistory(1) = R

end subroutine save_density_history

function Add_densityvector(R1, R2) result(R)
  !-----------------------------------------------------------------------------
  ! Add two density vectors together
  !-----------------------------------------------------------------------------
  type(DensityVector), intent(in) :: R1, R2
  type(DensityVector)             :: R
  real(KIND=dp) :: stor

$INITIALIZATION
$ADD

  allocate(R%chargedensity(nx,ny,nz))
  R%chargedensity = R1%chargedensity + R2%chargedensity

end function Add_densityvector

function multiply_densityvector(a, R1) result(R)
  !-----------------------------------------------------------------------------
  ! Add multiply a density vector with a real scalar.
  !-----------------------------------------------------------------------------
  type(DensityVector), intent(in) :: R1
  real(KIND=dp), intent(in)       :: a
  type(DensityVector)             :: R
  real(KIND=dp) :: stor

$INITIALIZATION
$MULTIPLY

  allocate(R%chargedensity(nx,ny,nz))
  R%chargedensity = a * R1%chargedensity 

end function multiply_densityvector
 
subroutine construct_canonical_basis(rho, kappa, rho_c, kappa_c)
  !-----------------------------------------------------------------------------
  ! High-level routine that constructs the canonical basis of spwfs.
  ! Does nothing for a HF or BCS calculation.
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Input :
  !  rho   : matrix elements of the density in the basis spanned by the spwfs 
  !          currently in memory
  !  kappa : matrix elements of the anomalous density in spanned by the spwfs 
  !          currently in memory
  ! Output:
  !  rho_can   : diagonal matrix elements of the density in the canonical basis
  !  kappa_can : "diagonal" matrix elements of kappa in the canonical basis
  !
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Side effects:
  !  *) The spwfs of canonical basis get constructed. Depending on whether 
  !     efficientHFB is active, these get stored in CanPsi or overwrite the 
  !     existing HFPsi array.
  !  *) Cantransfo is saved; this is the transformation matrix between 
  !     the HF-basis and the canonical basis. It is a trivial transformation if
  !     efficientHFB is active.
  !-----------------------------------------------------------------------------
  real(KIND=dp), intent(inout) :: rho(:,:), kappa(:,:)
  real(KIND=dp), intent(out)   :: rho_c(:), kappa_c(:)
  integer :: ifail, wave

  if(pairingtype.ne.2) return

  call start_timer(T_den_can)

    if(.not.allocated(Canenergies)) allocate(Canenergies(nwt)) 

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! a. construct the transformation CanTransfo that brings us into the 
    !    canonical basis by diagonalizing rho
    call Canonical(rho, kappa, rho_c, kappa_c, cantransfo,cancuttransfo,ifail)
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
      rho_pairing   = transform_mat(rho, cantransfo)
      ! Kappa transforms differently from rho, but in case of real 
      ! matrices this is largely irrelevant
      kappa = transform_mat(kappa, cantransfo)

      sphamil = transform_mat(sphamil, cantransfo)
      do wave=1,nwt
        canenergies(wave) = sphamil(wave,wave)
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
      canenergies = transform_diag(sphamil, cantransfo) 
    endif
    call stop_timer(T_den_can)

    ! If we are not doing HFB efficiently, we should rederive the spwfs
    if((.not. efficientHFB) .and. store_derivatives) call derivecan()

 end subroutine construct_canonical_basis
 
function densit(rho, kappa) result(R)
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
    external construct_charge_density
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
          ! We recalculate derivatives on the fly for this spwf
          do k=1,4
           call Derive_tot(denPsi(:,k,wave), sx(k,wave), sy(k,wave), sz(k,wave),&
           &                                           dendPsi(:,:,k,1),        &
           &                                           denddPsi(:,:,k,1))
          enddo
          der_index = 1
        endif
        !---------------------------------------------------------------------------

        do i=1,mv
$EXPRESSION
$DERIVATION_SUM_SPWF_PH
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

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Sum DivJ from the spwfs separately
    call sum_divJ_spwf(R)
    
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Construct the charge density
    call constructchargedensity(R)
    
    call stop_timer(T_densities)

end function densit

 subroutine ConstructChargeDensity(R)
    !---------------------------------------------------------------------------
    ! Construct the charge density from the proton and neutron densities,
    ! using various effective forms
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! TODO: document what this routine does precisely
    !---------------------------------------------------------------------------
    use Folding

    type(DensityVector),intent(inout) :: R
    real(KIND=dp)              :: temp(nx,ny,nz)
    integer                    :: i,j,k
#if(PASTA==1)
    real(KIND=dp)              :: rho_el, volume
#endif

    call start_timer(T_chargedensity)

    ! Deallocation such that rho_charge does not have the wrong dimensions
    if(allocated(R%chargedensity))      deallocate(R%chargedensity)
    if(.not.allocated(R%chargedensity)) allocate(R%chargedensity(nx,ny,nz))
    !---------------------------------------------------------------------------
    ! If we account for the finite extent of the charge of the nucleus, then
    ! we need to fold densities and potentials with gaussians. This sets up the
    ! required matrices.
    !
    ! Note: this little piece of code is duplicated, since in different
    !       runmodes of the code different Coulomb routines get called in
    !       different order; this makes sure we get no segfaults.
    !---------------------------------------------------------------------------
    if(any(protonsize .ne. 0.0_dp) .or. any(neutronsize.ne.0.0_dp)) then
      if(.not.allocated(Gaussx)) then
          allocate(Gaussx(nx,nx,2,2), Gaussy(ny,ny,2,2), Gaussz(nz,nz,2,2))
          Gaussx = 0.0 ;  Gaussy = 0.0 ; Gaussz = 0.0
      endif
      call ConstructFoldingMatrices(Gaussx,Gaussy,Gaussz,sx_rho, sy_rho, sz_rho)
    endif

    !---------------------------------------------------------------------------
    ! Proton contributions to the charge density.
    ! We start from the proton point density
    do k=1,nz
      do j=1,ny
        do i=1,nx
            temp(i,j,k) = R%D_I_I(meshindex(i,j,k),2)
        enddo
      enddo
    enddo

    if(protonsize(1).gt.0.0) then
        ! Fold the source with a Gaussian
        R%chargedensity = &
        & FoldGaussian(temp, GaussX(:,:,1,2), GaussY(:,:,1,2), GaussZ(:,:,1,2),&
        &                                                            nx, ny, nz)
    endif
    if(protonsize(2).gt.0.0) then
        ! Fold the source with another Gaussian, this time with minus sign.
        R%chargedensity = R%chargedensity + &
        & FoldGaussian(temp, GaussX(:,:,2,2), GaussY(:,:,2,2), GaussZ(:,:,2,2),&
        &                                                            nx, ny, nz)
    endif

    if(all(protonsize.eq.0.0)) then
        R%chargedensity = temp
    endif
    !---------------------------------------------------------------------------
    ! Neutron contributions to the charge density.
    if(any(neutronsize .gt. 0.0d0)) then 
      do k=1,nz
        do j=1,ny
          do i=1,nx
             temp(i,j,k) = R%D_I_I(meshindex(i,j,k),1)
          enddo
        enddo
      enddo

      if(neutronsize(1).gt.0.0) then
          ! Fold the source with a Gaussian
          R%chargedensity = R%chargedensity + &
          & FoldGaussian(temp, GaussX(:,:,1,1), GaussY(:,:,1,1), GaussZ(:,:,1,1),&
          &                                                            nx, ny, nz)
      endif
      if(neutronsize(2).gt.0.0) then
          ! Fold the source with a Gaussian, minus sign this time
          R%chargedensity = R%chargedensity - &
          & FoldGaussian(temp, GaussX(:,:,2,1), GaussY(:,:,2,1), GaussZ(:,:,2,1),&
          &                                                            nx, ny, nz)
      endif
    endif
    !---------------------------------------------------------------------------
    ! When performing simulations for nuclear pasta, one assumes the entire 
    ! volume is charge neutral: a constant background of electrons floods the 
    ! entire simulation volume. We subtract this backrgound here.
#if(PASTA==1)
    volume=nx*ny*nz*dv   ! simplification by WR: the physical volume simulated
                         ! can just be gotten by the volume element...
    rho_el=protons/volume
    !print *, protons,sum(rho_charge)*dv
    R%chargedensity = R%chargedensity -rho_el
#endif
    call stop_timer(T_chargedensity)
 end subroutine ConstructChargeDensity

subroutine sum_divJ_spwf(R)
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
    ! Output:
    !   R  : densityvector whose component divJ will be calculated
    !
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! The reason this thing is calculated separately is because the sum of 
    ! spwfs is represented much more accurately on the mesh than the sum of
    ! derivatives of J on the mesh. We use this summation for the calculation
    ! of multipole moments of J, but it should in W.R.'s opinion NOT be used
    ! in the calculation of any energy for consistency reasons.
    !---------------------------------------------------------------------------
    type(DensityVector), intent(inout) :: R
    !real(KIND=dp) :: temp(nx*ny*nz,4)
    !real(KIND=dp) :: weight
    !integer       :: wave,it

    R%divJ = 0.0d0
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
!    do wave=1,nwt
!        ! Isospin is neutron in the first half of blocks, proton in the rest
!        it = 2
!        if(wave.le.nwn) it = 1

!        ! For ordinary densities
!        weight  = rho_can(wave)

!        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!        ! I code this with calls to the Pauli and ImagMultiplySpinor functions 
!        ! to make no mistakes
!        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!        ! x y z
!        temp = Pauli(DenDPsi(:,2,:,wave), 3)
!        divJ(:,it) = divJ(:,it) &
!        &            + weight * ImagMultiplySpinor(DendPsi(:,1,:,wave), temp) 
!        ! y x z
!        temp = Pauli(DenDPsi(:,1,:,wave), 3)
!        divJ(:,it) = divJ(:,it) &
!        &            - weight * ImagMultiplySpinor(DendPsi(:,2,:,wave), temp) 
!        ! x z y 
!        temp = Pauli(DenDPsi(:,3,:,wave), 2)
!        divJ(:,it) = divJ(:,it) &
!        &            - weight * ImagMultiplySpinor(DendPsi(:,1,:,wave), temp) 
!        ! z x y 
!        temp = Pauli(DenDPsi(:,1,:,wave), 2)
!        divJ(:,it) = divJ(:,it) &
!        &            + weight * ImagMultiplySpinor(DendPsi(:,3,:,wave), temp) 
!        ! y z x 
!        temp = Pauli(DenDPsi(:,3,:,wave), 1)
!        divJ(:,it) = divJ(:,it) &
!        &            + weight * ImagMultiplySpinor(DendPsi(:,2,:,wave), temp) 
!        ! z y x 
!        temp = Pauli(DenDPsi(:,2,:,wave), 1)
!        divJ(:,it) = divJ(:,it) &
!        &            - weight * ImagMultiplySpinor(DendPsi(:,3,:,wave), temp) 
!    enddo
    ! Taking isospin combinations
    R%divJ(:,3) = R%divJ(:,1) + R%divJ(:,2)
    R%divJ(:,4) = R%divJ(:,1) - R%divJ(:,2)
  
end subroutine sum_divJ_spwf

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
    !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! TODO: this routine can be refactored by
    !       - adding a bunch of small functions for repeated instructions!
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
    integer       :: i, j, B, N, si, N2, N3, N4, wave, wave2
    integer       :: ranki, rankj, wave_global, wave2_global
    integer       :: designated_rank(2), calc_rank, T
    real(KIND=dp) :: NablaMElements(3,2,nwt,nwt), psi(mv,4)
    real(KIND=dp) :: derx(mv,4), dery(mv,4), derz(mv,4)
    logical       :: hidden_TR = .false.
#if(USE_MPI>0)
    integer       :: mpi_err
#endif

    call start_timer(T_NablaMElements)
    !---------------------------------------------------------------------------
    NablaMElements= 0.0_dp

    designated_rank(1) = 0
    designated_rank(2) = NPROCS - 1

    if(NPROCS .gt. 1) then
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
    endif

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
      do j=1,N3
        wave2_global = si+N+N2+j
        rankj        = rank_map(wave2_global)
        wave2        = spwf_inverse(wave2_global)

        ! sending, receing, copying etc to the array psi for calc_rank only
#if(USE_MPI>0)
        call transfer_derpsi(derz, wave2, 3,'DEN',.false., rankj, calc_rank)
#else
        call transfer_derpsi(derz, wave2, 3,'DEN',.false.)
#endif
        do i=1,N
          wave_global = si+i  ! this is the global index of the spwf
          ! .... and we have to find the MPI RANK and the index it is on
          ranki = rank_map(wave_global)
          wave  = spwf_inverse(wave_global)
          ! sending, receing, copying etc to the array psi for calc_rank only
#if(USE_MPI>0)
          call transfer_psi(psi, wave, 'DEN',ranki, calc_rank )
#else
          call transfer_psi(psi, wave, 'DEN')
#endif
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
$PBROKEN   do j=1,N
$PBROKEN     wave2_global = si+j
$PBROKEN     rankj        = rank_map(wave2_global)
$PBROKEN     wave2        = spwf_inverse(wave2_global)
$PBROKEN
$PBROKEN     ! sending, receing, copying etc to the array psi for calc_rank only
#if(USE_MPI>0)
$PBROKEN     call transfer_derpsi(derz,wave2, 3,'DEN',.false.,rankj,calc_rank)
#else
$PBROKEN     call transfer_derpsi(derz,wave2, 3,'DEN',.false.)
#endif
$PBROKEN     do i=1,N
$PBROKEN       wave_global = si+i  ! this is the global index of the spwf
$PBROKEN       ! .... and we have to find the MPI RANK and the index it is on
$PBROKEN       ranki = rank_map(wave_global)
$PBROKEN       wave  = spwf_inverse(wave_global)
#if(USE_MPI>0)
$PBROKEN       call transfer_psi(psi, wave, 'DEN', ranki, calc_rank)
#else
$PBROKEN       call transfer_psi(psi, wave, 'DEN')
#endif
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
      do j=1,N4
        wave2_global = si+N+N2+N3+j
        rankj        = rank_map(wave2_global)
        wave2        = spwf_inverse(wave2_global)

        ! sending, receing, copying etc to the array psi for calc_rank only
#if(USE_MPI>0)
        call transfer_derpsi(derz, wave2, 3, 'DEN',.false., rankj, calc_rank)
#else
        call transfer_derpsi(derz, wave2, 3, 'DEN',.false.)
#endif
        do i=1,N2
          wave_global  = si+N+i
          ranki        = rank_map(wave_global)
          wave         = spwf_inverse(wave_global)
          ! sending, receing, copying etc to the array psi for calc_rank only
#if(USE_MPI>0)
          call transfer_psi(psi, wave, 'DEN', ranki, calc_rank)
#else
          call transfer_psi(psi, wave, 'DEN')
#endif
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
$PBROKEN   do j=1,N2
$PBROKEN     wave2_global = si+N+j
$PBROKEN     rankj        = rank_map(wave2_global)
$PBROKEN     wave2        = spwf_inverse(wave2_global)
$PBROKEN
#if(USE_MPI>0)
$PBROKEN     call transfer_derpsi(derz,wave2,3,'DEN',.false.,rankj,calc_rank)
#else
$PBROKEN     call transfer_derpsi(derz,wave2,3,'DEN',.false.)
#endif
$PBROKEN     do i=1,N2
$PBROKEN       wave_global = si+N+i  ! this is the global index of the spwf
$PBROKEN       ! .... and we have to find the MPI RANK and the index it is on
$PBROKEN       ranki = rank_map(wave_global)
$PBROKEN       wave  = spwf_inverse(wave_global)
#if(USE_MPI>0)
$PBROKEN       call transfer_psi(psi, wave, 'DEN', ranki, calc_rank)
#else
$PBROKEN       call transfer_psi(psi, wave, 'DEN')
#endif
$PBROKEN 
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

      !----------------------------------------------------------------------------
      ! Re < p_x >, Im <p_y> 
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! hidden_TR = .true. => Include a time-reversal operation in the calculation
      !                       of matrix elements of these operators
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
$NTR  hidden_TR = .false.
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! hidden_TR = .true. => Do not include a time-reversal operation in the
      !                       calculation of matrix elements
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
$TR   hidden_TR = .true.
      !
      ! Note: I define these logicals here in order to cut down on the preprocessor
      !       directives needed just below.
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      ! Block 1 with block 4  (T-broken)
$NTR  do j=1, N4
$NTR     wave2_global = si+N+N2+N3+j
      ! Block 1 with block 3  (T-conserved)
$TR   do j=1, N3
$TR     wave2_global = si+N+j
        rankj        = rank_map(wave2_global)
        wave2        = spwf_inverse(wave2_global)
#if(USE_MPI>0)
        call transfer_derpsi(derx, wave2, 1,'DEN',hidden_TR,rankj,calc_rank)
        call transfer_derpsi(dery, wave2, 2,'DEN',hidden_TR,rankj,calc_rank)
#else
        call transfer_derpsi(derx, wave2, 1,'DEN',hidden_TR)
        call transfer_derpsi(dery, wave2, 2,'DEN',hidden_TR)
#endif
        do i=1,N
          wave_global  = si+i
          ranki        = rank_map(wave_global)
          wave         = spwf_inverse(wave_global)
          ! sending, receing, copying etc to the array psi for calc_rank only
#if(USE_MPI>0)
          call transfer_psi(psi, wave, 'DEN', ranki, calc_rank)
#else
          call transfer_psi(psi, wave, 'DEN')
#endif
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
      ! Block 1 with block 2  (T-broken, P-broken)
$PBROKEN $NTR do j=1, N2
$PBROKEN $NTR   wave2_global = si+N+j
$PBROKEN $NTR   rankj        = rank_map(wave2_global)
$PBROKEN $NTR   wave2        = spwf_inverse(wave2_global)
#if(USE_MPI>0)
$PBROKEN $NTR   call transfer_derpsi(derx,wave2,1,'DEN',hidden_TR,rankj,calc_rank)
$PBROKEN $NTR   call transfer_derpsi(dery,wave2,2,'DEN',hidden_TR,rankj,calc_rank)
#else
$PBROKEN $NTR   call transfer_derpsi(derx,wave2,1,'DEN',hidden_TR)
$PBROKEN $NTR   call transfer_derpsi(dery,wave2,2,'DEN',hidden_TR)
#endif

              ! Block 1 with block 1  (T-conserved, P-broken)
$PBROKEN $TR  do j=1, N
$PBROKEN $TR    wave2_global = si+j
$PBROKEN $TR    rankj        = rank_map(wave2_global)
$PBROKEN $TR    wave2        = spwf_inverse(wave2_global)
#if(USE_MPI>0)
$PBROKEN $TR    call transfer_derpsi(derx,wave2,1,'DEN',hidden_TR,rankj,calc_rank)
$PBROKEN $TR    call transfer_derpsi(dery,wave2,2,'DEN',hidden_TR,rankj,calc_rank)
#else
$PBROKEN $TR    call transfer_derpsi(derx,wave2,1,'DEN',hidden_TR)
$PBROKEN $TR    call transfer_derpsi(dery,wave2,2,'DEN',hidden_TR)
#endif
$PBROKEN        do i=1,N
$PBROKEN          wave_global  = si+i
$PBROKEN          ranki        = rank_map(wave_global)
$PBROKEN          wave         = spwf_inverse(wave_global)
$PBROKEN          ! sending, receing, copying etc to the array psi for calc_rank only
#if(USE_MPI>0)
$PBROKEN          call transfer_psi(psi, wave, 'DEN', ranki, calc_rank)
#else
$PBROKEN          call transfer_psi(psi, wave, 'DEN')
#endif

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

     ! Block 2 with block 3  (T-broken)
$NTR do j=1, N3
$NTR    wave2_global = si+N+N2+j
     ! Block 2 with block 4  (T-conserved)
$TR  do j=1, N4
$TR     wave2_global = si+N+N2+N3+j
        rankj        = rank_map(wave2_global)
        wave2        = spwf_inverse(wave2_global)
#if(USE_MPI>0)
        call transfer_derpsi(derx, wave2, 1,'DEN',hidden_TR,rankj,calc_rank)
        call transfer_derpsi(dery, wave2, 2,'DEN',hidden_TR,rankj,calc_rank)
#else
        call transfer_derpsi(derx, wave2, 1,'DEN',hidden_TR)
        call transfer_derpsi(dery, wave2, 2,'DEN',hidden_TR)
#endif
        do i=1,N2
          wave_global  = si+N+i
          ranki        = rank_map(wave_global)
          wave         = spwf_inverse(wave_global)
          ! sending, receing, copying etc to the array psi for calc_rank only
#if(USE_MPI>0)
          call transfer_psi(psi, wave, 'DEN', ranki, calc_rank)
#else
          call transfer_psi(psi, wave, 'DEN')
#endif
          if(MPI_RANK.eq.calc_rank) then
            NablaMElements(1,1,wave_global,wave2_global) = dv*                 &
            & sum(                 derx(:,1) * psi(:,1) + derx(:,2) * psi(:,2) &
            &                   +  derx(:,3) * psi(:,3) + derx(:,4) * psi(:,4))!

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

              ! Block 2 with block 2 (P-broken, T-conserved)
              ! Note: there is no "time-reversal broken" version of the following
              ! code block since "block 1 with 2" is already taken into account
              ! above..
$PBROKEN $TR  do j=1, N2
$PBROKEN $TR   wave2_global = si+N+j
$PBROKEN $TR   rankj        = rank_map(wave2_global)
$PBROKEN $TR   wave2        = spwf_inverse(wave2_global)
#if(USE_MPI>0)
$PBROKEN $TR   call transfer_derpsi(derx,wave2,1,'DEN',hidden_TR,rankj,calc_rank)
$PBROKEN $TR   call transfer_derpsi(dery,wave2,2,'DEN',hidden_TR,rankj,calc_rank)
#else
$PBROKEN $TR   call transfer_derpsi(derx,wave2,1,'DEN',hidden_TR)
$PBROKEN $TR   call transfer_derpsi(dery,wave2,2,'DEN',hidden_TR)
#endif

$PBROKEN $TR   do i=1,N2
$PBROKEN $TR      wave_global  = si+N+i
$PBROKEN $TR      ranki        = rank_map(wave_global)
$PBROKEN $TR      wave         = spwf_inverse(wave_global)
#if(USE_MPI>0)
$PBROKEN $TR      call transfer_psi(psi, wave, 'DEN', ranki, calc_rank)
#else
$PBROKEN $TR      call transfer_psi(psi, wave, 'DEN')
#endif
$PBROKEN $TR      if(MPI_RANK.eq.calc_rank) then
$PBROKEN $TR        NablaMElements(1,1,wave_global,wave2_global) = dv*           &
$PBROKEN $TR        & sum(         derx(:,1) * psi(:,1) + derx(:,2) * psi(:,2)   &
$PBROKEN $TR        &           +  derx(:,3) * psi(:,3) + derx(:,4) * psi(:,4))

$PBROKEN $TR        NablaMElements(2,2,wave_global,wave2_global) = dv*           &
$PBROKEN $TR        & sum(          dery(:,2) * psi(:,1) - dery(:,1) * psi(:,2)  &
$PBROKEN $TR        &            -  dery(:,3) * psi(:,4) + dery(:,4) * psi(:,3))

$PBROKEN $TR        NablaMElements(1,1,wave2_global,wave_global) = &
$PBROKEN $TR        &               - NablaMElements(1,1,wave_global,wave2_global)
$PBROKEN $TR        NablaMElements(2,2,wave2_global,wave_global) = &
$PBROKEN $TR        &                 NablaMElements(2,2,wave_global,wave2_global)
$PBROKEN $TR      endif
$PBROKEN $TR   enddo
$PBROKEN $TR  enddo
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

      ! Broadcast all matrix elements from the calculating rank to all the rest
      call MPI_BCAST(NablaMElements(:,:,si+1:si+T, si+1:si+T), 6*T**2, &
      &              MPI_REAL8, calc_rank, MPI_COMM_WORLD, mpi_err)
      si = si + T 
    enddo
#endif
    !---------------------------------------------------------------------------
    call stop_timer(T_NablaMElements)

end function CompNablaMelements

subroutine print_boxsize_check(R)
  !-----------------------------------------------------------------------------
  ! Print the maximum values of the densities D_I_I and DP_I_I at the edges of 
  ! the box, to check that we are not dealing with a too small box.
  !-----------------------------------------------------------------------------
  type(DensityVector), intent(in), target :: R
  real(KIND=dp), pointer                  :: rho3D(:,:,:,:)
  real(KIND=dp), pointer                  :: rhoP_3D(:,:,:,:)

  1 format ('------------------------- Box Size Check ---------------------------')
 11 format (' Normal  density        rho  =  D_I_I')
 12 format (' Pairing density \tilde{rho} = DP_I_I')
  2 format (' Xmax = (nx+0.5)dx = ', f10.3, ' fm,  max(rho(X=Xmax)) = ', es12.3 )
! 21 format (' Xmin =-(nx+0.5)dx = ', f6.3, 'fm,  max(rho(X=Xmin)) = ', e12.3 )

  3 format (' Ymax = (ny+0.5)dx = ', f10.3, ' fm,  max(rho(Y=Ymax)) = ', es12.3 )
! 31 format (' Ymin =-(ny+0.5)dx = ', f5.3, 'fm,  max(rho(Y=Ymax)) = ', e12.3 )

  4 format (' Zmax = (nz+0.5)dx = ', f10.3, ' fm,  max(rho(Z=Zmax)) = ', es12.3 )
$PBROKEN 41 format (' Zmin =-(nz+0.5)dx = ', f10.3, ' fm,  max(rho(Z=Zmin)) = ', es12.3 )
  
  rho3D(1:nx,1:ny,1:nz,1:2)   => R%D_I_I
  
  print 1
  print 11
  print 2, meshX(nx) , maxval(sum(rho3D(nx,:,:,:),3))
  print 3 , meshY(ny), maxval(sum(rho3D(:,ny,:,:),3))
  print 4 , meshZ(nz), maxval(sum(rho3D(:,:,nz,:),3))
$PBROKEN  print 41, meshZ(1) , maxval(sum(rho3D(:,:,1,:),3))

  if(pairingtype .ne. 0) then
    rhoP_3D(1:nx,1:ny,1:nz,1:2) => R%DP_I_I

    print 12
    print 2, meshX(nx) , maxval(abs(sum(rhoP_3D(nx,:,:,:),3)))
    print 3 , meshY(ny), maxval(abs(sum(rhoP_3D(:,ny,:,:),3)))
    print 4 , meshZ(nz), maxval(abs(sum(rhoP_3D(:,:,nz,:),3)))
$PBROKEN  print 41, meshZ(1) , maxval(abs(sum(rhoP_3D(:,:,1,:),3)))
  endif
  
end subroutine print_boxsize_check

end module densities
