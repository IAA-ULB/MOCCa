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
! CLEANING        : [WAY too long to include here]
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
!===============================================================================
use compilation
use geninfo
use wavefunctions
use pairing
use derivatives 
use preconditioning 
use timing

implicit none

    !---------------------------------------------------------------------------
    ! Type declaration of the various densities
$DECLARATION   

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
    ! Pointer to which basis is supposed to be used to calculate the densities
    ! Based on pairingtype
    !  (0) HF  => use the HF basis
    !  (1) BCS => use the HF basis
    !  (2) HFB => Use the canonical basis
    real(KIND=dp), pointer ::      DenPsi(:,:,:)
    real(KIND=dp), pointer ::   DendPsi(:,:,:,:)
    real(KIND=dp), pointer ::  DenddPsi(:,:,:,:)
    real(KIND=dp), pointer :: DendddPsi(:,:,:,:)
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

subroutine densit(ifail, SaveRho)
    !---------------------------------------------------------------------------
    ! Calculate all of the densities. 
    ! If SaveRho=.false., do not save the previous values to history!
    !---------------------------------------------------------------------------
    integer, intent(out) :: ifail

    integer      :: i, it, wave, wave2, B, N, si, N2, T
    real(KIND=dp):: weight
    logical      :: SaveRho
    real(KIND=dp), allocatable :: kappa_cut(:,:)

    call start_timer(T_densities)

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Allocation and initialization
$INITIALIZATION

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

    if(PairingType.eq.2) then
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      ! Construct the transformations to 
      ! a) Canonical basis, where rho_pairing is diagonal
      ! b) Cut-canonical basis, where kappa_pairing with cutoffs is diagonal
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      call start_timer(T_den_can)
      call Canonical(rho_pairing, kappa_pairing, rho_can, kappa_can,           &
      &               cantransfo,cancuttransfo,ifail)
     
      ! Apply this transformation
      call ConstructCanonicalBasis(cantransfo)
      call stop_timer(T_den_can)

      call derivecan()
    endif

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
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    call start_timer(T_den_ph)
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! PARTICLE-HOLE DENSITIES
    do wave=1,nwt
        ! Isospin is neutron in the first half of blocks, proton in the rest
        it = 2
        if(wave.le.nwn) it = 1
        
        ! For ordinary densities
        weight  = rho_can(wave) 

        do i=1,mv
$EXPRESSION
        enddo
    enddo
    call stop_timer(T_den_ph)

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Construct the basis where kappa (with cutoffs) is canonical 

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
      do wave=1,nwt
          ! Isospin is neutron in the first half of blocks, proton in the rest
          it = 2
          if(wave.le.nwn) it = 1
          weight  = 2 * kappa_can(wave) * Pcutoffs(wave)**2
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
      !------------------------------------
      ! HFB calculations: full summations.
      !------------------------------------
      
      ! For now, sum the pairing densities in the HFbasis, not the canonical
      ! basis.
      DenPsi    => HFPsi   ; DenDPsi   => HFdPsi 
      DenddPsi  => HFddPsi ; DendddPsi => HFdddpsi
      
      si = 0
      do B=1,8,2
        N = HFBlocks(B) ;  if (N.eq.0) cycle
        N2= HFBlocks(B+1)
        T = N+N2
        it = 2          
        if( B.le. 4) it = 1

        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        ! Calculation of the pairing cutoffs * kappa
        kappa_cut = kappa_pairing(si+1:si+T,si+1:si+T)
        if ((.not. diagsphamil)) then  
          ! Transform to the HF-basis
          kappa_cut =matmul(transpose(HFtransfo(si+1:si+T,si+1:si+T)),kappa_cut)
          kappa_cut =matmul(          kappa_cut, HFtransfo(si+1:si+T,si+1:si+T))
        endif
        ! Multiply with the cutoffs
        do wave=1,N
$TR          do wave2=wave,N      
$NTR          do wave2=N+1,N+N2      
              kappa_cut(wave,wave2) = kappa_cut(wave,wave2) &
              &                     *Pcutoffs(si+wave)*Pcutoffs(si+wave2)
$NTR              kappa_cut(wave2,wave) = kappa_cut(wave2,wave) &
$NTR              &                     *Pcutoffs(si+wave)*Pcutoffs(si+wave2)              
          enddo
        enddo
        if((.not. diagsphamil) ) then 
          ! Transform back to the basis in memory
          kappa_cut=matmul(HFtransfo(si+1:si+T,si+1:si+T), kappa_cut)
          kappa_cut=matmul( kappa_cut,transpose(HFtransfo(si+1:si+T,si+1:si+T)))
        endif

        do wave=1,N
$TR          do wave2=wave,N      
$NTR          do wave2=N+1,N+N2      

            weight=2*kappa_cut(wave,wave2)
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
    end select
    call stop_timer(T_den_pp)
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

    call stop_timer(T_densities)

end subroutine densit

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
      print *, 'Invalid iso argument to couple_iso. iso = ', iso
      stop    
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
    !
    ! Some notes:
    ! 1) These are, in general, complex numbers!
    ! 2) Symmetries restrict them in weird ways:
    !    Same p             => <Nx> = <Ny> = <Nz> = 0
    !    Same signature     => <Nx> = <Ny> = 0
    !    Different signature=> <Nz> = 0
    !    Different isospin  => <Nx> = <Ny> = <Nz> = 0
    ! 3) We must not forget the time-reversed states when time reversal is
    !    conserved!
    !---------------------------------------------------------------------------
    integer       :: i,j, B, N, si, N2, N3, N4, wave, wave2
    real(KIND=dp) :: NablaMElements(3,2,nwt,nwt), psi(mv,4)
    real(KIND=dp) :: derx(mv,4), dery(mv,4), derz(mv,4)

    !---------------------------------------------------------------------------
    NablaMElements= 0.0_dp
    
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

    si = 0 
    do B=1,8,4 ! This loop is essentially over protons vs neutrons
      N = HFblocks(B)  ;  N2= HFBlocks(B+1)
      N3= HFBlocks(B+2);  N4= HFBlocks(B+3)

      !-------------------------------------------------------------------------
      ! Re < p_z > : only non-zero matrix elements when equal signature
      !              and opposite parity (if it is conserved)
      
      ! Block 1 with block 3
      do i=1,N
        wave  = si+i
        psi   = DenPsi(:,:,wave)
        do j=1,N3
          wave2 = si+N+N2+j
          Derz  =  DendPsi(:,3,:,wave2)

          NablaMElements(3,1,wave,wave2) = dv*                                 &
          & sum(                    derz(:,1) * psi(:,1) + derz(:,2) * psi(:,2)&
          &                      +  derz(:,3) * psi(:,3) + derz(:,4) * psi(:,4))

          NablaMElements(3,1,wave2,wave) =   NablaMElements(3,1,wave,wave2)
        enddo
      enddo
      
      ! Block 1 with block 1 (parity broken)
$PBROKEN      do i=1,N
$PBROKEN        wave  = si+i
$PBROKEN        psi   = DenPsi(:,:,wave)
$PBROKEN        do j=1,N
$PBROKEN          wave2 = si+j
$PBROKEN          Derz  =  DendPsi(:,3,:,wave2)
$PBROKEN
$PBROKEN          NablaMElements(3,1,wave,wave2) = dv*                                 &
$PBROKEN          & sum(                    derz(:,1) * psi(:,1) + derz(:,2) * psi(:,2)&
$PBROKEN          &                      +  derz(:,3) * psi(:,3) + derz(:,4) * psi(:,4))
$PBROKEN
$PBROKEN          NablaMElements(3,1,wave2,wave) =   NablaMElements(3,1,wave,wave2)
$PBROKEN        enddo
$PBROKEN      enddo
  
      ! Block 2 with block 4
      do i=1,N2
        wave  = si+N+i
        psi   = DenPsi(:,:,wave)
        do j=1,N4
          wave2 = si+N+N2+N3+j
          Derz  =  DendPsi(:,3,:,wave2)
          ! Re < p_z > 
          NablaMElements(3,1,wave,wave2) = dv*                                 &
          & sum(                    derz(:,1) * psi(:,1) + derz(:,2) * psi(:,2)&
          &                      +  derz(:,3) * psi(:,3) + derz(:,4) * psi(:,4))

          NablaMElements(3,1,wave2,wave) =   NablaMElements(3,1,wave,wave2)
        enddo
      enddo
      
      ! Block 2 with block 2 (parity broken)
$PBROKEN      do i=1,N2
$PBROKEN        wave  = si+N+i
$PBROKEN        psi   = DenPsi(:,:,wave)
$PBROKEN        do j=1,N2
$PBROKEN          wave2 = si+N+j
$PBROKEN          Derz  =  DendPsi(:,3,:,wave2)
$PBROKEN          ! Re < p_z > 
$PBROKEN          NablaMElements(3,1,wave,wave2) = dv*                                 &
$PBROKEN          & sum(                    derz(:,1) * psi(:,1) + derz(:,2) * psi(:,2)&
$PBROKEN          &                      +  derz(:,3) * psi(:,3) + derz(:,4) * psi(:,4))
$PBROKEN
$PBROKEN          NablaMElements(3,1,wave2,wave) =   NablaMElements(3,1,wave,wave2)
$PBROKEN        enddo
$PBROKEN      enddo

      !-------------------------------------------------------------------------
      ! Re < p_x >, Im <p_y> 

      do i=1,N
        wave  = si+i
        psi   = DenPsi(:,:,wave)
           ! Block 1 with block 4  (T-broken)    
$NTR       do j=1, N4
$NTR          wave2 = si+N+N2+N3+j
$NTR          Derx  =  DendPsi(:,1,:,wave2)
$NTR          Dery  =  DendPsi(:,2,:,wave2)
           ! Block 1 with block 3  (T-conserved)    
$TR        do j=1, N3
$TR           wave2 = si+N+j
$TR           Derx  =  TimeReverse(DendPsi(:,1,:,wave2))
$TR           Dery  =  TimeReverse(DendPsi(:,2,:,wave2))

          NablaMElements(1,1,wave,wave2) = dv*                                 &
          & sum(                 derx(:,1) * psi(:,1) + derx(:,2) * psi(:,2)   &
          &                   +  derx(:,3) * psi(:,3) + derx(:,4) * psi(:,4))

          NablaMElements(2,2,wave,wave2) = dv*                                 &
          & sum(                  dery(:,2) * psi(:,1) - dery(:,1) * psi(:,2)  &
          &                    -  dery(:,3) * psi(:,4) + dery(:,4) * psi(:,3))

          NablaMElements(1,1,wave2,wave) = - NablaMElements(1,1,wave,wave2)
          NablaMElements(2,2,wave2,wave) =   NablaMElements(2,2,wave,wave2)
        enddo
      enddo
      
$PBROKEN      do i=1,N
$PBROKEN        wave  = si+i
$PBROKEN        psi   = DenPsi(:,:,wave)
                    ! Block 1 with block 2  (T-broken, P-broken)    
$PBROKEN $NTR       do j=1, N2
$PBROKEN $NTR          wave2 = si+N+j
$PBROKEN $NTR          Derx  =  DendPsi(:,1,:,wave2)
$PBROKEN $NTR          Dery  =  DendPsi(:,2,:,wave2)
                    ! Block 1 with block 1  (T-conserved, P-broken)    
$PBROKEN $TR        do j=1, N
$PBROKEN $TR           wave2 = si+j
$PBROKEN $TR           Derx  =  TimeReverse(DendPsi(:,1,:,wave2))
$PBROKEN $TR           Dery  =  TimeReverse(DendPsi(:,2,:,wave2))

$PBROKEN          NablaMElements(1,1,wave,wave2) = dv*                                 &
$PBROKEN          & sum(                 derx(:,1) * psi(:,1) + derx(:,2) * psi(:,2)   &
$PBROKEN          &                   +  derx(:,3) * psi(:,3) + derx(:,4) * psi(:,4))

$PBROKEN          NablaMElements(2,2,wave,wave2) = dv*                                 &
$PBROKEN          & sum(                  dery(:,2) * psi(:,1) - dery(:,1) * psi(:,2)  &
$PBROKEN          &                    -  dery(:,3) * psi(:,4) + dery(:,4) * psi(:,3))

$PBROKEN          NablaMElements(1,1,wave2,wave) = - NablaMElements(1,1,wave,wave2)
$PBROKEN          NablaMElements(2,2,wave2,wave) =   NablaMElements(2,2,wave,wave2)
$PBROKEN        enddo
$PBROKEN      enddo

      do i=1,N2
        wave  = si+N+i
        psi   = DenPsi(:,:,wave)
           ! Block 2 with block 3  (T-broken)    
$NTR       do j=1, N3
$NTR          wave2 = si+N+N2+j
$NTR          Derx  =  DendPsi(:,1,:,wave2)
$NTR          Dery  =  DendPsi(:,2,:,wave2)
           ! Block 2 with block 4  (T-conserved)    
$TR        do j=1, N4
$TR           wave2 = si+N+N2+N3+j
$TR           Derx  =  TimeReverse(DendPsi(:,1,:,wave2))
$TR           Dery  =  TimeReverse(DendPsi(:,2,:,wave2))

          NablaMElements(1,1,wave,wave2) = dv*                                 &
          & sum(                 derx(:,1) * psi(:,1) + derx(:,2) * psi(:,2)   &
          &                   +  derx(:,3) * psi(:,3) + derx(:,4) * psi(:,4))

          NablaMElements(2,2,wave,wave2) = dv*                                 &
          & sum(                  dery(:,2) * psi(:,1) - dery(:,1) * psi(:,2)  &
          &                    -  dery(:,3) * psi(:,4) + dery(:,4) * psi(:,3))

          NablaMElements(1,1,wave2,wave) = - NablaMElements(1,1,wave,wave2)
          NablaMElements(2,2,wave2,wave) =   NablaMElements(2,2,wave,wave2)
        enddo
      enddo
      
$PBROKEN      do i=1,N2
$PBROKEN        wave  = si+N+i
$PBROKEN        psi   = DenPsi(:,:,wave)

                    ! Block 2 with block 1 (P-broken, T-broken)
                    ! This one is likely superfluous, as we have this 
                    ! combination already once above
$PBROKEN $NTR       do j=1, N
$PBROKEN $NTR          wave2 = si+j
$PBROKEN $NTR          Derx  =  DendPsi(:,1,:,wave2)
$PBROKEN $NTR          Dery  =  DendPsi(:,2,:,wave2)

                    ! Block 2 with block 2 (P-broken, T-conserved)
$PBROKEN $TR        do j=1, N2
$PBROKEN $TR           wave2 = si+N+j
$PBROKEN $TR           Derx  =  TimeReverse(DendPsi(:,1,:,wave2))
$PBROKEN $TR           Dery  =  TimeReverse(DendPsi(:,2,:,wave2))

$PBROKEN          NablaMElements(1,1,wave,wave2) = dv*                                 &
$PBROKEN          & sum(                 derx(:,1) * psi(:,1) + derx(:,2) * psi(:,2)   &
$PBROKEN          &                   +  derx(:,3) * psi(:,3) + derx(:,4) * psi(:,4))

$PBROKEN          NablaMElements(2,2,wave,wave2) = dv*                                 &
$PBROKEN          & sum(                  dery(:,2) * psi(:,1) - dery(:,1) * psi(:,2)  &
$PBROKEN          &                    -  dery(:,3) * psi(:,4) + dery(:,4) * psi(:,3))

$PBROKEN          NablaMElements(1,1,wave2,wave) = - NablaMElements(1,1,wave,wave2)
$PBROKEN          NablaMElements(2,2,wave2,wave) =   NablaMElements(2,2,wave,wave2)
$PBROKEN        enddo
$PBROKEN      enddo

      si = si + N + N2 + N3 + N4 
    enddo    
    !---------------------------------------------------------------------------
end function CompNablaMelements

subroutine print_boxsize_check()
  !-----------------------------------------------------------------------------
  ! Print the maximum values of the density at the edges of the box, to
  ! check that we are not dealing with a too small box.
  !-----------------------------------------------------------------------------

  real(KIND=dp), pointer :: rho3D(:,:,:,:)

  1 format ('----------------------- Box Size Check -------------------------')
  2 format (' Xmax = (nx+0.5)dx = ', f6.3, ' fm,  max(rho(X=Xmax)) = ', es12.3 )
! 21 format (' Xmin =-(nx+0.5)dx = ', f6.3, 'fm,  max(rho(X=Xmin)) = ', e12.3 )

  3 format (' Ymax = (ny+0.5)dx = ', f6.3, ' fm,  max(rho(Y=Ymax)) = ', es12.3 )
! 31 format (' Ymin =-(ny+0.5)dx = ', f5.3, 'fm,  max(rho(Y=Ymax)) = ', e12.3 )

  4 format (' Zmax = (nz+0.5)dx = ', f6.3, ' fm,  max(rho(Z=Zmax)) = ', es12.3 )
$PBROKEN 41 format (' Zmin =-(nz+0.5)dx = ', f6.3, 'fm,  max(rho(Z=Zmin)) = ', es12.3 )
  
  rho3D(1:nx,1:ny,1:nz,1:2) => D_I_I
  
  print 1
  print 2, meshX(nx) , maxval(sum(rho3D(nx,:,:,:),3))
!  print 21, meshX(1) , maxval(sum(rho3D(1,:,:),3))

  print 3 , meshY(ny), maxval(sum(rho3D(:,ny,:,:),3))
!  print 31, meshY(1) , maxval(sum(rho3D(:,1,:,:),3)

  print 4 , meshZ(nz), maxval(sum(rho3D(:,:,nz,:),3))
$PBROKEN  print 41, meshZ(1) , maxval(sum(rho3D(:,:,1,:),3))

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
