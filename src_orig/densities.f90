module densities
!===============================================================================
!  #######   ##   #    # #####   ##   #      #    #  ####
!     #     #  #  ##   #   #    #  #  #      #    # #
!     #    #    # # #  #   #   #    # #      #    #  ####
!     #    ###### #  # #   #   ###### #      #    #      #
!     #    #    # #   ##   #   #    # #      #    # #    #
!     #    #    # #    #   #   #    # ######  ####   ####
!
!  Copyright W. Ryssens & M. Bender
!
!===============================================================================
!
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
    ! Previous value(s) of the density rho.
    real(KIND=dp), allocatable :: D_I_I_hist(:,:,:)
    !---------------------------------------------------------------------------
    ! Charge density of the protons, possibly including the correction for the 
    ! finite size of the proton. It is stored here, as both the moments module 
    ! and the coulomb module need it, even though Coulomb depends on the moments
    ! module.
    real(KIND=dp), allocatable :: chargedensity(:,:,:)
    !---------------------------------------------------------------------------
    ! The amount of iterations to keep in memory for the density mixing
    integer           :: memory = 1
    
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
    
contains

subroutine densit(SaveRho)
    !---------------------------------------------------------------------------
    ! Calculate all of the densities. 
    ! If SaveRho=.false., do not save the previous values to history!
    !---------------------------------------------------------------------------
    integer      :: i, it, wave, wave2, B, N, si, N2
    real(KIND=dp):: weight
    logical      :: SaveRho

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
    endif   
    if(SaveRho) then
      do i=1,memory-1
          D_I_I_hist(:,:,memory-i+1) = D_I_I_hist(:,:,memory-i)
      enddo
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
      ! Saving the input density for mixing    
      D_I_I_hist(:,:,1) = D_I_I
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
      &               cantransfo,cancuttransfo)
     
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
        it = 2          
        if( B.le. 4) it = 1
        do wave=1,N
$TR          do wave2=wave,N      
$NTR          do wave2=N+1,N+N2      

            weight  =  2 *kappa_pairing(si+wave,si+wave2)*                     &
            &                               Pcutoffs(si+wave)*Pcutoffs(si+wave2)
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
    
    if(all(D_I_I_hist.eq.0.0)) return
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Compute the residual
    resid = D_I_I - D_I_I_hist(:,:,1)

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
    end select
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Safeguard
    where(D_I_I.lt.1d-10) D_I_I = 0 
end subroutine MassageDensity

function CompNablaMelements() result(NablaMelements)
    !---------------------------------------------------------------------------
    ! Computes the matrix elements of Nabla
    !
    !   < Psi_i | \nabla | \Psi_j >
    !
    ! Some notes:
    ! 1) They are, in general, complex!
    ! 2) Symmetries restrict them in weird ways:
    !    Same p             => <Nx> = <Ny> = <Nz> = 0
    !    Same signature     => <Nx> = <Ny> = 0
    !    Different signature=> <Nz> = 0
    !    Different isospin  => <Nx> = <Ny> = <Nz> = 0
    ! 3) We must not forget the time-reversed states when time reversal is
    !    conserved!
    !---------------------------------------------------------------------------
    ! IMPORTANT NOTE
    !---------------------------------------------------------------------------
    ! The symmetries under (i <=> j) have been simply set to 1 and not been 
    ! checked. I implemented this lazily, since the 2-body COMcorrection
    ! only depends on the matrix elements squared.
    !---------------------------------------------------------------------------
    ! Why do we need them?
    !  1) Calculation of 2-body COM correction in Energy module.
    !---------------------------------------------------------------------------
    ! This is one of the most complicated routines in MOCCa, due to the
    ! complications of the symmetries here. Any suggestions to make this
    ! less complicated are absolutely welcome!
    ! P.S. I nominate the subroutine pipj in CR8/EV8/EV4 as the place in the
    ! respective codes that needs clarification the most, as several
    ! transformation are totally implicit.
    !---------------------------------------------------------------------------
  
    integer       :: i,j
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
    
    do i=1, HFBlocks(1)
     do j=HFBlocks(1)+HFBlocks(2)+1, HFBlocks(1)+HFBlocks(2)+HFBlocks(3)
	  !-------------------------------------------------------------------------
      ! Positive parity neutrons with s=+i
      Derx  =  TimeReverse(DendPsi(:,1,:,j))
      Dery  =  TimeReverse(DendPsi(:,2,:,j))
      Derz  =  DendPsi(:,3,:,j)
      psi   =  DenPsi(:,:,i)
      !-------------------------------------------------------------------------
      NablaMElements(1,1,i,j) = dv*                                            &
      & sum(                     derx(:,1) * psi(:,1) + derx(:,2) * psi(:,2)   &
      &                       +  derx(:,3) * psi(:,3) + derx(:,4) * psi(:,4))
      NablaMElements(2,2,i,j) = dv*                                            &
      & sum(                     dery(:,2) * psi(:,1) - dery(:,1) * psi(:,2)   &
      &                       -  dery(:,3) * psi(:,4) + dery(:,4) * psi(:,3))
      NablaMElements(3,1,i,j) = dv*                                            &
      & sum(                     derz(:,1) * psi(:,1) + derz(:,2) * psi(:,2)   &
      &                       +  derz(:,3) * psi(:,3) + derz(:,4) * psi(:,4))

      !-------------------------------------------------------------------------     
      ! Negative parity neutrons with s=+i can be obtained with symmetry
      !-------------------------------------------------------------------------
      NablaMElements(1,1,j,i) = - NablaMElements(1,1,i,j)
      NablaMElements(2,2,j,i) =   NablaMElements(2,2,i,j)
      NablaMElements(3,1,j,i) =   NablaMElements(3,1,i,j)
     enddo
    enddo

    do i=sum(HFblocks(1:4))+1,sum(HFblocks(1:5))
     ! Positive parity protons with s=+i
     do j=sum(HFBlocks(1:6))+1,sum(HFblocks(1:7))
      Derx  =  TimeReverse(DendPsi(:,1,:,j))
      Dery  =  TimeReverse(DendPsi(:,2,:,j))
      Derz  =  DendPsi(:,3,:,j)
      psi   =  DenPsi(:,:,i)

      !-------------------------------------------------------------------------
      NablaMElements(1,1,i,j) = dv*                                            &
      & sum(                     derx(:,1) * psi(:,1) + derx(:,2) * psi(:,2)   &
      &                        + derx(:,3) * psi(:,3) + derx(:,4) * psi(:,4))
      NablaMElements(2,2,i,j) = dv*                                            &
      & sum(                     dery(:,2) * psi(:,1) - dery(:,1) * psi(:,2)   &
      &                        - dery(:,3) * psi(:,4) + dery(:,4) * psi(:,3))
      NablaMElements(3,1,i,j) = dv*                                            &
      & sum(                     derz(:,1) * psi(:,1) + derz(:,2) * psi(:,2)   &
      &                        + derz(:,3) * psi(:,3) + derz(:,4) * psi(:,4))
      !-------------------------------------------------------------------------     
      ! Negative parity neutrons with s=+i can be obtained with symmetry
      !-------------------------------------------------------------------------
      NablaMElements(1,1,j,i) = - NablaMElements(1,1,i,j)
      NablaMElements(2,2,j,i) =   NablaMElements(2,2,i,j)
      NablaMElements(3,1,j,i) =   NablaMElements(3,1,i,j)
     enddo
    enddo

    !---------------------------------------------------------------------------
end function CompNablaMelements

subroutine clean_densities
$CLEANING
    if(allocated(D_I_I_hist))    deallocate(D_I_I_hist)
    if(allocated(chargedensity)) deallocate(chargedensity)
    nullify(DenPsi)
    nullify(DendPsi)
    nullify(DenddPsi)
    nullify(DendddPsi)

end subroutine clean_densities

!subroutine writedensity(den ,fname)
!  !-----------------------------------------------------------------------------
!  ! Write the density to file for plotting afterwards.
!  !-----------------------------------------------------------------------------

!  real(KIND=dp),intent(in), target :: den(:,:)
!  real(KIND=dp), pointer           :: dn(:,:,:), dp(:,:,:)
!  character(len=*), intent(in)     :: fname
!  integer                          :: io, i,j,k
!  

!  open(1,file=fname, iostat=io)
!  if(io.ne.0) then    
!    print *, 'Something went wrong with writing a density to file.'
!    print *, 'filename = ', fname
!    stop
!  endif

!  dn(1:nx,1:ny,1:nz)  => den(:,1)
!  dp(1:nx,1:ny,1:nz)  => den(:,2)
!  
!  write(1, fmt='(3i3, f15.3)') nx,ny, nz, dx
!  do k=1,nz
!    do j=1,ny
!      do i=1,nx
!        write(1, fmt='(2f8.3, 3es25.12)')                                      & 
!        &               meshx(i), meshx(j), meshz(k), dn(i,j,k), dp(i,j,k) 
!      enddo
!    enddo
!  enddo

!  close(1)
!end subroutine writedensity

end module densities
