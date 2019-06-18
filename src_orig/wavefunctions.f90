module wavefunctions
 !==============================================================================
 !  #######   ##   #    # #####   ##   #      #    #  ####
 !     #     #  #  ##   #   #    #  #  #      #    # #
 !     #    #    # # #  #   #   #    # #      #    #  ####
 !     #    ###### #  # #   #   ###### #      #    #      #
 !     #    #    # #   ##   #   #    # #      #    # #    #
 !     #    #    # #    #   #   #    # ######  ####   ####
 !
 !  Copyright W. Ryssens & M. Bender
 !
 !==============================================================================
 !
 ! Module containing the single-particle wave-functions (spwfs for short)
 ! for the Tantalus program.
 !
 !
 !
 !
 !==============================================================================
 ! Hephaestos:
 ! BLOCKS $BLOCKS
 !==============================================================================
 use compilation
 use derivatives
 use geninfo
 use nil8
 
 implicit none
 
 !------------------------------------------------------------------------------
 ! Array containing the values of the spwfs in the Hartree-Fock basis
 ! and their derivatives
 !
 ! Note that higher-order derivative tensors are stored in lexicographical order
 ! in order to cut down on the number of indices and wasted computation.
 !            1    2    3    4    5    6    7    8    9    10
 ! 1st order: Dx   Dy   Dz
 ! 2nd order: Dxx  Dxy  Dxz  Dyy  Dyz  Dzz
 ! 3rd order: Dxxx Dxxy Dxxz Dxyy Dxyz Dxzz Dyyy Dyyz Dyzz Dzzz
 real(KIND=dp), allocatable, target ::      HFPsi(:,:,:)
 real(KIND=dp), allocatable, target ::   HFdPsi(:,:,:,:)!First order derivatives
 real(KIND=dp), allocatable, target ::  HFddPsi(:,:,:,:)!Second order derivatives
 real(KIND=dp), allocatable, target :: HFdddPsi(:,:,:,:)!Third order derivatives
 
 !------------------------------------------------------------------------------
 ! Array containing the values of the spwfs in the Canonical basis
 ! and their derivatives.
 real(KIND=dp), allocatable, target ::     CANPsi(:,:,:)
 real(KIND=dp), allocatable, target ::  CANdPsi(:,:,:,:)!First order derivatives
 real(KIND=dp), allocatable, target ::CANddPsi(:,:,:,:)!Second order derivatives
 real(KIND=dp), allocatable, target ::CANdddPsi(:,:,:,:)!Third order derivatives
 !------------------------------------------------------------------------------
 ! Density matrix rho and anomalous density matrix kappa
 ! Dimensions (nwt, nwt) (although many are zero when symmetries are conserved)
 !real(KIND=dp), allocatable :: rho(:,:), kappa(:,:)
 !------------------------------------------------------------------------------
 ! Single-particle energies, diagonal elements of the single-particle
 ! hamiltonian
 ! \langle psi_i | h | psi_i \rangle
 real(KIND=dp), allocatable :: spenergies(:) 
 real(KIND=dp), allocatable :: dispersions(:)
 real(KIND=dp), allocatable :: canenergies(:)
 !------------------------------------------------------------------------------
 ! Number of the blocks with the same quantum numbers that divide up the 
 ! HFBasis. Any possibility has a maximum of two spatial operators that 
 ! introduce a quantum number, while proton-neutron symmetry adds another one. 
 ! The number of blocks thus needs to be decided on compile time by a
 ! replacement script. HFBlocks contains the sizes of the various blocks.  
 ! For ease of reference, we also store the number of neutron and proton 
 ! spwfs independently.
 !------------------------------------------------------------------------------
 ! Examples:
 !    * EV8-like calculation:     8 blocks (P,Rz,T3)
 !    * EV4-like calculation:     4 blocks (  Rz,T3)
 !    * Rx broken, Sx conserved : 4 blocks (  Sx,T3)
 !------------------------------------------------------------------------------
 integer, parameter   :: Blocks          =$BLOCKS
 integer              :: HFBlocks(Blocks)=0
 integer              :: nwn, nwp
 !------------------------------------------------------------------------------
 ! Properties of the single-particle wave-functions with regard to reflections
 ! of the axes.
 integer, allocatable :: sx(:,:), sy(:,:), sz(:,:)
 !------------------------------------------------------------------------------
 ! Oscillator frequencies to use for the initialization with a Nilsson  
 ! hamiltonian.
 real(KIND=dp) :: osc_freq(3) = 0.2

contains 

  subroutine ReadWFdata(file_number)
    !---------------------------------------------------------------------------
    ! Read the number of single-particle neutron and proton wave-functions.
    !
    !---------------------------------------------------------------------------

    integer(dp), intent(in), optional   :: file_number   

    namelist /wfs/ nwn, nwp, osc_freq

    if(present(file_number)) then
      read(unit=file_number, nml = wfs)
    else
      read(unit=*, nml = wfs)
    endif
    nwt = nwn + nwp
  end subroutine ReadWFdata

  subroutine iniwavefunctions()   
    !---------------------------------------------------------------------------
    ! Build harmonic oscillator eigenfunctions
    ! a) in an EV8-like box
    ! b) expanding to the full box
    ! c) restricting again to the box desired 
    !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Also initialized:
    !  *) Diagonal matrix elements of <h>
    ! 
    ! Not initialized here:
    !  *) Delta for the gaps. Since this module can not know what kind of 
    !     pairing is needed, it cannot correctly guess a structure. 
    !---------------------------------------------------------------------------
    
    real(KIND=dp)             :: homegax, homegay,homegaz, alpha,qqq
    !real(KIND=dp),allocatable :: fullbox(:,:,:)
    !integer                   :: j,k, p
    integer                   :: i
    integer, allocatable      :: kparz(:)
        
    alpha = 0.2    
    qqq   = 1.0    
    homegaz  = alpha*qqq**(-2.0/3.0)
    homegax  = alpha*qqq**(-2*cos(-2*pi/3)/3)
    homegay  = alpha*qqq**(-2*cos(+2*pi/3)/3)
    
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
    ! a) Generating the nilsson wave-functions in an EV8-box   
    call nilsson (HFPsi,kparz,spenergies,6,5,nwt,nwp,nwn,                      &
    &           floor(neutrons),floor(protons),nx,ny,nz,dx,osc_freq)
    allocate(dispersions(nwt))
    allocate(sx(4,nwt), sy(4,nwt), sz(4,nwt))
    
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! b) blow-up into the full box
!    allocate(fullbox(2*nx, 2*ny, 2*nz, 4, nwt))
!    do wave=1,nwt
!        ! Copy the original
!        fullbox(nx+1: 2*nx,ny+1: 2*ny, nz+1: 2*nz,:,wave) = HFPsi   (:,:,:,:,wave)     
!        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!        ! Use y-time-simplex to fill in the y-axis.
!        do j=1,ny
!                ! S^t_y Psi (x,y,z,sigma) = Psi^*(x,-y,z,sigma)
!                fullbox(nx+1: 2*nx,j, nz+1:2*nz,1,wave) =   fullbox(nx+1: 2*nx, 2*ny - j +1, nz+1:2*nz,1,wave)
!                fullbox(nx+1: 2*nx,j, nz+1:2*nz,2,wave) = - fullbox(nx+1: 2*nx, 2*ny - j +1, nz+1:2*nz,2,wave)
!                fullbox(nx+1: 2*nx,j, nz+1:2*nz,3,wave) =   fullbox(nx+1: 2*nx, 2*ny - j +1, nz+1:2*nz,3,wave)
!                fullbox(nx+1: 2*nx,j, nz+1:2*nz,4,wave) = - fullbox(nx+1: 2*nx, 2*ny - j +1, nz+1:2*nz,4,wave)
!        enddo
!        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!        ! Use z-signature to fill in the x-axis
!        do j=1,2*ny
!            do i=1,nx
!                ! R_z Psi (x,y,z,sigma) = -i sigma Psi(-x,-y,z, sigma)
!                fullbox(i,j,nz+1:2*nz,1,wave) =   fullbox(2*nx -i +1, 2*ny - j +1, nz+1:2*nz,2,wave)
!                fullbox(i,j,nz+1:2*nz,2,wave) = - fullbox(2*nx -i +1, 2*ny - j +1, nz+1:2*nz,1,wave)
!                fullbox(i,j,nz+1:2*nz,3,wave) = - fullbox(2*nx -i +1, 2*ny - j +1, nz+1:2*nz,4,wave)
!                fullbox(i,j,nz+1:2*nz,4,wave) =   fullbox(2*nx -i +1, 2*ny - j +1, nz+1:2*nz,3,wave)
!            enddo
!        enddo
!        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!        ! Use parity to fill in the z-axis
!        p = kparz(wave)
!        do k=1,nz
!                do j=1,2*ny
!                        do i=1,2*nx
!                                fullbox(i,j,k,1,wave) = p*fullbox(2*nx -i +1, 2*ny - j +1, 2*nz-k+1,1,wave)
!                                fullbox(i,j,k,2,wave) = p*fullbox(2*nx -i +1, 2*ny - j +1, 2*nz-k+1,2,wave)
!                                fullbox(i,j,k,3,wave) = p*fullbox(2*nx -i +1, 2*ny - j +1, 2*nz-k+1,3,wave)
!                                fullbox(i,j,k,4,wave) = p*fullbox(2*nx -i +1, 2*ny - j +1, 2*nz-k+1,4,wave)
!                        enddo
!                enddo
!        enddo
!    enddo
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    do i=1,nwn
        if(kparz(i) .gt. 0) HFBlocks(1) = HFBlocks(1) +1
        if(kparz(i) .lt. 0) HFBlocks(3) = HFBlocks(3) +1
    enddo
    do i=nwn+1,nwt
        if(kparz(i) .gt. 0) HFBlocks(5) = HFBlocks(5) +1
        if(kparz(i) .lt. 0) HFBlocks(7) = HFBlocks(7) +1
    enddo
    
    do i=1, HFBlocks(1)
        sx(1,i) =  1 ; sy(1,i) = +1 ; sz(1,i) = +1
        sx(2,i) = -1 ; sy(2,i) = -1 ; sz(2,i) = +1 
        sx(3,i) = -1 ; sy(3,i) = +1 ; sz(3,i) = -1
        sx(4,i) =  1 ; sy(4,i) = -1 ; sz(4,i) = -1
    enddo
    do i=HFBlocks(1) + 1,HFBlocks(1) + HFBlocks(3)
        sx(1,i) =  1 ; sy(1,i) = +1 ; sz(1,i) = -1
        sx(2,i) = -1 ; sy(2,i) = -1 ; sz(2,i) = -1 
        sx(3,i) = -1 ; sy(3,i) = +1 ; sz(3,i) = +1
        sx(4,i) =  1 ; sy(4,i) = -1 ; sz(4,i) = +1
    enddo
    do i=HFBlocks(1) + HFBlocks(3)+1,HFBlocks(1) + HFBlocks(3) +HFBlocks(5)
        sx(1,i) =  1 ; sy(1,i) = +1 ; sz(1,i) = +1
        sx(2,i) = -1 ; sy(2,i) = -1 ; sz(2,i) = +1 
        sx(3,i) = -1 ; sy(3,i) = +1 ; sz(3,i) = -1
        sx(4,i) =  1 ; sy(4,i) = -1 ; sz(4,i) = -1
    enddo
    do i=HFBlocks(1)+HFBlocks(3)+HFBlocks(5) + 1,                      &
    &       HFBlocks(1)+HFBlocks(3)+HFBlocks(5) + HFBLocks(7)
        sx(1,i) =  1 ; sy(1,i) = +1 ; sz(1,i) = -1
        sx(2,i) = -1 ; sy(2,i) = -1 ; sz(2,i) = -1 
        sx(3,i) = -1 ; sy(3,i) = +1 ; sz(3,i) = +1
        sx(4,i) =  1 ; sy(4,i) = -1 ; sz(4,i) = +1
    enddo

!    HFPsi = HFPsi*sqrt(2.)
    ! Simply because I distrust the nilsson routine
!    call GramSchmidt
    
  end subroutine iniwavefunctions
  
  subroutine deriveHF()
    !---------------------------------------------------------------------------
    ! Derives all of the single-particle wave-functions. 
    ! b) In the canonical basis
    !---------------------------------------------------------------------------
    integer :: wave,k
    
    if(.not.allocated(HFdPsi)) then
        allocate(HFdPsi(nx*ny*nz,3,4,nwt))
        allocate(HFddPsi(nx*ny*nz,6,4,nwt))
    endif
    
$N3    if(.not.allocated(HFdddpsi)) then
$N3        allocate(HFdddPsi(nx*ny*nz,10,4,nwt))
$N3    endif

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Currently EV8 symmetries are hardcoded.
    do wave=1,nwt
        do k=1,4

$N2        call Derive_tot(HFPsi(:,k,wave), sx(k,wave), sy(k,wave), sz(k,wave),&
$N2        &                                           HFdPsi(:,:,k,wave),     &
$N2        &                                           HFddPsi(:,:,k,wave))

$N3        call Derive_tot(HFPsi(:,k,wave), sx(k,wave), sy(k,wave), sz(k,wave),&
$N3        &                                           HFdPsi(:,:,k,wave),     &
$N3        &                                           HFddPsi(:,:,k,wave),    &
$N3        &                                           HFdddPsi(:,:,k,wave))

        enddo
    enddo
  end subroutine DeriveHF
  
  subroutine deriveCan()
    !---------------------------------------------------------------------------
    ! Derives all of the single-particle wave-functions. 
    ! b) In the canonical basis
    !---------------------------------------------------------------------------
    integer :: wave,k
      
    if(allocated(CanPsi)) then
      if(.not.allocated(CANdPsi)) then
          allocate( CANdPsi(nx*ny*nz,3,4,nwt))
          allocate(CANddPsi(nx*ny*nz,6,4,nwt))
      endif
    endif

$N3    if(allocated(CanPsi)) then
$N3       if(.not.allocated(CANdddpsi)) then
$N3         llocate(CandddPsi(nx*ny*nz,10,4,nwt))
$N3       endif
$N3    endif
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Currently EV8 symmetries are hardcoded.
    if(allocated(CanPsi)) then
      do wave=1,nwt
        do k=1,4

$N2        call Derive_tot(CANPsi(:,k,wave), sx(k,wave),sy(k,wave), sz(k,wave),&
$N2        &                                           CANdPsi(:,:,k,wave),    &
$N2        &                                           CANddPsi(:,:,k,wave))

$N3        call Derive_tot(CANPsi(:,k,wave), sx(k,wave),sy(k,wave), sz(k,wave),&
$N3        &                                           CANdPsi(:,:,k,wave),    &
$N3        &                                           CANddPsi(:,:,k,wave),   &
$N3        &                                           CANdddPsi(:,:,k,wave))

        enddo
      enddo
    endif
    
  end subroutine DeriveCan
  
  function OrderSpwfsISO(Isospin) result(Indices)
    !---------------------------------------------------------------------------
    ! Orders the wavefunctions within an isospin block. 
    !---------------------------------------------------------------------------
    integer, intent(in)        :: Isospin
    
    integer, allocatable       :: Indices(:)
    real(Kind=dp), allocatable :: Energies(:)
    integer                    :: i, nwf,  HolePos, ToInsertIndex
    real(Kind=dp)              :: ToInsert
    
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !Count the number of proton and neutron wavefunctions
    if (Isospin .eq. -1) then
        nwf = nwn
    else
        nwf = nwp
    endif
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !Filling Energies & Indices
    allocate(Indices(nwf), Energies(nwf))
    do i=1,nwf
       Indices(i) = i 
    enddo
    if(Isospin.eq.-1) then
        Energies = spenergies(1:nwn)
    else
        Indices  = Indices + nwn
        Energies = spenergies(nwn+1:nwn+nwp)
    endif
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !Sort the energies
    do i=2,nwf
      !Make a hole at index i
      ToInsert = Energies(i)
      HolePos  = i
      ToInsertIndex = Indices(i)
      do while(ToInsert.lt.Energies(HolePos-1))
        !Move the hole one place down
        Energies(HolePos) = Energies(HolePos-1)
        Indices(HolePos) = Indices(HolePos-1)
        HolePos = HolePos - 1
        if(HolePos.eq.1.0_dp) exit
      enddo
      !Insert the energy at the correct place
      Energies(HolePos) = ToInsert
      Indices(HolePos)  = ToInsertIndex
    enddo
  end function OrderSpwfsISO
  
  function OrderSpwfsSym(block) result(indices)
    !---------------------------------------------------------------------------
    ! Sort the single-particle wave-functions in the given symmetry-block 
    ! by single-particle energy.
    !---------------------------------------------------------------------------
    
    integer, intent(in)        :: block
    integer, allocatable       :: Indices(:)
    real(KIND=dp), allocatable :: Energies(:)
    integer                    :: nwf, startind, HolePos, ToInsertIndex, i
    real(Kind=dp)              :: ToInsert
    
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !Count the number of relevant wavefunctions
    nwf      = HFBlocks(block) 
    startind = sum(HFBlocks(1:block-1))
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !Filling Energies & Indices
    allocate(Indices(nwf), Energies(nwf))
    do i=1,nwf
       Indices(i) = startind + i 
    enddo
    Energies = spenergies(startind+1:startind+nwf)
    
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !Sort the energies
    do i=2,nwf
      !Make a hole at index i
      ToInsert = Energies(i)
      HolePos  = i
      ToInsertIndex = Indices(i)
      do while(ToInsert.lt.Energies(HolePos-1))
        !Move the hole one place down
        Energies(HolePos) = Energies(HolePos-1)
        Indices(HolePos) = Indices(HolePos-1)
        HolePos = HolePos - 1
        if(HolePos.eq.1.0_dp) exit
      enddo
      !Insert the energy at the correct place
      Energies(HolePos) = ToInsert
      Indices(HolePos)  = ToInsertIndex
    enddo
  end function OrderSpwfsSym

  subroutine GramSchmidt
    !---------------------------------------------------------------------------
    ! This subroutine uses a Gram-Schmidt scheme to orthonormalise the Spwfs in
    ! the HF basis. Small point of interest: the orthonormalisation is done in 
    ! order of ascending energy within each symmetry block, in order to avoid
    ! wasting CPU cycles reordering levels.
    !---------------------------------------------------------------------------
    integer  :: b, i,j,nw, mw,l
    integer  :: indices(maxval(HFBlocks))
    real(KIND=dp) ::  norm
    
    do b = 1, Blocks 
        indices = 0
        indices(1:HFblocks(b)) = OrderSpwfsSym(b)
        do i = 1, HFBlocks(b)
            !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ! Normalize wave-function nw
            nw = indices(i)
            norm = sum(HFpsi(:,:,nw)**2) * dv
            HFPsi(:,:,nw) = (sqrt(1.0/norm)) * HFPsi(:,:,nw) 
            
            !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ! Then subtract the projection on \Psi_{nw} from all the following
            ! Spwf.
            ! Re(\Psi(\sigma)_{mw}) = Re(\Psi(\sigma)_{mw})
            !                - Re(< \Psi_{nw}|\Psi_{mw} >) Re(\Psi(\sigma)_{nw})
            !                + Im(< \Psi_{nw}|\Psi_{mw} >) Im(\Psi(\sigma)_{nw})
            ! Im(\Psi(\sigma)_{mw}) = Im(\Psi(\sigma)_{mw})
            !                - Re(< \Psi_{nw}|\Psi_{mw} >) Im(\Psi(\sigma)_{nw})
            !                - Im(< \Psi_{nw}|\Psi_{mw} >) Re(\Psi(\sigma)_{nw})
            !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ! Note that the imaginary part of the inproduct only needs to be 
            ! taken into account when there is no antilinear, hermitian 
            ! symmetry that is conserved.
            !
            ! THIS IS NOT IMPLEMENTED YET HOWEVER!
            !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ! What is also missing is an orthonormalisation versus the spwf
            ! that are assumed to be present but not represented numerically.
            ! The MOCCa example is of course conserved time-reversal but broken
            ! signature.
            !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            do j= i+1, HFBlocks(b)
                mw = indices(j)    
                ! Real part of the inproduct
                norm = sum(HFpsi(:,:,nw)*HFpsi(:,:,mw)) * dv
                do l=1,4*nx*ny*nz
                    HFPsi(l,1,mw) = HFPsi(l,1,mw) - norm * HFPsi(l,1,nw)
                enddo
            enddo
        enddo
    enddo
   
  end subroutine GramSchmidt
  
  function TimeReverse(psi) result(Tpsi)
    !---------------------------------------------------------------------------
    ! Perform a time-reversal on the input spinor.
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in) :: Psi(mv,4)
    real(KIND=dp)             :: Tpsi(mv,4)
    
    TPsi(:,1) =   Psi(:,3)
    TPsi(:,2) = - Psi(:,4)
    TPsi(:,3) = - Psi(:,1)
    TPsi(:,4) =   Psi(:,2)
    
  end function TimeReverse
  
  function angmom_x_real(wf2, wf1, dwf1) result(angmom)
    !---------------------------------------------------------------------------
    ! calculate the matrix element
    !
    !     Re < wf2 | j_x | wf1 >
    !
    ! Note that this routine assumes that all checks regarding symmetries have 
    ! been performed.
    !
    !       J_x = 1/2*( 0  1 ) + i z \partial_y - i y\partial_z
    !                 ( 1  0 )
    !---------------------------------------------------------------------------

    real(KIND=dp), intent(in) :: wf1(mv,4), wf2(mv,4), dwf1(mv,3,4)
    integer                   :: i  
    real(KIND=dp)             :: angmom

    angmom = 0    

    do i=1,mv
      ! Spin part
      angmom = angmom       + 0.5*(  wf2(i,1) * wf1(i,3)                       &
      &                            + wf2(i,2) * wf1(i,4)                       & 
      &                            + wf2(i,3) * wf1(i,1)                       &
      &                            + wf2(i,4) * wf1(i,2))                        
      ! Orbital part
      angmom = angmom &
      &           + wf2(i,2) * meshgrid(i,3) * dwf1(i,2,1)                     &
      &           - wf2(i,2) * meshgrid(i,2) * dwf1(i,3,1)                     &
      !
      &           - wf2(i,1) * meshgrid(i,3) * dwf1(i,2,2)                     &
      &           + wf2(i,1) * meshgrid(i,2) * dwf1(i,3,2)                     &
      !
      &           + wf2(i,4) * meshgrid(i,3) * dwf1(i,2,3)                     &
      &           - wf2(i,4) * meshgrid(i,2) * dwf1(i,3,3)                     &
      !
      &           - wf2(i,3) * meshgrid(i,3) * dwf1(i,2,4)                     &
      &           + wf2(i,3) * meshgrid(i,2) * dwf1(i,3,4)
    enddo
    angmom = angmom * dv

  end function angmom_x_real

  function angmom_x_imag(wf2, wf1, dwf1) result(angmom)
    !---------------------------------------------------------------------------
    ! calculate the matrix element
    !
    !    Im < wf2 | j_x | wf1 >
    !
    ! Note that this routine assumes that all checks regarding symmetries have 
    ! been performed.
    !
    !       J_x = 1/2*( 0  1 ) + i z \partial_y - i y\partial_z
    !                 ( 1  0 )
    !---------------------------------------------------------------------------

    real(KIND=dp), intent(in) :: wf1(mv,4), wf2(mv,4), dwf1(mv,3,4)
    integer                   :: i
    real(KIND=dp)             :: angmom

    angmom = 0    

    do i=1,mv
      ! Spin part
      angmom = angmom       + 0.5*(  wf2(i,1) * wf1(i,4)                       &
      &                            - wf2(i,2) * wf1(i,3)                       & 
      &                            + wf2(i,3) * wf1(i,1)                       &
      &                            - wf2(i,4) * wf1(i,2))                        
      ! Orbital part
      angmom = angmom &
      &           + wf2(i,1) * meshgrid(i,3) * dwf1(i,2,1)                     &
      &           - wf2(i,1) * meshgrid(i,2) * dwf1(i,3,1)                     &
      !
      &           + wf2(i,2) * meshgrid(i,3) * dwf1(i,2,2)                     &
      &           - wf2(i,2) * meshgrid(i,2) * dwf1(i,3,2)                     &
      !
      &           + wf2(i,3) * meshgrid(i,3) * dwf1(i,2,3)                     &
      &           - wf2(i,3) * meshgrid(i,2) * dwf1(i,3,3)                     &
      !
      &           + wf2(i,4) * meshgrid(i,3) * dwf1(i,2,4)                     &
      &           - wf2(i,4) * meshgrid(i,2) * dwf1(i,3,4)
    enddo
    angmom = angmom * dv

  end function angmom_x_imag

  function angmom_xt_real(wf2, wf1, dwf1) result(angmom)
    !---------------------------------------------------------------------------
    ! calculate the matrix element
    !
    !     < wf2 | j_x T | wf1 >
    !
    ! Note that this routine assumes that all checks regarding symmetries have 
    ! been performed.
    !
    !       J_x = 1/2*( 0  1 ) + i z \partial_y - i y\partial_z
    !                 ( 1  0 )
    !---------------------------------------------------------------------------

    real(KIND=dp), intent(in) :: wf1(:,:), wf2(:,:), dwf1(:,:,:)
    integer                   :: i
    real(KIND=dp)             :: angmom

    angmom = 0    
    do i=1,mv
      ! Spin part
      angmom = angmom       + 0.5*(- wf2(i,1) * wf1(i,1)                       &
      &                            + wf2(i,2) * wf1(i,2)                       & 
      &                            + wf2(i,3) * wf1(i,3)                       &
      &                            - wf2(i,4) * wf1(i,4))                        
      ! Orbital part
      angmom = angmom &
      &           + wf2(i,2) * meshgrid(i,3) * dwf1(i,2,3)                     &
      &           - wf2(i,2) * meshgrid(i,2) * dwf1(i,3,3)                     &
      !
      &           + wf2(i,1) * meshgrid(i,3) * dwf1(i,2,4)                     &
      &           - wf2(i,1) * meshgrid(i,2) * dwf1(i,3,4)                     &
      !
      &           - wf2(i,4) * meshgrid(i,3) * dwf1(i,2,1)                     &
      &           + wf2(i,4) * meshgrid(i,2) * dwf1(i,3,1)                     &
      !
      &           - wf2(i,3) * meshgrid(i,3) * dwf1(i,2,2)                     &
      &           + wf2(i,3) * meshgrid(i,2) * dwf1(i,3,2)
    enddo
    angmom = angmom * dv
  end function angmom_xt_real

  function angmom_x_quad(wf2, dwf2, wf1, dwf1) result(angmom)
    !---------------------------------------------------------------------------
    ! calculate the matrix element
    !
    !     < wf2 | j^\dagger_x j_x | wf1 >
    !
    ! Note that this routine assumes that all checks regarding symmetries have 
    ! been performed.
    !
    !       J_x = 1/2*( 0  1 ) + i z \partial_y - i y\partial_z
    !                 ( 1  0 )
    !---------------------------------------------------------------------------
    
    real(KIND=dp), intent(in) :: wf1(:,:), wf2(:,:), dwf1(:,:,:), dwf2(:,:,:)
    integer                   :: i
    real(KIND=dp)             :: angmom, l1,l2,l3,l4, r1,r2,r3,r4

    angmom = 0    
    do i=1,mv

       ! Action of J_x to the right
       r1 = - meshgrid(i,3)* dwf1(i,2,2) + meshgrid(i,2)* dwf1(i,3,2)          & 
       &    + 0.5_dp       * wf1(i,3)
       r2 = + meshgrid(i,3)* dwf1(i,2,1) - meshgrid(i,2)* dwf1(i,3,1)          & 
       &    + 0.5_dp       * wf1(i,4)
       r3 = - meshgrid(i,3)* dwf1(i,2,4) + meshgrid(i,2)* dwf1(i,3,4)          & 
       &    + 0.5_dp       * wf1(i,1)
       r4 = + meshgrid(i,3)* dwf1(i,2,3) - meshgrid(i,2)* dwf1(i,3,3)          & 
       &    + 0.5_dp       * wf1(i,2)
            
       ! Action of J_x to the left
       l1 = - meshgrid(i,3)* dwf2(i,2,2) + meshgrid(i,2)* dwf2(i,3,2)          & 
       &    + 0.5_dp       * wf2(i,3)
       l2 = + meshgrid(i,3)* dwf2(i,2,1) - meshgrid(i,2)* dwf2(i,3,1)          & 
       &    + 0.5_dp       * wf2(i,4)
       l3 = - meshgrid(i,3)* dwf2(i,2,4) + meshgrid(i,2)* dwf2(i,3,4)          & 
       &    + 0.5_dp       * wf2(i,1)
       l4 = + meshgrid(i,3)* dwf2(i,2,3) - meshgrid(i,2)* dwf2(i,3,3)          & 
       &    + 0.5_dp       * wf2(i,2)
            
        angmom = angmom + l1*r1 + l2*r2 + l3*r3 + l4*r4
    enddo
    angmom = angmom * dv
  end function angmom_x_quad

  function angmom_y_real(wf2, wf1, dwf1) result(angmom)
    !---------------------------------------------------------------------------
    ! calculate the matrix element
    !
    !     Re < wf2 | j_y | wf1 >
    !
    ! Note that this routine assumes that all checks regarding symmetries have 
    ! been performed.
    !
    !       J_z = 1/2*( 0 -i ) + i x \partial_z - i z\partial_x
    !                 ( i  0 )
    !---------------------------------------------------------------------------

    real(KIND=dp), intent(in) :: wf1(mv,4), wf2(mv,4), dwf1(mv,3,4)
    integer                   :: i 
    real(KIND=dp)             :: angmom

    angmom = 0
    do i=1,mv
      ! Spin part
      angmom = angmom       + 0.5*(  wf2(i,1) * wf1(i,4)                       &
      &                            - wf2(i,2) * wf1(i,3)                       & 
      &                            - wf2(i,3) * wf1(i,2)                       &
      &                            + wf2(i,4) * wf1(i,1))                        
      ! Orbital part
      angmom = angmom                                                          &
      &           + wf2(i,2) * meshgrid(i,1) * dwf1(i,3,1)                     &
      &           - wf2(i,2) * meshgrid(i,3) * dwf1(i,1,1)                     &
      !
      &           - wf2(i,1) * meshgrid(i,1) * dwf1(i,3,2)                     &
      &           + wf2(i,1) * meshgrid(i,3) * dwf1(i,1,2)                     &
      !
      &           + wf2(i,4) * meshgrid(i,1) * dwf1(i,3,3)                     &
      &           - wf2(i,4) * meshgrid(i,3) * dwf1(i,1,3)                     &
      !
      &           - wf2(i,3) * meshgrid(i,1) * dwf1(i,3,4)                     &
      &           + wf2(i,3) * meshgrid(i,3) * dwf1(i,1,4) 
     !-------------------------------------------------------------------------
    enddo
    angmom = angmom * dv

  end function angmom_y_real

  function angmom_y_imag(wf2, wf1, dwf1) result(angmom)
    !---------------------------------------------------------------------------
    ! calculate the matrix element
    !
    !     Im < wf2 | j_y | wf1 >
    !
    ! Note that this routine assumes that all checks regarding symmetries have 
    ! been performed.
    !
    !       J_z = 1/2*( 0 -i ) + i x \partial_z - i z\partial_x
    !                 ( i  0 )
    !---------------------------------------------------------------------------

    real(KIND=dp), intent(in) :: wf1(mv,4), wf2(mv,4), dwf1(mv,3,4)
    integer                   :: i
    real(KIND=dp)             :: angmom

    angmom = 0
    do i=1,mv
      ! Spin part
      angmom = angmom       + 0.5*(  wf2(i,1) * wf1(i,4)                       &
      &                            - wf2(i,2) * wf1(i,3)                       & 
      &                            + wf2(i,3) * wf1(i,1)                       &
      &                            - wf2(i,4) * wf1(i,2))                        
      ! Orbital part
      angmom = angmom &
      &           + wf2(i,1) * meshgrid(i,3) * dwf1(i,2,1)                     &
      &           - wf2(i,1) * meshgrid(i,2) * dwf1(i,3,1)                     &
      !
      &           + wf2(i,2) * meshgrid(i,3) * dwf1(i,2,2)                     &
      &           - wf2(i,2) * meshgrid(i,2) * dwf1(i,3,2)                     &
      !
      &           + wf2(i,3) * meshgrid(i,3) * dwf1(i,2,3)                     &
      &           - wf2(i,3) * meshgrid(i,2) * dwf1(i,3,3)                     &
      !
      &           + wf2(i,4) * meshgrid(i,3) * dwf1(i,2,4)                     &
      &           - wf2(i,4) * meshgrid(i,2) * dwf1(i,3,4)
    enddo
    angmom = angmom * dv

  end function angmom_y_imag

  function angmom_yt_imag(wf2, wf1, dwf1) result(angmom)
    !---------------------------------------------------------------------------
    ! calculate the matrix element
    !
    !     < wf2 | j_yT | wf1 >
    !
    ! Note that this routine assumes that all checks regarding symmetries have 
    ! been performed.
    !
    !       J_y = 1/2*( 0 -i ) + i x \partial_z - i z\partial_x
    !                 ( i  0 )
    !---------------------------------------------------------------------------

    real(KIND=dp), intent(in) :: wf1(mv,4), wf2(mv,4), dwf1(mv,3,4)
    integer                   :: i
    real(KIND=dp)             :: angmom

    angmom = 0
    do i=1,mv
      ! Spin part
      angmom = angmom       + 0.5*(- wf2(i,1) * wf1(i,1)                       &
      &                            + wf2(i,2) * wf1(i,2)                       & 
      &                            - wf2(i,3) * wf1(i,3)                       &
      &                            + wf2(i,4) * wf1(i,4))                        
      ! Orbital part
      angmom = angmom                                                          &
      &           + wf2(i,1) * meshgrid(i,3) * dwf1(i,1,3)                     &
      &           - wf2(i,1) * meshgrid(i,1) * dwf1(i,3,3)                     &
      !
      &           - wf2(i,2) * meshgrid(i,3) * dwf1(i,1,4)                     &
      &           + wf2(i,2) * meshgrid(i,1) * dwf1(i,3,4)                     &
      !
      &           - wf2(i,3) * meshgrid(i,3) * dwf1(i,1,1)                     &
      &           + wf2(i,3) * meshgrid(i,1) * dwf1(i,3,1)                     &
      !
      &           + wf2(i,4) * meshgrid(i,3) * dwf1(i,1,2)                     &
      &           - wf2(i,4) * meshgrid(i,1) * dwf1(i,3,2) 
     !-------------------------------------------------------------------------
    enddo
    angmom = angmom * dv

  end function angmom_yt_imag

  function angmom_y_quad(wf2, dwf2, wf1, dwf1) result(angmom)
    !---------------------------------------------------------------------------
    ! calculate the matrix element
    !
    !     < wf2 | j^\dagger_y j_y | wf1 >
    !
    ! Note that this routine assumes that all checks regarding symmetries have 
    ! been performed.
    !
    !       J_y = 1/2*( 0 -i ) + i x \partial_z - i z\partial_x
    !                 ( i  0 )    
    !---------------------------------------------------------------------------
    
    real(KIND=dp), intent(in) :: wf1(:,:), wf2(:,:), dwf1(:,:,:), dwf2(:,:,:)
    integer                   :: i
    real(KIND=dp)             :: angmom, l1,l2,l3,l4, r1,r2,r3,r4

    angmom = 0    
    do i=1,mv
       ! Action of J_y to the right
       r1 = - meshgrid(i,1)* dwf1(i,3,2) + meshgrid(i,3)* dwf1(i,1,2)          & 
       &    + 0.5_dp       * wf1(i,4)
       r2 = + meshgrid(i,1)* dwf1(i,3,1) - meshgrid(i,3)* dwf1(i,1,1)          & 
       &    - 0.5_dp       * wf1(i,3)
       r3 = - meshgrid(i,1)* dwf1(i,3,4) + meshgrid(i,3)* dwf1(i,1,4)          & 
       &    - 0.5_dp       * wf1(i,2)
       r4 = + meshgrid(i,1)* dwf1(i,3,3) - meshgrid(i,3)* dwf1(i,1,3)          & 
       &    + 0.5_dp       * wf1(i,1)
            
       ! Action of J_y to the left
       l1 = - meshgrid(i,1)* dwf2(i,3,2) + meshgrid(i,3)* dwf2(i,1,2)          & 
       &    + 0.5_dp       * wf2(i,4)
       l2 = + meshgrid(i,1)* dwf2(i,3,1) - meshgrid(i,3)* dwf2(i,1,1)          & 
       &    - 0.5_dp       * wf2(i,3)
       l3 = - meshgrid(i,1)* dwf2(i,3,4) + meshgrid(i,3)* dwf2(i,1,4)          & 
       &    - 0.5_dp       * wf2(i,2)
       l4 = + meshgrid(i,1)* dwf2(i,3,3) - meshgrid(i,3)* dwf2(i,1,3)          & 
       &    + 0.5_dp       * wf2(i,1)
        angmom = angmom + l1*r1 + l2*r2 + l3*r3 + l4*r4
    enddo
    angmom = angmom * dv
  end function angmom_y_quad

  function angmom_z_real(wf2, wf1, dwf1) result(angmom)
    !---------------------------------------------------------------------------
    ! calculate the matrix element
    !
    !     < wf2 | j_z | wf1 >
    !
    ! Note that this routine assumes that all checks regarding symmetries have 
    ! been performed.
    !
    !       J_z = 1/2*( 1  0 ) + i y \partial_x - i x\partial_y
    !                 ( 0 -1 ) 
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in) :: wf1(mv,4), wf2(mv,4), dwf1(mv,3,4)
    real(KIND=dp)             :: angmom
    integer                   :: i

    angmom = 0

    do i=1,nx*ny*nz
      !-------------------------------------------------------------------------
      ! Real part  
      ! Spin part
      angmom = angmom    + 0.5*(         wf2(i,1) * wf1(i,1)                   &
      &                                + wf2(i,2) * wf1(i,2)                   & 
      &                                - wf2(i,3) * wf1(i,3)                   &
      &                                - wf2(i,4) * wf1(i,4))                        
      ! Orbital part
      angmom = angmom              +                                           &
      &                              wf2(i,2) * meshgrid(i,2) * dwf1(i,1,1)    &
      &                            - wf2(i,2) * meshgrid(i,1) * dwf1(i,2,1)    &
      !                 
      &                            - wf2(i,1) * meshgrid(i,2) * dwf1(i,1,2)    &
      &                            + wf2(i,1) * meshgrid(i,1) * dwf1(i,2,2)    &
      !
      &                            + wf2(i,4) * meshgrid(i,2) * dwf1(i,1,3)    &
      &                            - wf2(i,4) * meshgrid(i,1) * dwf1(i,2,3)    &
      !                 
      &                            - wf2(i,3) * meshgrid(i,2) * dwf1(i,1,4)    &
      &                            + wf2(i,3) * meshgrid(i,1) * dwf1(i,2,4)
    enddo

    angmom = angmom * dv

  end function angmom_z_real

  function angmom_z_imag(wf2, wf1, dwf1) result(angmom)
    !---------------------------------------------------------------------------
    ! calculate the matrix element
    !
    !    Im < wf2 | j_z | wf1 >
    !
    ! Note that this routine assumes that all checks regarding symmetries have 
    ! been performed.
    !
    !       J_z = 1/2*( 1  0 ) + i y \partial_x - i x\partial_y
    !                 ( 0 -1 ) 
    !---------------------------------------------------------------------------

    real(KIND=dp), intent(in) :: wf1(mv,4), wf2(mv,4), dwf1(mv,3,4)
    integer                   :: i
    real(KIND=dp)             :: angmom

    angmom = 0    

    do i=1,mv
      !-------------------------------------------------------------------------
      ! Real part  

      ! Spin part
      angmom = angmom       + 0.5*(- wf2(i,2) * wf1(i,1)                       &
      &                            + wf2(i,1) * wf1(i,2)                       & 
      &                            + wf2(i,4) * wf1(i,3)                       &
      &                            - wf2(i,3) * wf1(i,4))                        
      ! Orbital part
      angmom = angmom &
      &           + wf2(i,1) * meshgrid(i,2) * dwf1(i,1,1)                     &
      &           - wf2(i,1) * meshgrid(i,1) * dwf1(i,2,1)                     &
      !
      &           + wf2(i,2) * meshgrid(i,2) * dwf1(i,1,2)                     &
      &           - wf2(i,2) * meshgrid(i,1) * dwf1(i,2,2)                     &
      !
      &           + wf2(i,3) * meshgrid(i,2) * dwf1(i,1,3)                     &
      &           - wf2(i,3) * meshgrid(i,1) * dwf1(i,2,3)                     &
      !
      &           + wf2(i,4) * meshgrid(i,2) * dwf1(i,1,4)                     &
      &           - wf2(i,4) * meshgrid(i,1) * dwf1(i,2,4)
    enddo
    angmom = angmom * dv

  end function angmom_z_imag

  function angmom_zt_real(wf2, wf1, dwf1) result(angmom)
    !---------------------------------------------------------------------------
    ! calculate the matrix element
    !
    !     < wf2 | j_z T | wf1 >
    !
    ! Note that this routine assumes that all checks regarding symmetries have 
    ! been performed.
    !
    !       J_z = 1/2*( 1  0 ) + i y \partial_x - i x\partial_y
    !                 ( 0 -1 ) 
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in) :: wf1(mv,4), wf2(mv,4), dwf1(mv,3,4)
    integer                   :: i 
    real(KIND=dp)             :: angmom

    angmom = 0
    do i=1,mv
      !-------------------------------------------------------------------------
      ! Real part  

      ! Spin part
      angmom = angmom    + 0.5*(         wf2(i,1) * wf1(i,3)                   &
      &                                - wf2(i,2) * wf1(i,4)                   & 
      &                                + wf2(i,3) * wf1(i,1)                   &
      &                                - wf2(i,4) * wf1(i,2))                        
      ! Orbital part
      angmom = angmom                                                          &
      &                            + wf2(i,2) * meshgrid(i,2) * dwf1(i,1,3)    &
      &                            - wf2(i,2) * meshgrid(i,1) * dwf1(i,2,3)    &
      !                 
      &                            + wf2(i,1) * meshgrid(i,2) * dwf1(i,1,4)    &
      &                            - wf2(i,1) * meshgrid(i,1) * dwf1(i,2,4)    &
      !
      &                            - wf2(i,4) * meshgrid(i,2) * dwf1(i,1,1)    &
      &                            + wf2(i,4) * meshgrid(i,1) * dwf1(i,2,1)    &
      !                 
      &                            - wf2(i,3) * meshgrid(i,2) * dwf1(i,1,2)    &
      &                            + wf2(i,3) * meshgrid(i,1) * dwf1(i,2,2)
    enddo
    angmom = dv * angmom
  end function angmom_zt_real

  function angmom_z_quad(wf2, dwf2, wf1, dwf1) result(angmom)
    !---------------------------------------------------------------------------
    ! calculate the matrix element
    !
    !     < wf2 | j^\dagger_z j_z | wf1 >
    !
    ! Note that this routine assumes that all checks regarding symmetries have 
    ! been performed.
    !
    !       J_z = 1/2*( 1  0 ) + i y \partial_x - i x\partial_y
    !                 ( 0 -1 ) 
    !---------------------------------------------------------------------------
    
    real(KIND=dp), intent(in) :: wf1(:,:), wf2(:,:), dwf1(:,:,:), dwf2(:,:,:)
    integer                   :: i
    real(KIND=dp)             :: angmom, l1,l2,l3,l4, r1,r2,r3,r4

    angmom = 0    
    do i=1,mv
       ! Action of J_z to the right
       r1 = - meshgrid(i,2)* dwf1(i,1,2) + meshgrid(i,1)* dwf1(i,2,2)          & 
       &    + 0.5_dp       * wf1(i,1)
       r2 = + meshgrid(i,2)* dwf1(i,1,1) - meshgrid(i,1)* dwf1(i,2,1)          & 
       &    + 0.5_dp       * wf1(i,2)
       r3 = - meshgrid(i,2)* dwf1(i,1,4) + meshgrid(i,1)* dwf1(i,2,4)          & 
       &    - 0.5_dp       * wf1(i,3)
       r4 = + meshgrid(i,2)* dwf1(i,1,3) - meshgrid(i,1)* dwf1(i,2,3)          & 
       &    - 0.5_dp       * wf1(i,4)
            
       ! Action of J_z to the left
       l1 = - meshgrid(i,2)* dwf2(i,1,2) + meshgrid(i,1)* dwf2(i,2,2)          & 
       &    + 0.5_dp       * wf2(i,1)
       l2 = + meshgrid(i,2)* dwf2(i,1,1) - meshgrid(i,1)* dwf2(i,2,1)          & 
       &    + 0.5_dp       * wf2(i,2)
       l3 = - meshgrid(i,2)* dwf2(i,1,4) + meshgrid(i,1)* dwf2(i,2,4)          & 
       &    - 0.5_dp       * wf2(i,3)
       l4 = + meshgrid(i,2)* dwf2(i,1,3) - meshgrid(i,1)* dwf2(i,2,3)          & 
       &    - 0.5_dp       * wf2(i,4)
       
       angmom = angmom + l1*r1 + l2*r2 + l3*r3 + l4*r4
    enddo
    angmom = angmom * dv
  end function angmom_z_quad

  subroutine clean_wavefunctions()

    if(allocated(HFPsi))    deallocate(HFPsi)
    if(allocated(HFdPsi))   deallocate(HFdPsi)
    if(allocated(HFddPsi))  deallocate(HFddPsi)
    if(allocated(HFdddPsi)) deallocate(HFdddPsi)

    if(allocated(CANPsi))    deallocate(CANPsi)
    if(allocated(CANdPsi))   deallocate(CANdPsi)
    if(allocated(CANddPsi))  deallocate(CANddPsi)
    if(allocated(CANdddPsi)) deallocate(CANdddPsi)

    if(allocated(spenergies))  deallocate(spenergies)
    if(allocated(dispersions)) deallocate(dispersions)
    if(allocated(canenergies)) deallocate(canenergies)

    if(allocated(sx)) deallocate(sx,sy,sz)

  end subroutine clean_wavefunctions

end module wavefunctions
