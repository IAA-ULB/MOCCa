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
 !==============================================================================
 ! Hephaestos:
 ! 
 !    N2  : $N2
 !    N3  : $N3
 !
 ! ININX  : $ININX
 ! ININY  : $ININY
 ! ININZ  : $ININZ
 !
 ! ININWN : $ININWN
 ! ININWP : $ININWP
 ! ININWT : $ININWT
 ! 
 !==============================================================================
 use compilation
 use derivatives
 use geninfo
 use nil8
 use timing 

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
 ! Single-particle energies, diagonal elements of the single-particle
 ! hamiltonian: \langle psi_i | h | psi_i \rangle
 real(KIND=dp), allocatable :: spenergies(:) 
 real(KIND=dp), allocatable :: current_sph(:,:)
 ! Dispersions of the spwfs with respect to h
 real(KIND=dp), allocatable :: dispersions(:)
 ! expectation values of the single-particle hamiltonian in the canonical basis
 real(KIND=dp), allocatable :: canenergies(:)
 !------------------------------------------------------------------------------
 ! Angular momentum properties of the spwfs in both the HF and canonical basis
 ! "Ordinary" <Jx>, <Jy>, <Jz> in the HF and canonical basis
 real(KIND=dp), allocatable :: spwf_J(:,:), can_J(:,:)
 ! Squared   <Jx^2>, <Jy^2>, <Jz^2>
 real(KIND=dp), allocatable :: spwf_J2(:,:), can_J2(:,:)
 ! With an extra time-reversal operator < Jx T >, < Jy T >, < Jz T >
 ! Both real and imaginary parts
 real(KIND=dp), allocatable :: spwf_JTR(:,:), can_JTR(:,:)
 real(KIND=dp), allocatable :: spwf_JTI(:,:), can_JTI(:,:)
 ! Total angular momentum "quantum number", i.e. the number J such that 
 !  J (J+1) = <J^2_x> +  <J^2_y> + <J^2_z>
 real(KIND=dp), allocatable :: spwf_JJ(:), can_JJ(:)
 !------------------------------------------------------------------------------
 ! Number of the blocks with the same quantum numbers that divide up the 
 ! the single-particle wavefunctions.
 ! 
 ! Any of the possible symmetry combinations give rise to at most eight 
 ! different symmetry blocks.
 ! 
 !   2 for protons <-> neutrons
 !   2 for a hermitian, linear operator      ( parity      in EV8) 
 !   2 for an antihermitian, linear operator ( z-signature in EV8)
 ! 
 ! If symmetries are not conserved, we can simply eliminate blocks by setting 
 ! their size to zero.
 !------------------------------------------------------------------------------
 integer, parameter   :: Blocks          =8 ! This can always be fixed
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
 !------------------------------------------------------------------------------

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
    ! Build harmonic oscillator eigenfunctions in an EV8-like box
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Also initialized:
    !  *) Diagonal matrix elements of <h>
    ! 
    ! Not initialized here:
    !  *) Delta for the gaps. Since this module can not know what kind of 
    !     pairing is needed, it cannot correctly guess a structure. 
    !---------------------------------------------------------------------------
    
!    real(KIND=dp)             :: homegax, homegay,homegaz, alpha,qqq
    integer                   :: i
    integer, allocatable      :: kparz(:)
    
!    alpha = 0.5    
!    qqq   = 2.0  
!    homegaz  = alpha*qqq**(-2.0/3.0)
!    homegax  = alpha*qqq**(-2*cos(-2*pi/3)/3)
!    homegay  = alpha*qqq**(-2*cos(+2*pi/3)/3)

    allocate(hfpsi($ININX*$ININY*$ININZ,4,$ININWT)) ; hfpsi = 0.0d0
    if (allocated(kparz))  deallocate(kparz)       
   
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
    ! a) Generating the nilsson wave-functions in an EV8-box   
    call nilsson (HFPsi,kparz,spenergies,8,7,$ININWT,$ININWP,$ININWN,          &
    &           floor(neutrons),floor(protons),$ININX,$ININY,$ININZ,dx,osc_freq)


    allocate(dispersions($ININWT)) ; dispersions  = 0
    allocate(sx(4,$ININWT), sy(4,$ININWT), sz(4,$ININWT))

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
    ! b) fill in the right symmetry properties for the wavefunctions
    hfblocks = 0
    do i=1,$ININWN
        if(kparz(i) .gt. 0) HFBlocks(1) = HFBlocks(1) +1
        if(kparz(i) .lt. 0) HFBlocks(3) = HFBlocks(3) +1
    enddo
    do i=$ININWN+1,$ININWT
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
    deallocate(kparz)

    call GramSchmidt
  end subroutine iniwavefunctions
  
  subroutine deriveHF()
    !---------------------------------------------------------------------------
    ! Derives all of the single-particle wave-functions. 
    ! a) In the HF basis
    !---------------------------------------------------------------------------
    integer :: wave,k
    
    call start_timer(T_derivatives)

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
    call stop_timer(T_derivatives)
  end subroutine DeriveHF
  
  subroutine deriveCan()
    !---------------------------------------------------------------------------
    ! Derives all of the single-particle wave-functions. 
    ! b) In the canonical basis
    !---------------------------------------------------------------------------
    integer :: wave,k
      
    call start_timer(T_derivatives_can)

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

    call stop_timer(T_derivatives_can)
    
  end subroutine DeriveCan
  
  function OrderSpwfsISO(Isospin, canonical) result(Indices)
    !---------------------------------------------------------------------------
    ! Orders the wavefunctions within an isospin block. 
    !---------------------------------------------------------------------------
    integer, intent(in)        :: Isospin
    
    integer, allocatable       :: Indices(:)
    real(Kind=dp), allocatable :: Energies(:)
    integer                    :: i, nwf,  HolePos, ToInsertIndex
    real(Kind=dp)              :: ToInsert
    logical, intent(in), optional :: canonical
    
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !Count the number of proton and neutron wavefunctions
    if (Isospin .eq. -1) then
        nwf = nwn
    else
        nwf = nwp
    endif
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !Filling Energies & Indices
    if(allocated(indices))  deallocate(indices)
    if(allocated(Energies)) deallocate(Energies)
    allocate(Indices(nwf), Energies(nwf))
    do i=1,nwf
       Indices(i) = i 
    enddo

    if(.not. present(canonical)) then
      if(Isospin.eq.-1) then
          Energies = spenergies(1:nwn)
      else
          Indices  = Indices + nwn
          Energies = spenergies(nwn+1:nwn+nwp)
      endif
    elseif(canonical) then
      if(Isospin.eq.-1) then
          Energies = canenergies(1:nwn)
      else
          Indices  = Indices + nwn
          Energies = canenergies(nwn+1:nwn+nwp)
      endif
    else
      if(Isospin.eq.-1) then
          Energies = spenergies(1:nwn)
      else
          Indices  = Indices + nwn
          Energies = spenergies(nwn+1:nwn+nwp)
      endif
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

    deallocate(Energies)
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
    if(allocated(Indices))   deallocate(indices)
    if(allocated(Energies))  deallocate(energies)
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
    deallocate(energies)
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
    
    call start_timer(T_ortho)

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
            !$$OMP PARALLEL private(j, mw, norm, l)
            !$$OMP DO
            do j= i+1, HFBlocks(b)
                mw = indices(j)    
                ! Real part of the inproduct
                norm = sum(HFpsi(:,:,nw)*HFpsi(:,:,mw)) * dv
                do l=1,4*nx*ny*nz
                    HFPsi(l,1,mw) = HFPsi(l,1,mw) - norm * HFPsi(l,1,nw)
                enddo
            enddo
           !$$OMP END DO
           !$$OMP END PARALLEL 

        enddo
    enddo

    call stop_timer(T_ortho)

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

  subroutine update_spwf_angmom()
    !---------------------------------------------------------------------------
    ! Calculate all the angular momentum properties of the spwfs.
    ! For now, we only calculate properties that would be accessible in a 
    ! CR8-geometry.
    !---------------------------------------------------------------------------
    integer       :: wave 

    if(.not.allocated(spwf_J)) then
      allocate(spwf_J(3,nwt))   ; spwf_J = 0.0
      allocate(spwf_JTR(3,nwt)) ; spwf_JTR= 0.0
      allocate(spwf_JTI(3,nwt)) ; spwf_JTI= 0.0
      allocate(spwf_J2(3,nwt))  ; spwf_J2= 0.0
      allocate(spwf_JJ(nwt))    ; spwf_JJ= 0.0
    endif

    do wave=1,nwt
      spwf_JTR(1,wave) = & 
            & angmom_xt_real(HFPsi(:,:,wave),HFPsi(:,:,wave),HFdPsi(:,:,:,wave))
      spwf_JTI(2,wave) = &
            & angmom_yt_imag(HFPsi(:,:,wave),HFPsi(:,:,wave),HFdPsi(:,:,:,wave))
      spwf_J(3,wave)   = & 
            & angmom_z_real (HFPsi(:,:,wave),HFPsi(:,:,wave),HFdPsi(:,:,:,wave))

      spwf_J2(1,wave)  = &
        &   angmom_x_quad(HFPsi(:,:,wave),HFdPsi(:,:,:,wave), &
        &                 HFPsi(:,:,wave),HFdPsi(:,:,:,wave)) 
      spwf_J2(2,wave)  = &
        &   angmom_y_quad(HFPsi(:,:,wave),HFdPsi(:,:,:,wave), &
        &                 HFPsi(:,:,wave),HFdPsi(:,:,:,wave)) 
      spwf_J2(3,wave)  = &
        &   angmom_z_quad(HFPsi(:,:,wave),HFdPsi(:,:,:,wave), &
        &                 HFPsi(:,:,wave),HFdPsi(:,:,:,wave)) 
  
      spwf_JJ(wave) = (-1. + sqrt(1. + 4*sum(spwf_J2(:,wave))))/2.
    enddo

    if(allocated(CANPSI)) then
      if(.not.allocated(can_J)) then
        allocate(can_J(3,nwt))   ; can_J  = 0.0
        allocate(can_JTR(3,nwt)) ; can_JTR= 0.0
        allocate(can_JTI(3,nwt)) ; can_JTI= 0.0
        allocate(can_J2(3,nwt))  ; can_J2 = 0.0
        allocate(can_JJ(nwt))    ; can_JJ = 0.0
      endif

      do wave=1,nwt
        can_JTR(1,wave) = & 
        & angmom_xt_real(CanPsi(:,:,wave),CanPsi(:,:,wave),CanDPsi(:,:,:,wave))
        can_JTI(2,wave) = &
        & angmom_yt_imag(CanPsi(:,:,wave),CanPsi(:,:,wave),CanDPsi(:,:,:,wave))
        can_J(3,wave)   = & 
        & angmom_z_real (CanPsi(:,:,wave),CanPsi(:,:,wave),CanDPsi(:,:,:,wave))

        can_J2(1,wave)  = &
          &   angmom_x_quad(CanPsi(:,:,wave),CandPsi(:,:,:,wave), &
          &                 CanPsi(:,:,wave),CandPsi(:,:,:,wave)) 
        can_J2(2,wave)  = &
          &   angmom_y_quad(CanPsi(:,:,wave),CandPsi(:,:,:,wave), &
          &                 CanPsi(:,:,wave),CandPsi(:,:,:,wave)) 
        can_J2(3,wave)  = &
          &   angmom_z_quad(CanPsi(:,:,wave),CandPsi(:,:,:,wave), &
          &                 CanPsi(:,:,wave),CandPsi(:,:,:,wave)) 
    
        can_JJ(wave) = (-1. + sqrt(1. + 4*sum(can_J2(:,wave))))/2.
      enddo
    endif

  end subroutine update_spwf_angmom
  
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
      angmom = angmom       + 0.5*(- wf2(i,1) * wf1(i,3)                       &
      &                            - wf2(i,2) * wf1(i,4)                       & 
      &                            + wf2(i,3) * wf1(i,1)                       &
      &                            + wf2(i,4) * wf1(i,2))                        
      ! Orbital part
      angmom = angmom &
      &           - wf2(i,1) * meshgrid(i,3) * dwf1(i,1,1)                     &
      &           + wf2(i,1) * meshgrid(i,1) * dwf1(i,3,1)                     &
      !
      &           - wf2(i,2) * meshgrid(i,3) * dwf1(i,1,2)                     &
      &           + wf2(i,2) * meshgrid(i,1) * dwf1(i,3,2)                     &
      !
      &           - wf2(i,3) * meshgrid(i,3) * dwf1(i,1,3)                     &
      &           + wf2(i,3) * meshgrid(i,1) * dwf1(i,3,3)                     &
      !
      &           - wf2(i,4) * meshgrid(i,3) * dwf1(i,1,4)                     &
      &           + wf2(i,4) * meshgrid(i,1) * dwf1(i,3,4)
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
      angmom = angmom       + 0.5*(+ wf2(i,1) * wf1(i,1)                       &
      &                            - wf2(i,2) * wf1(i,2)                       & 
      &                            + wf2(i,3) * wf1(i,3)                       &
      &                            - wf2(i,4) * wf1(i,4))                        
      ! Orbital part
      angmom = angmom                                                          &
      &           - wf2(i,1) * meshgrid(i,3) * dwf1(i,1,3)                     &
      &           + wf2(i,1) * meshgrid(i,1) * dwf1(i,3,3)                     &
      !
      &           + wf2(i,2) * meshgrid(i,3) * dwf1(i,1,4)                     &
      &           - wf2(i,2) * meshgrid(i,1) * dwf1(i,3,4)                     &
      !
      &           + wf2(i,3) * meshgrid(i,3) * dwf1(i,1,1)                     &
      &           - wf2(i,3) * meshgrid(i,1) * dwf1(i,3,1)                     &
      !
      &           - wf2(i,4) * meshgrid(i,3) * dwf1(i,1,2)                     &
      &           + wf2(i,4) * meshgrid(i,1) * dwf1(i,3,2) 
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

    if(allocated(sx)) deallocate(sx)
    if(allocated(sy)) deallocate(sy)
    if(allocated(sz)) deallocate(sz)

  end subroutine clean_wavefunctions

end module wavefunctions
