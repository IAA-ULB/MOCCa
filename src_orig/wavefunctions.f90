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
 use nil8
 
 implicit none
 
 !------------------------------------------------------------------------------
 ! Array containing the values of the spwfs in the Hartree-Fock basis
 ! and their derivatives
 ! Dimensions (nx,ny,nz,4,nwt)
 real(KIND=dp), allocatable ::   HFPsi(:,:,:,:,:)
 real(KIND=dp), allocatable ::  HFdPsi(:,:,:,:,:,:)   ! First order derivatives
 real(KIND=dp), allocatable :: HFddPsi(:,:,:,:,:,:,:) ! Second order derivatives
 real(KIND=dp), allocatable ::HFlapPsi(:,:,:,:,:)     ! Second order derivatives
 
 !------------------------------------------------------------------------------
 ! Density matrix rho and anomalous density matrix kappa
 ! Dimensions (nwt, nwt) (although many are zero when symmetries are conserved)
 !real(KIND=dp), allocatable :: rho(:,:), kappa(:,:)
 !------------------------------------------------------------------------------
 ! Occupations of the single-particle wave-functions, i.e. the eigenvalues
 ! of rho. 
 real(KIND=dp), allocatable :: occupations(:)
 !------------------------------------------------------------------------------
 ! Single-particle energies, diagonal elements of the single-particle
 ! hamiltonian
 ! \langle psi_i | h | psi_i \rangle
 real(KIND=dp), allocatable :: spenergies(:)
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
 !
 integer, parameter   :: Blocks          =$BLOCKS
 integer              :: HFBlocks(Blocks)=0
 integer              :: nwn, nwp

contains 

  subroutine iniwavefunctions()   
    !--------------------------------------------------------------------
    ! Build harmonic oscillator eigenfunctions
    ! a) in an EV8-like box
    ! b) expanding to the full box
    ! c) restricting again to the box desired 
    !--------------------------------------------------------------------
    
    real(KIND=dp)             :: homegax, homegay,homegaz, alpha,qqq
    real(KIND=dp),allocatable :: fullbox(:,:,:,:,:)
    integer                   :: i,j,k, wave, p
    integer, allocatable      :: kparz(:)
        
    alpha = 0.2    
    qqq   = 1.0    
    homegaz  = alpha*qqq**(-2.0/3.0)
    homegax  = alpha*qqq**(-2*cos(-2*pi/3)/3)
    homegay  = alpha*qqq**(-2*cos(+2*pi/3)/3)
    
    nwn = 10
    nwp = 10

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! a) Generating the nilsson wave-functions in an EV8-box   
    call nilsson (HFPsi,kparz,spenergies,2,1,nwt,nwn,nwp,                      &
    &           floor(neutrons),floor(protons),nx,ny,nz,0.8d0,0.2d0,0.2d0,0.2d0)
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! b) blow-up into the full box
    allocate(fullbox(2*nx, 2*ny, 2*nz, 4, nwt))
    do wave=1,nwt
        ! Copy the original
        fullbox(nx+1: 2*nx,ny+1: 2*ny, nz+1: 2*nz,:,wave) = HFPsi   (:,:,:,:,wave)     
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        ! Use y-time-simplex to fill in the y-axis.
        do j=1,ny
                ! S^t_y Psi (x,y,z,sigma) = Psi^*(x,-y,z,sigma)
                fullbox(nx+1: 2*nx,j, nz+1:2*nz,1,wave) =   fullbox(nx+1: 2*nx, 2*ny - j +1, nz+1:2*nz,1,wave)
                fullbox(nx+1: 2*nx,j, nz+1:2*nz,2,wave) = - fullbox(nx+1: 2*nx, 2*ny - j +1, nz+1:2*nz,2,wave)
                fullbox(nx+1: 2*nx,j, nz+1:2*nz,3,wave) =   fullbox(nx+1: 2*nx, 2*ny - j +1, nz+1:2*nz,3,wave)
                fullbox(nx+1: 2*nx,j, nz+1:2*nz,4,wave) = - fullbox(nx+1: 2*nx, 2*ny - j +1, nz+1:2*nz,4,wave)
        enddo
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        ! Use z-signature to fill in the x-axis
        do j=1,2*ny
            do i=1,nx
                ! R_z Psi (x,y,z,sigma) = -i sigma Psi(-x,-y,z, sigma)
                fullbox(i,j,nz+1:2*nz,1,wave) =   fullbox(2*nx -i +1, 2*ny - j +1, nz+1:2*nz,2,wave)
                fullbox(i,j,nz+1:2*nz,2,wave) = - fullbox(2*nx -i +1, 2*ny - j +1, nz+1:2*nz,1,wave)
                fullbox(i,j,nz+1:2*nz,3,wave) = - fullbox(2*nx -i +1, 2*ny - j +1, nz+1:2*nz,4,wave)
                fullbox(i,j,nz+1:2*nz,4,wave) =   fullbox(2*nx -i +1, 2*ny - j +1, nz+1:2*nz,3,wave)
            enddo
        enddo
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        ! Use parity to fill in the z-axis
        p = kparz(wave)
        do k=1,nz
                do j=1,2*ny
                        do i=1,2*nx
                                fullbox(i,j,k,1,wave) = p*fullbox(2*nx -i +1, 2*ny - j +1, 2*nz-k+1,1,wave)
                                fullbox(i,j,k,2,wave) = p*fullbox(2*nx -i +1, 2*ny - j +1, 2*nz-k+1,2,wave)
                                fullbox(i,j,k,3,wave) = p*fullbox(2*nx -i +1, 2*ny - j +1, 2*nz-k+1,3,wave)
                                fullbox(i,j,k,4,wave) = p*fullbox(2*nx -i +1, 2*ny - j +1, 2*nz-k+1,4,wave)
                        enddo
                enddo
        enddo
    enddo
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    do i=1,10
        if(kparz(i) .gt. 0) HFBlocks(1) = HFBlocks(1) +1
        if(kparz(i) .lt. 0) HFBlocks(3) = HFBlocks(3) +1
    enddo
    do i=11,20
        if(kparz(i) .gt. 0) HFBlocks(5) = HFBlocks(5) +1
        if(kparz(i) .lt. 0) HFBlocks(7) = HFBlocks(7) +1
    enddo
    

  end subroutine iniwavefunctions
  
  subroutine deriveall()
    !---------------------------------------------------------------------------
    ! Derives all of the single-particle wave-functions. (For now in the HFbasis)
    !---------------------------------------------------------------------------
    integer :: wave
    
    if(.not.allocated(HFdPsi)) then
        allocate(HFdPsi(nx,ny,nz,4,3,nwt))
        allocate(HFddPsi(nx,ny,nz,4,3,3,nwt))
    endif
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Currently EV8 symmetries are hardcoded, as well as the
    do wave=1,HFBlocks(1)
        call Derive_tot(HFPsi(:,:,:,1,wave), +1, +1, +1,                           &
        &                                              HFdPsi(:,:,:,1,:,wave),     &
        &                                           HFddPsi(:,:,:,1,:,:,wave), 1)
        call Derive_tot(HFPsi(:,:,:,2,wave), -1, -1, +1,                           &
        &                                              HFdPsi(:,:,:,2,:,wave),     &
        &                                           HFddPsi(:,:,:,2,:,:,wave), 1)
        call Derive_tot(HFPsi(:,:,:,3,wave), -1, +1, -1,                           &
        &                                              HFdPsi(:,:,:,3,:,wave),     &
        &                                           HFddPsi(:,:,:,3,:,:,wave), 1)
        call Derive_tot(HFPsi(:,:,:,4,wave), +1, -1, -1,                           &
        &                                              HFdPsi(:,:,:,4,:,wave),     &
        &                                           HFddPsi(:,:,:,4,:,:,wave), 1)
    enddo
    
    do wave=HFBlocks(1) + 1,HFBlocks(1) + HFBlocks(3)
        call Derive_tot(HFPsi(:,:,:,1,wave), +1, +1, -1,                           &
        &                                              HFdPsi(:,:,:,1,:,wave),     &
        &                                           HFddPsi(:,:,:,1,:,:,wave), 1)
        call Derive_tot(HFPsi(:,:,:,2,wave), -1, -1, -1,                           &
        &                                              HFdPsi(:,:,:,2,:,wave),     &
        &                                           HFddPsi(:,:,:,2,:,:,wave), 1)
        call Derive_tot(HFPsi(:,:,:,3,wave), -1, +1, +1,                           &
        &                                              HFdPsi(:,:,:,3,:,wave),     &
        &                                           HFddPsi(:,:,:,3,:,:,wave), 1)
        call Derive_tot(HFPsi(:,:,:,4,wave), +1, -1, +1,                           &
        &                                              HFdPsi(:,:,:,4,:,wave),     &
        &                                           HFddPsi(:,:,:,4,:,:,wave), 1)
    enddo

    do wave=HFBlocks(1) + HFBlocks(3)+1,HFBlocks(1) + HFBlocks(3)+HFBlocks(5)
        call Derive_tot(HFPsi(:,:,:,1,wave), +1, +1, +1,                           &
        &                                              HFdPsi(:,:,:,1,:,wave),     &
        &                                           HFddPsi(:,:,:,1,:,:,wave), 1)
        call Derive_tot(HFPsi(:,:,:,2,wave), -1, -1, +1,                           &
        &                                              HFdPsi(:,:,:,2,:,wave),     &
        &                                           HFddPsi(:,:,:,2,:,:,wave), 1)
        call Derive_tot(HFPsi(:,:,:,3,wave), -1, +1, -1,                           &
        &                                              HFdPsi(:,:,:,3,:,wave),     &
        &                                           HFddPsi(:,:,:,3,:,:,wave), 1)
        call Derive_tot(HFPsi(:,:,:,4,wave), +1, -1, -1,                           &
        &                                              HFdPsi(:,:,:,4,:,wave),     &
        &                                           HFddPsi(:,:,:,4,:,:,wave), 1)
    enddo
!   
    do wave=HFBlocks(1)+HFBlocks(3)+HFBlocks(5) + 1,                           &
    &       HFBlocks(1)+HFBlocks(3)+HFBlocks(5) + HFBLocks(7)
        call Derive_tot(HFPsi(:,:,:,1,wave), +1, +1, -1,                           &
        &                                              HFdPsi(:,:,:,1,:,wave),     &
        &                                           HFddPsi(:,:,:,1,:,:,wave), 1)
        call Derive_tot(HFPsi(:,:,:,2,wave), -1, -1, -1,                           &
        &                                              HFdPsi(:,:,:,2,:,wave),     &
        &                                           HFddPsi(:,:,:,2,:,:,wave), 1)
        call Derive_tot(HFPsi(:,:,:,3,wave), -1, +1, +1,                           &
        &                                              HFdPsi(:,:,:,3,:,wave),     &
        &                                           HFddPsi(:,:,:,3,:,:,wave), 1)
        call Derive_tot(HFPsi(:,:,:,4,wave), +1, -1, +1,                           &
        &                                              HFdPsi(:,:,:,4,:,wave),     &
        &                                           HFddPsi(:,:,:,4,:,:,wave), 1)
    enddo
  end subroutine DeriveAll
  
  subroutine PrintSpwfs
    !---------------------------------------------------------------------------
    ! Print the info on the single-particle wave-functions in the HFBasis.
    !---------------------------------------------------------------------------
    
    10 format (21 ('-'), ' Sp wavefunctions ', 41('-'))
    20 format (94 ('-'))
    30 format (94 ('_'),/,3x , 'Neutron wavefunctions')
    40 format (94 ('_'),/,3x , 'Proton  wavefunctions')
    50 format (94 ('_'),/,3x , 'HF Basis')
    
    integer :: wave,k
    integer :: ProtonOrder(nwp), NeutronOrder(nwn)
    
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Order the spwfs according to growing energy
    ProtonOrder = OrderSpwfsISO(+1)
    NeutronOrder= OrderSpwfsISO(-1)
    
    print 10
    print 50
    print 30
    
    do k=1,nwn 
        wave = NeutronOrder(k)
        print *, wave, occupations(wave), spenergies(wave)
    enddo
    
    print 40  
    
    do k=1,nwp
        wave = ProtonOrder(k)
        print *, wave, occupations(wave), spenergies(wave)
    enddo
    print 20
  end subroutine PrintSpwfs
  
  function OrderSpwfsISO(Isospin) result(Indices)
    !---------------------------------------------------------------------------
    ! Orders the wavefunctions, but keeps the neutrons and protons separate.
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
  
end module wavefunctions
