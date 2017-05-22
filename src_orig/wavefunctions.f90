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
    !
    !
    !--------------------------------------------------------------------   
    
    real(KIND=dp)        :: homegax, homegay,homegaz, alpha,qqq
    integer              :: i
    integer, allocatable :: kparz(:)
    !--------------------------------------------------------------------
    ! Build harmonic oscillator eigenfunctions by constructing them in  
    ! an EV8-like box and then expanding them to the entire box. 
    !--------------------------------------------------------------------
        
    alpha = 0.2    
    qqq   = 1.0    
    homegaz  = alpha*qqq**(-2.0/3.0)
    homegax  = alpha*qqq**(-2*cos(-2*pi/3)/3)
    homegay  = alpha*qqq**(-2*cos(+2*pi/3)/3)
    
    nwn = 10
    nwp = 10
    
    call nilsson (HFPsi,kparz,spenergies,2,1,nwt,nwn,nwp,                      &
    &           floor(neutrons),floor(protons),nx,ny,nz,0.8d0,0.2d0,0.2d0,0.2d0)
  
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
        allocate(HFlapPsi(nx,ny,nz,4,nwt))
    endif
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Currently EV8 symmetries are hardcoded.
    do wave=1,HFBlocks(1)
        call Derive(HFPsi(:,:,:,1,wave), +1, +1, +1,                           &
        &                                              HFdPsi(:,:,:,1,1,wave), &
        &                                              HFdPsi(:,:,:,1,2,wave), &
        &                                              HFdPsi(:,:,:,1,3,wave), &
        &                                              HFlapPsi(:,:,:,1,wave))
        call Derive(HFPsi(:,:,:,2,wave), -1, -1, +1,                           &
        &                                              HFdPsi(:,:,:,2,1,wave), &
        &                                              HFdPsi(:,:,:,2,2,wave), &
        &                                              HFdPsi(:,:,:,2,3,wave), &
        &                                              HFlapPsi(:,:,:,2,wave))
        call Derive(HFPsi(:,:,:,3,wave), -1, +1, -1,                           &
        &                                              HFdPsi(:,:,:,3,1,wave), &
        &                                              HFdPsi(:,:,:,3,2,wave), &
        &                                              HFdPsi(:,:,:,3,3,wave), &
        &                                              HFlapPsi(:,:,:,3,wave))
        call Derive(HFPsi(:,:,:,4,wave), +1, -1, -1,                           &
        &                                              HFdPsi(:,:,:,4,1,wave), &
        &                                              HFdPsi(:,:,:,4,2,wave), &
        &                                              HFdPsi(:,:,:,4,3,wave), &
        &                                              HFlapPsi(:,:,:,4,wave))
    enddo
    
!    do wave=1,4
!        print *, HFPsi(1:nx,1,1,wave,1)
!        print *
!        print *, HFdPsi(1:nx,1,1,wave,1,1)
!        print *, HFdPsi(1,1:ny,1,wave,2,1)
!        print *, HFdPsi(1,1,1:nz,wave,3,1)
!        print *
!        print * ,'  - - - - - - - - - - - - - - - - - '
!    enddo
!    stop
    
    do wave=HFBlocks(1) + 1,HFBlocks(1) + HFBlocks(3)
        call Derive(HFPsi(:,:,:,1,wave), +1, +1, -1,                           &
        &                                              HFdPsi(:,:,:,1,1,wave), &
        &                                              HFdPsi(:,:,:,1,2,wave), &
        &                                              HFdPsi(:,:,:,1,3,wave), &
        &                                              HFlapPsi(:,:,:,1,wave))
        call Derive(HFPsi(:,:,:,2,wave), -1, -1, -1,                           &
        &                                              HFdPsi(:,:,:,2,1,wave), &
        &                                              HFdPsi(:,:,:,2,2,wave), &
        &                                              HFdPsi(:,:,:,2,3,wave), &
        &                                              HFlapPsi(:,:,:,2,wave))
        call Derive(HFPsi(:,:,:,3,wave), -1, +1, +1,                           &
        &                                              HFdPsi(:,:,:,3,1,wave), &
        &                                              HFdPsi(:,:,:,3,2,wave), &
        &                                              HFdPsi(:,:,:,3,3,wave), &
        &                                              HFlapPsi(:,:,:,3,wave))
        call Derive(HFPsi(:,:,:,4,wave), +1, -1, +1,                           &
        &                                              HFdPsi(:,:,:,4,1,wave), &
        &                                              HFdPsi(:,:,:,4,2,wave), &
        &                                              HFdPsi(:,:,:,4,3,wave), &
        &                                              HFlapPsi(:,:,:,4,wave))
    enddo

    do wave=HFBlocks(1) + HFBlocks(3)+1,HFBlocks(1) + HFBlocks(3)+HFBlocks(5)
        call Derive(HFPsi(:,:,:,1,wave), +1, +1, +1,                           &
        &                                              HFdPsi(:,:,:,1,1,wave), &
        &                                              HFdPsi(:,:,:,1,2,wave), &
        &                                              HFdPsi(:,:,:,1,3,wave), &
        &                                              HFlapPsi(:,:,:,1,wave))
        call Derive(HFPsi(:,:,:,2,wave), -1, -1, +1,                           &
        &                                              HFdPsi(:,:,:,2,1,wave), &
        &                                              HFdPsi(:,:,:,2,2,wave), &
        &                                              HFdPsi(:,:,:,2,3,wave), &
        &                                              HFlapPsi(:,:,:,2,wave))
        call Derive(HFPsi(:,:,:,3,wave), -1, +1, -1,                           &
        &                                              HFdPsi(:,:,:,3,1,wave), &
        &                                              HFdPsi(:,:,:,3,2,wave), &
        &                                              HFdPsi(:,:,:,3,3,wave), &
        &                                              HFlapPsi(:,:,:,3,wave))
        call Derive(HFPsi(:,:,:,4,wave), +1, -1, -1,                           &
        &                                              HFdPsi(:,:,:,4,1,wave), &
        &                                              HFdPsi(:,:,:,4,2,wave), &
        &                                              HFdPsi(:,:,:,4,3,wave), &
        &                                              HFlapPsi(:,:,:,4,wave))
    enddo
   
    do wave=HFBlocks(1)+HFBlocks(3)+HFBlocks(5) + 1,                           &
    &       HFBlocks(1)+HFBlocks(3)+HFBlocks(5) + HFBLocks(7)
        call Derive(HFPsi(:,:,:,1,wave), +1, +1, -1,                           &
        &                                              HFdPsi(:,:,:,1,1,wave), &
        &                                              HFdPsi(:,:,:,1,2,wave), &
        &                                              HFdPsi(:,:,:,1,3,wave), &
        &                                              HFlapPsi(:,:,:,1,wave))
        call Derive(HFPsi(:,:,:,2,wave), -1, -1, -1,                           &
        &                                              HFdPsi(:,:,:,2,1,wave), &
        &                                              HFdPsi(:,:,:,2,2,wave), &
        &                                              HFdPsi(:,:,:,2,3,wave), &
        &                                              HFlapPsi(:,:,:,2,wave))
        call Derive(HFPsi(:,:,:,3,wave), -1, +1, +1,                           &
        &                                              HFdPsi(:,:,:,3,1,wave), &
        &                                              HFdPsi(:,:,:,3,2,wave), &
        &                                              HFdPsi(:,:,:,3,3,wave), &
        &                                              HFlapPsi(:,:,:,3,wave))
        call Derive(HFPsi(:,:,:,4,wave), +1, -1, +1,                           &
        &                                              HFdPsi(:,:,:,4,1,wave), &
        &                                              HFdPsi(:,:,:,4,2,wave), &
        &                                              HFdPsi(:,:,:,4,3,wave), &
        &                                              HFlapPsi(:,:,:,4,wave))
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
