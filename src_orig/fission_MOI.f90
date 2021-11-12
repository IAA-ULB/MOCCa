module fission_MOI
 !==============================================================================
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
 !==============================================================================
 ! Hephaestos keywords
 ! 
 ! PBROKEN : $PBROKEN
 !
 !
 !
 !------------------------------------------------------------------------------

  use geninfo
  use densities
  use parameterization
  use moments

  implicit none

  !-----------------------------------------------------------------------------
  !
  !
  real(KIND = dp), allocatable :: collective_inertia(:,:,:)
  
contains 

  subroutine calc_collective_inertia()
    !---------------------------------------------------------------------------
    !
    !
    !---------------------------------------------------------------------------
    real(KIND=dp), allocatable :: coll_me(:,:)
    real(KIND=dp)              :: s
    integer :: i
    
    
    
    allocate(coll_me(nwt,nwt))
    
    coll_me = Qlm_spme(2,0,.false.)
    
    s = 0 
    do i=1, nwt
      s = s + rho_can(i) * coll_me(i,i)
      print *, rho_can(i) , coll_me(i,i)
    enddo
    print *, 'Collective Q20', s
    
    coll_me = Qlm_spme(3,0,.false.)
    
    s = 0 
    do i=1, nwt
      s = s + rho_can(i) * coll_me(i,i)
      print *, rho_can(i) , coll_me(i,i)
    enddo
    print *, 'Collective Q30', s
    
    
    stop
  end subroutine calc_collective_inertia

  function Qlm_spme(l,m,Imaginary) result(me)
    !---------------------------------------------------------------------------
    ! Routine to calculate single-particle matrix elements of a given multipole
    ! moment in the CANONICAL basis. 
    ! 
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !     l,m       : integers, 
    !     Imaginary : to calculate < i | Re Q_lm | j > or < i |Im Q_lm | j >
    !                 logical
    ! Output:
    !     me        : single-particle matrix elements of the multipole moment
    !                 in the canonical basis, real(nwt,nwt)
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Notes:
    ! *) does not allow for the calculation of matrix elements with odd values 
    !    of m yet; those connect single-particle states with different signature,
    !    for the traditional definition of spherical harmonics at least.
    !---------------------------------------------------------------------------

    integer, intent(in)        :: l, m
    logical, intent(in)        :: Imaginary
    real(KIND=dp)              :: me(nwt,nwt)
    
    real(KIND=dp), allocatable         :: SpherHarmMesh(:,:,:,:,:,:)
    real(KIND=dp), allocatable, target :: Qlm(:)
    real(KIND=dp), pointer             :: harm_3D(:,:,:)
    
    integer :: i, j, Bi, Bj, Ni, Nj, si, sj, im, k

    if(mod(m,2) .eq. 1) then
      print *, 'Qlm_spme does not allow for the calculation of multipole matrix elements with odd m YET.'
      stop
    endif
    if(Imaginary) then 
      print *, 'Qlm_spme does not allow for the calculation of imaginary multipole matrix elements YET.'
      stop
    endif
    if((m .lt. 0) .or. (m .gt. l)) then
      print *, 'Invalid value of m in Qlm_spe.'      
      stop
    endif
    
    ! Initialize to zero
    me = 0.0d0

    ! Generate the necessary spherical harmonic
    allocate(SpherHarmMesh(nx,ny,nz,0:l,0:l,2), Qlm(nx*ny*nz))
    call GenSphericalHarmonics(l,nx,ny,nz, &
    &                          meshx_shifted,meshy_shifted,meshz_shifted,    & 
    &                          SpherHarmMesh,quantisationaxis,secondaryaxis) 

    if(Imaginary) then
      im  = 2
    else
      im  = 1
    endif
    
    harm_3D(1:nx,1:ny,1:nz) => Qlm
    harm_3D                 = SpherHarmMesh(:,:,:,l,m,im)
    deallocate(SpherHarmMesh)
    
    print *, Qlm(1), Qlm(2), Qlm(3)
   
    !--------------------------------------------------------------------------- 
    ! Loop over the neutron single-particle states
    !---------------------------------------------------------------------------
    si = 0
    do Bi = 1, 4
      Ni =  HFBlocks(Bi) ; if(Ni.eq.0) cycle
      do Bj = 1, 4
        Nj = HFBlocks(Bj) ; if(Nj.eq.0) cycle

        if(bj.eq.1) then
          sj = 0
        else
          sj = sum(HFblocks(1:Bj-1))
        endif

        ! Parity selection rule: if l = even, then only single-particle states
        !                        of identical parity (and signature) contribute    
$PCONSERVED        if(mod(l,2) .eq. 0 .and.  Bi .ne. Bj)        cycle
        ! Parity selection rule: if l = odd, then only single-particle states
        !                        of different parity contribute
$PCONSERVED        if(mod(l,2) .eq. 1 .and. abs(Bi - Bj).ne.2 ) cycle
        ! (THIS is only needed if parity conserved of course)
$PBROKEN if( Bi .ne. Bj ) cycle

        do i=1,NI    
         do j=1,NJ
          ! me = Int d^3r Sum_sigma psi^*_i(r,sigma) psi_j(r,sigma) Qlm(r)
          me(si+i,sj+j) = 0
          select case(pairingtype)
          case(0,1)
            do k=1,4
              me(si+i,sj+j) = me(si+i,sj+j) & 
              &                    + sum(HFpsi(:,k,si+i)*HFpsi(:,k,sj+j)*Qlm(:))
            enddo
          case(2)
            do k=1,4
              me(si+i,sj+j) = me(si+i,sj+j) &
              &                  + sum(Canpsi(:,k,si+i)*Canpsi(:,k,sj+j)*Qlm(:))
            enddo
          end select
          ! All these matrix elements are real if 
          ! (i)  we consider only real multipole moments
          ! (ii) time simplex is conserved
          ! 
          ! which means the matrix we store them in is symmetric
          me(si+i,sj+j) = me(si+i,sj+j) * dv 
          me(sj+j,si+i) = me(si+i,sj+j)
         enddo
        enddo
      enddo
      si = si + NI
    enddo
    
    !--------------------------------------------------------------------------- 
    ! Loop over the proton single-particle states
    !---------------------------------------------------------------------------
    si = nwn
    do Bi = 5, 8
      Ni =  HFBlocks(Bi) ; if(Ni.eq.0) cycle
      do Bj = 5, 8
        Nj = HFBlocks(Bj) ; if(Nj.eq.0) cycle
        sj = sum(HFblocks(1:Bj-1))

        ! Parity selection rule: if l = even, then only single-particle states
        !                        of identical parity (and signature) contribute    
$PCONSERVED        if(mod(l,2) .eq. 0 .and.  Bi .ne. Bj)        cycle
        ! Parity selection rule: if l = odd, then only single-particle states
        !                        of different parity contribute
$PCONSERVED        if(mod(l,2) .eq. 1 .and. abs(Bi - Bj).ne.2 ) cycle
        ! (THIS is only needed if parity conserved of course)
$PBROKEN if( Bi .ne. Bj ) cycle

        do i=1,NI    
         do j=1,NJ
          ! me = Int d^3r Sum_sigma psi^*_i(r,sigma) psi_j(r,sigma) Qlm(r)
          me(si+i,sj+j) = 0
          select case(pairingtype)
          case(0,1)
            do k=1,4
              me(si+i,sj+j) = me(si+i,sj+j) & 
              &                    + sum(HFpsi(:,k,si+i)*HFpsi(:,k,sj+j)*Qlm(:))
            enddo
          case(2)
            do k=1,4
              me(si+i,sj+j) = me(si+i,sj+j) &
              &                  + sum(Canpsi(:,k,si+i)*Canpsi(:,k,sj+j)*Qlm(:))
            enddo
          end select
          ! All these matrix elements are real if 
          ! (i)  we consider only real multipole moments
          ! (ii) time simplex is conserved
          ! 
          ! which means the matrix we store them in is symmetric
          me(si+i,sj+j) = me(si+i,sj+j) * dv 
          me(sj+j,si+i) = me(si+i,sj+j)
         enddo
        enddo
      enddo
      si = si + NI
    enddo
      
  end function Qlm_spme


end module fission_MOI
