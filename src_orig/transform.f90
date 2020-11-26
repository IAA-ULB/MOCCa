module transform
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
 ! Module containing all the code necessary to transform calculations: 
 !   *) Change the number of mesh points in any direction 
 !   *) Change the mesh constant dx  (not yet implemented)
 !   *) Break a set of symmetries
 !         => transformation of the spwfs
 !         => transformation of the densities (NOT IMPLEMENTED YET)
 !         => transformation of the pairing gaps
 ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 ! Hephaestos keywords
 ! 
 !  NONSPATIAL : $NONSPATIAL
 !
 !==============================================================================
  use geninfo
  use wavefunctions
  use pairing

  implicit none
  
  ! Indicates whether we need to transform the spwfs read on input.
  logical :: symtransfo_needed = .false.

contains

  subroutine Transformspwfs( wfs , blocks, oldnx, oldny, oldnz)
    !---------------------------------------------------------------------------
    ! Transform a set of spwfs, in a blockstructure dictated by blocks, into
    ! a set of spwfs with less symmetries.
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(inout), allocatable :: wfs(:,:,:)
    integer, intent(inout)                    :: blocks(8)
    integer, intent(in)                       :: oldnx, oldny, oldnz
    real(KIND=dp), allocatable                :: temp(:,:,:), tempe(:)
    real(KIND=dp), allocatable                :: tempd(:), tempr(:)
    real(KIND=dp), allocatable                :: tempgaps(:,:)
    integer, allocatable                :: tempsx(:,:), tempsy(:,:), tempsz(:,:)

    integer  :: wave, N, B, si, sb,i, wave2

    if(.not. allocated(rho_can)) then
        ! There is one case where this array might not be allocated upon entry
        ! in this routine: when initializing from scratch
        allocate(rho_can(sum(blocks))) ; rho_can = 0.0d0
    endif

    temp = wfs ; tempe = spenergies ; tempd = dispersions ; tempr = rho_can
    deallocate(wfs)         ; allocate(wfs(nx*ny*nz,4,nwt))    
    deallocate(dispersions) ; allocate(dispersions(nwt))
    deallocate(spenergies)  ; allocate(spenergies(nwt))   
    deallocate(rho_can)     ; allocate(rho_can(nwt))

    if( $NONSPATIAL ) then
      ! Use an antilinear, antihermitian symmetry operator 
      ! (usually time-reversal) to construct partners of the old wfs.

      si = 0
      sb = 0
      do B = 1,8
        !-----------------------------------------------------------------------
        ! Loop over the blocks. 
        !- - - - - - - - - - - - 
        ! Note that blocks = 2,4,6,8 are always of zero size in this loop, as we 
        ! are breaking the antilinear, anithermitian conserved symmetry
        N = blocks(B); if(N .eq. 0) cycle
   
        ! Copy the wavefunctions that were already in storage
        wfs(:,:,sb+1:sb+N)    = temp(:,:,si+1:si+N)
        dispersions(sb+1:sb+N)= tempd(si+1:si+N)
        spenergies(sb+1:sb+N) = tempe(si+1:si+N)
        rho_can(sb+1:sb+N)    = tempr(si+1:si+N) /2.0 ! Note the factor 1/2

        ! Use the antilinear symmetry to obtain the transformed spwfs
        do wave = 1, N        
         do i=1, oldnx*oldny*oldnz
           wfs(i,1,sb+N+wave) = $TRANSFO_NONSPATIAL_1 
           wfs(i,2,sb+N+wave) = $TRANSFO_NONSPATIAL_2
           wfs(i,3,sb+N+wave) = $TRANSFO_NONSPATIAL_3
           wfs(i,4,sb+N+wave) = $TRANSFO_NONSPATIAL_4
         enddo
        enddo
        dispersions(sb+N+1:sb+2*N)  = tempd(si+1:si+N)
        spenergies(sb+N+1 :sb+2*N)  = tempe(si+1:si+N)
        rho_can(sb+N+1:sb+2*N)      = tempr(si+1:si+N)/2.0 ! Note the factor 1/2

        si = si +   N
        sb = sb + 2*N
      enddo
      !-------------------------------------------------------------------------
      ! Create the second block in every pair
      HFBlocks(1) = blocks(1) ; HFBlocks(2) = blocks(1)
      HFBlocks(3) = blocks(3) ; HFBlocks(4) = blocks(3)
      HFBlocks(5) = blocks(5) ; HFBlocks(6) = blocks(5)
      HFBlocks(7) = blocks(7) ; HFBlocks(8) = blocks(7)

      !-------------------------------------------------------------------------
      ! And now we determine the signs of the reflections of the spwfs
      tempsx = sx ; tempsy = sy ; tempsz = sz
      deallocate(sx, sy, sz)
      allocate(sx(4,nwt), sy(4,nwt), sz(4,nwt))

      si = 0
      sb = 0
      do B = 1,8
        N = blocks(B) ; if(N .eq. 0) cycle
        do wave=1,N
          sx(:,sb + wave) = tempsx(:,si+wave)
          sy(:,sb + wave) = tempsy(:,si+wave)
          sz(:,sb + wave) = tempsz(:,si+wave)

          sx(:,sb + N + wave) = - tempsx(:,si+wave)
          sy(:,sb + N + wave) =   tempsy(:,si+wave)
          sz(:,sb + N + wave) = - tempsz(:,si+wave)
        enddo
        si = si +     N
        sb = sb + 2 * N
      enddo
      !-------------------------------------------------------------------------
      ! Transformation of the pairing gaps
      tempgaps = HFBgaps

      deallocate(HFBgaps) ; allocate(HFBgaps(nwt, nwt)) ; HFBgaps = 0

      si = 0  ; sb = 0
      do B = 1,8
        N = blocks(B) ; if(N .eq. 0) cycle
        do wave=1,N
          do wave2=1,N
            HFBgaps(sb + wave    , sb + wave2 + N)  = tempgaps(si+wave,si+wave2)
            HFBgaps(sb + wave + N, sb + wave2    )  =-tempgaps(si+wave,si+wave2)
          enddo
        enddo

!        print *, 'Gaps on file ', B, N, si
!        do wave=1,N
!          print ('(99f10.3)'), tempgaps(si+wave, si+1:si+N)
!        enddo
!        print *
!        print *, 'Gaps in block ', B, N, si
!        do wave=1,2*N
!          print ('(99f10.3)'), HFBgaps(sb+wave, sb+1:sb+2*N)
!        enddo
!        print *
        si = si +   N
        sb = sb + 2*N
      enddo

      ! Clean up
      deallocate(temp, tempsx, tempsy, tempsz)    
    endif
  end subroutine Transformspwfs

  subroutine TransformInput(filenx,fileny,filenz,filenwn,filenwp, filedx,      &
  &                         fileblocks, extraspwfs)
    !---------------------------------------------------------------------------
    ! Transform the input from file to the parameters of the new calculation.
    !---------------------------------------------------------------------------

 1 format &
     &  (/,8x,' ___________________________________________________________', &
     &   /,8x,'| Transformation of the input                              |', &
     &   /,8x,'| On file:                                                 |', &
     &   /,8x,'|    nx, ny, nz    = ', 3i5, '                       |', &
     &   /,8x '|    dx            = ', f10.7, ' (fm)                       |', &
     &   /,8x,'|    nwn, nwn      = ', 2i7 , '                        |', &
     &   /,8x,'|    (n+,n-,p+,p-) = (', 4i5, ')                |', &
     &   /,8x,'|This calculation:                                         |' , &
     &   /,8x,'|    nx, ny, nz    = ', 3i5, '                       |'       , &
     &   /,8x '|    dx            = ', f10.7, '(fm)                        |', &
     &   /,8x,'|    nwn, nwn      = ', 2i7, '                        |'      , &
     &   /,8x,'|    (n+,n-,p+,p-) = (', 4i5, ')                |'            , &
     &   /,8x,'|__________________________________________________________|')

  2 format (' Interpolation (changing of dx) not yet allowed.')

    integer, intent(in)        :: filenx,fileny,filenz,filenwn, filenwp
    integer, intent(in)        :: fileblocks(blocks),extraspwfs(blocks)
    real(KIND=dp), intent(in)  :: filedx
    real(KIND=dp), allocatable :: extended(:,:,:), newenergy(:), temp(:)
    real(KIND=dp), allocatable :: temp2(:,:)

    logical :: ChangeBoxSize = .false., transformed= .false.
    integer :: wave,i,j, b, sb, sf

    transformed = .false.
    !---------------------------------------------------------------------------
    ! Change of box size (unchanged dx but changed nx/ny/nz)
    ChangeBoxSize = .false.
    if ( dx .eq. filedx ) then
      if(nx.ne.filenx .or. ny .ne.fileny .or. nz .ne. filenz) then
        ChangeBoxSize = .true.
      endif
    else
      print 2
      stop
    endif    

    if ( ChangeBoxSize ) then
        !-----------------------------------------------------------------------
        ! If we are just changing the box size, we only need to attack the 
        ! SPWFS, as the potentials are taken care of in their own ReadPotential
        ! routine in functional.f90
        !-----------------------------------------------------------------------
        call ChangeBoxSizeSpwf(HFPsi,filenx, fileny, filenz, filenwn+filenwp)
        transformed = .true.
    endif

    !---------------------------------------------------------------------------
    ! Add in extra wavefunctions
    if(filenwn .ne. nwn .or. filenwp .ne. nwp) then
          !---------------------------------------------------------------------
          !  First a bunch of sanity checks
          if(nwn .ne. sum(fileblocks(1:4)) + sum(extraspwfs(1:4)) ) then
            print *, 'Inconsistent number of neutron wavefunctions.'
            stop
          endif
          if(nwp .ne. sum(fileblocks(5:8)) + sum(extraspwfs(5:8)) ) then
            print *, 'Inconsistent number of proton wavefunctions.'
            stop
          endif
          if(nwn .lt. filenwn) then
            print *, ' Nwn lower than nwn on file.'
            stop
          endif  
          if(nwp .lt. filenwp) then
            print *, ' Nwp lower than nwp on file.'
            stop
          endif 

          if((extraspwfs(2).ne.0) .or. &
          &  (extraspwfs(4).ne.0) .or. & 
          &  (extraspwfs(6).ne.0) .or. &
          &  (extraspwfs(8).ne.0) ) then
            print *, 'Blocks 2,4,6,8 not allowed in this version of Tantalus.'
            stop
          endif
          !---------------------------------------------------------------------
          allocate(extended(nx*ny*nz,4,nwt)) ; allocate(newenergy(nwt))
          extended = 0.0
          hfblocks = fileblocks + extraspwfs ; newenergy = 1000.0

          sb = 0 ; sf = 0  
          do b = 1, blocks
            do i=1, fileblocks(b)
                extended(:,:,sb+i) = hfpsi(:,:,sf+i)
                newenergy(sb+i)    = spenergies(sf+i)
            enddo
            do i=sb+fileblocks(b)+1, sb+ fileblocks(b) + extraspwfs(b)
                call random_number(extended(:,:,i))
                newenergy(i)    = +1000.0
            enddo
            sb  = sb + hfblocks(b)
            sf  = sf + fileblocks(b)
          enddo
          hfpsi      = extended
          spenergies = newenergy
          deallocate(dispersions) ; allocate(dispersions(nwt)) ; dispersions=0.0
          deallocate(rho_can)     ; allocate(rho_can(nwt))     ; rho_can    =0.0
          !---------------------------------------------------------------------
          ! Dealing with the gaps
          select case (pairingtype)
          case(0)
            ! HF, nothing to do
          case(1)
            ! BCS, need to move the gaps into the correct position
            temp = BCSgaps ; deallocate(BCSgaps)
            allocate(BCSgaps(nwt)) ; BCSgaps = 0.0
    
            sb = 0 ; sf = 0  
            do b = 1, blocks
              do i=1, fileblocks(b)
                 BCSgaps(sb+i) = temp(sf+i)
              enddo
              sb  = sb + hfblocks(b)
              sf  = sf + fileblocks(b)
            enddo            
          case(2)
            !HFB, need to move the gaps into the correct position
            temp2 = HFBgaps ; deallocate(HFBgaps)
            allocate(HFBgaps(nwt,nwt)) ; HFBgaps = 0.0 
    
            sb = 0 ; sf = 0  
            do b = 1, blocks
              do i=1, fileblocks(b)
                do j=1, fileblocks(b)
                 HFBgaps(sb+i, sb+j) = temp2(sf+i, sf+j)
                enddo
              enddo
              sb  = sb + hfblocks(b)
              sf  = sf + fileblocks(b)
            enddo

            temp2 = kappa_pairing; deallocate(kappa_pairing)
            allocate(kappa_pairing(nwt,nwt)) ; kappa_pairing = 0.0
            
            sb = 0 ; sf = 0  
            do b = 1, blocks
              do i=1, fileblocks(b)
                do j=1, fileblocks(b)
                 kappa_pairing(sb+i, sb+j) = temp2(sf+i, sf+j)
                enddo
              enddo
              sb  = sb + hfblocks(b)
              sf  = sf + fileblocks(b)
            enddo
          end select
    else
      ! Copy this information
      hfblocks = fileblocks
    endif 

    print 1, filenx, fileny, filenz, filedx, filenwn, filenwp,                &
    &        fileblocks(1:8:2), nx,ny,nz,dx,nwn,nwp, hfblocks(1:8:2)

  end subroutine TransformInput

  subroutine ChangeBoxSizeSpwf(Phi,filenx,fileny, filenz, filenwt) 
    !---------------------------------------------------------------------------
    ! Put a given wave-function from a mesh with dx to the new mesh 
    ! (nx,ny,nz,dx) with same dx.
    !---------------------------------------------------------------------------
    
    real(KIND=dp), intent(inout), allocatable, target :: phi(:,:,:)
    integer, intent(in)                               :: filenx, fileny, filenz
    integer, intent(in)                               :: filenwt
    real(KIND=dp), allocatable, target                :: temp(:,:,:)
    integer                                           :: wave, l

    temp = phi ;  deallocate(phi)
    allocate(phi(nx*ny*nz,4,filenwt)) ; HFPSI = 0

    do wave=1, filenwt
      do l=1,4
         phi(:,l,wave) = changeboxsize_function(temp(:,l,wave),        &
         &                                               filenx, fileny, filenz)         
      enddo
    enddo

  end subroutine ChangeBoxSizeSpwf

  function changeboxsize_function(f, filenx, fileny, filenz) result(ft)
    !---------------------------------------------------------------------------
    ! Put a given function from a mesh with dx to the new mesh 
    ! (nx,ny,nz,dx) with same dx.
    !---------------------------------------------------------------------------
    integer, intent(in)                :: filenx, fileny, filenz
    real(KIND=dp), intent(in), target  :: f(:)
    real(KIND=dp), allocatable, target :: ft(:)
    real(KIND=dp), pointer             :: f3 (:,:,:), ft3(:,:,:)
    integer                            :: endx, endy, endz

    allocate(ft(nx*ny*nz)) ; ft = 0
  
    ft3(1:nx, 1:ny, 1:nz)             => ft(1:nx*ny*nz)
    f3(1:filenx, 1:fileny, 1:filenz)  => f (1:filenx*fileny*filenz)
       
    endx = min(nx, filenx)
    endy = min(ny, fileny)
    endz = min(nz, filenz)

    ft3(1:endx,1:endy,1:endz) = f3(1:endx,1:endy,1:endz)

  end function changeboxsize_function

end module transform
