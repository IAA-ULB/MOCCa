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
 !  SPATIAL    : $SPATIAL
 !  EXPANDX    : $EXPANDX
 !  EXPANDY    : $EXPANDY
 !  EXPANDZ    : $EXPANDZ
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
  1 format("------------------------------------------------------------------")
  2 format(" Number of spwfs on STDIN does not match what is required.")
  3 format(" Symmetry-broken calculation should have (nwn,nwp) = ", 2i5)
  4 format(" STDIN says                              (nwn,nwp) = ", 2i5)

    real(KIND=dp), intent(inout), allocatable, target :: wfs(:,:,:)
    integer, intent(inout)                    :: blocks(8)
    integer, intent(in)                       :: oldnx, oldny, oldnz

    real(KIND=dp), allocatable                :: temp(:,:,:), tempe(:)
    real(KIND=dp), allocatable                :: tempd(:), tempr(:)
    real(KIND=dp), allocatable, target        :: wftarget(:,:)
    real(KIND=dp), allocatable                :: tempgaps(:,:), tempkap(:,:)

    real(KIND=dp),  pointer                   :: right3D(:,:,:,:), left3D(:,:,:,:) 

    integer  :: wave, N, B, si, sb,i, wave2, offset, j, k

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

      ! A more elegant way of stopping if the number of spwfs does not match
      if((2*sum(blocks(1:4)) .ne. nwn)  .or. (2*sum(blocks(5:8)) .ne. nwp)) then
        print 1
        print 2
        print 3, 2*sum(blocks(1:4)), 2*sum(blocks(5:8))
        print 4, nwn, nwp
        print 1
        stop
      endif
  
      si = 0
      sb = 0
      do B = 1,8
        !-----------------------------------------------------------------------
        ! Loop over the blocks. 
        !- - - - - - - - - - - - 
        ! Note that blocks = 2,4,6,8 are always of zero size in this loop, as we 
        ! are breaking the antilinear, antihermitian conserved symmetry
        N = blocks(B); if(N .eq. 0) cycle
   
        ! Copy the wavefunctions that were already in storage
        wfs(:,:,sb+1:sb+N)    = temp(:,:,si+1:si+N)
        dispersions(sb+1:sb+N)= tempd(si+1:si+N)
        spenergies(sb+1:sb+N) = tempe(si+1:si+N)
        rho_can(sb+1:sb+N)    = tempr(si+1:si+N) /2.0 ! Note the factor 1/2

        ! Use the antilinear symmetry to obtain the transformed spwfs
        do wave = 1, N  
         wftarget = temp(:,:,si+wave)      
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
      ! Transformation of the pairing gaps and pairing tensor kappa
      if(pairingtype.eq.2) then
        ! Only do this if HFB gaps have been read from file, otherwise we rely
        ! on the initialization routine for gaps
        if(allocated(HFBgaps)) then 
          tempgaps = HFBgaps
          tempkap  = kappa_pairing 
        
          deallocate(HFBgaps)       ; allocate(HFBgaps(nwt, nwt))  
          deallocate(kappa_pairing) ; allocate(kappa_pairing(nwt,nwt))

          HFBgaps = 0; kappa_pairing = 0

          si = 0  ; sb = 0
          do B = 1,8
            N = blocks(B) ; if(N .eq. 0) cycle
            do wave=1,N
              do wave2=1,N
                HFBgaps(sb + wave    , sb + wave2 + N)  = &
                &                                     tempgaps(si+wave,si+wave2)
                HFBgaps(sb + wave + N, sb + wave2    )  = &
                &                                    -tempgaps(si+wave,si+wave2)

                kappa_pairing(sb + wave    , sb + wave2 + N) = & 
                &                                     tempgaps(si+wave,si+wave2)
                kappa_pairing(sb + wave + N, sb + wave2    ) = &
                &                                    -tempgaps(si+wave,si+wave2)
              enddo
            enddo

            si = si +   N
            sb = sb + 2*N
          enddo
        endif
      endif
    endif
    
    if( $SPATIAL ) then
      if($EXPANDX) then
          print *, 'Extending to the full X-axis not implemented yet'
          stop
      elseif($EXPANDY) then   
          print *, 'Extending to the full Y-axis not implemented yet'
          stop
      elseif($EXPANDZ) then
          ! Expand the wavefunctions in the first relevant block
          si = 0
          sb = 0
          offset = 0
          do B=1,8,4 ! This is essentially an isospin loop now
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
            ! First block is not modified
            do wave=1,HFBlocks(B)
              wftarget = temp(:,:,sb+wave)
              right3D(1:oldnx,1:oldny, 1:oldnz,1:4) => wftarget(:,:)
              left3D(1:nx, 1:ny, 1:nz,1:4)          => wfs(:,:, sb+HFBlocks(B))
            
              ! Copy the positive z-axis
              wfs(oldnx*oldny*oldnz+1:nx*ny*nz,1:4,sb+wave) = wftarget
              do i=1,oldnx*oldny*oldnz
!                wfs(i,1,sb+wave) = $TRANSFO_Z_1_B1
!                wfs(i,2,sb+wave) = $TRANSFO_Z_2_B1
!                wfs(i,3,sb+wave) = $TRANSFO_Z_3_B1
!                wfs(i,4,sb+wave) = $TRANSFO_Z_4_B1
              enddo
            enddo
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
            ! What was the third block on file, now becomes part of the first
            ! block
            offset = HFBlocks(B)
            do wave=1,HFBlocks(B+2)
              ! Copy the positive z-axis
!              wfs(1:oldnx*oldny*oldnz,1:4,sb+offset+wave) = &
!              &      temp(:,:,sb+sum(HFBlocks(B:B+1))+1)
              do i=1,oldnx*oldny*oldnz

              enddo
            enddo
            sb = sb + sum(HFBlocks(1:B+3))            
          enddo
      endif
      !-------------------------------------------------------------------------
      ! Merge the blocks associated with a linear, hermitian symmetry
      HFBlocks(1) = blocks(1) + HFBlocks(3)
      HFBlocks(2) = blocks(2) + HFBlocks(4)
      HFBlocks(5) = blocks(5) + HFBlocks(7)
      HFBlocks(6) = blocks(6) + HFBlocks(8)
      ! Remove the old blocks
      HFBlocks(3:4) = 0
      HFBlocks(7:8) = 0
    endif
    
  end subroutine Transformspwfs

  subroutine TransformInput(filenx,fileny,filenz,filenwn,filenwp, filedx,      &
  &                         fileblocks, extraspwfs)
    !---------------------------------------------------------------------------
    ! Transform the input from file to the parameters of the new calculation.
    ! Note that this means either 
    !   (*) A change of mesh parameters
    !   (*) A change of number of wavefunctions
    !
    ! but all of this with the same conserved/broken symmetries. These 
    ! manipulations are not compatible in the same run with the breaking 
    ! additional symmetries, achieved by the Transformspwfs routine.
    !---------------------------------------------------------------------------

 1 format &
     &  (/,8x,' ___________________________________________________________', &
     &   /,8x,'| Transformation of the input                              |', &
     &   /,8x,'| On file:                                                 |', &
     &   /,8x,'|    nx, ny, nz    = ', 3i5, '                       |', &
     &   /,8x '|    dx            = ', f10.7, ' (fm)                       |', &
     &   /,8x,'|    nwn, nwn      = ', 2i7 , '                        |', &
     &   /,8x,'|    (n+,n-)       = (', 4i5, ')                |', &
     &   /,8x,'|    (p+,p-)       = (', 4i5, ')                |', &
     &   /,8x,'| This calculation:                                        |' , &
     &   /,8x,'|    nx, ny, nz    = ', 3i5, '                       |'       , &
     &   /,8x '|    dx            = ', f10.7, '(fm)                        |', &
     &   /,8x,'|    nwn, nwn      = ', 2i7, '                        |'      , &
     &   /,8x,'|    (n+,n-)       = (', 4i5, ')                |', &
     &   /,8x,'|    (p+,p-)       = (', 4i5, ')                |', &
     &   /,8x,'|__________________________________________________________|')

  2 format (' Interpolation (changing of dx) not yet allowed.')

    integer, intent(in)        :: filenx,fileny,filenz,filenwn, filenwp
    integer, intent(in)        :: fileblocks(blocks),extraspwfs(blocks)
    real(KIND=dp), intent(in)  :: filedx
    real(KIND=dp), allocatable :: extended(:,:,:), newenergy(:), temp(:)
    real(KIND=dp), allocatable :: temp2(:,:)

    logical :: ChangeBoxSize = .false., transformed= .false.
    integer :: i,j, b, sb, sf

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

$TR       if((extraspwfs(2).ne.0) .or. &
$TR          &  (extraspwfs(4).ne.0) .or. & 
$TR          &  (extraspwfs(6).ne.0) .or. &
$TR          &  (extraspwfs(8).ne.0) ) then
$TR           print *, 'Blocks 2,4,6,8 not allowed with T conserved.'
$TR           stop
$TR       endif
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
          ! Cleaning up some stuff
          if(allocated(dispersions)) deallocate(dispersions) 
          allocate(dispersions(nwt)) ; dispersions=0.0
          if(allocated(rho_can))     deallocate(rho_can)     
          allocate(rho_can(nwt))     ; rho_can    =0.0
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
            if(allocated(HFBgaps)) then
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
              else
                 ! We do nothing if no gaps were initialised yet.
              endif 
          end select
          !---------------------------------------------------------------------
    else
      ! Copy this information
      hfblocks = fileblocks
    endif 

    print 1, filenx, fileny, filenz, filedx, filenwn, filenwp,                &
    &        fileblocks(1:8), nx,ny,nz,dx,nwn,nwp, hfblocks(1:8)

  end subroutine TransformInput

  subroutine ChangeBoxSizeSpwf(Phi,filenx,fileny, filenz, filenwt) 
    !---------------------------------------------------------------------------
    ! Put a given wave-function from a mesh with dx to the new mesh 
    ! (nx,ny,nz,dx) with same dx.
    !---------------------------------------------------------------------------
    
    real(KIND=dp), intent(inout), allocatable, target :: phi(:,:,:)
    integer, intent(in)                               :: filenx, fileny, filenz
    integer, intent(in)                               :: filenwt
    real(KIND=dp), target                :: temp(filenx*fileny*filenz,4,filenwt)
    integer                              :: wave, l

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
