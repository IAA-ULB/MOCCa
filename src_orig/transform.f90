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

  use geninfo
  use wavefunctions
  use pairing

  implicit none

contains

  subroutine TransformInput(filenx,fileny,filenz,filenwn,filenwp, filedx,      &
  &                         fileblocks, extraspwfs)
    !---------------------------------------------------------------------------
    ! Transform the input from file to the parameters of the new calculation.
    !
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
