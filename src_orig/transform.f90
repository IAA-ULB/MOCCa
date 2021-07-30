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
 !
 ! PBROKEN     : $PBROKEN
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
    real(KIND=dp), allocatable                :: tempbogo(:,:), temprho(:,:)
    real(KIND=dp), allocatable                :: temptransfo(:,:), tempsph(:,:)
    real(KIND=dp), allocatable                :: tempconfig(:)
    
    real(KIND=dp),  pointer                   :: right3D(:,:,:,:) 
    real(KIND=dp),  pointer                   :: left3D(:,:,:,:) 


    integer  :: wave, N, B, si, sb,i, wave2, offset_left, offset_right, j, k, sc

    if(.not. allocated(rho_can)) then
        ! There is one case where this array might not be allocated upon entry
        ! in this routine: when initializing from scratch
        allocate(rho_can(sum(blocks))) ; rho_can = 0.0d0
    endif

    temp        = wfs         ; tempe   = spenergies 
    tempd       = dispersions ; tempr   = rho_can
    temptransfo = hftransfo   ; tempsph = current_sph
    tempconfig  = configmatrix
    
    deallocate(wfs)         ; allocate(wfs(nx*ny*nz,4,nwt))    
    deallocate(dispersions) ; allocate(dispersions(nwt))    ; dispersions  = 0
    deallocate(spenergies)  ; allocate(spenergies(nwt))     ; spenergies   = 0
    deallocate(rho_can)     ; allocate(rho_can(nwt))        ; rho_can      = 0
    deallocate(hftransfo)   ; allocate(hftransfo(nwt,nwt))  ; hftransfo    = 0
    deallocate(current_sph) ; allocate(current_sph(nwt,nwt)); current_sph  = 0
    deallocate(configmatrix); allocate(configmatrix(2*nwt)) ; configmatrix = 0

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

        ! Use the antilinear, antihermitian symmetry to obtain the transformed
        ! spwfs. The remapping into 3D functions is superfluous here, but 
        ! it makes the Hephaestos coding more flexible. 
        do wave = 1, N  
         wftarget = temp(:,:,si+wave)
         right3D(1:oldnx,1:oldny, 1:oldnz,1:4) => wftarget(:,:)
         left3D(1:nx, 1:ny, 1:nz,1:4)          => wfs(:,:,sb+N+wave)      
         do k=1, oldnz
          do j=1,oldny
           do i=1, oldnx
             left3D(i,j,k,1) = $TRANSFO_NONSPATIAL_1 
             left3D(i,j,k,2) = $TRANSFO_NONSPATIAL_2
             left3D(i,j,k,3) = $TRANSFO_NONSPATIAL_3
             left3D(i,j,k,4) = $TRANSFO_NONSPATIAL_4
           enddo
          enddo
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
      ! Transformation of pairing quantities
      ! (1)  the pairing gaps
      ! (2)  the density matrix rho
      ! (3)  the anomalous density matrix
      ! (4)  the Bogoliubov transformation
      ! (5)  the configuration matrix
      ! (6)  the HF transformation
      ! (7)  the current single-particle hamiltonian
      if(pairingtype.eq.2) then
        ! Only do this if HFB gaps have been read from file, otherwise we rely
        ! on the initialization routine for gaps
        if(allocated(HFBgaps)) then 
          tempgaps = HFBgaps
          tempkap  = kappa_pairing 
          tempbogo = Bogoliubov    
          temprho  = rho_pairing
          
          deallocate(HFBgaps)       ; allocate(HFBgaps(nwt, nwt))  
          deallocate(kappa_pairing) ; allocate(kappa_pairing(nwt,nwt))
          deallocate(Bogoliubov)    ; allocate(Bogoliubov(2*nwt,2*nwt))
          deallocate(rho_pairing)   ; allocate(rho_pairing(nwt,nwt))
          
          HFBgaps = 0; rho_pairing = 0 ; kappa_pairing = 0 ; Bogoliubov = 0

          si = 0  ; sb = 0 ; sc = 0
          do B = 1,8
            N = blocks(B) ; if(N .eq. 0) cycle
            !-------------------------------------------------------------------  
            ! (1)  the pairing gaps
            ! (2)  the density matrix rho
            ! (3)  the anomalous density matrix
            do wave=1,N
              do wave2=1,N
                HFBgaps(sb + wave     , sb + wave2 + N)  = &
                &                                     tempgaps(si+wave,si+wave2)
                HFBgaps(sb + wave  + N, sb + wave2    )  = &
                &                                    -tempgaps(si+wave,si+wave2)
                
                rho_pairing(sb+wave  , sb+wave2)    = &
                &                                      temprho(si+wave,si+wave2)
                rho_pairing(sb+wave+N, sb+wave2+N)  = &
                &                                      temprho(si+wave,si+wave2)

                kappa_pairing(sb + wave    , sb + wave2 + N) = & 
                &                                      tempkap(si+wave,si+wave2)
                kappa_pairing(sb + wave+ N , sb + wave2     ) = &
                &                                     -tempkap(si+wave,si+wave2)
              enddo 
            enddo
            !-------------------------------------------------------------------  
            ! (4) The current Bogoliubov transformation
            do wave=1,N
                ! On file, the Bogoliubov transformation has the following form
                !
                !  W =   ( V^T,+  U^+)     => time-reversal invariant, i.e.
                !        ( U^T,+  V^+)        half of all columns 
                !
                ! but we need to produce a Bogoliubov transform that reads
                ! (in every pair of blocks linked by an antihermitian, linear 
                !  symmetry)
                ! 
                !
                !       (  V^*+   0     U+  0   )
                !  W =  (  0      V^*-  0   U-  )      
                !       (  0      U^*-  0   V-  )
                !       (  U^*+   0     V+  0   )
                !
                ! with U^+ = U^- and V^- = - V^+.

                ! We start by getting the r.h.s. columns correct
                ! - - - - - - - - - - - - - - - - - - - - - - - -
                ! U^+
                Bogoliubov(sc    +1:sc+  N, sc+2*N+wave) = &
                &                              tempbogo(sb  +1:sb+  N,sb+N+wave)
                ! V^+ (note the minus sign!)
                Bogoliubov(sc+3*N+1:sc+4*N, sc+2*N+wave) = &
                &                            - tempbogo(sb+N+1:sb+2*N,sb+N+wave)
                ! U^-
                Bogoliubov(sc+  N+1:sc+2*N, sc+3*N+wave) = &
                &                              tempbogo(sb  +1:sb+  N,sb+N+wave)
                ! V^- 
                Bogoliubov(sc+2*N+1:sc+3*N, sc+3*N+wave) = &
                &                              tempbogo(sb+N+1:sb+2*N,sb+N+wave)
                ! - - - - - - - - - - - - - - - - - - - - - - - -
                ! and then we could construct the l.h.s. columns by symmetry ,
                ! but this is never used by the code.        
            enddo
            
            !-------------------------------------------------------------------
            ! (5) The configuration matrix
            configmatrix(sc    +1:sc  +N) = tempconfig(sb+1:sb+N)
            configmatrix(sc+  N+1:sc+2*N) = tempconfig(sb+1:sb+N)
            configmatrix(sc+2*N+1:sc+3*N) = tempconfig(sb+N+1:sb+2*N)
            configmatrix(sc+3*N+1:sc+4*N) = tempconfig(sb+N+1:sb+2*N)
            !-------------------------------------------------------------------
            ! (6) The current HF transformation
            hftransfo(sb  +1:sb+  N,sb  +1:sb+  N) = &
            &                                   temptransfo(si+1:si+N,si+1:si+N)
            hftransfo(sb+N+1:sb+2*N,sb+N+1:sb+2*N) = &
            &                                   temptransfo(si+1:si+N,si+1:si+N)
            !-------------------------------------------------------------------
            ! (7) The current single-particle hamiltonian
            current_sph(sb  +1:sb+  N,sb  +1:sb  +N)         = &
            &                                       tempsph(si+1:si+N,si+1:si+N)
            current_sph(sb+N+1:sb+2*N,sb+N+1:sb+2*N) = &
            &                                       tempsph(si+1:si+N,si+1:si+N)

            si = si +   N
            sb = sb + 2*N
            sc = sc + 4*N
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
          sb = 0
          do B=1,8,4 ! This is essentially an isospin loop now
            
            !-------------------------------------------------------------------
            ! First take are of all auxiliary matrices

            ! First block does not get modified
            offset_left  = 0
            offset_right = 0
            dispersions(sb+offset_left +1:sb+offset_left +blocks(B)) &
            &   = tempd(sb+offset_right+1:sb+offset_right+blocks(B))
            spenergies (sb+offset_left +1:sb+offset_left +blocks(B)) &
            &   = tempe(sb+offset_right+1:sb+offset_right+blocks(B))
            rho_can    (sb+offset_left +1:sb+offset_left +blocks(B)) &
            &   = tempr(sb+offset_right+1:sb+offset_right+blocks(B))

            ! Third block on file becomes part of the first block
            offset_left  = blocks(B)
            offset_right = blocks(B) + blocks(B+1)
            dispersions(sb+offset_left +1:sb+offset_left +blocks(B+2)) &
            &   = tempd(sb+offset_right+1:sb+offset_right+blocks(B+2))
            spenergies (sb+offset_left +1:sb+offset_left +blocks(B+2)) &
            &   = tempe(sb+offset_right+1:sb+offset_right+blocks(B+2))
            rho_can    (sb+offset_left +1:sb+offset_left +blocks(B+2)) &
            &   = tempr(sb+offset_right+1:sb+offset_right+blocks(B+2))
            
            ! Second block on file becomes first part of second block
            offset_left  = blocks(B) + blocks(B+2)
            offset_right = blocks(B)
            dispersions(sb+offset_left +1:sb+offset_left +blocks(B+1)) &
            &   = tempd(sb+offset_right+1:sb+offset_right+blocks(B+1))
            spenergies (sb+offset_left +1:sb+offset_left +blocks(B+1)) &
            &   = tempe(sb+offset_right+1:sb+offset_right+blocks(B+1))
            rho_can    (sb+offset_left +1:sb+offset_left +blocks(B+1)) &
            &   = tempr(sb+offset_right+1:sb+offset_right+blocks(B+1))

            ! Fourth  block on file becomes second part of second block
            offset_left  = sum(blocks(B:B+2))
            offset_right = sum(blocks(B:B+2))
            dispersions(sb+offset_left +1:sb+offset_left +blocks(B+3)) &
            &   = tempd(sb+offset_right+1:sb+offset_right+blocks(B+3))
            spenergies (sb+offset_left +1:sb+offset_left +blocks(B+3)) &
            &   = tempe(sb+offset_right+1:sb+offset_right+blocks(B+3))
            rho_can    (sb+offset_left +1:sb+offset_left +blocks(B+3)) &
            &   = tempr(sb+offset_right+1:sb+offset_right+blocks(B+3))

            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
            ! First block is not modified
            offset_left = 0
            offset_right= 0
            do wave=1,Blocks(B)
              wftarget = temp(:,:,sb+offset_right+wave)
              right3D(1:oldnx,1:oldny, 1:oldnz,1:4) &
              &                                 => wftarget(:,:)
              left3D(1:nx, 1:ny, 1:nz,1:4) &
              &                                 => wfs(:,:, sb+offset_left+wave)
            
              !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
              ! Actual spatial transformations
              do j=1,oldny
               do i=1, oldnx
                  ! Copy the positive z-axis
                  left3D(i,j,oldnz+1:2*oldnz,:) = right3D(i,j,:,:)
               enddo
              enddo
              
              do k=1, oldnz
                do j=1,oldny
                 do i=1, oldnx
                  left3D(i,j,k,1) = $TRANSFO_Z_1_B1
                  left3D(i,j,k,2) = $TRANSFO_Z_2_B1
                  left3D(i,j,k,3) = $TRANSFO_Z_3_B1
                  left3D(i,j,k,4) = $TRANSFO_Z_4_B1
                 enddo
                enddo
              enddo
            enddo
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
            ! What was the 3rd block on file, now becomes part of the 1st block
            offset_right = Blocks(B) + Blocks(B+1) ! Block B+2 from file
            offset_left  = Blocks(B) 
            do wave=1,Blocks(B+2)
              wftarget = temp(:,:,sb+offset_right+wave)
              right3D(1:oldnx,1:oldny, 1:oldnz,1:4) &
              &                                 => wftarget(:,:)
              left3D(1:nx, 1:ny, 1:nz,1:4) &
              &                                 => wfs(:,:, sb+offset_left+wave)
              
              !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
              ! Actual spatial transformations
              do j=1,oldny
               do i=1, oldnx
                  ! Copy the positive z-axis
                  left3D(i,j,oldnz+1:2*oldnz,:) = right3D(i,j,:,:)
               enddo
              enddo
            
              do k=1, oldnz
                do j=1,oldny
                 do i=1, oldnx
                  left3D(i,j,k,1) = $TRANSFO_Z_1_B3
                  left3D(i,j,k,2) = $TRANSFO_Z_2_B3
                  left3D(i,j,k,3) = $TRANSFO_Z_3_B3
                  left3D(i,j,k,4) = $TRANSFO_Z_4_B3
                 enddo
                enddo
              enddo
              !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            enddo
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
            ! Second block on file stays the second block
            offset_right = Blocks(B)                    ! Block B+1 from file
            offset_left  = Blocks(B) + Blocks(B+2) 
            do wave=1,Blocks(B+1)
              wftarget = temp(:,:,sb+offset_right+wave)
              right3D(1:oldnx,1:oldny, 1:oldnz,1:4) &
              &                                 => wftarget(:,:)
              left3D(1:nx, 1:ny, 1:nz,1:4) &
              &                                 => wfs(:,:, sb+offset_left+wave)
            
              !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
              ! Actual spatial transformations
              do j=1,oldny
               do i=1, oldnx
                  ! Copy the positive z-axis
                  left3D(i,j,oldnz+1:2*oldnz,:) = right3D(i,j,:,:)
               enddo
              enddo
              
              do k=1,oldnz
                do j=1,oldny
                 do i=1, oldnx
                  left3D(i,j,k,1) = $TRANSFO_Z_1_B2
                  left3D(i,j,k,2) = $TRANSFO_Z_2_B2
                  left3D(i,j,k,3) = $TRANSFO_Z_3_B2
                  left3D(i,j,k,4) = $TRANSFO_Z_4_B2
                 enddo
                enddo
              enddo
              !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -             
            enddo
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
            ! Fourth block on file becomes part of the second block
            offset_right = sum(Blocks(B:B+2))         ! Block B+3 from file
            offset_left  = sum(Blocks(B:B+2)) 
            do wave=1,Blocks(B+3)
              wftarget = temp(:,:,sb+offset_right+wave)
              right3D(1:oldnx,1:oldny, 1:oldnz,1:4) &
              &                                 => wftarget(:,:)
              left3D(1:nx, 1:ny, 1:nz,1:4) &
              &                                 => wfs(:,:, sb+offset_left+wave)
              
              !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
              ! Actual spatial transformations
              do j=1,oldny
               do i=1, oldnx
                  ! Copy the positive z-axis
                  left3D(i,j,oldnz+1:2*oldnz,:) = right3D(i,j,:,:)
               enddo
              enddo
               do k=1, oldnz
                do j=1,oldny
                 do i=1, oldnx
                  left3D(i,j,k,1) = $TRANSFO_Z_1_B4
                  left3D(i,j,k,2) = $TRANSFO_Z_2_B4
                  left3D(i,j,k,3) = $TRANSFO_Z_3_B4
                  left3D(i,j,k,4) = $TRANSFO_Z_4_B4                 
                 enddo
                enddo
               enddo
            enddo            
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
            sb = sb + sum(Blocks(B:B+3))              
          enddo
      endif
      !-------------------------------------------------------------------------
      ! Merge the blocks associated with a linear, hermitian symmetry
      HFBlocks(1) = blocks(1) + blocks(3)
      HFBlocks(2) = blocks(2) + blocks(4)
      HFBlocks(5) = blocks(5) + blocks(7)
      HFBlocks(6) = blocks(6) + blocks(8)
      ! Remove the old blocks
      HFBlocks(3:4) = 0
      HFBlocks(7:8) = 0
    endif
    
  end subroutine Transformspwfs

  subroutine TransformInput(filenx,fileny,filenz,filenwn,filenwp, filedx,      &
  &                         fileblocks, file_HFB_blocks, extraspwfs)
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
    integer, intent(in)        :: fileblocks(8),extraspwfs(8)
    integer, intent(in)        :: file_HFB_blocks(8)
    integer                    :: bogo_blocks(8)

    real(KIND=dp), intent(in)  :: filedx
    real(KIND=dp), allocatable :: extended(:,:,:), newenergy(:), temp(:)
    real(KIND=dp), allocatable :: temp2(:,:), U(:,:), V(:,:)
    
    real(KIND=dp), pointer :: Unew(:,:), Vnew(:,:)

    logical :: ChangeBoxSize = .false., transformed= .false.
    logical :: gradient_detected
    integer :: i,j, b, sb, sf
    integer :: N1HF   , N2HF   , N1F   , N2F   , TF, THF
    integer :: N1HF_sp, N2HF_sp, N1F_sp, N2F_sp, TF_sp, THF_sp
    

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
$TR          print *, 'Blocks 2,4,6,8 not allowed with time-reversal conserved.'
$TR          stop
$TR       endif

$NTR      do b=1,8,2
$NTR        if(extraspwfs(b) .ne. extraspwfs(b+1)) then
$NTR         print *, 'Spwf number with Rz = +i needs to match the number with Rz = -i.'
$NTR         print *, 'Block = ', B, ' extraspwfs = ', extraspwfs(b), extraspwfs(b+1)
$NTR         stop
$NTR        endif
$NTR      enddo

$PBROKEN  do b=1,8,4 ! Essentially isospin loop
$PBROKEN    if(extraspwfs(b+2).ne.0 .or. extraspwfs(b+3).ne. 0) then
$PBROKEN      print *, 'Parity is broken, so spwfs can only be added in the first few blocks.'
$PBROKEN      print *, 'Valid input is thus of the form'
$PBROKEN      print *, ' extraspwfs = a, b, 0, 0 , c, d, 0 ,0'
$PBROKEN      stop
$PBROKEN    endif
$PBROKEN  enddo

          !---------------------------------------------------------------------
          allocate(extended(nx*ny*nz,4,nwt)) ; allocate(newenergy(nwt))
          extended = 0.0
          hfblocks = fileblocks + extraspwfs ; newenergy = 1000.0
          
          gradient_detected =.false.
          do b=1,8
            if(fileblocks(B).ne.file_HFB_blocks(B)) then
              gradient_detected = .true.
            endif
          enddo
          
          if(gradient_detected) then
            bogo_blocks = file_HFB_blocks + extraspwfs
          else
            bogo_blocks = hfblocks
          endif
          
          sb = 0 ; sf = 0  
          do b = 1, 8
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
          ! Dealing with the quantities related to the pairing subproblem
          ! (1)  the pairing gaps
          ! (2)  the density matrix rho
          ! (3)  the anomalous density matrix
          ! (4)  the Bogoliubov transformation
          ! (5)  the configuration matrix
          ! (6)  the HF transformation
          ! (7)  the current single-particle hamiltonian
          select case (pairingtype)
          case(0)
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
            ! HF, nothing to do
          case(1)
            ! BCS, need to move the gaps into the correct position
            temp = BCSgaps ; deallocate(BCSgaps)
            allocate(BCSgaps(nwt)) ; BCSgaps = 0.0
    
            sb = 0 ; sf = 0  
            do b = 1,8,2
              do i=1, fileblocks(b)
                 BCSgaps(sb+i) = temp(sf+i)
              enddo
              sb  = sb + hfblocks(b)
              sf  = sf + fileblocks(b)
            enddo            
          case(2)
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ! HFB pairing 
            if(allocated(HFBgaps)) then
              !-----------------------------------------------------------------
              ! (1) Pairing gaps
              temp2 = HFBgaps ; deallocate(HFBgaps)
              allocate(HFBgaps(nwt,nwt)) ; HFBgaps = 0.0 
      
              sb = 0 ; sf = 0  
              do b = 1,8,2
                do i=1, fileblocks(b) 
                  do j=1, fileblocks(b+1) 
                   HFBgaps(sb+i, sb+hfblocks(b)+j) = &
                   &                             temp2(sf+i, sf+fileblocks(b)+j)
                   HFBgaps(sb+hfblocks(b)+j, sb+i) = &
                   &                             temp2(sf+fileblocks(b)+j, sf+i)
                  enddo
                enddo

                sb  = sb +    hfblocks(b) +   hfblocks(b+1)
                sf  = sf +  fileblocks(b) + fileblocks(b+1)
              enddo
              !-----------------------------------------------------------------
              ! (2) density matrix rho
              temp2 = rho_pairing; deallocate(rho_pairing)
              allocate(rho_pairing(nwt,nwt)) ; rho_pairing = 0.0
              
              sb = 0 ; sf = 0  
              do b = 1, 8
                do i=1, fileblocks(b)
                  do j=1, fileblocks(b) 
                   rho_pairing(sb+i, sb+j) = temp2(sf+i, sf+j)
                  enddo
                enddo
                sb  = sb +    hfblocks(b) 
                sf  = sf +  fileblocks(b) 
              enddo
              !-----------------------------------------------------------------
              ! (3) anomalous density matrix kappa
              temp2 = kappa_pairing; deallocate(kappa_pairing)
              allocate(kappa_pairing(nwt,nwt)) ; kappa_pairing = 0.0
              
              sb = 0 ; sf = 0  
              do b = 1, 8, 2
                do i=1, fileblocks(b) 
                  do j=1, fileblocks(b+1) 
                   kappa_pairing(sb+i, sb+hfblocks(b)+j) = &
                   &                             temp2(sf+i, sf+fileblocks(b)+j)
                   kappa_pairing(sb+hfblocks(b)+j, sb+i) = &
                   &                             temp2(sf+fileblocks(b)+j, sf+i)
                  enddo
                enddo
                sb  = sb +    hfblocks(b) +   hfblocks(b+1)
                sf  = sf +  fileblocks(b) + fileblocks(b+1)
              enddo
              !-----------------------------------------------------------------
              ! (4) the Bogoliubov transformation
              temp2 = Bogoliubov; deallocate(Bogoliubov)
              allocate(Bogoliubov(2*nwt,2*nwt)) ; Bogoliubov = 0.0
              
              sb = 0 ; sf = 0  
              do b = 1, 8, 2 
                !
                ! We need to keep track of a whole bunch of indexes:
                !
                ! *) single-particle indices (_sp)
                ! *) quasiparticle indices  (no suffix)
                ! *) file indices (F)
                ! *) memory indices (HF)
                !
                ! We need to take both single-particle and quasi-particle
                ! indices, as the Bogoliubov transformation on the file might
                ! have been generated by a gradient solver, and blocking 
                ! might have changed the effective qp-block size.                
                N1HF = bogo_Blocks(B)   ; N1F = file_HFB_blocks(b)    
                N2HF = bogo_Blocks(B+1) ; N2F = file_HFB_blocks(b+1)
                
                N1HF_sp = HFBlocks(B)   ; N1F_sp = fileblocks(b)    
                N2HF_sp = HFBlocks(B+1) ; N2F_sp = fileblocks(b+1)
                
                THF   = N1HF   +N2HF    ; TF   = N1F    + N2F
                THF_sp= N1HF_sp+N2HF_sp ; TF_sp= N1F_sp + N2F_sp

                ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                ! Copy the Bogoliubov transformation on file into the new one
                !
                !       (  V^*+   0     U+  0   )
                !  W =  (  0      V^*-  0   U-  )      
                !       (  0      U^*-  0   V-  )
                !       (  U^*+   0     V+  0   )
                !
                ! Note, we only fill in the r.h.s., as that is the only part
                ! that gets propagated by the code. 
                ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                ! U^+ and V^+
                U = temp2(sf             +1:sf+N1F_sp,sf+TF+1:sf+TF+N1F)  
                V = temp2(sf+TF_sp+N2F_sp+1:         ,sf+TF+1:sf+TF+N1F)  
                Unew => Bogoliubov(sb               +1:sb+  THF_sp,sb+THF+1:sb+THF+N1HF)
                Vnew => Bogoliubov(sb+THF_sp+N2HF_sp+1:sb+2*THF_sp,sb+THF+1:sb+THF+N1HF) 

                Unew(1:N1F_sp, 1:N1F) = U
                Vnew(1:N1F_sp, 1:N1F) = V
                
                !---------------------------------------------------------------
                ! I initially tried to set up some random numbers in the
                ! new U and V columns, but this made for a quite unpredictable
                ! start of the iterative process for a gradient solver; either
                !
                ! (i)  I obtained a gas-like solution with all levels
                !      fractionally occupied
                ! (ii) The energy moved by several tens of MeV for just a few
                !      iterations. 
                !---------------------------------------------------------------
                
!                ! We have the old information in the right places, but now 
!                ! we need to fill something in the new columns
!                call random_number(Unew(N1F_sp+1:N1HF_sp,N1F+1:N1HF))
!                call random_number(Vnew(N1F_sp+1:N1HF_sp,N1F+1:N1HF))

!                ! We have the old information in the right places, but now 
!                ! we need to fill something in the new columns
!                Unew(N1F_sp+1:N1HF_sp,N1F+1:N1HF) = &
!                &                    Unew(N1F_sp+1:N1HF_sp,N1F+1:N1HF) * 0.001
!                Vnew(N1F_sp+1:N1HF_sp,N1F+1:N1HF) = &
!                &                    Vnew(N1F_sp+1:N1HF_sp,N1F+1:N1HF) * 0.001

                ! And keep this levels primarily as "unoccupied"
                do i=1,extraspwfs(b)
                    Unew(N1F_sp+i,N1F+i) =  Unew(N1F_sp+i,N1F+i)+ 1.0
                enddo

                ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                ! U^- and V^-
                U = temp2(sf+N1F_sp+1:sf+N1F_sp+N2F_sp,sf+TF+N1F+1:)  
                V = temp2(sf+ TF_sp+1:sf+ TF_sp+N1F_sp,sf+TF+N1F+1:)  

                Unew => Bogoliubov(sb+N1HF_sp+1:sb+THF_sp        ,sb+THF+N1HF+1:)
                Vnew => Bogoliubov(sb+ THF_sp+1:sb+THF_sp+N1HF_sp,sb+THF+N1HF+1:) 

                Unew(1:N2F_sp, 1:N2F) = U
                Vnew(1:N2F_sp, 1:N2F) = V

                !---------------------------------------------------------------
                ! Random part inactive
                !---------------------------------------------------------------
!                ! We have the old information in the right places, but now 
!                ! we need to fill something in the new columns
!                call random_number(Unew(N2F_sp+1:N2HF_sp,N2F+1:N2HF))
!                call random_number(Vnew(N2F_sp+1:N2HF_sp,N2F+1:N2HF))

!                ! But make these random numbers not mess up the calculation 
!                ! TOO much
!                Unew(N2F_sp+1:N2HF_sp,N2F+1:N2HF) = &
!                &                    Unew(N2F_sp+1:N2HF_sp,N2F+1:N2HF) * 0.001
!                Vnew(N2F_sp+1:N2HF_sp,N2F+1:N2HF) = &
!                &                    Vnew(N2F_sp+1:N2HF_sp,N2F+1:N2HF) * 0.001

                ! And keep this levels primarily as "unoccupied"
                do i=1,extraspwfs(b+1)
                    Unew(N2F_sp+i,N2F+i) = Unew(N2F_sp+i,N2F+i) + 1.0
                enddo

                sb  = sb + 2*bogo_blocks(b)     + 2*bogo_blocks(b+1)
                sf  = sf + 2*file_HFB_blocks(b) + 2*file_HFB_blocks(b+1)
              enddo

!             ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!             ! Debugging print statements
!             ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!              sb = 0 ; sf = 0  
!              do b = 1, 8, 2 
!                N1HF = bogo_Blocks(B)   ; N1F = file_HFB_blocks(b)    
!                N2HF = bogo_Blocks(B+1) ; N2F = file_HFB_blocks(b+1)
!                
!                N1HF_sp = HFBlocks(B)   ; N1F_sp = fileblocks(b)    
!                N2HF_sp = HFBlocks(B+1) ; N2F_sp = fileblocks(b+1)
!                
!                THF   = N1HF   +N2HF    ; TF   = N1F    + N2F
!                THF_sp= N1HF_sp+N2HF_sp ; TF_sp= N1F_sp + N2F_sp
!              
!                print *, 'BOGO'
!                do i=1,2*(N1HF_sp+N2HF_sp)
!                  print ('(99f10.3)'), Bogoliubov(sb+i, sb+THF+1:sb+2*THF)
!                enddo
!                print *
!                print *, 'TEMP2'
!                do i=1,2*(N1F_sp+N2F_sp)
!                  print ('(99f10.3)'), temp2(sf+i, sf+TF+1:sf+2*TF)
!                enddo
!                print *
!                sb  = sb + 2*bogo_blocks(b)   + 2*bogo_blocks(b+1)
!                sf  = sf + 2*file_HFB_blocks(b) + 2*file_HFB_blocks(b+1)
!              enddo
!             ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
              
              ! And now make sure all columns are orthogonal
              ! (NOTE: this routine only works on the right-most half of the 
              !        transformation, but that is sufficient)
              call ortho_bogo(Bogoliubov, HFBlocks)
              
              !-----------------------------------------------------------------
              ! (5) the configuration matrix
              temp = configmatrix 
              deallocate(configmatrix) ; allocate(configmatrix(2*nwt))
              sb = 0 ; sf = 0
              do b=1,8,2
                N1HF = bogo_Blocks(B)   ; N1F = file_HFB_blocks(b)    
                N2HF = bogo_Blocks(B+1) ; N2F = file_HFB_blocks(b+1)
                THF  = N1HF+N2HF     ; TF  = N1F + N2F

                configmatrix(sb     +1:sb+N1F)       = temp(sf    +1:sf+N1F)
                configmatrix(sb+N1HF+1:sb+N1HF+N2F ) = temp(sf+N1F+1:sf+TF)

                configmatrix(sb+THF     +1:sb+THF+N1F) = &
                &                                  temp(sf+TF    +1:sf+  TF+N1F)
                configmatrix(sb+THF+N1HF+1:sb+THF+N1HF+N2F) = &
                &                                  temp(sf+TF+N1F+1:sf+2*TF)

                ! Filling in the occupation numbers for the new qps
                configmatrix(sb+N1F+1:sb+N1HF)        = 0
                configmatrix(sb+N1HF+N2F+1:sb+THF)    = 0

                configmatrix(sb+THF+N1F     +1:sb+THF+N2HF)= 1
                configmatrix(sb+THF+N1HF+N2F+1:sb+2*THF)   = 1
              
                sb  = sb + 2*bogo_blocks(b)   + 2*bogo_blocks(b+1)
                sf  = sf + 2*file_HFB_blocks(b) + 2*file_HFB_blocks(b+1)
              enddo
              
              !-----------------------------------------------------------------
              ! (6) The HF transformation
              temp2 = HFtransfo
              deallocate(hftransfo) ; allocate(hftransfo(nwt,nwt))
              hftransfo = 0
              
              sb = 0; sf = 0
              do b=1,8
                  hftransfo(sb+1:sb+fileblocks(b),sb+1:sb+fileblocks(b)) &
                  & = temp2(sf+1:sf+fileblocks(b),sf+1:sf+fileblocks(b))
                  do i=1,extraspwfs(b)
                      hftransfo(sb+fileblocks(b)+i,sb+fileblocks(b)+i) = 1.0d0
                  enddo
                  
                  sb  = sb + hfblocks(b)
                  sf  = sf + fileblocks(b)
              enddo
              !-----------------------------------------------------------------
              ! (7) The current single-particle hamiltonian 
              temp2 = current_sph
              deallocate(current_sph); allocate(current_sph(nwt,nwt))
              current_sph = 0.0d0
                
              sb = 0; sf = 0
              do b=1,8
                  current_sph(sb+1:sb+fileblocks(b),sb+1:sb+fileblocks(b)) &
                  & = temp2(sf+1:sf+fileblocks(b),sf+1:sf+fileblocks(b))
                  do i=1,extraspwfs(b)
                      current_sph(sb+fileblocks(b)+i,sb+fileblocks(b)+i) = &
                      &                    spenergies(sb+fileblocks(b)+i)
                  enddo
                  
                  sb  = sb + hfblocks(b)
                  sf  = sf + fileblocks(b)
              enddo
              !-----------------------------------------------------------------
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
