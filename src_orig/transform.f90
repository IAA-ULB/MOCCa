!===============================================================================
!     __  __  ___   ____ ____
!    |  \/  |/ _ \ / ___/ ___|__ _
!    | |\/| | | | | |  | |   / _` |
!    | |  | | |_| | |__| |__| (_| |
!    |_|  |_|\___/ \____\____\__,_|
!
! Written mainly by W. Ryssens & M. Bender
!
! Opensource software distributed under the GNU AGPLv3 licence, see the
!  LICENCE file in the root of this project.
!===============================================================================
module transform
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
 ! TR          : $TR
 ! NTR         : $NTR
 !==============================================================================
  use geninfo
  use wavefunctions
  use pairing

  implicit none

  ! Indicates whether we need to transform the spwfs read on input.
  ! TODO: refactor things such that symtransfo_needed is no longer a global setting!
  logical :: sym_transfo_needed = .false.

  interface changeboxsize_function
     module procedure changeboxsize_function_complex
     module procedure changeboxsize_function_real
  end interface

contains

  subroutine Transformspwfs(wfs,oldnx, oldny, oldnz,                          &
  &                       blocks, blocks_local, old_rank_map, old_spwf_inverse)
    !---------------------------------------------------------------------------
    ! Transform a set of spwfs (wfs) with a given set of symmetries as 
    ! characterized by a set of symmetry classes (blocks) that are defined on 
    ! a mesh with (oldnx,oldny,oldnz) points into a set of spwfs with the 
    ! symmetry choices of the calculation on the new mesh.
    !
    ! This routine performs EXACTLY ONE of the following actions:
    !  (1) Breaks a single non-spatial symmetry (like T-reversal)
    !  (2) Breaks a single spatial symmetry (like P)
    ! but it cannot do more in a single go.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    !
    ! Input:
    !   wfs          : set of old spwfs to be transformed
    !   oldnx/y/z    : parameters of the original mesh
    !   blocks       : symmetry blocks of the old spwfs. These need to be GLOBAL
    !                  values, i.e. the total number of spwfs in a block across
    !                  all MPI ranks. 
    !   blocks_local : symmetry blocks of the old spwfs. 
    !   old_rank_map     : MPI rank mapping for the old spwfs
    !   old_spwf_inverse : inverse mapping for the old spwfs
    !
    ! Output:
    !   wfs       : set of transformed spwfs
    !
    ! Side-effects:
    !   The code effects a new load balancing of the spwfs among the available
    !   MPI ranks. At the start of this routine, there is a set of wfs that are
    !   dispersed in a certain way among the available MPI ranks. The new set
    !   of symmetry blocks is determined, prompting another look at the 
    !   load balancing and a reshuffling of the spwfs among the ranks.
    !---------------------------------------------------------------------------
  1 format("------------------------------------------------------------------")
  2 format(" Number of spwfs on STDIN does not match what is required.")
  3 format(" Symmetry-broken calculation should have (nwn,nwp) = ", 2i5)
  4 format(" STDIN says                              (nwn,nwp) = ", 2i5)
  5 format &
     &  (/,8x,' ___________________________________________________________', &
     &   /,8x,'| Symmetry-transformation of the input                     |', &
     &   /,8x,'| On file:                                                 |', &
     &   /,8x,'|    nx, ny, nz    = ', 3i5, '                       |', &
     &   /,8x,'|    dx            = ', f10.7, ' (fm)                       |', &
     &   /,8x,'|    nwn, nwn      = ', 2i7 , '                        |', &
     &   /,8x,'|    (n+,n-)       = (', 4i5, ')                |', &
     &   /,8x,'|    (p+,p-)       = (', 4i5, ')                |', &
     &   /,8x,'| This calculation:                                        |' , &
     &   /,8x,'|    nx, ny, nz    = ', 3i5, '                       |'       , &
     &   /,8x,'|    dx            = ', f10.7, '(fm)                        |', &
     &   /,8x,'|    nwn, nwn      = ', 2i7, '                        |'      , &
     &   /,8x,'|    (n+,n-)       = (', 4i5, ')                |', &
     &   /,8x,'|    (p+,p-)       = (', 4i5, ')                |', &
     &   /,8x,'|__________________________________________________________|')


    real(KIND=dp), intent(inout), allocatable, target :: wfs(:,:,:)
    integer, intent(inout)                    :: blocks(8), blocks_local(8)
    integer, intent(in)                       :: oldnx, oldny, oldnz
    integer, intent(in)                       :: old_rank_map(:)
    integer, intent(in)                       :: old_spwf_inverse(:)

    integer, allocatable                      :: prev_rank_map(:)
    integer, allocatable                      :: prev_spwf_inverse(:)
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
    integer  :: recv_rank, send_rank, tempind, wfind
#if(USE_MPI>0)
    integer :: mpi_err
#endif

    if(.not. allocated(rho_can)) then
        ! There is one case where this array might not be allocated upon entry
        ! in this routine: when initializing from scratch
        allocate(rho_can(sum(blocks))) ; rho_can = 0.0d0
    endif

    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Decide on what the new block structure will be
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if( $NONSPATIAL ) then
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      ! We will be using an antilinear, antihermitian symmetry operator to 
      ! double the total number of spwfs and construct partner states.
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      ! A sanity check
      if((2*sum(blocks(1:4)) .ne. nwn)  .or. (2*sum(blocks(5:8)) .ne. nwp)) then
        print 1
        print 2
        print 3, 2*sum(blocks(1:4)), 2*sum(blocks(5:8))
        print 4, nwn, nwp
        print 1
        call stp('')
      endif

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      ! Create the second block in every pair
      HFBlocks_global(1) = blocks(1) ; HFBlocks_global(2) = blocks(1)
      HFBlocks_global(3) = blocks(3) ; HFBlocks_global(4) = blocks(3)
      HFBlocks_global(5) = blocks(5) ; HFBlocks_global(6) = blocks(5)
      HFBlocks_global(7) = blocks(7) ; HFBlocks_global(8) = blocks(7)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    endif 
    if( $SPATIAL ) then 
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      ! Merge the blocks associated with a linear, hermitian symmetry
      HFBlocks_global(1) = blocks(1) + blocks(3)
      HFBlocks_global(2) = blocks(2) + blocks(4)
      HFBlocks_global(5) = blocks(5) + blocks(7)
      HFBlocks_global(6) = blocks(6) + blocks(8)
      ! Remove the old blocks
      HFBlocks_global(3:4) = 0
      HFBlocks_global(7:8) = 0
    endif

    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! .. and now we have each rank decide what (transformed) spwfs to take 
    !    from file
    call loadbalance(HFBlocks_global,                 &               ! inputs
    &       HFblocks, spwf_map, rank_map, spwf_inverse)               ! outputs
    ! this particular will hold nwt_local spwfs at the end of the transformation
    nwt_local = sum(HFblocks)

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Copying old information into temporary arrays
    temp        = wfs         ; tempe   = spenergies
    tempd       = dispersions ; tempr   = rho_can
    temptransfo = hftransfo   ; tempsph = sphamil
    !... and making space for the new set of spwfs
    deallocate(wfs)         ; allocate(wfs(nx*ny*nz,4,nwt_local))
    deallocate(dispersions) ; allocate(dispersions(nwt))    ; dispersions  = 0
    deallocate(spenergies)  ; allocate(spenergies(nwt))     ; spenergies   = 0
    deallocate(rho_can)     ; allocate(rho_can(nwt))        ; rho_can      = 0
    deallocate(hftransfo)   ; allocate(hftransfo(nwt,nwt))  ; hftransfo    = 0
    deallocate(sphamil) ; allocate(sphamil(nwt,nwt)); sphamil  = 0

    if(pairingtype.eq.2) then
      ! The configmatrix is not necessarily initialised, hence a few more lines
      ! of code to deal with it.
      if(allocated(configmatrix)) then
        tempconfig  = configmatrix
        deallocate(configmatrix)
      else
        allocate(tempconfig(2*nwt)) ; tempconfig = 0.0d0
      endif
      allocate(configmatrix(2*nwt)) ; configmatrix = 0
    endif

    !---------------------------------------------------------------------------
    ! Use an antilinear, antihermitian symmetry operator (usually time-reversal)
    ! to double the total number of spwfs and construct the partner states.
    !---------------------------------------------------------------------------
    if( $NONSPATIAL ) then
     N = size(old_rank_map)
     allocate(prev_rank_map(2*N), prev_spwf_inverse(2*N))

     ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     ! First make sure that the code knows how the spwfs will be distributed
     ! after the transformation
     sb = 0 ; si = 0
     do B=1,8
          N = blocks(B) ; if (N.eq.0) cycle ! -> global index

          ! Both the original |psi> and T|psi> are stored on the same MPI rank
          prev_rank_map(sb  +1:sb+  N) = old_rank_map(si+1:si+N)
          prev_rank_map(sb+N+1:sb+2*N) = old_rank_map(si+1:si+N)

          ! .... but their relative indices in the total calculation change
          prev_spwf_inverse(sb  +1:sb+  N) = &
          &               old_spwf_inverse(si+1:si+N) + sum(blocks_local(1:B-1))
          prev_spwf_inverse(sb+N+1:sb+2*N) = &
          &               old_spwf_inverse(si+1:si+N) + sum(blocks_local(1:B))
          
          si = si +     N
          sb = sb + 2 * N
     enddo
     ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

     si = 0 ; sb = 0
     do B = 1,8
        !-----------------------------------------------------------------------
        ! Transformation of the spwfs that are dispersed over MPI ranks
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ! Note that blocks = 2,4,6,8 are always of zero size in this loop, as we
        ! are breaking the antilinear, antihermitian conserved symmetry
        N = blocks_local(B); if(N .eq. 0) cycle   ! -> local index

        ! Copy the wavefunctions that were already in storage
        wfs(:,:,sb+1:sb+N)    = temp(:,:,si+1:si+N)

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
        si = si +   N
        sb = sb + 2*N
      enddo

      si = 0 ; sb = 0
      do B = 1,8 
        !-----------------------------------------------------------------------
        ! Transformation of quantities that are shared across all processors
        !-----------------------------------------------------------------------
        N = blocks(B); if(N .eq. 0) cycle   ! -> global index

        dispersions(sb+1:sb+N)= tempd(si+1:si+N)
        spenergies(sb+1:sb+N) = tempe(si+1:si+N)
        rho_can(sb+1:sb+N)    = tempr(si+1:si+N) /2.0 ! Note the factor 1/2

        dispersions(sb+N+1:sb+2*N) = tempd(si+1:si+N)
        spenergies(sb+N+1 :sb+2*N) = tempe(si+1:si+N)
        rho_can(sb+N+1:sb+2*N)     = tempr(si+1:si+N)/2.0 ! Note the factor 1/2

        !-------------------------------------------------------------------
        ! The current HF transformation
        hftransfo(sb  +1:sb+  N,sb  +1:sb+  N) = &
        &                                   temptransfo(si+1:si+N,si+1:si+N)
        hftransfo(sb+N+1:sb+2*N,sb+N+1:sb+2*N) = &
        &                                   temptransfo(si+1:si+N,si+1:si+N)
        !-------------------------------------------------------------------
        ! The current single-particle hamiltonian
        sphamil(sb  +1:sb+  N,sb  +1:sb  +N) = &
        &                                       tempsph(si+1:si+N,si+1:si+N)
        sphamil(sb+N+1:sb+2*N,sb+N+1:sb+2*N) = &
        &                                       tempsph(si+1:si+N,si+1:si+N)

        si = si +   N 
        sb = sb + 2*N
      enddo

      !-------------------------------------------------------------------------
      ! Transformation of pairing quantities (shared across MPI ranks)
      ! (1)  the pairing gaps
      ! (2)  the density matrix rho
      ! (3)  the anomalous density matrix
      ! (4)  the Bogoliubov transformation
      ! (5)  the configuration matrix
      ! (6)  the HF transformation
      ! (7)  the current single-particle hamiltonian
      !
      ! Note that (6) and (7) have been moved above, since they need to always
      ! be performed, even if we are having new gaps initialized.
      !-------------------------------------------------------------------------
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
            N = blocks(B) ; if(N .eq. 0) cycle ! -> global index
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

            si = si +   N
            sb = sb + 2*N
            sc = sc + 4*N
          enddo
        endif
      endif
    endif
    !---------------------------END OF NONSPATIAL TRANSFORMATION ---------------

    !---------------------------------------------------------------------------
    ! Extend the calculation to a larger part of the simulation volume
    !
    if( $SPATIAL ) then
      ! Simple copies: breaking spatial symmetries does not create additional
      !                spwfs.
      prev_rank_map     = old_rank_map
      prev_spwf_inverse = old_spwf_inverse

      if($EXPANDX) then
          call stp('Extending to the full X-axis not implemented yet')
      elseif($EXPANDY) then
          call stp('Extending to the full Y-axis not implemented yet')
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

            hftransfo        (sb+offset_left +1:sb+offset_left +blocks(B), &
            &                 sb+offset_left +1:sb+offset_left +blocks(B)) &
            &   = temptransfo(sb+offset_left +1:sb+offset_left +blocks(B), &
            &                 sb+offset_left +1:sb+offset_left +blocks(B))

            sphamil      (sb+offset_left +1:sb+offset_left +blocks(B),  &
            &                 sb+offset_left +1:sb+offset_left +blocks(B))  &
            &   = tempsph    (sb+offset_left +1:sb+offset_left +blocks(B),  &
            &                 sb+offset_left +1:sb+offset_left +blocks(B))

            ! Third block on file becomes part of the first block
            offset_left  = blocks(B)
            offset_right = blocks(B) + blocks(B+1)
            dispersions(sb+offset_left +1:sb+offset_left +blocks(B+2)) &
            &   = tempd(sb+offset_right+1:sb+offset_right+blocks(B+2))
            spenergies (sb+offset_left +1:sb+offset_left +blocks(B+2)) &
            &   = tempe(sb+offset_right+1:sb+offset_right+blocks(B+2))
            rho_can    (sb+offset_left +1:sb+offset_left +blocks(B+2)) &
            &   = tempr(sb+offset_right+1:sb+offset_right+blocks(B+2))

            hftransfo        (sb+offset_left +1:sb+offset_left +blocks(B+2), &
            &                 sb+offset_left +1:sb+offset_left +blocks(B+2)) &
            &   = temptransfo(sb+offset_right+1:sb+offset_right+blocks(B+2), &
            &                 sb+offset_right+1:sb+offset_right+blocks(B+2))

            sphamil      (sb+offset_left +1:sb+offset_left +blocks(B+2), &
            &                 sb+offset_left +1:sb+offset_left +blocks(B+2)) &
            &   = tempsph    (sb+offset_right+1:sb+offset_right+blocks(B+2), &
            &                 sb+offset_left +1:sb+offset_left +blocks(B+2))

            ! Second block on file becomes first part of second block
            offset_left  = blocks(B) + blocks(B+2)
            offset_right = blocks(B)
            dispersions(sb+offset_left +1:sb+offset_left +blocks(B+1)) &
            &   = tempd(sb+offset_right+1:sb+offset_right+blocks(B+1))
            spenergies (sb+offset_left +1:sb+offset_left +blocks(B+1)) &
            &   = tempe(sb+offset_right+1:sb+offset_right+blocks(B+1))
            rho_can    (sb+offset_left +1:sb+offset_left +blocks(B+1)) &
            &   = tempr(sb+offset_right+1:sb+offset_right+blocks(B+1))

            hftransfo        (sb+offset_left +1:sb+offset_left +blocks(B+1), &
            &                 sb+offset_left +1:sb+offset_left +blocks(B+1)) &
            &   = temptransfo(sb+offset_right+1:sb+offset_right+blocks(B+1), &
            &                 sb+offset_left +1:sb+offset_left +blocks(B+1))

            sphamil      (sb+offset_left +1:sb+offset_left +blocks(B+1), &
            &                 sb+offset_left +1:sb+offset_left +blocks(B+1)) &
            &   = tempsph    (sb+offset_right+1:sb+offset_right+blocks(B+1), &
            &                 sb+offset_right+1:sb+offset_right+blocks(B+1))

            ! Fourth  block on file becomes second part of second block
            offset_left  = sum(blocks(B:B+2))
            offset_right = sum(blocks(B:B+2))
            dispersions(sb+offset_left +1:sb+offset_left +blocks(B+3)) &
            &   = tempd(sb+offset_right+1:sb+offset_right+blocks(B+3))
            spenergies (sb+offset_left +1:sb+offset_left +blocks(B+3)) &
            &   = tempe(sb+offset_right+1:sb+offset_right+blocks(B+3))
            rho_can    (sb+offset_left +1:sb+offset_left +blocks(B+3)) &
            &   = tempr(sb+offset_right+1:sb+offset_right+blocks(B+3))

            hftransfo        (sb+offset_left +1:sb+offset_left +blocks(b+3), &
            &                 sb+offset_left +1:sb+offset_left +blocks(b+3)) &
            &   = temptransfo(sb+offset_right+1:sb+offset_right+blocks(b+3), &
            &                 sb+offset_right+1:sb+offset_right+blocks(b+3))

            sphamil      (sb+offset_left +1:sb+offset_left +blocks(b+3), &
            &                 sb+offset_right+1:sb+offset_right+blocks(b+3)) &
            &   = tempsph    (sb+offset_right+1:sb+offset_right+blocks(b+3), &
            &                 sb+offset_right+1:sb+offset_right+blocks(b+3))

            !-------------------------------------------------------------------
            ! ... and then start operations on the spwfs
            ! Step 1 = do all the symmetry transformations locally, i.e. on the 
            !          MPI rank containing the spwf.
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ! First block is not modified
            offset_left = 0
            offset_right= 0
            do wave=1,Blocks_local(B)  !  Local index
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
            offset_right= Blocks_local(B)+Blocks_local(B+1) ! Block B+2 from file
            offset_left = Blocks_local(B)
            do wave=1,Blocks_local(B+2) !  Local index
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
            offset_right = Blocks_local(B)                 ! Block B+1 from file
            offset_left  = Blocks_local(B) + Blocks_local(B+2)
            do wave=1,Blocks_local(B+1) !  Local index
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
            offset_right = sum(Blocks_local(B:B+2))        ! Block B+3 from file
            offset_left  = sum(Blocks_local(B:B+2))
            do wave=1,Blocks(B+3) !  Local index
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
            sb = sb + sum(Blocks_local(B:B+3))
          enddo
      endif
    endif
    !---------------------------END OF SPATIAL TRANSFORMATION ---------------

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Step 2 = load rebalancing: redivide the spwfs among all the MPI ranks.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    deallocate(temp) ; temp = wfs ! copy again to the temp array
    deallocate(wfs)  ; allocate(wfs(mv,4,nwt_local))
    do i=1,nwt ! Loop over all wavefunctions
      send_rank = prev_rank_map(i)
      recv_rank =      rank_map(i)

      wfind     =      spwf_inverse(i)
      tempind   = prev_spwf_inverse(i)

#if(USE_MPI>0)
      if((MPI_RANK .eq. recv_rank) .AND. (send_rank .eq. recv_rank)) then
          ! nothing to send or receive; a simple copy will work
          wfs(:,:,wfind) = temp(:,:,tempind)
      elseif(MPI_RANK.eq.recv_rank) then
          ! recv_rank receives a message into the wfs array
          call MPI_RECV(wfs(:,:,wfind), 4*mv, MPI_REAL8, send_rank,  &
          &                        1,MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_err)
      else if(MPI_RANK .eq. send_rank)  then
          ! send_rank sends an spwf from the temp array
          call MPI_SEND(temp(:,:,tempind), 4*mv, MPI_REAL8, recv_rank, &
          &                                         1, MPI_COMM_WORLD, mpi_err)
      endif
#else
      ! No MPI shenanigans
      wfs(:,:,wfind) = temp(:,:,tempind)
#endif
    enddo

    !---------------------------------------------------------------------------
    if(MPI_RANK.eq.0) then
      print 5, oldnx, oldny, oldnz, dx, sum(blocks(1:4)), sum(blocks(5:8)),    &
      &       blocks, nx,ny,nz,dx,nwn,nwp, hfblocks
    endif

    deallocate(temp) ! just for safety, because this is a huge array
  end subroutine Transformspwfs

  subroutine TransformInput(filenx,fileny,filenz,filenwn,filenwp, filedx,      &
  &                         fileblocks, file_HFB_blocks, file_spwf_map,        &
  &                          file_rank_map, file_spwf_inverse,extraspwfs)
    !---------------------------------------------------------------------------
    ! Transform the input from file to the parameters of the new calculation.
    ! Note that this means either
    !   (*) A change of mesh parameters
    !   (*) A change of number of wavefunctions
    !
    ! but all of this with the same conserved/broken symmetries. These
    ! manipulations are not compatible in the same run with the breaking
    ! additional symmetries, achieved by the Transformspwfs routine.
    !
    ! - - - - - - - - -- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! TODO: 
    !  - Document this routine 
    !  - Fix this routine to be MPI-compatible
    !---------------------------------------------------------------------------

 1 format &
     &  (/,8x,' ___________________________________________________________', &
     &   /,8x,'| Transformation of the input                              |', &
     &   /,8x,'| On file:                                                 |', &
     &   /,8x,'|    nx, ny, nz    = ', 3i5, '                       |', &
     &   /,8x,'|    dx            = ', f10.7, ' (fm)                       |', &
     &   /,8x,'|    nwn, nwn      = ', 2i7 , '                        |', &
     &   /,8x,'|    (n+,n-)       = (', 4i5, ')                |', &
     &   /,8x,'|    (p+,p-)       = (', 4i5, ')                |', &
     &   /,8x,'| This calculation:                                        |' , &
     &   /,8x,'|    nx, ny, nz    = ', 3i5, '                       |'       , &
     &   /,8x,'|    dx            = ', f10.7, '(fm)                        |', &
     &   /,8x,'|    nwn, nwn      = ', 2i7, '                        |'      , &
     &   /,8x,'|    (n+,n-)       = (', 4i5, ')                |', &
     &   /,8x,'|    (p+,p-)       = (', 4i5, ')                |', &
     &   /,8x,'|__________________________________________________________|')

    integer, intent(in)        :: filenx,fileny,filenz,filenwn, filenwp
    integer, intent(in)        :: fileblocks(8),extraspwfs(8)
    integer, intent(in)        :: file_HFB_blocks(8)
    integer, intent(in)        :: file_spwf_map(:), file_rank_map(:)
    integer, intent(in)        :: file_spwf_inverse(:)
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
      call stp(' Interpolation (changing of dx) not yet allowed.')
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
#if(USE_MPI>0)
    call stp('Changing spwf number is not allowed for MPI calculations.')
#endif
          !---------------------------------------------------------------------
          !  Sanity checks on the number of spwfs
          if(nwn .ne. sum(fileblocks(1:4)) + sum(extraspwfs(1:4)) ) then
            print *, ' nwn in this calculation = ', nwn
            print *, ' nwn on file             = ', sum(fileblocks(1:4))
            print *, ' extra neutron spwfs     = ', sum(extraspwfs(1:4))
            call stp('Inconsistent number of neutron wavefunction.')
          endif
          if(nwp .ne. sum(fileblocks(5:8)) + sum(extraspwfs(5:8)) ) then
            print *, ' nwp in this calculation = ', nwp
            print *, ' nwp on file             = ', sum(fileblocks(5:8))
            print *, ' extra proton spwfs      = ', sum(extraspwfs(5:8))
            call stp('Inconsistent number of proton wavefunctions.')
          endif

$TR       if((extraspwfs(2).ne.0) .or. &
$TR         &  (extraspwfs(4).ne.0) .or. &
$TR         &  (extraspwfs(6).ne.0) .or. &
$TR         &  (extraspwfs(8).ne.0) ) then
$TR         call stp('Blocks 2,4,6,8 not allowed with time-reversal conserved.')
$TR       endif

$NTR      do b=1,8,2
$NTR        if(extraspwfs(b) .ne. extraspwfs(b+1)) then
$NTR         call stp('Spwf number with Rz = +i needs to match the number with Rz = -i.')
$NTR        endif
$NTR      enddo

$PBROKEN  do b=1,8,4 ! Essentially isospin loop
$PBROKEN    if(extraspwfs(b+2).ne.0 .or. extraspwfs(b+3).ne. 0) then
$PBROKEN      print *, 'Parity is broken, so spwfs can only be added in the first few blocks.'
$PBROKEN      print *, 'Valid input is thus of the form'
$PBROKEN      print *, ' extraspwfs = a, b, 0, 0 , c, d, 0 ,0'
$PBROKEN      print *, extraspwfs
$PBROKEN      call stp('')
$PBROKEN    endif
$PBROKEN  enddo

          !---------------------------------------------------------------------
          allocate(extended(nx*ny*nz,4,nwt)) ; allocate(newenergy(nwt))
          extended = 0.0
          hfblocks = fileblocks + extraspwfs ; newenergy = 1000.0
          hfblocks_global = hfblocks

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
          ! have each rank decide what (transformed) spwfs to take from file
          call loadbalance(HFBlocks_global,                    &   ! inputs
          &       HFblocks, spwf_map, rank_map, spwf_inverse)      ! outputs
          ! this particular will hold nwt_local spwfs at the end of the transformation
          nwt_local = sum(HFblocks)
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
                U = temp2(sf             +1:sf+  N1F_sp,sf+TF+1:sf+TF+N1F)
                V = temp2(sf+TF_sp+N2F_sp+1:sf+2*TF_sp ,sf+TF+1:sf+TF+N1F)
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
                U = temp2(sf+N1F_sp+1:sf+N1F_sp+N2F_sp,sf+TF+N1F+1:sf+2*TF)
                V = temp2(sf+ TF_sp+1:sf+ TF_sp+N2F_sp,sf+TF+N1F+1:sf+2*TF)

                Unew => Bogoliubov(sb+N1HF_sp+1:sb+THF_sp        ,sb+THF+N1HF+1:sb+2*THF)
                Vnew => Bogoliubov(sb+ THF_sp+1:sb+THF_sp+N2HF_sp,sb+THF+N1HF+1:sb+2*THF)

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
              temp2 = sphamil
              deallocate(sphamil); allocate(sphamil(nwt,nwt))
              sphamil = 0.0d0

              sb = 0; sf = 0
              do b=1,8
                  sphamil(sb+1:sb+fileblocks(b),sb+1:sb+fileblocks(b)) &
                  & = temp2(sf+1:sf+fileblocks(b),sf+1:sf+fileblocks(b))
                  do i=1,extraspwfs(b)
                      sphamil(sb+fileblocks(b)+i,sb+fileblocks(b)+i) = &
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
      hfblocks     = fileblocks
      spwf_map     = file_spwf_map
      rank_map     = file_rank_map
      spwf_inverse = file_spwf_inverse
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

  function changeboxsize_function_real(f, filenx, fileny, filenz) result(ft)
    !---------------------------------------------------------------------------
    ! Transform a real function defined on a mesh characterized by
    !
    !         (filenx, fileny, filenz, dx)
    !
    ! to a mesh defined by
    !         (    nx,     ny,      nz,dx)
    !
    ! with same mesh size dx.
    !
    ! Input:
    !   filenx, fileny, filenz : dimensions of the input function
    !   f                      : input function
    !
    ! Output:
    !   ft                     : output function
    !
    !---------------------------------------------------------------------------
    integer, intent(in)                :: filenx, fileny, filenz
    real(KIND=dp), intent(in), target  :: f(:)
    real(KIND=dp), allocatable, target :: ft(:)
    real(KIND=dp), pointer             :: f3 (:,:,:), ft3(:,:,:)
    integer                            :: endx, endy
    integer                            :: endz_left, endz_right
    integer                            :: startz_right, startz_left

    allocate(ft(nx*ny*nz)) ; ft = 0

    ft3(1:nx, 1:ny, 1:nz)             => ft(1:nx*ny*nz)
    f3(1:filenx, 1:fileny, 1:filenz)  => f (1:filenx*fileny*filenz)

    endx = min(nx, filenx)
    endy = min(ny, fileny)

    if (filenz .gt. nz) then
        ! Removing points along the z-axis
        startz_right  = 1
        endz_right    = nz

        startz_left   = 1
        endz_left     = nz

$PBROKEN  startz_right = 1  + (filenz - nz)/2
$PBROKEN  endz_right   = nz + (filenz - nz)/2
    elseif(nz .gt. filenz) then
        ! Adding points along the z-axis
        startz_right  = 1
        endz_right    = filenz

        startz_left   = 1
        endz_left     = filenz

$PBROKEN  startz_left  = 1         +   (nz - filenz)/2
$PBROKEN  endz_left    = filenz    +   (nz - filenz)/2

    else
       startz_right = 1 ; endz_right = nz
       startz_left  = 1 ; endz_left  = nz
    endif

    ft3(1:endx,1:endy,startz_left:endz_left) =&
    &                                  f3(1:endx,1:endy,startz_right:endz_right)

  end function changeboxsize_function_real

  function changeboxsize_function_complex(f, filenx, fileny, filenz) result(ft)
    !---------------------------------------------------------------------------
    ! Transform a complex function defined on a mesh characterized by
    !
    !         (filenx, fileny, filenz, dx)
    !
    ! to a mesh defined by
    !         (    nx,     ny,      nz,dx)
    !
    ! with same mesh size dx.
    !
    !
    ! Input:
    !   filenx, fileny, filenz : dimensions of the input function
    !   f                      : input function
    !
    ! Output:
    !   ft                     : output function
    !
    !---------------------------------------------------------------------------
    integer, intent(in)                   :: filenx, fileny, filenz
    complex(KIND=dp), intent(in), target  :: f(:)
    complex(KIND=dp), allocatable, target :: ft(:)
    complex(KIND=dp), pointer             :: f3 (:,:,:), ft3(:,:,:)
    integer                               :: endx, endy
    integer                               :: endz_left, endz_right
    integer                               :: startz_right, startz_left

    allocate(ft(nx*ny*nz)) ; ft = 0

    ft3(1:nx, 1:ny, 1:nz)             => ft(1:nx*ny*nz)
    f3(1:filenx, 1:fileny, 1:filenz)  => f (1:filenx*fileny*filenz)

    endx = min(nx, filenx)
    endy = min(ny, fileny)

    if (filenz .gt. nz) then
        ! Removing points along the z-axis
        startz_right  = 1
        endz_right    = nz

        startz_left   = 1
        endz_left     = nz

$PBROKEN  startz_right = 1  + (filenz - nz)/2
$PBROKEN  endz_right   = nz + (filenz - nz)/2
    elseif(nz .gt. filenz) then
        ! Adding points along the z-axis
        startz_right  = 1
        endz_right    = filenz

        startz_left   = 1
        endz_left     = filenz

$PBROKEN  startz_left  = 1         +   (nz - filenz)/2
$PBROKEN  endz_left    = filenz    +   (nz - filenz)/2

    else
       startz_right = 1 ; endz_right = nz
       startz_left  = 1 ; endz_left  = nz
    endif

    ft3(1:endx,1:endy,startz_left:endz_left) =&
    &                                  f3(1:endx,1:endy,startz_right:endz_right)

  end function changeboxsize_function_complex

end module transform
