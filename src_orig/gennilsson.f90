!==============================================================================
!  #######   ##   #    # #####   ##   #      #    #  ####
!     #     #  #  ##   #   #    #  #  #      #    # #
!     #    #    # # #  #   #   #    # #      #    #  ####
!     #    ###### #  # #   #   ###### #      #    #      #
!     #    #    # #   ##   #   #    # #      #    # #    #
!     #    #    # #    #   #   #    # ######  ####   ####
!
!  Copyright W. Ryssens & M. Bender
!------------------------------------------------------------------------------
! Small main routine to generate nilsson wavefunctions and output them
! separately.
!==============================================================================
module nilsson_extra

  implicit none

contains 

  subroutine nilsson_print(selected)
    !---------------------------------------------------------------------------
    ! 
    !---------------------------------------------------------------------------
    use wavefunctions
  
    1 format (22 ('-'), ' Sp wavefunctions ', 30('-'))
    2 format (65 ('-'))
    6 format (3x,'i',3x,'P', 1x, 'iso', 11x,'E',8x, 'JxT',7x, 'JyT', 7x ,'Jz', 9x, 'J')    

    integer, intent(in) :: selected
    integer       :: wave,k, iso
    !integer       :: ProtonOrder(nwp), NeutronOrder(nwn)
    real(KIND=dp) :: p, Jx, Jy, Jz, JJ  
    character(len=3) :: sel


    print 1
    print 6
    print 2
    do k=1,nwt 
        wave = k
        
        iso = +1
        if(wave .le. nwn) iso = -1
        
        p = +1
        if(k .gt. HFBlocks(1) .and. k .lt. sum(HFBlocks(1:3))) p = -1
        if(k .gt. sum(HFBlocks(1:5))) p = -1

        !-----------------------------------------------------------------------
        ! Depending on the symmetries, select different quantities to print 
        ! for Jx, Jy, Jz. Currently configured for EV8 mode!
        Jx = angmom_xt_real(HFPsi(:,:,wave),HFPsi(:,:,wave),HFdPsi(:,:,:,wave))
        Jy = angmom_yt_imag(HFPsi(:,:,wave),HFPsi(:,:,wave),HFdPsi(:,:,:,wave))
        Jz = angmom_z_real (HFPsi(:,:,wave),HFPsi(:,:,wave),HFdPsi(:,:,:,wave))

        JJ = & 
        &   angmom_x_quad(HFPsi(:,:,wave),HFdPsi(:,:,:,wave), &
        &                 HFPsi(:,:,wave),HFdPsi(:,:,:,wave)) &
        & + angmom_y_quad(HFPsi(:,:,wave),HFdPsi(:,:,:,wave), &
        &                 HFPsi(:,:,wave),HFdPsi(:,:,:,wave)) &
        & + angmom_z_quad(HFPsi(:,:,wave),HFdPsi(:,:,:,wave), &
        &                 HFPsi(:,:,wave),HFdPsi(:,:,:,wave)) 
        JJ = (-1. + sqrt(1. + 4*JJ))/2.

        sel = ''
        if(wave .eq. selected) sel = '(*)'
        print ('(3i4,3x, a3, 99f10.3)'), wave, int(p), iso, sel,spenergies(wave),   & 
        &                                                        Jx, Jy, Jz, JJ
    enddo
    print 2

    wave = selected
       
    p = +1
    if(k .gt. HFBlocks(1) .and. k .lt. sum(HFBlocks(1:3))) p = -1
    if(k .gt. sum(HFBlocks(1:5))) p = -1

    !-----------------------------------------------------------------------
    ! Depending on the symmetries, select different quantities to print 
    ! for Jx, Jy, Jz. Currently configured for EV8 mode!
    Jx = angmom_xt_real(HFPsi(:,:,wave),HFPsi(:,:,wave),HFdPsi(:,:,:,wave))
    Jy = angmom_yt_imag(HFPsi(:,:,wave),HFPsi(:,:,wave),HFdPsi(:,:,:,wave))
    Jz = angmom_z_real (HFPsi(:,:,wave),HFPsi(:,:,wave),HFdPsi(:,:,:,wave))
  

    JJ = & 
    &   angmom_x_quad(HFPsi(:,:,wave),HFdPsi(:,:,:,wave), &
    &                 HFPsi(:,:,wave),HFdPsi(:,:,:,wave)) &
    & + angmom_y_quad(HFPsi(:,:,wave),HFdPsi(:,:,:,wave), &
    &                 HFPsi(:,:,wave),HFdPsi(:,:,:,wave)) &
    & + angmom_z_quad(HFPsi(:,:,wave),HFdPsi(:,:,:,wave), &
    &                 HFPsi(:,:,wave),HFdPsi(:,:,:,wave)) 
    JJ = (-1. + sqrt(1. + 4*JJ))/2.

    sel = '*'
    print *, 'Selected:'
    print ('(2i4,3x, a3, 99f10.3)'), wave, int(p), sel,spenergies(wave),   & 
    &                                                        Jx, Jy, Jz, JJ
    print 2
  end subroutine nilsson_print


end module nilsson_extra

program generate_nilson

  use geninfo
  use wavefunctions
  use derivatives
  use nil8
  use nilsson_extra

  implicit none

  1 format (3i3, f8.3, 2i3)
  2 format (99f18.15)

  10 format('-----------------------------------------------------------------')
  11 format(' Nilsson generating program')
  12 format(' (nx,ny,nz,dx)    = ', 3i3, f8.3, ' fm')
  13 format(' oscillator freq. = ', 3f8.3)
  14 format(' nwn, nwp         = ', 2i3)
  15 format(' N, Z             = ', 2i3)
  16 format(' filename         = ', a40)
  17 format(' spwf selected    = ', i3)
  integer                    ::  npp, npn, selection = 1
  real(KIND=dp)              ::  osc_x, osc_y, osc_z

  character(len=40)          ::  fname = 'model.spwf'

  ! Practical redefinition
  real(KIND=dp), pointer             :: wf3d(:,:,:)
  integer :: i,j,k,l, wave, par, it

  !-----------------------------------------------------------------------------
  ! Input phase
  namelist /nil/  nx, ny, nz, dx, neutrons, protons, nwn, nwp,  &
  &                   osc_x, osc_y, osc_z, fname, selection

  read(unit=*, NML=nil)
  print 10
  print 11
  print 12,nx,ny,nz,dx
  print 13, osc_x, osc_y, osc_z
  print 14, nwn, nwp
  print 15, int(neutrons), int(protons)
  print 16, fname
  print 17, selection
  print 10

  ! Dealing with input in a better way
  osc_freq(1) =osc_x ;  osc_freq(2) =osc_y ; osc_freq(3) =osc_z
  npp      = int(protons)  ;  npn      = int(neutrons)  
  nwt = nwn + nwp ; nwt_local = nwt ; HFBlocks_global = HFBLocks

  !-----------------------------------------------------------------------------
  ! Initializing everything in the code
  dv = (dx**3)*8
  mv = nx*ny*nz
  call inimesh(meshx, meshy, meshz, nx, ny,nz, meshgrid, 0.0d0,0.0d0, 0.0d0)
  call inimesh(meshx_shifted, meshy_shifted, meshz_shifted, nx, ny,nz, meshgrid_shifted, 0.0d0,0.0d0, 0.0d0)
  call iniwavefunctions(nx,ny,nz,nwn,nwp)
  call inilag()                               ! Initialize derivative matrices. 
  call add_timer('HF-basis Derivatives'       , T_derivatives)  
  call deriveHF()
  !-----------------------------------------------------------------------------

  wave = selection
  if(wave.gt.sum(hfblocks(1:3))) then
    it = 2 
    par=+1
    if(wave.gt.sum(hfblocks(1:5))) then
       par = -1
    endif
  else
    it  =  1
    par = +1
    if(wave.gt.(hfblocks(1))) then
       par = -1
    endif
  endif


  call nilsson_print(wave)
  !-----------------------------------------------------------------------------
  ! Writing the selected spwf to file
  open(unit=6, file = fname)


  write(unit=6, fmt=1) nx,ny,nz,dx,it, par
  do l=1,4
    wf3d(1:nx,1:ny,1:nz) => hfpsi(1:nx*ny*nz,l,wave)
    do k=1,nz
      do j=1,ny
        do i=1, nx
          write(unit=6, fmt=2) wf3d(i,j,k)
        enddo
      enddo
    enddo 
  enddo
  close(unit=6)

end program generate_nilson

