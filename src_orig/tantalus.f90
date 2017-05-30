program Tantalus
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

 use compilation
 use geninfo
 use wavefunctions
 use IO
 
 implicit none

 100 format (/,' ___________________________________________________________', &
     &       /,'|                                                          |', &
     &       /,'|                                    version VERY UNSTABLE |', &
     &       /,'|                                                          |', &
     &       /,'|  #######   ##   #    # #####   ##   #      #    #  ####  |', &
     &       /,'|     #     #  #  ##   #   #    #  #  #      #    # #      |', &
     &       /,'|     #    #    # # #  #   #   #    # #      #    #  ####  |', &
     &       /,'|     #    ###### #  # #   #   ###### #      #    #      # |', &
     &       /,'|     #    #    # #   ##   #   #    # #      #    # #    # |', &
     &       /,'|     #    #    # #    #   #   #    # ######  ####   ####  |', &
     &       /,'|                                                          |', &
     &       /,'|  Copyright  P.-H. Heenen, M.Bender & W. Ryssens          |', &
     &       /,'|__________________________________________________________|')

 print 100
 
 call PrintInput
 call Test()
 
end program Tantalus


subroutine Test
     
    use compilation
    use derivatives
    use wavefunctions
    use constants
    use densities
    use hartreefock
    use functional
    
    implicit none
   
    integer :: i,lol,k, par,iso
 
    open (12,form='unformatted',file='MOCCa.test')
    read(12)
    read(12)
    read(12)
    read(12)
    read(12)
    allocate(Occupations(nwt))
    allocate(spenergies(nwt))
    allocate(HFPsi(nx,ny,nz,4,nwt))
    
    do i=1,nwt
        read(12) HFPsi(:,:,:,:,i)
        read(12), occupations(i), spenergies(i), lol, lol, lol, lol, lol, iso, k,par,k,k
    enddo
    
    nwn=10
    nwp=10
    
    HFBlocks(1) = 7
    HFBlocks(3) = 3
    HFBlocks(5) = 7
    HFBlocks(7) = 3
    
    !call iniwavefunctions()
    call inilag()
    call NaiveFill(occupations)
    
    call printSpwfs
    call deriveall()
    
    call calcedfcoefs()
    call printedfcoefs()     
    call densit
    Kinetic = CompKinetic()
    
    call CompSkyrme()
    call PrintSkyrme
    print *, Kinetic
!    print *
!    print *, 'Kinetic'
!    print *
!    print ('(10f12.3)'), Kinetic
!    print *
!    print *, 'NLO'
!    print *
!    print ('(10f12.3)'),  Skyrme_LO()
!    print *
!    print *, 'NLO'
!    print *
!    print ('(10f12.3)'),  Skyrme_NLO()
!    print *
!    print *, 'N2LO'
!    print *
!    print ('(10f12.3)'),  Skyrme_N2LO()
!    print *, hfbasis(:,1,1,1, 1)
!    print *, hfbasis(:,1,1,1, 11)

end subroutine Test
