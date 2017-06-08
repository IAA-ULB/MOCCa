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
    !
    ! Testing subroutine for various uses.
    !
    use compilation
    use derivatives
    use wavefunctions
    use constants
    use densities
    use hartreefock
    use functional
    use evolution
    
    implicit none
   
    integer :: i,lol,k, par,iso, iter
    
 
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
    
    allocate(sx(4,nwt), sy(4,nwt), sz(4,nwt))
    
    do i=1, HFBlocks(1)
        sx(1,i) =  1 ; sy(1,i) = +1 ; sz(1,i) = +1
        sx(2,i) = -1 ; sy(2,i) = -1 ; sz(2,i) = +1 
        sx(3,i) = -1 ; sy(3,i) = +1 ; sz(3,i) = -1
        sx(4,i) =  1 ; sy(4,i) = -1 ; sz(4,i) = -1
    enddo
    do i=HFBlocks(1) + 1,HFBlocks(1) + HFBlocks(3)
        sx(1,i) =  1 ; sy(1,i) = +1 ; sz(1,i) = -1
        sx(2,i) = -1 ; sy(2,i) = -1 ; sz(2,i) = -1 
        sx(3,i) = -1 ; sy(3,i) = +1 ; sz(3,i) = +1
        sx(4,i) =  1 ; sy(4,i) = -1 ; sz(4,i) = +1
    enddo
    do i=HFBlocks(1) + HFBlocks(3)+1,HFBlocks(1) + HFBlocks(3) +HFBlocks(5)
        sx(1,i) =  1 ; sy(1,i) = +1 ; sz(1,i) = +1
        sx(2,i) = -1 ; sy(2,i) = -1 ; sz(2,i) = +1 
        sx(3,i) = -1 ; sy(3,i) = +1 ; sz(3,i) = -1
        sx(4,i) =  1 ; sy(4,i) = -1 ; sz(4,i) = -1
    enddo
    do i=HFBlocks(1)+HFBlocks(3)+HFBlocks(5) + 1,                      &
    &       HFBlocks(1)+HFBlocks(3)+HFBlocks(5) + HFBLocks(7)
        sx(1,i) =  1 ; sy(1,i) = +1 ; sz(1,i) = -1
        sx(2,i) = -1 ; sy(2,i) = -1 ; sz(2,i) = -1 
        sx(3,i) = -1 ; sy(3,i) = +1 ; sz(3,i) = +1
        sx(4,i) =  1 ; sy(4,i) = -1 ; sz(4,i) = +1
    enddo
    
    
    
    !call iniwavefunctions()
    call inilag()
    call NaiveFill(occupations)
    
    call printSpwfs
    call deriveall()
    
    call calcedfcoefs()
    call printedfcoefs()     
    call densit(0)
    
    Kinetic = CompKinetic()
    call CompSkyrme()
    
    call PrintSkyrme

    call calcFields()

    do iter=1,1000 
        print *,  '*************************************'
        print *,  ' Iteration ', iter
        print *,  '*************************************'
        
        call Evolve_graddesc(iter)
        call deriveall()
        call densit(iter)
        call calcFields()
        if(mod(iter,100).eq.0) then
            call CompSkyrme()
            Kinetic = CompKinetic()
            call PrintSkyrme
        endif
    enddo
    

 
end subroutine Test
