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
 
 call ReadInput
 call PrintInput
 call ReachForWaterAndFood
 
end program Tantalus


subroutine ReachForWaterAndFood
    use compilation
    use derivatives
    use wavefunctions
    use constants
    use densities
    use hartreefock
    use functional
    use evolution
    
    implicit none
   
    integer :: iter, lol,i
   
    call iniwavefunctions()
    call NaiveFill(occupations)
 
!-------------------------------------------------------------------------------   
! Use this to  quick-and-dirty read wavefunctions from a MOCCa file until I 
! get around to a proper IO module.
!-------------------------------------------------------------------------------
!    open (12,form='unformatted',file='wf/MOCCa.C12.SLy4.dx=0.53.wf')
!    read(12)
!    read(12)
!    read(12)
!    read(12)
!    read(12)
!    nwt = nwp + nwn
!    
!    
!    do i=1,nwt
!        read(12) HFPsi(:,:,i)
!        read(12), occupations(i), spenergies(i), lol, lol, lol, lol, lol, lol, lol,lol,lol,lol
!    enddo
    
    call inilag()
    
    call printSpwfs
    call deriveall()
    
    call calcedfcoefs()
    call printedfcoefs()     
    call densit(0)
    call calcFields()
    call CalcEnergy()
    call PrintEnergy 
    
    do iter=1,maxiter
        print *,  '*************************************'
        print *,  ' Iteration ', iter
        print *,  ' Energy =  ', totalE
        print *,  ' Spwfs  =  ', spwfenergy
        print *,  ' GradNorm =  ', gradientnorm
        print *,  '*************************************'
        
        call Evolve(iter)
        call deriveall()
        call NaiveFill(occupations)

        call densit(iter)
       
        call calcFields()
        call CalcEnergy()
        
        if(mod(iter,PrintIter).eq.0) then
            call PrintSpwfs
            call PrintEnergy            
        endif
    enddo

!    open (12,file='SST.dat')
!    do i=1,nx
!      j = i
!      k = i
!      m = i + (j-1)*nx + (k-1)*nx*nx
!!      print *, i, m
!      r = ( dx/2 + (i-1)*dx) 
!      write(12, '(5f10.5)') r,  D_NmNkNq_NmNkNq(i,1), D_NmNkNq_NmNkNq(i,2) &
!      &                    , sum(D_NmNkNq_NmNkNq(i,:))
!      
!    enddo

end subroutine ReachForWaterAndFood

subroutine ReadInput
    !--------------------------------------------
    ! Subroutine to read all the data from STDIN.
    !
    !--------------------------------------------

    use GenInfo,       only : ReadGenInfo
    use Evolution,     only : ReadEvolution
    use wavefunctions, only : ReadWFdata
    use densities,     only : ReadDensit
    
    call ReadGenInfo
    call ReadEvolution
    call ReadDensit
    call ReadWFdata

end subroutine ReadInput


!    open (12,form='unformatted',file='MOCCa.test')
!    read(12)
!    read(12)
!    read(12)
!    read(12)
!    read(12)
!    allocate(Occupations(nwt))
!    allocate(spenergies(nwt)); allocate(dispersions(nwt))
!    allocate(HFPsi(nx,ny,nz,4,nwt))
!    
!    do i=1,nwt
!        read(12) HFPsi(:,:,:,:,i)
!        read(12), occupations(i), spenergies(i), lol, lol, lol, lol, lol, iso, k,par,k,k
!    enddo
!    
!    nwn=10
!    nwp=10
!    
!    HFBlocks(1) = 7
!    HFBlocks(3) = 3
!    HFBlocks(5) = 7
!    HFBlocks(7) = 3
!    
!    allocate(sx(4,nwt), sy(4,nwt), sz(4,nwt))

