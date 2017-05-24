module constants
!=======================================================================
!  #######   ##   #    # #####   ##   #      #    #  ####
!     #     #  #  ##   #   #    #  #  #      #    # #
!     #    #    # # #  #   #   #    # #      #    #  ####
!     #    ###### #  # #   #   ###### #      #    #      #
!     #    #    # #   ##   #   #    # #      #    # #    #
!     #    #    # #    #   #   #    # ######  ####   ####
!
!  Copyright W. Ryssens & M. Bender
!
!=======================================================================
! Module governing the constant that go into a mean-field program. 
! Static (values of SLy4) for the moment, but will get changeable at 
! runtime rather quickly. 
!=======================================================================
    use compilation

    implicit none

    real(KIND=dp) :: t0=-2488.913 
    real(KIND=dp) :: x0=0.834
    real(KIND=dp) :: t1=486.818
    real(KIND=dp) :: x1=-0.344
    real(KIND=dp) :: t2=-546.395
    real(KIND=dp) :: x2=-1.0
    real(KIND=dp) :: t3a=13777.0
    real(KIND=dp) :: x3a=1.354
    real(KIND=dp) :: yt3a=0.166666666666666666667 
    real(KIND=dp) :: t3b=0.0 
    real(KIND=dp) :: x3b=0.0
    real(KIND=dp) :: yt3b=0
    real(KIND=dp) :: te=0.0
    real(KIND=dp) :: to=0.0
    real(KIND=dp) :: wso=123.0
    real(KIND=dp) :: wsoq=123.0
    
    real(KIND=dp) :: t1n2=24.3409
    real(KIND=dp) :: t2n2=-27.31975
    real(KIND=dp) :: x1n2=-0.344
    real(KIND=dp) :: x2n2=-1.0   
    
    !Functional parameters in the BFH representation
    real(KIND=dp) :: B1,B2,B3,B4,B5,B6,B7,B8,Byt3,B9,B9q
    real(KIND=dp) :: B10,B11,B12,B13,B14
    real(KIND=dp) :: B15,B16, B17, B18,B19,B20,B21
    !------------------------------------------------------------------------------
    ! Physical constants.
    real(KIND=dp):: e2             =  1.43996446_dp
    real(KIND=dp):: hbar           =  6.58211928_dp
    real(KIND=dp):: clum           = 29.9792458_dp
    real(KIND=dp):: nucleonmass(2) = (/939.565379_dp , 938.272046_dp /)
    real(KIND=dp):: hbm(2)         = 20.73551910_dp
    
    
    !Functional parameters in the C-representation
    !Note: first component is the isoscalar constant, second the isovector part
    real(KIND=dp), public :: Crho(2),Crhosat(2), Ctau(2), Cdrho(2), CnablaJ(2)
    real(KIND=dp), public :: Cs(2), Cssat(2), Ct(2), Cf(2), Cds(2), Cnablas(2)
    !Recouplings of the tensor terms
    real(KIND=dp), public :: CJ0(2), CJ1(2), CJ2(2)
    ! J^2 terms and whether or not to take average nucleon masses
    logical               :: J2terms=.false., averagemass=.true.
    
    !---------------------------------------------------------------------------
    ! N2LO terms
    ! Ordering 
    ! C(1) => Delta Rho Delta Rho_0
    ! C(2) => Delta Rho Delta Rho_1
    ! C(3) => M(rho)_0
    ! C(4) => M(rho)_1
    ! C(5) => Delta s Delta s_0
    ! C(6) => Delta s Delta s_1
    ! C(7) => M(s)_0
    ! C(8) => M(s)_1
    real(KIND=dp) :: CN2LO(8), BN2LO(8)
    !---------------------------------------------------------------------------
    ! Coupling constants that are in principle calculated from the BN2LO, but 
    ! allow for independently determining the contribution of all the terms.
    !               (Delta rho)^2 (Delta s)^2 rho Q      tau^2     (Re tmn)^2
    real(KIND=dp) :: N2D2rho(2),    N2D2s(2), N2rhoQ(2), N2tau(2), N2rtaumn(2)
    !                (Im tmn)^2  tmn DmDnrho   (Dm jm)^2   (j Pi)   (Dm Jmn)^2
    real(KIND=dp) :: N2itaumn(2),  N2tddr(2),   N2Dvecj(2), N2jpi(2), N2DJ(2)
    !                J V       s S       (T_k)^2  (Re Tmn)^2   (Im Tmn)^2 
    real(KIND=dp) :: N2JV(2), N2sS(2), N2vecT(2), N2ReTmn(2), N2ImTmn(2)
    !                Tmn DmDns
    real(KIND=dp) :: N2TmnD2s(2) 
    !---------------------------------------------------------------------------
    ! Coupling constants that are in principle calculated from the BN3LO, but 
    ! allow for independently determining the contribution of all the terms.
    !               rho Delta Delta Delta rho
    real(KIND=dp) :: N3D3rho(2), N3tauDtau(2), N3taumnDtaumn(2),N3DrhoDtau(2)
    real(KIND=dp) :: N3tauDmntau(2), N3tauDDrho(2)
    real(KIND=dp) :: C3D3rho(2), C3tauDtau(2), C3taumnDtaumn(2),C3DrhoDtau(2)
    real(KIND=dp) :: C3tauDmntau(2), C3tauDDrho(2)

contains 

    subroutine calcedfcoefs()
        real(KIND=dp), parameter :: rhosat=0.16_dp
        !---------------------------------------------------------------------------
        !Calculating the B-s
        !---------------------------------------------------------------------------
        B1   = t0*(1.0_dp   + 1/2.0_dp*x0)/2.0_dp
        B2   =-t0*(1/2.0_dp + x0         )/2.0_dp

        B3   = (t1*( 1      + 1/2.0_dp*x1) + t2*(1      + 1/2.0_dp*x2))/4.0_dp
        B4   =-(t1*( 1/2.0_dp + x1       ) - t2*(1/2.0_dp +        x2))/4.0_dp

        B5   =-(3.0_dp*t1*(1  + x1/2.0_dp)-t2*(1  + x2/2.0_dp))/16.0_dp
        B6   = (3.0_dp*t1*(x1 +  1/2.0_dp)+t2*(x2 +  1/2.0_dp))/16.0_dp

        B7  = t3a*( 1   + x3a/2.0_dp)/12.0_dp
        B8  =-t3a*( x3a +   1/2.0_dp)/12.0_dp
        !---------------------------------------------------------------------------
        ! B9 and B9q are the constants related to the spin-orbit interaction.
        ! Starting from a Skyrme Force, they should in principle be equal for a
        ! true force, but apparently there are some physics cases where a difference
        ! between the isoscalar and isovector coupling is wanted.
        !---------------------------------------------------------------------------
        B9   = -wso /2.0_dp
        B9q  = -wsoq/2.0_dp

        Byt3= yt3a
        
        !---------------------------------------------------------------------------
        ! N2LO terms
        if(t1n2 .ne. 0.0d0 .or. t2n2.ne.0.0d0) then
            CN2LO(1) = 9/128.0d0 * t1n2   - 0.0d0    * t1n2 *x1n2                &
            &                             + t2n2 * (-5/128.0  - 4 /128.0d0 * x2n2)
            CN2LO(2) =-3/128.0d0 * t1n2   - 3/64.0d0 * t1n2 *x1n2                &
            &                             + t2n2 * (-1/128.0  - 1 / 64.0d0 * x2n2)
            CN2LO(3) = 3/ 32.0d0 * t1n2   - 0.0d0    * t1n2 *x1n2                &
            &                             + t2n2 * ( 5/ 32.0  + 1 /  8.0d0 * x2n2)
            CN2LO(4) =-1/ 32.0d0 * t1n2   - 1/16.0d0 * t1n2 *x1n2                &
            &                             + t2n2 * ( 1/ 32.0  + 1 / 16.0d0 * x2n2)
            CN2LO(5) =-3/128.0d0 * t1n2   + 3/64.0d0 * t1n2 *x1n2                &
            &                             + t2n2 * (-1/128.0  - 1 / 64.0d0 * x2n2)
            CN2LO(6) =-3/128.0d0 * t1n2   - 0.0d0    * t1n2 *x1n2                &
            &                             + t2n2 * (-1/128.0  + 0.0d0      * x2n2)
            CN2LO(7) =-1/ 32.0d0 * t1n2   + 1/16.0d0 * t1n2 *x1n2                &
            &                             + t2n2 * ( 1/ 32.0  + 1 / 16.0d0 * x2n2)
            CN2LO(8) =-1/ 32.0d0 * t1n2   - 0.0d0    * t1n2 *x1n2                &
            &                             + t2n2 * ( 1/ 32.0  + 0.0d0      * x2n2)
        
            BN2LO(1) =   CN2LO(1) -  CN2LO(2)
            BN2LO(2) =             2*CN2LO(2)
            BN2LO(3) =   CN2LO(3) -  CN2LO(4)
            BN2LO(4) =             2*CN2LO(4)
            BN2LO(5) =   CN2LO(5) -  CN2LO(6)
            BN2LO(6) =             2*CN2LO(6)
            BN2LO(7) =   CN2LO(7) -  CN2LO(8)
            BN2LO(8) =             2*CN2LO(8)  
            
            N2D2rho(1) = BN2LO(1) ; N2D2rho(2) = BN2LO(2)
            
            N2D2s(1)   = BN2LO(5) ; N2D2s(2)   = BN2LO(6)
            
            N2rhoQ(1)  = BN2LO(3) ; N2rhoQ(2)  = BN2LO(4)
            N2tau(1)   = BN2LO(3) ; N2tau(2)   = BN2LO(4)
            N2rtaumn(1)= BN2LO(3) ; N2rtaumn(2)= BN2LO(4)
            N2itaumn(1)= BN2LO(3) ; N2itaumn(2)= BN2LO(4)
            N2tddr(1)  = BN2LO(3) ; N2tddr(2)  = BN2LO(4)
            N2Dvecj(1) = BN2LO(3) ; N2Dvecj(2) = BN2LO(4)
            N2jpi(1)   = BN2LO(3) ; N2jpi(2)   = BN2LO(4)
            
            N2DJ(1)    = BN2LO(7) ; N2DJ(2)    = BN2LO(8)
            N2JV(1)    = BN2LO(7) ; N2JV(2)    = BN2LO(8)
            N2sS(1)    = BN2LO(7) ; N2sS(2)    = BN2LO(8)
            N2vecT(1)  = BN2LO(7) ; N2vecT(2)  = BN2LO(8)
            N2ReTmn(1) = BN2LO(7) ; N2reTmn(2) = BN2LO(8)
            N2ImTmn(1) = BN2LO(7) ; N2ImTmn(2) = BN2LO(8)
            N2TmnD2s(1)= BN2LO(7) ; N2TmnD2s(2)= BN2LO(8)
        endif
        
        
        return
    end subroutine calcedfcoefs
    
    

end module 
