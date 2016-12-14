module functional
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
 !
 ! Temporary module holding the definition of the functional. 
 !
 !
 !=======================================================================
 
 use compilation
 use geninfo
 
 implicit none
 
 !------------------------------------------------------------------------------
 !Name of the force
 character(len=200), public   :: force
 
 !------------------------------------------------------------------------------
 !Skyrme Force parameters
 real(KIND=dp) :: t0,x0,t1,x1,t2,x2,t3a,x3a,yt3a,t3b,x3b,yt3b,te,to
 real(KIND=dp) :: wso,wsoq
 !Functional parameters in the BFH representation
 real(KIND=dp) :: B1,B2,B3,B4,B5,B6,B7a,B7b,B8a,B8b,Byt3a,Byt3b,B9,B9q
 real(KIND=dp) :: B10,B11,B12a,B12b,B13a,B13b,B14
 real(KIND=dp) :: B15,B16, B17, B18,B19,B20,B21
 
 ! J^2 terms and whether or not to take average nucleon masses
 logical               :: J2terms=.false., averagemass=.true.
 
 real(KIND=dp):: e2=1.43996446_dp, hbar = 6.58211928_dp, clum  =  29.9792458_dp
 real(KIND=dp):: nucleonmass(2) = (/939.565379_dp , 938.272046_dp /)
 real(KIND=dp):: hbm(2)      = 20.73551910_dp
 
 integer :: COM1Body=2, COM2Body=0
 
 contains

 subroutine ReadForceInfo()
    !---------------------------------------------------------------------------
    ! Subroutine governing the user input of parameters regarding this module.
    !---------------------------------------------------------------------------
    NameList /Functional/ force

    !Read the correct name from input
    read(unit=*, NML=Functional)

    ! Assigning the correct variables to the Force coefficients.
    call ReadForce
    !Calculating all EDF coefficients
    call CalcEdfCoef
  end subroutine ReadForceInfo
  
  subroutine ReadForce
  !-----------------------------------------------------------------------------
  ! Subroutine that reads a force from the forces.param file that should be
  ! provided.
  !-----------------------------------------------------------------------------

    integer            :: inunit, io
    character(len=200) :: Name, UpName, UpAfor
    logical            :: exists

    NameList /skf/ Name,t0,x0,t1,x1,t2,x2,t3a,x3a,yt3a,t3b,x3b,yt3b,te,to, &
    &                    wso,wsoq,                                         &
    &                    hbm,e2,                                           &
    &                    COM1body, COM2body,                               &
    &                    J2Terms   ,                                       &
    &                    averagemass, hbar, nucleonmass

    inquire(file='forces.param', exist=exists)
    if (.not. exists) stop

    !Setting some default values
    t0=0.0_dp;      x0=0.0_dp
    t1=0.0_dp;      t2=0.0_dp;       x2=0.0_dp
    t3a=0.0_dp;     yt3a=0.0_dp;     t3b=0.0_dp;
    x3b=0.0_dp;     yt3b=0.0_dp
    te=0.0_dp;      to=0.0_dp
    wso=0.0_dp;     wsoq=0.0_dp
    COM1body=0;     COM2body=0
    J2Terms=.false.
    !call ResetConstants

    inunit=1
    open(unit=inunit, file='forces.param',iostat=io)
    if (io.ne.0) stop
    do
        read(NML=skf, unit=inunit)
        !Creating uppercase strings
        call to_upper(force, UpAfor)
        call to_upper(Name , UpName)

        if(trim(upname).eq.trim(upafor)) then
            exit
        elseif(trim(upname).eq.trim('END')) then
            stop
        else
             !Resetting some default values
            t0=0.0_dp;      x0=0.0_dp
            t1=0.0_dp;      t2=0.0_dp;       x2=0.0_dp
            t3a=0.0_dp;     yt3a=0.0_dp;     t3b=0.0_dp;
            x3a=0.0_dp
            x3b=0.0_dp;     yt3b=0.0_dp
            te=0.0_dp;      to=0.0_dp
            wso=0.0_dp;     wsoq=0.0_dp
            b14=0.0_dp;     b15=0.0_dp;     b16=0.0_dp;      b17=0.0_dp
            COM1body=0;     COM2body=0
            J2Terms=.false.
            !call ResetConstants
        endif
    enddo
    close(inunit)
    !Averaging masses if needed:
    if(Averagemass) nucleonmass = (nucleonmass(1) + nucleonmass(2))/2
!     !Multiply hbm by two
     hbm = hbm*2.0

  end subroutine ReadForce
  
  subroutine CalcEDFCoef()
  !-----------------------------------------------------------------------------
  ! This subroutine calculates the Bi, the coefficients in the EDF, as a
  ! function of the Skyrme parameters. See the formulas in the 24Mg paper and
  ! the tensor notes.
  !-----------------------------------------------------------------------------
  ! Note for the calculation of the C's: we use the formulas in terms of B's,
  ! because these work whether or  not we are dealing with a strict force.
  !-----------------------------------------------------------------------------
  ! Note again for the C's: they are BUGGED!
  !-----------------------------------------------------------------------------
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

    B7a  = t3a*( 1   + x3a/2.0_dp)/12.0_dp
    B8a  =-t3a*( x3a +   1/2.0_dp)/12.0_dp
    B7b  = t3b*( 1   + x3b/2.0_dp)/12.0_dp
    B8b  =-t3b*( x3b +   1/2.0_dp)/12.0_dp
    !---------------------------------------------------------------------------
    ! B9 and B9q are the constants related to the spin-orbit interaction.
    ! Starting from a Skyrme Force, they should in principle be equal for a
    ! true force, but apparently there are some physics cases where a difference
    ! between the isoscalar and isovector coupling is wanted.
    !---------------------------------------------------------------------------
    B9   = -wso /2.0_dp
    B9q  = -wsoq/2.0_dp

    Byt3a= yt3a
    Byt3b= yt3b

!    if(.not.TRC) then
!       B10  = t0*x0/4.0_dp
!       B11  =-t0   /4.0_dp
!       B12a = t3a*x3a/24.0_dp
!       B13a =-t3a    /24.0_dp
!       B12b = t3b*x3b/24.0_dp
!       B13b =-t3b    /24.0_dp
!    endif

!    if(J2Terms) then
!      B14  = -(1.0_dp/8.0_dp) * (t1*x1 + t2*x2)
!      B15  =  (1.0_dp/8.0_dp) * (t1 - t2)
!      ! Tensor contribution
!      B14 = B14 + (1.0_dp/4.0_dp) * (te + to)
!      B15 = B15 - (1.0_dp/4.0_dp) * (te - to)
!    endif
!    !
!    B16  =-(3.0_dp/8.0_dp) * (te + to)
!    B17  = (3.0_dp/8.0_dp) * (te - to)
!    !
!    if(J2Terms .and. (.not. TRC)) then
!      B18 =-(1.0_dp/32.0_dp)* (3 * t1 * x1 - t2 *x2)
!      B19 = (1.0_dp/32.0_dp)* (3 * t1      + t2    )
!      !Tensor contribution
!      B18 = B18 + (1.0_dp/16.0_dp) * (3*te - to)
!      B19 = B19 - (1.0_dp/16.0_dp) * (3*te + to)
!    endif
!    if(.not.TRC) then
!      B20 = (3.0_dp/16.0_dp) * (3*te - to)
!      B21 =-(3.0_dp/16.0_dp) * (3*te + to)
!    endif
  end subroutine CalcEDFCoef
end module functional
