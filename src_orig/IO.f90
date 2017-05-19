module IO
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
 ! Module governing the in- and output of Tantalus. 
 !
 !==============================================================================

use geninfo
use wavefunctions

implicit none

contains
    
  subroutine PrintInput
  !-----------------------------------------------------------------------------
  ! This subroutine prints all relevant information of the input, both from the
  ! user and from the wavefunction file.
  !-----------------------------------------------------------------------------
   
    use wavefunctions
   
    1 format ( 20('-'), 'General Information ', 20('-'))
    2 format ( 'Mesh parameters' )
    3 format ( '   nx = ', i5 , ' ny = ' , i5 , ' nz = ' , i5)
    4 format ( '   dx = ', f20.10,' (fm  ) ')
    5 format ( '   dv = ', f5.2,' (fm^3) ')
    6 format ( 'Nucleus')
    7 format ( '    N = ', f10.5  ,'  Z = ', f10.5)
    8 format ( 'Wavefunctions')
    9 format ( '  nwt = ', i5, / &
    &          '  nwn = ', i5, / &
    &          '  nwp = ', i5 )

    print 1
    print 2
    print 3 , nx, ny, nz
    print 4 , dx
    print 5 , dv
    print 6
    print 7 , neutrons, protons
    print 8
    print 9 , nwt,nwn,nwp
    
  end subroutine PrintInput



end module IO
