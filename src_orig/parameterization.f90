module parameterization
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
 
 use compilation
 
 implicit none
 
contains

  !=============================================================================
  ! Various functions that might be useful to define coupling constants in 
  ! the .func files.
  !=============================================================================
  
  real(KIND=dp) function Cc(t,x,p,s,it) result(c)
    !---------------------------------------------------------------------------
    ! Michaels convention for coupling coefficients of central Skyrme terms.
    ! 
    !
    ! C^+/-_cST       +              -
    !              t      tx      t    tx
    !       
    !    c00     +3/8     0     +5/8  +1/2 
    !    c01     -1/8   -1/4    +1/8  +1/4
    !    c10     -1/8   +1/4    +1/8  +1/4
    !    c11     -1/8     0     +1/8   0
    !---------------------------------------------------------------------------
    
    integer, intent(in)      :: p,s,it
    real(KIND=dp),intent(in) :: t,x
    
    c = 0
    
    select case(p)
    !---------------------------------------------------------------------------
    case(+1) ! C^+_c
      select case(s)          ! C^+_cS  
        case(0)               ! C^+_c0 
          select case(it)     ! C^+_c0T
          case(0)
            c = +3.0/8.0 * t                     ! C^+_c00
          case(1)
            c = -1.0/8.0 * t  - 1.0/4.0 * t * x  ! C^+_c01 
          end select
        case(1)               ! C^+_c1 
          select case(it)     ! C^+_c1T
          case(0)
            c = -1.0/8.0 * t  - 1.0/4.0 * t * x  ! C^+_c10
          case(+1)
            c = -1.0/8.0 * t                     ! C^+_c11
          end select
      end select
    !---------------------------------------------------------------------------
    case(-1) ! C^-_c
      select case(s)          ! C^-_cS  
        case(0)               ! C^-_c0 
          select case(it)     ! C^-_c0T
          case(0)
            c = +5.0/8.0 * t  + 1.0/2.0 * t * x  ! C^-_c00
          case(1)
            c = +1.0/8.0 * t  + 1.0/4.0 * t * x  ! C^-_c01 
          end select
        case(1)               ! C^-_c1 
          select case(it)     ! C^-_c1T
          case(0)
            c = +1.0/8.0 * t  + 1.0/4.0 * t * x  ! C^-_c10
          case(+1)
            c = +1.0/8.0 * t                     ! C^-_c11
          end select
      end select
    end select
    !---------------------------------------------------------------------------
  end function Cc


end module parameterization
