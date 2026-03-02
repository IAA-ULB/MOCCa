subroutine calc_$NAME(coupl_constant, sigma, sigma_b, DENLIST, ngrid, $NAME)
    !--------------------------------------------------------------------------------------------
    ! Calculate the mean-field potential $NAME from the densities
    !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !    coupl_constant: an real*8 array of all the coupling constants
    !    ngrid         : total number of grid points
    !
    ! Output:
    !    $NAME    : the calculated mean-field potential $NAME at every grid point
    !               dimension ( ngrid, 2 )
    !               with the final index indicating the isospin (1 = neutrons, 2 = protons)
    !--------------------------------------------------------------------------------------------

    implicit none

DECLARATION
    real*8, intent(in)  :: coupl_constant(:)
    integer, intent(in) :: ngrid
    real*8, intent(out) :: $NAME(ngrid,2)

    real*8 :: temp

$EXPRESSION

end subroutine calc_$NAME

pure function pow(f,alpha) result(pf)
    !---------------------------------------------------------------------------
    ! Safely take powers of a REAL density f, avoiding the raising of negative
    ! numbers to powers that are 0 or negative. This is achieved by adding a
    ! small (positive value) to the density.
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in) :: f(mv), alpha
    real(KIND=dp)             :: pf(mv)

    if(alpha .lt. 0) then
      pf = (f + eps)**(alpha)
    else
      pf = (f)**(alpha)
    endif
end function pow
