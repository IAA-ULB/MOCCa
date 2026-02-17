subroutine calc_$NAME(rho_1b, psi, dpsi, ddpsi, hfblocks, ngrid, npsi, $NAME)
    !--------------------------------------------------------------------------------------------
    ! A stand-alone function capable of calculating a local density $NAME from the one-body
    ! density matrix in the canonical basis and the single-particle wavefunctions and their 
    ! derivatives. 
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !    rho_1b   : one-body density matrix in the canonical basis, dimension (npsi)
    !    psi      : single-particle wavefunctions, dimension ( ngrid, 4, npsi )
    !    dpsi     : derivatives of single-particle wavefunctions, dimension ( ngrid, 3, 4, npsi )
    !    ddpsi    : derivatives of single-particle wavefunctions, dimension ( ngrid, 6, 4, npsi )
    !               Note: storage of the second derivatives is 'reduced', meaning that 
    !               Component   1   2   3    4  5  6
    !               Derivative  xx  xy  xz   yy yz zz
    !               The missing derivatives can be constructed from the ones that are stored.
    !               For example, d2 \Psi /dy dx = d2 \Psi /dx dy = component 2    
    !    hfblocks : block division of the single-particle wavefunctions
    !    ngrid    : total number of grid points
    !    npsi     : total number of single-particle wavefunctions
    !
    ! Output:
    !    $NAME    : the calculated local density $NAME at every grid point
    !               dimension ( ngrid$DIM, 2 )
    !               with the final index indicating the isospin (1 = neutrons, 2 = protons)
    !--------------------------------------------------------------------------------------------

    implicit none

    real*8, intent(in)  :: psi(ngrid,4,npsi), dpsi(ngrid,3,4,npsi), ddpsi(ngrid,6,4,npsi)
    real*8, intent(in)  :: rho_1b(npsi)
    integer, intent(in) :: hfblocks(8), ngrid, npsi
    real*8, intent(out) :: $NAME(ngrid$DIM,2)

    real*8 :: temp

    real*8 :: weight
    integer :: B, si, Nblock, it, wave, i

    ! Zero the first part
    $NAME = 0.0d0 

    si = 0
    do B=1,8
        Nblock = hfblocks(B) ; if(Nblock.eq.0) cycle
        it = 1 ; if(B.gt.4) it = 2
        do wave=si+1, si+Nblock 
            weight = rho_1b(wave)
            do i=1,ngrid 
$Expression
            enddo
        enddo
        si = si + Nblock
    enddo
end subroutine calc_$NAME
