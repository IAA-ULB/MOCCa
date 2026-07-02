subroutine calc_G_I_NS(Der_D_I_I, C_I_NS, coupl_constant, ngrid, G_I_NS)
  !--------------------------------------------------------------------------------------------
  ! Calculate the mean-field potential G_I_NS from the densities
  !
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Input:
  !    Der_D_I_I     : real*8 array, dimension (ngrid, 3, 2)
  !                    the gradient of the mean-field density D_I_I at every mesh point
  !                    the final index indicating ISOSPIN, i.e. 1-> t=0 & 2-> t=1.
  !
  !    C_I_NS        : real*8 array, dimension (ngrid,3,3,2)
  !                    the spin-current mean-field density C_I_NS at every mesh point
  !                     1st index -> space
  !                     2nd index -> Nabla
  !                     3rd index -> Sigma
  !                     4th index -> isospin, t=0/1
  !
  !    coupl_constant: real*8 array of length 14
  !                    listing all the coupling constants in the ordering of the BXL.func file
  !
  !    ngrid         : integer
  !                    total number of grid points
  !
  ! Output:
  !    G_I_NS        : real*8 array dimension (ngrid,3,3,2)
  !                    the calculated mean-field potential G_I_NS at every grid point
  !                     1st index -> space
  !                     2nd index -> Nabla
  !                     3rd index -> Sigma
  !                     4th index -> isospin, t=0/1
  !--------------------------------------------------------------------------------------------

  integer, intent(in) :: ngrid
  real*8, intent(in)  :: Der_D_I_I(ngrid,3,2), C_I_NS(ngrid,3,3,2), coupl_constant(14)
  real*8, intent(out) :: G_I_NS(ngrid,3,3,2)

  !---------------------------------------------------------------------------
  ! Calculation of G_I_NS
  G_I_NS = 0.0
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Isospin = 0
  G_I_NS(:,1,1,1) = G_I_NS(:,1,1,1)  &
       & + coupl_constant(11)  * C_I_NS(:,1,1,1)

  G_I_NS(:,1,2,1) = G_I_NS(:,1,2,1)  &
       & + coupl_constant(11)  * C_I_NS(:,1,2,1)

  G_I_NS(:,1,3,1) = G_I_NS(:,1,3,1)  &
       & + coupl_constant(11)  * C_I_NS(:,1,3,1)

  G_I_NS(:,2,1,1) = G_I_NS(:,2,1,1)  &
       & + coupl_constant(11)  * C_I_NS(:,2,1,1)

  G_I_NS(:,2,2,1) = G_I_NS(:,2,2,1)  &
       & + coupl_constant(11)  * C_I_NS(:,2,2,1)

  G_I_NS(:,2,3,1) = G_I_NS(:,2,3,1)  &
       & + coupl_constant(11)  * C_I_NS(:,2,3,1)

  G_I_NS(:,3,1,1) = G_I_NS(:,3,1,1)  &
       & + coupl_constant(11)  * C_I_NS(:,3,1,1)

  G_I_NS(:,3,2,1) = G_I_NS(:,3,2,1)  &
       & + coupl_constant(11)  * C_I_NS(:,3,2,1)

  G_I_NS(:,3,3,1) = G_I_NS(:,3,3,1)  &
       & + coupl_constant(11)  * C_I_NS(:,3,3,1)

  G_I_NS(:,1,1,1) = G_I_NS(:,1,1,1)  &
       & + coupl_constant(11)  * C_I_NS(:,1,1,1)

  G_I_NS(:,1,2,1) = G_I_NS(:,1,2,1)  &
       & + coupl_constant(11)  * C_I_NS(:,1,2,1)

  G_I_NS(:,1,3,1) = G_I_NS(:,1,3,1)  &
       & + coupl_constant(11)  * C_I_NS(:,1,3,1)

  G_I_NS(:,2,1,1) = G_I_NS(:,2,1,1)  &
       & + coupl_constant(11)  * C_I_NS(:,2,1,1)

  G_I_NS(:,2,2,1) = G_I_NS(:,2,2,1)  &
       & + coupl_constant(11)  * C_I_NS(:,2,2,1)

  G_I_NS(:,2,3,1) = G_I_NS(:,2,3,1)  &
       & + coupl_constant(11)  * C_I_NS(:,2,3,1)

  G_I_NS(:,3,1,1) = G_I_NS(:,3,1,1)  &
       & + coupl_constant(11)  * C_I_NS(:,3,1,1)

  G_I_NS(:,3,2,1) = G_I_NS(:,3,2,1)  &
       & + coupl_constant(11)  * C_I_NS(:,3,2,1)

  G_I_NS(:,3,3,1) = G_I_NS(:,3,3,1)  &
       & + coupl_constant(11)  * C_I_NS(:,3,3,1)

  G_I_NS(:,2,3,1) = G_I_NS(:,2,3,1)  &
       & - coupl_constant(13)  * Der_D_I_I(:,1,1)

  G_I_NS(:,3,2,1) = G_I_NS(:,3,2,1)  &
       & + coupl_constant(13)  * Der_D_I_I(:,1,1)

  G_I_NS(:,3,1,1) = G_I_NS(:,3,1,1)  &
       & - coupl_constant(13)  * Der_D_I_I(:,2,1)

  G_I_NS(:,1,3,1) = G_I_NS(:,1,3,1)  &
       & + coupl_constant(13)  * Der_D_I_I(:,2,1)

  G_I_NS(:,1,2,1) = G_I_NS(:,1,2,1)  &
       & - coupl_constant(13)  * Der_D_I_I(:,3,1)

  G_I_NS(:,2,1,1) = G_I_NS(:,2,1,1)  &
       & + coupl_constant(13)  * Der_D_I_I(:,3,1)

  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Isospin = 1
  G_I_NS(:,1,1,2) = G_I_NS(:,1,1,2)  &
       & + coupl_constant(12)  * C_I_NS(:,1,1,2)

  G_I_NS(:,1,2,2) = G_I_NS(:,1,2,2)  &
       & + coupl_constant(12)  * C_I_NS(:,1,2,2)

  G_I_NS(:,1,3,2) = G_I_NS(:,1,3,2)  &
       & + coupl_constant(12)  * C_I_NS(:,1,3,2)

  G_I_NS(:,2,1,2) = G_I_NS(:,2,1,2)  &
       & + coupl_constant(12)  * C_I_NS(:,2,1,2)

  G_I_NS(:,2,2,2) = G_I_NS(:,2,2,2)  &
       & + coupl_constant(12)  * C_I_NS(:,2,2,2)

  G_I_NS(:,2,3,2) = G_I_NS(:,2,3,2)  &
       & + coupl_constant(12)  * C_I_NS(:,2,3,2)

  G_I_NS(:,3,1,2) = G_I_NS(:,3,1,2)  &
       & + coupl_constant(12)  * C_I_NS(:,3,1,2)

  G_I_NS(:,3,2,2) = G_I_NS(:,3,2,2)  &
       & + coupl_constant(12)  * C_I_NS(:,3,2,2)

  G_I_NS(:,3,3,2) = G_I_NS(:,3,3,2)  &
       & + coupl_constant(12)  * C_I_NS(:,3,3,2)

  G_I_NS(:,1,1,2) = G_I_NS(:,1,1,2)  &
       & + coupl_constant(12)  * C_I_NS(:,1,1,2)

  G_I_NS(:,1,2,2) = G_I_NS(:,1,2,2)  &
       & + coupl_constant(12)  * C_I_NS(:,1,2,2)

  G_I_NS(:,1,3,2) = G_I_NS(:,1,3,2)  &
       & + coupl_constant(12)  * C_I_NS(:,1,3,2)

  G_I_NS(:,2,1,2) = G_I_NS(:,2,1,2)  &
       & + coupl_constant(12)  * C_I_NS(:,2,1,2)

  G_I_NS(:,2,2,2) = G_I_NS(:,2,2,2)  &
       & + coupl_constant(12)  * C_I_NS(:,2,2,2)

  G_I_NS(:,2,3,2) = G_I_NS(:,2,3,2)  &
       & + coupl_constant(12)  * C_I_NS(:,2,3,2)

  G_I_NS(:,3,1,2) = G_I_NS(:,3,1,2)  &
       & + coupl_constant(12)  * C_I_NS(:,3,1,2)

  G_I_NS(:,3,2,2) = G_I_NS(:,3,2,2)  &
       & + coupl_constant(12)  * C_I_NS(:,3,2,2)

  G_I_NS(:,3,3,2) = G_I_NS(:,3,3,2)  &
       & + coupl_constant(12)  * C_I_NS(:,3,3,2)

  G_I_NS(:,2,3,2) = G_I_NS(:,2,3,2)  &
       & - coupl_constant(14)  * Der_D_I_I(:,1,2)

  G_I_NS(:,3,2,2) = G_I_NS(:,3,2,2)  &
       & + coupl_constant(14)  * Der_D_I_I(:,1,2)

  G_I_NS(:,3,1,2) = G_I_NS(:,3,1,2)  &
       & - coupl_constant(14)  * Der_D_I_I(:,2,2)

  G_I_NS(:,1,3,2) = G_I_NS(:,1,3,2)  &
       & + coupl_constant(14)  * Der_D_I_I(:,2,2)

  G_I_NS(:,1,2,2) = G_I_NS(:,1,2,2)  &
       & - coupl_constant(14)  * Der_D_I_I(:,3,2)

  G_I_NS(:,2,1,2) = G_I_NS(:,2,1,2)  &
       & + coupl_constant(14)  * Der_D_I_I(:,3,2)

end subroutine calc_G_I_NS
