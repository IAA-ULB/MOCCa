!===============================================================================
!     __  __  ___   ____ ____
!    |  \/  |/ _ \ / ___/ ___|__ _
!    | |\/| | | | | |  | |   / _` |
!    | |  | | |_| | |__| |__| (_| |
!    |_|  |_|\___/ \____\____\__,_|
!
! Written mainly by W. Ryssens & M. Bender
!
! Opensource software distributed under the GNU AGPLv3 licence, see the
!  LICENCE file in the root of this project.
!===============================================================================
module compilation
  !-----------------------------------------------------------------------------
  ! Module to propagate the definition of a double and single precision real
  !-----------------------------------------------------------------------------

  implicit none (external)

  public

  integer, parameter :: dp = selected_real_kind(15,307)
  integer, parameter :: sp = selected_real_kind(6,37)

end module compilation

