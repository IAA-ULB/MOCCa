!===============================================================================
!     __  __  ___   ____ ____
!    |  \/  |/ _ \ / ___/ ___|__ _
!    | |\/| | | | | |  | |   / _` |
!    | |  | | |_| | |__| |__| (_| |
!    |_|  |_|\___/ \____\____\__,_|
!
!    Copyright (C) 2026 W. Ryssens and M. Bender
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU Affero General Public License as published
!    by the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU Affero General Public License for more details.
!
!    You should have received a copy of the GNU Affero General Public License
!    along with this program.  If not, see <https://www.gnu.org/licenses/>.
!
!===============================================================================
module preconditioning
!===============================================================================
!
! Module that governs all of the possible preconditioning that can be applied.
!
!===============================================================================
!  W.R.-style preconditioning for the potentials
!===============================================================================
!
!     Invert the following operator on the mesh
!         Pf(r) = b*f(r) + a*Delta[f(r)]
!
!     Not by constructing its inverse, but by repeatedly applying that
!     operator in a conjugate gradient scheme.
!===============================================================================
    use compilation, only : dp
    use geninfo, only: nx,ny,nz, mv, dv

    implicit none (external)

    public

contains

 function PreconditionPotential(pot,A,B,sx,sy,sz) result(invpot)
    !---------------------------------------------------------------------------
    ! Precondition a potential with the matrix
    !
    !   ( B  + A * \Delta)^{-1}
    !
    ! Calculated through the repeated application of its inverse in a
    ! conjugate gradient algorithm.
    !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !    potential : potential to be preconditioned (all components!)
    !    A, B      : numbers characterizing the preconditioning matrix
    !    sx/y/z    : integers characterizing the reflection symmetries of the
    !                potential
    ! Output:
    !    invpot    : preconditioned potential (all components!)
    !---------------------------------------------------------------------------

    real(KIND=dp), intent(in) :: a, b
    real(KIND=dp), intent(in) :: pot(nx*ny*nz,4)
    integer, intent(in) :: sx,sy,sz

    real(KIND=dp) :: residual(nx*ny*nz), update(nx*ny*nz), aCG
    real(KIND=dp) :: invpot(nx*ny*nz,4), direction(nx*ny*nz), bCG
    real(KIND=dp) :: newresnorm, oldresnorm

    integer:: it, iter

    !---------------------------------------------------------------------------
    invpot  = 0.0
    do it=1,2
        Residual         = pot(:,it)
        Direction        = Residual
        newresnorm       = sum(direction**2)*dv
        if(abs(newresnorm)<1e-30) then
            ! Don't iterate if we are already good enough. This also takes care
            ! of the possible explosion of this subroutine if pot is zero everywhere.
            !
            ! Note: this tolerance is tighter than the one inside the iteration
            ! loop to ensure we get an update of the potentials more often than not.
           invpot(:,it) = pot(:,it)
           cycle
        endif
        !-----------------------------------------------------------------------
        do iter=1,300
          update   = preconoperator(direction,a,b,sx,sy,sz)

          aCG    = NewResNorm/(sum(Direction*update)*dv)

          invpot(:,it)   = invpot(:,it)  + aCG * Direction
          residual       = residual      - aCG * update

          oldresnorm = newresnorm
          newresnorm = sum(residual**2)*dv

          BCG        = NewResNorm/OldResNorm
          Direction  = Residual + bCG * Direction
          if(newresnorm<1d-16) exit
        enddo
        !-----------------------------------------------------------------------
    enddo
    invpot(:,3) = invpot(:,1) + invpot(:,2)
    invpot(:,4) = invpot(:,1) - invpot(:,2)

    !---------------------------------------------------------------------------
!    print *, 'Inv. Pot., iter = ', iter, newresnorm, sum(invpot(:,1))*dv,  &
!    &                     sum(invpot(:,2))*dv

 end function Preconditionpotential

 function KerkerPreconditionPotential(pot,a,k0,sx,sy,sz) result(invpot)
    !---------------------------------------------------------------------------
    ! Precondition a potential with the matrix
    !
    !      - \Delta
    !   --------------------
    !   ( k0**2  - a^{-1} \Delta)
    !
    ! Calculated through
    !   (a) applying function Preconditionpotential to obtain the multiplication
    !       of (k0**2 - \Delta)^-1 to the potential
    !   (b) applying \Delta to the result
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !    potential : potential to be preconditioned (all components!)
    !    k0        : wavevector associated with the Kerker screening
    !    a         : other parameter, equivalent to 1/mixingstepsize if k0 = 0
    !    sx/y/z    : integers characterizing the reflection symmetries of the
    !                potential
    ! Output:
    !    invpot    : preconditioned potential (all components!)
    !---------------------------------------------------------------------------
    use Derivatives, only : Derive_lap

    real(KIND=dp), intent(in) :: k0, a
    real(KIND=dp), intent(in) :: pot(nx*ny*nz,4)
    integer, intent(in)       :: sx,sy,sz
    integer                   :: it
    real(KIND=dp)             :: invpot(nx*ny*nz,4), temp(nx*ny*nz,4)

    ! Calculate - (k0**2 -a \Delta)^-1 [!note minus sign!]
    temp = - PreconditionPotential(pot,-1/a,k0**2,sx,sy,sz)
    ! Calculate \Delta temp
    do it=1,2
      call Derive_lap(temp(:,it),sx,sy,sz,invpot(:,it))
    enddo
    invpot(:,3) = invpot(:,1) + invpot(:,2)
    invpot(:,4) = invpot(:,1) - invpot(:,2)
 end function KerkerPreconditionpotential

 function preconoperator(f,a,b,sx,sy,sz) result(Pf)
    !---------------------------------------------------------------------------
    ! Implement the preconditioning operator
    !
    ! Pf(r) = B(r)*f(r) + a(r) * Delta(f(r))
    !---------------------------------------------------------------------------
    use Derivatives, only : Derive_lap

    real(KIND=dp), intent(in) :: f(nx*ny*nz)
    real(KIND=dp)             :: Pf(nx*ny*nz), df(nx*ny*nz)
    real(KIND=dp),intent(in)  :: a, b
    integer, intent(in)       :: sx,sy,sz

    call Derive_lap(f,sx,sy,sz,df)
    Pf = b*f + a*df

  end function preconoperator

end module preconditioning
