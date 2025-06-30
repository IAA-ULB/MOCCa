!>
!!
!! This module is based on the GMRES implementation of Choral of C. Pierre
!! https://plmlab.math.cnrs.fr/cpierre1/choral
!!
!!
!! <B> GMRES LINEAR SOLVER  </B>
!!
!! Source = Youssef SAAD, 'Iterative methods for sparse linear system'
!! https://www-users.cs.umn.edu/~saad/IterMethBook_2ndEd.pdf
!>

module gmres

  use evolution


  implicit none
  private

  public :: do_gmres   !! To be tested

  contains


  !> GMRES (no preconditioning)
  !> 
  !>
  !> for the linear system \f$ Ax = b\f$
  !>
  !> INPUT/OUTPUT :
  !> \li            x = initial guess / solution
  !>
  !> OUTPUT :
  !> \li            res  = final residual
  !> \li            iter = number of performed iterations
  !> \li            iter = -1 = resolution failure     
  !>
  !> INPUT :
  !> \li            b      = RHS
  !> \li            A      = \f$ A~:~~~x \mapsto A x \f$
  !>                         matrix/vector product (procedural)
  !> \li            norm_2 = \f$ f~:~~~x \mapsto real \f$
  !>                         norm (procedural)
  !> \li            norm_2 = \f$ f~:~~~(x1, x2) \mapsto real \f$
  !>                         Scalar product (procedural)
  !> \li            tol    = tolerance
  !> \li            itMax  = maximal iteration-number
  !> \li            rst    = restart number
  !> \li            verb   = verbosity
  !>
  subroutine do_gmres(x, iter, nbPrd, res, &
       & b, A, norm_2, ScalProd, tol, itmax, rst, verb)

    real(KIND=dp), dimension(:), intent(inout) :: x
    integer                    , intent(out)   :: iter, nbPrd
    real(KIND=dp)              , intent(out)   :: res
    procedure(dHtodH)                          :: A
    procedure(dHtoreal)                        :: norm_2
    procedure(dHdHtoreal)                      :: ScalProd
    real(KIND=dp), dimension(:), intent(in)    :: b
    real(KIND=dp)              , intent(in)    :: tol
    integer                    , intent(in)    :: itmax, rst, verb

    real(KIND=dp), dimension(rst+1, size(x,1) ) :: V
    real(KIND=dp), dimension(rst+1, rst       ) :: H

    real(KIND=dp), dimension(size(x,1)) :: r, w
    real(KIND=dp), dimension(rst      ) :: sn, cs, y
    real(KIND=dp), dimension(rst +1   ) :: s

    real(KIND=dp) :: nb2, nr2, temp
    integer       :: nn, ii, kk

    if (verb>1) write(*,*) 'gmres       : gmres'

    iter  = 0
    nbPrd = 0


    nb2 = norm_2(b)
    nn = size(x,1)
    if (nb2/real(nn, dp)<1E-8_dp) nb2=sqrt(real(nn, dp))
    nb2 = 1._dp/nb2

    call A(r, x)
    nbPrd = 1
    r   = b-r
    nr2 = norm_2(r)
    res = nr2 * nb2

    if (verb>2) write(*,*)'  iter', iter,' residual', res
    if (res<tol) return

    V  = 0._dp
    H  = 0._dp
    cs = 0._dp
    sn = 0._dp

    do iter=1, itmax

       if (verb>2) write(*,*)'  iter', iter,' residual', res

       s     = 0._dp
       s(1)  = nr2

       V(1,:)  = r / nr2
      
       do ii = 1, rst  ! orthonormal basis using Gram-Schmidt

          call A(w, V(ii,:))
          nbPrd = nbPrd + 1

          do kk = 1, ii
             H(kk,ii)= ScalProd( w, V(kk,:))
             w = w - H(kk,ii)*V(kk,:)
          end do

          H(ii+1,ii) = norm_2(w)

          V(ii+1, : ) = w / H(ii+1,ii)

          ! apply Givens rotation
          !
          do kk = 1, ii-1 
             temp       =  cs(kk)*H(kk,ii) + sn(kk)*H(kk+1,ii)
             H(kk+1,ii) = -sn(kk)*H(kk,ii) + cs(kk)*H(kk+1,ii)
             H(kk,ii)   = temp
          end do

          call grotmat(cs(ii), sn(ii), H(ii,ii), H(ii+1,ii) )

          temp    =  cs(ii)*s(ii)             ! approximate residual norm
          s(ii+1) = -sn(ii)*s(ii)
          s(ii)   = temp

          H(ii,ii)   = cs(ii)*H(ii,ii) + sn(ii)*H(ii+1,ii)
          H(ii+1,ii) = 0._dp

          res  = abs(s(ii+1)) * nb2

          if (verb>2) write(*,*)'  iter', iter,' residual', res

          if ( res < tol ) then

             call invtrisup(y(1:ii), H(1:ii,1:ii), s(1:ii))
             do kk=1, ii
                x    = x + V(kk, : )*y(kk)
             end do

             return

          end if

       end do

       ! update approximation
       !
       call invtrisup(y(1:rst), H(1:rst,1:rst), s(1:rst))
       do kk=1, rst
          x    = x + V(kk, : )*y(kk)
       end do

       ! compute residual
       !
       call  A(r, x)                  
       nbPrd = nbPrd + 1
       r     = b-r
       nr2   = norm_2(r)
       res   = nr2 * nb2
       if (verb>2) write(*,*)'  iter', iter,' residual', res

       if ( res < tol ) return

    end do

    ! Convergence failure flag
    print *, "gmres: gmres: not converged"
    iter = -1

  end subroutine do_gmres

  !> Matrice de rotation de Givens
  subroutine grotmat(cs, sn, a, b)

    real(KIND=dp), intent(in)  :: a,b
    real(KIND=dp), intent(out) :: cs,sn

    real(KIND=dp)              :: tmp

    if (abs(b)<1E-12_dp) then ! tolerance might be needed to adjusted 
       cs = 1._dp
       sn = 0._dp

    else if (abs(b)>abs(a)) then
       tmp = a/b
       sn  = 1._dp / sqrt( 1._dp + tmp**2)
       cs  = tmp*sn

    else
       tmp = b/a
       cs  = 1._dp / sqrt( 1._dp + tmp**2)
       sn  = tmp*cs

    end if

  end subroutine grotmat

   subroutine invtrisup(res,M,vec)

    real(KIND=dp), dimension(:,:), intent(in)  :: M
    real(KIND=dp), dimension(:)  , intent(in)  :: vec
    real(KIND=dp), dimension(:)  , intent(out) :: res

    integer :: ii,ji,ni

    ni=size(M,1)
    res(ni)=vec(ni)/M(ni,ni)

    do ii=ni-1,1,-1
       res(ii)=vec(ii)
       do ji=ii+1,ni
          res(ii)=res(ii)-M(ii,ji)*res(ji)
       end do
       res(ii)=res(ii)/M(ii,ii)
    end do

  end subroutine invtrisup


  subroutine dHtodH(MdH, dH)
    ! abstract preocedure dH -> dH required to pass procedure to gmres
    ! to be updated to the objects of the dimensions of the perturbed
    ! sp hamiltonian dh and ddelta (in HF basis)
    real(KIND=dp), dimension(:), intent(out) :: MdH
    real(KIND=dp), dimension(:), intent(in)  :: dH

  end subroutine

  function dHtoreal(dH) result(res)
    ! abstract preocedure dH -> real required to pass procedure to gmres
    ! to be updated to the objects of the dimensions of the perturbed
    ! sp hamiltonian dh and ddelta (in HF basis)
    real(KIND=dp), dimension(:), intent(in)  :: dH
    real(KIND=dp)                            :: res

  end function

  function dHdHtoreal(dHl, dHr) result(res)
    ! abstract preocedure (dH,dH) -> complex required to pass procedure to gmres
    ! to be updated to the objects of the dimensions of the perturbed
    ! sp hamiltonian dh and ddelta (in HF basis)
    real(KIND=dp), dimension(:), intent(in)  :: dHl, dHr
    real(KIND=dp)                            :: res

  end function

end module gmres


