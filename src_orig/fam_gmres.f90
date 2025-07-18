module gmres

  use evolution

  implicit none


  abstract interface
    subroutine vectovec(vec_in, vec_out)
      import :: dp
      complex(KIND=dp), dimension(:), target, intent(in)   :: vec_in
      complex(KIND=dp), dimension(:), target, intent(out)  :: vec_out
    end subroutine

    function vectoreal(vec_in) result(scalar)
      import :: dp
      complex(KIND=dp), dimension(:), intent(in)  :: vec_in
      real(KIND=dp)                               :: scalar
    end function

    function vecvectocmplx(vec_l, vec_r) result(scalar)
      import :: dp
      complex(KIND=dp), dimension(:), intent(in)  :: vec_l, vec_r
      complex(KIND=dp)                            :: scalar
    end function

  end interface

  subroutine do_gmres(x, iter, nbPrd, res, &
       & b, A, norm_2, ScalProd, tol, itmax, rst, verb)
  ! 
  ! GMRES LINEAR SOLVER 
  ! 
  ! This module is based on the GMRES implementation of Choral of C. Pierre
  ! https://plmlab.math.cnrs.fr/cpierre1/choral, which itself references the 
  ! book of Youssef SAAD, 'Iterative methods for sparse linear system'
  ! https://www-users.cs.umn.edu/~saad/IterMethBook_2ndEd.pdf
  !
  ! GMRES (no preconditioning)
  ! 
  !
  ! for the linear system Ax = b
  !
  ! INPUT/OUTPUT :
  !     x = initial guess / solution
  !
  ! OUTPUT :
  !     res  = final residual
  !     iter = number of performed iterations
  !     iter = -1 = resolution failure     
  !
  ! INPUT :
  !     b        = RHS
  !     A        = A : x -> A x 
  !                           matrix/vector product (procedural)
  !     norm_2   = f : x ->  real 
  !                           norm (procedural)
  !     ScalProd = f : (x1, x2) ->  real 
  !                           Scalar product (procedural)
  !     tol      = tolerance
  !     itMax    = maximal iteration-number
  !     rst      = restart number
  !     verb     = verbosity
  ! 
  ! Remark :
  !   This implementation differs from the one of C. Pierre in the ordering of 
  !   of the arguments in the procedure A : x -> A x. Here, A is expected to be 
  !   subroutine A(x_in, x_out). 

    complex(KIND=dp), dimension(:), intent(inout) :: x
    integer                    , intent(out)   :: iter, nbPrd
    real(KIND=dp)              , intent(out)   :: res
    procedure(vectovec)                        :: A
    procedure(vectoreal)                       :: norm_2
    procedure(vecvectocmplx)                 :: ScalProd
    complex(KIND=dp), dimension(:), intent(in)    :: b
    real(KIND=dp)              , intent(in)    :: tol
    integer                    , intent(in)    :: itmax, rst, verb

    complex(KIND=dp), dimension(rst+1, size(x,1) ) :: V
    complex(KIND=dp), dimension(rst+1, rst       ) :: H

    complex(KIND=dp), dimension(size(x,1)) :: r, w
    complex(KIND=dp), dimension(rst      ) :: sn, cs, y
    complex(KIND=dp), dimension(rst +1   ) :: s

    real(KIND=dp) :: nb2, nr2, temp
    integer       :: nn, ii, kk

    if (verb>1) write(*,*) 'gmres       : gmres'

    iter  = 0
    nbPrd = 0


    ! compute nb2 = 1./||b|| 
    nb2 = norm_2(b)
    nn = size(x,1) 
    ! if ||b|| is small, use sqrt of the length of the array instead
    if (nb2/real(nn, dp)<1E-8_dp) nb2=sqrt(real(nn, dp))
    nb2 = 1._dp/nb2


    ! compute the initial residual vector r0 = b - A x0 
    call A(x, r)
    nbPrd = 1
    r   = b-r
    nr2 = norm_2(r) ! corresponds to beta in Y.SAAD
    res = nr2 * nb2 ! residual is ||r|| / ||b||

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

          call A(V(ii,:), w)
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
       call A(x, r)                  
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

    complex(KIND=dp), intent(in)  :: a,b
    complex(KIND=dp), intent(out) :: cs,sn

    complex(KIND=dp)              :: tmp

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

    complex(KIND=dp), dimension(:,:), intent(in)  :: M
    complex(KIND=dp), dimension(:)  , intent(in)  :: vec
    complex(KIND=dp), dimension(:)  , intent(out) :: res

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
  end subroutine invtrisu

end module gmres


