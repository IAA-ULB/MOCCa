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


  public

    integer           :: gmres_itmax = 100
    real(kind=dp)     :: gmres_tol = 1.0e-6_dp
    integer           :: gmres_histmax = 10
    integer           :: gmres_iter = 0
    real(kind=dp)     :: gmres_res = 100.0
    logical           :: verbose = .false.

    procedure(vectovec), pointer          :: apply_A
    complex(KIND=dp), allocatable         :: b(:), r0(:), beta(:), x_guess(:), x_gmres(:), y_minres(:)
    procedure(vectoreal), pointer         :: norm
    procedure(vecvectocmplx), pointer     :: dotprod
    ! NOTE : Y. SAAD assumes the maths convention for the innner product on a 
    !        complex vector space, i.e. linearity in the first component. 


    complex(KIND=dp), allocatable :: H(:,:) ! Hessenberg matrix
    !                                  | '-> 1:histmax
    !                                  '---> 1:histmax+1
    complex(KIND=dp), allocatable :: Q(:,:) ! history of the GMRES approximants v
    !                                  | '-> iter index 1:histmax+1
    !                                  '---> 1:xsize
    ! columns of Q contain the various orthonormalised Arnoldi vectors, Q(:,i) = v_i  
    ! where i ranges from 1 to m+1 cfr. Y.SAAD


  contains


  subroutine alloc_gmres(A_proc, b_rhs, itmax, histmax, tol, xsize, norm_proc, dotprod_proc)

    procedure(vectovec)          :: A_proc
    complex(KIND=dp), intent(in) :: b_rhs(:)
    integer         , intent(in) :: itmax, histmax, xsize
    real(kind=dp)   , intent(in) :: tol
    procedure(vectoreal)         :: norm_proc
    procedure(vecvectocmplx)     :: dotprod_proc


    ! set GMRES params
    gmres_itmax = itmax
    gmres_histmax = histmax
    gmres_tol = tol
    gmres_iter = 0

    ! set the procedure pointers
    apply_A => A_proc
    norm => norm_proc
    dotprod => dotprod_proc


    ! allocate the GMRES work space
    allocate(H(histmax+1,histmax))
    allocate(Q(xsize,histmax+1))
    allocate(b(xsize))
    allocate(r0(xsize))
    allocate(x_guess(xsize))
    allocate(x_gmres(xsize))
    allocate(y_minres(histmax+1))
    allocate(beta(histmax+1))

    H(:,:)      = 0
    beta(:)     = 0
    Q(:,:)      = 0
    x_guess(:)  = 0
    x_gmres(:)  = 0
    y_minres(:) = 0


    b = b_rhs

    print *, 'Set up GMRES work space : '
    print *, '    histmax = ', gmres_histmax
    print *, '    itmax   = ', gmres_itmax
    print *, '    tol     = ', gmres_tol
    print *, '    x dim   = ', xsize
    print *, '   ||b||    = ', norm(b)


  end subroutine alloc_gmres

  subroutine init_gmres(x0)
    complex(KIND=dp), intent(in)   :: x0(:)
    integer                        :: i

    ! Initialize the GMRES solver and the first Arnoldi vector Q(:,1).
    print *, "GMRES iter", gmres_iter

    x_guess = x0

    ! compute the initial residual vector r0 = b - A x0 
    call apply_A(x0, r0)

    r0 = b - r0
    beta(1) = norm(r0)
    Q(:,1)  = r0(:) / beta(1)

    gmres_res =  beta(1) / norm(b)
    gmres_iter = 1

    print * , '    res = ', gmres_res

    if (verbose) then 

      print *, 'r0 = b - A*x0 = '
      do i=1,size(r0,1)
          print "(*('(', F8.5, ',', F8.5, ') ', :))", r0(i)
      enddo

      print *, 'H = '
      do i=1,size(H,1)
          print "(*('(', F8.5, ',', F8.5, ') ', :))",  H(i, :)
      enddo
      print *, 'Q = '
      do i=1,size(Q,1)
          print "(*('(', F8.5, ',', F8.5, ') ', :))",  Q(i, :)
      enddo

    endif


  end subroutine init_gmres

  subroutine iterate_gmres()
    complex(KIND=dp) :: wj(size(Q,1))
    integer          :: i

    print *, "GMRES iter", gmres_iter ! idx j in Y. Saad

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! perform one iteration wj =  A(vj)

    call apply_A(Q(:,gmres_iter), wj)

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! preform modified Gram-Schmidt orthogonalisation

    do i = 1, gmres_iter
        H(i,gmres_iter) = dotprod( wj , Q(:,i) );
        ! note that DIRQFAM uses dot_product(Q, w) since they follow the physics 
        ! convention dot_product(x, y) = conjg(x) * y, while here the math convention 
        ! of Y SAAD is used. 
        wj = wj - H(i,gmres_iter) * Q(:,i);
    end do

    H(gmres_iter+1,gmres_iter) = norm(wj);

    if(abs(H(gmres_iter+1,gmres_iter)) < 1E-8_dp) then
      print *, '||wj|| < 1e-8. Stopping GMRES'
    else
      Q(:,gmres_iter+1) = wj / norm(wj);
    endif
    
    if (verbose) then 
      print *, 'H = '
      do i=1,size(H,1)
          print "(*('(', F8.5, ',', F8.5, ') ', :))",  H(i, :)
      enddo
      print *, 'Q = '
      do i=1,size(Q,1)
          print "(*('(', F8.5, ',', F8.5, ') ', :))",  Q(i, :)
      enddo

    endif

    call extract_x_gmres()

    gmres_iter = gmres_iter + 1

  end subroutine iterate_gmres

  subroutine extract_x_gmres()
    complex(KIND=dp) :: H_tmp(gmres_iter+1,gmres_iter)
    complex(KIND=dp) :: beta_tmp(gmres_iter+1)
    complex(KIND=dp), allocatable :: work(:)
    complex(KIND=dp)              :: workquery(1)
    integer           :: info_zgels
    complex(KIND=dp) :: Ax(size(Q,1))
    integer          :: i, j

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! solve the min res problem, y = argmin( || beta - H y|| )
    ! the lapack routine zgels

    ! take temporary copies since lapack modifies the input matrices

    H_tmp   ( 1:gmres_iter+1 , 1:gmres_iter ) = H   ( 1:gmres_iter+1 , 1:gmres_iter )
    beta_tmp( 1:gmres_iter+1                ) = beta( 1:gmres_iter+1                )

    ! perform a so-called workquery to obtain the optimal work dimension
    call zgels('N', gmres_iter+1 , gmres_iter, 1, H_tmp, size(H_tmp,1),   & 
        &       beta_tmp, size(beta_tmp), workquery , -1, info_zgels)

    print *, 'zgels exit : ', info_zgels

    allocate( work( 1:int(real(workquery(1) + 0.5d0)) ) );

    ! solve the overdetermined minimisation problem with zgels
    call zgels('N', gmres_iter+1 , gmres_iter, 1, H_tmp, size(H_tmp,1),   &
        &       beta_tmp, size(beta_tmp), work, size(work), info_zgels)

    print *, 'zgels exit : ', info_zgels


    y_minres = 0
    y_minres(1:gmres_iter) = beta_tmp( 1:gmres_iter)

    gmres_res = abs(beta_tmp(gmres_iter+1)) / norm(b)

    print *, '   res : min( || beta - H y|| ) / ||b|| = ', abs(beta_tmp(gmres_iter+1)) / norm(b)
    print *, '   res :    || beta - H y_min|| / ||b|| = ', norm( beta - matmul(H,y_minres) ) / norm(b)

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! calculate GMRES approximant Xm = X0 + Vm * Ym 

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! by hand
    x_gmres = x_guess

    do j=1,gmres_iter
      x_gmres = x_gmres + Q(:, j) * y_minres(j)
    enddo 

    if (verbose) then 
      print *, 'Ym = '
      do i=1,size(y_minres)
        print "(*('(', F8.5, ',', F8.5, ') ', :))", y_minres(i)
      enddo

      print *, 'Xm = '
      do i=1,size(x_gmres)
        print "(*('(', F8.5, ',', F8.5, ') ', :))", x_gmres(i)
      enddo

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! by LAPACK zgemv

      x_gmres = x_guess

      call zgemv( 'N' , size(x_gmres) , gmres_iter , dcmplx(1.0_dp, 0.0_dp), Q , size(Q,1) , &
                & y_minres , 1 , dcmplx(1.0_dp, 0.0_dp), x_gmres , 1 )

      print *, 'Xm lapack = '
      do i=1,size(x_gmres)
          print "(*('(', F8.5, ',', F8.5, ') ', :))", x_gmres(i)
      enddo

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! eval res in different ways 

      call apply_A(x_gmres, Ax)

      print *, 'A*Xm = '
      do i=1,size(x_gmres)
        print "(*('(', F8.5, ',', F8.5, ') ', :))", Ax(i)
      enddo

      print *, ' rm = b - A*xm = '
      do i=1,size(x_gmres)
        print "(*('(', F8.5, ',', F8.5, ') ', :))", b(i) - Ax(i)
      enddo

      print *, 'res : ', norm(b - Ax) / norm(b)
    endif

  end subroutine extract_x_gmres

  

  subroutine do_gmres_choral(x, iter, nbPrd, res, &
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
  ! 
  ! /!\ NOTE : this implementation is NOT operational for complex A !

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
    integer       :: nn, ii, kk, i

    if (verb>1) write(*,*) 'gmres       : gmres'

    iter  = 0
    nbPrd = 0


    ! compute nb2 = 1./||b|| 
    nb2 = norm_2(b)
    if (verb>2) print*, '||b||', nb2
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
             H(kk,ii)= conjg(ScalProd( w, V(kk,:)))
             w = w - H(kk,ii)*V(kk,:)
          end do

          H(ii+1,ii) = norm_2(w)

          V(ii+1, : ) = w / H(ii+1,ii)

          print *, 'H = '
          do i=1,rst+1
              print "(*('(', F8.5, ',', F8.5, ') ', :))",  H(i, :)
          enddo
          print *, 'Q = '
          do i=1,6
              print "(*('(', F8.5, ',', F8.5, ') ', :))",  V(i, :)
          enddo

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

  end subroutine do_gmres_choral

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
  end subroutine invtrisup

end module gmres


