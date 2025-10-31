module gmres

   use geninfo, only : dp

  implicit none


  abstract interface
    subroutine vectovec(vec_in, vec_out)
      import :: dp
      complex(KIND=dp), dimension(:), intent(in)   :: vec_in
      complex(KIND=dp), dimension(:), intent(out)  :: vec_out
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

    integer           :: gmres_itermax = 100 ! Max number of iterations, i.e. # apply_A(v) calls. 
    integer           :: gmres_iter = 0      ! Current number of iterations, i.e. # apply_A(v) calls. 
    integer           :: gmres_histmax = 10  ! Max history size, i.e. # of vectors in the GMRES workspace
    integer           :: gmres_hist = 0      ! current number of vectors in GMRES workspace
    real(kind=dp)     :: gmres_res = 100.0
    real(kind=dp)     :: gmres_precision = 1.0e-6_dp
    logical           :: gmres_verbose = .false.

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


  subroutine alloc_gmres(A_proc, b_rhs, itmax, histmax, prec, norm_proc, dotprod_proc)
    !---------------------------------------------------------------------------------------
    ! Allocate the GMRES workspace and set up all routines
    !
    ! INPUT : 
    !    A_proc (procedural) = the update routine A : x -> A(x), must be of the form vectovec()
    !    b_rhs (complex(:))  = right-hand-side of Ax=b one tries to solve
    !    itmax (integer)     = maximal number of GMRES iterations, i.e. of A_proc calls
    !    histmax (integer)   = maximal size of GMRES work space, i.e. number of vectors stored
    !    prec (real)         = precision of the GMRES solver, stopped when residual < prec
    !    norm_proc (procedural) 
    !                        = routine to compute the norm of a vector x, 
    !                          must be of the form vectoreal()
    !    dotprod_proc (procedural) 
    !                        = routine to compute the dotproduct of two vectors, 
    !                          must be of the form vecvectocmplx()
    ! 
    !---------------------------------------------------------------------------------------

    procedure(vectovec)          :: A_proc
    complex(KIND=dp), intent(in) :: b_rhs(:)
    integer         , intent(in) :: itmax, histmax
    integer                      :: xsize
    real(kind=dp)   , intent(in) :: prec
    procedure(vectoreal)         :: norm_proc
    procedure(vecvectocmplx)     :: dotprod_proc

    1 format(50('='), /, "Set up GMRES work space : ")
    2 format(6x, "histmax = ", i8)
    3 format(6x, "itmax   = ", i8)
    4 format(6x, "prec     = ", es8.1)
    5 format(6x, "x dim   = ", i8)
    6 format(6x, "||b||   = ", f8.2)

    xsize = size(b_rhs, 1)

    ! set GMRES params
    gmres_itermax = itmax
    gmres_histmax = histmax
    gmres_precision = prec
    gmres_hist = 0
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

    print 1
    print 2, gmres_histmax
    print 3, gmres_itermax
    print 4, gmres_precision
    print 5, xsize
    print 6, norm(b)


  end subroutine alloc_gmres

  subroutine dealloc_gmres()

    deallocate(H)
    deallocate(Q)
    deallocate(b)
    deallocate(r0)
    deallocate(x_guess)
    deallocate(x_gmres)
    deallocate(y_minres)
    deallocate(beta)


  end subroutine dealloc_gmres

  subroutine init_gmres(x0)
    !-------------------------------------------------------------------------------
    ! Initialize the GMRES solver and the first Arnoldi vector Q(:,1)
    !
    ! INPUT : 
    !    x0 (complex(:)) = initial guess of a solution of Ax = b
    !-------------------------------------------------------------------------------
    complex(KIND=dp), intent(in)   :: x0(:)
    integer                        :: i


    1 format(50('-'))
    2 format('GMRES iteration : ', i5)
    3 format('  residual = ', es10.3)

    print 1
    print 2, gmres_iter

    x_guess = x0

    ! compute the initial residual vector r0 = b - A x0 
    call apply_A(x0, r0)
    gmres_iter = gmres_iter + 1

    r0 = b - r0
    beta(1) = norm(r0)
    Q(:,1)  = r0(:) / beta(1)

    gmres_res =  beta(1) / norm(b)
    gmres_hist = 1

    print 3, gmres_res

    if (gmres_verbose) then 

      print *, "||r0||", norm(r0)
      print *, "||x0||", norm(x0)

      print *, 'x0 = '
      do i=1,10
          print "(*('(', E10.2, ',', E10.2, ') ', :))", x0(i)
      enddo


      print *, 'r0 = b - A*x0 = '
      do i=1,10
          print "(*('(', E10.2, ',', E10.2, ') ', :))", r0(i)
      enddo

      print *, 'H = '
      do i=1,size(H,1)
          print "(*('(', F8.5, ',', F8.5, ') ', :))",  H(i, :)
      enddo
      print *, 'Q = '
      do i=1,10
          print "(*('(', F8.5, ',', F8.5, ') ', :))",  Q(i, :)
      enddo

    endif


  end subroutine init_gmres

  subroutine iterate_gmres()
    !-------------------------------------------------------------------------------
    ! Perform one GMRES iteration. 
    ! Note that the approximate solution x_gmres is NOT computed. This is only to 
    ! to be done once GMRES is converged. 
    !-------------------------------------------------------------------------------

    complex(KIND=dp) :: wj(size(Q,1))
    integer          :: i

    1 format(50('-'))
    2 format('GMRES iteration : ', i5)

    print 1
    print 2, gmres_iter

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! perform one Arnoldi iteration wj =  A(vj)

    call apply_A(Q(:,gmres_hist), wj)
    gmres_iter = gmres_iter + 1

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! preform modified Gram-Schmidt orthogonalisation

    do i = 1, gmres_hist
        H(i,gmres_hist) = dotprod( wj , Q(:,i) );
        ! note that DIRQFAM uses dot_product(Q, w) since they follow the physics 
        ! convention dot_product(x, y) = conjg(x) * y, while here the math convention 
        ! of Y SAAD is used. 
        wj = wj - H(i,gmres_hist) * Q(:,i);
    end do

    H(gmres_hist+1,gmres_hist) = norm(wj);

    if(abs(H(gmres_hist+1,gmres_hist)) < 1E-8_dp) then
      print *, '||wj|| < 1e-8. Stopping GMRES'
      gmres_itermax = gmres_iter
    else
      Q(:,gmres_hist+1) = wj / norm(wj);
    endif
    
    if (gmres_verbose) then 
      print *, 'H = '
      do i=1,size(H,1)
          print "(*('(', F8.5, ',', F8.5, ') ', :))",  H(i, :)
      enddo
      print *, 'Q = '
      do i=1,10
          print "(*('(', F8.5, ',', F8.5, ') ', :))",  Q(i, :)
      enddo

    endif

    call minimize_residual()

    gmres_hist = gmres_hist + 1

    if (gmres_hist > gmres_histmax) then
      print *, 'restart gmres solver'

      ! compute the current best GMRES approximation to X
      call extract_x_gmres()

      ! reinitialise GMRES solver with this approximant as the initial guess
      H(:,:)      = 0
      Q(:,:)      = 0
      call init_gmres(x_gmres)

    endif

  end subroutine iterate_gmres

  subroutine minimize_residual()
    !-------------------------------------------------------------------------------
    ! Solve the min res problem, y = argmin( || beta - H y|| ) using 
    ! the lapack routine zgels
    !-------------------------------------------------------------------------------

    complex(KIND=dp) :: H_tmp(gmres_hist+1,gmres_hist)
    complex(KIND=dp) :: beta_tmp(gmres_hist+1)
    complex(KIND=dp), allocatable :: work(:)
    complex(KIND=dp)              :: workquery(1)
    integer          :: info_zgels

    1 format('  residual = ', es10.3)
    
    ! take temporary copies since lapack modifies the input matrices
    H_tmp   ( 1:gmres_hist+1 , 1:gmres_hist ) = H   ( 1:gmres_hist+1 , 1:gmres_hist )
    beta_tmp( 1:gmres_hist+1                ) = beta( 1:gmres_hist+1                )

    ! perform a so-called workquery to obtain the optimal work dimension
    call zgels('N', gmres_hist+1 , gmres_hist, 1, H_tmp, size(H_tmp,1),   & 
        &       beta_tmp, size(beta_tmp), workquery , -1, info_zgels)

    if (info_zgels .ne. 0) then
      print *, 'LAPACK zgels workquery failed'
    endif

    allocate( work( 1:int(real(workquery(1) + 0.5d0)) ) );

    ! solve the overdetermined minimisation problem with zgels
    call zgels('N', gmres_hist+1 , gmres_hist, 1, H_tmp, size(H_tmp,1),   &
        &       beta_tmp, size(beta_tmp), work, size(work), info_zgels)

    if (info_zgels .ne. 0) then
      print *, 'LAPACK zgels failed'
    endif

    y_minres = 0
    y_minres(1:gmres_hist) = beta_tmp( 1:gmres_hist)

    gmres_res = abs(beta_tmp(gmres_hist+1)) / norm(b)

    if (abs(gmres_res - norm( beta - matmul(H,y_minres) ) / norm(b)) > 1e-10_dp) then
      print *, 'GMRES residuals differ, i.e. min( || beta - H y|| ) /= || beta - H y_min|| '
    endif

    print 1, gmres_res

  end subroutine minimize_residual


  subroutine extract_x_gmres()
    !-------------------------------------------------------------------------------
    ! calculate GMRES approximant Xm = X0 + Vm * Ym 
    !-------------------------------------------------------------------------------

    complex(KIND=dp) :: Ax(size(Q,1))
    integer          :: i

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! by hand
    x_gmres = x_guess

    do i=1,gmres_hist
      x_gmres = x_gmres + Q(:, i) * y_minres(i)
    enddo 

    if (gmres_verbose) then 

      print *, "||xm||", norm(x_gmres)
      print *, "||x0||", norm(x_guess)
      print *, "||Ym||", norm(y_minres)

      do i=1,gmres_hist
        print *, "||vi||", norm(Q(:, i))
      enddo

      print *, 'Ym = '
      do i=1,size(y_minres)
        print "(*('(',E10.2, ',',E10.2, ') ', :))", y_minres(i)
      enddo

      print *, 'Xm = '
      do i=1,10
        print "(*('(',E10.2, ',',E10.2, ') ', :))", x_gmres(i)
      enddo

      ! ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! ! by LAPACK zgemv

      ! x_gmres = x_guess

      ! call zgemv( 'N' , size(x_gmres) , gmres_iter , dcmplx(1.0_dp, 0.0_dp), Q , size(Q,1) , &
      !           & y_minres , 1 , dcmplx(1.0_dp, 0.0_dp), x_gmres , 1 )

      ! print *, 'Xm lapack = '
      ! do i=1,10
      !     print "(*('(', F8.5, ',', F8.5, ') ', :))", x_gmres(i)
      ! enddo

      ! ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! ! eval res in different ways 

      ! call apply_A(x_gmres, Ax)

      ! print *, 'A*Xm = '
      ! do i=1,10
      !   print "(*('(', F8.5, ',', F8.5, ') ', :))", Ax(i)
      ! enddo

      ! print *, ' rm = b - A*xm = '
      ! do i=1,10
      !   print "(*('(', F8.5, ',', F8.5, ') ', :))", b(i) - Ax(i)
      ! enddo

      ! print *, 'res : ', norm(b - Ax) / norm(b)
    
    endif

  end subroutine extract_x_gmres


end module gmres


