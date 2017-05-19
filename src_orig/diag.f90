module diag
 !=======================================================================
 !  #######   ##   #    # #####   ##   #      #    #  ####
 !     #     #  #  ##   #   #    #  #  #      #    # #
 !     #    #    # # #  #   #   #    # #      #    #  ####
 !     #    ###### #  # #   #   ###### #      #    #      #
 !     #    #    # #   ##   #   #    # #      #    # #    #
 !     #    #    # #    #   #   #    # ######  ####   ####
 !
 !  Copyright W. Ryssens & M. Bender
 !
 !=======================================================================

 use compilation

contains
      subroutine diagon (a,ndim,n,v,d,wd)
!c..............................................................................
!c  diagonalization of a real symmetric matrix (a(i,j)                         .
!c     input : a  n*n matrix           (with declared dimensions ndim*ndim)    .
!c             a(i,k) * v(k,j) = v(i,k) * d(k)                                 .
!c     output: v block of eigenvectors (with declared dimensions ndim*ndim)    .
!c             d eigenvalues in ascending order                                .
!c             wd working array                                                .
!c             the content of a is lost (actualy, a(i,i) = d(i)                .
!c             a(i,j) = v(i,k) * d(k) * v(j,k)                                 .
!c                      v is an orthogonal matrix : v(i,k)*v(j,k) = delta_ij   .
!c..............................................................................

      implicit real*8 (a-h,o-z)

      parameter (zero=0.0d0,one=1.0d0,two=2.0d0)
      parameter (eps=9.0d-12,epsd=1.0d-16,tol=1.0d-36,jstop=30)
      dimension a(ndim,ndim),v(ndim,ndim),d(ndim),wd(ndim)

!c........................................................................
      v(1,1) = one
      d(1)   = a(1,1)
      if (n.le.1) return

      do i=1,ndim*ndim
        v(i,1) = a(i,1)
      enddo

      do 9 ii=2,n
        i     = n+2-ii
        l     = i-1
        h     = zero
        scale = zero
        if (l.eq.1) go to 100
        do k=1,l
          scale = scale + abs(v(i,k))
        enddo
        if (scale.gt.tol) go to 3
  100   continue
        wd(i) = v(i,l)
        d(i) = h
        go to 9
    3   do k=1,l
          v(i,k) = v(i,k)/scale
          h      = h + v(i,k)**2
        enddo
        f      = v(i,l)
        g      =-sign(sqrt(h),f)
        wd(i)  = g*scale
        h      = h-f*g
        v(i,l) = f-g
        f      = zero
        do j=1,l
          v(j,i) = v(i,j)/(h*scale)
          g = zero
          do k=1,j
            g = g + v(j,k)*v(i,k)
          enddo
          j1 = j + 1
          if (j1.le.l) then
            do k=j1,l
              g = g + v(k,j)*v(i,k)
            enddo
          endif
          wd(j) = g/h
          f     = f + wd(j)*v(i,j)
        enddo
        hh = f/(h+h)
        do j=1,l
          f     = v(i,j)
          g     = wd(j) - hh*f
          wd(j) = g
          do k=1,j
            v(j,k) = v(j,k) - f*wd(k) - g*v(i,k)
          enddo
        enddo
        do k=1,l
          v(i,k) = scale * v(i,k)
        enddo
        d(i)=h
    9 continue

      d(1)  = zero
      wd(1) = zero
      do i=1,n
        l=i-1
        if ((abs(d(i)).ge.epsd).and.(l.ne.0)) then
          do j=1,l
            g = zero
            do k=1,l
              g = g + v(i,k)*v(k,j)
            enddo
            do k=1,l
              v(k,j) = v(k,j) - g*v(k,i)
            enddo
          enddo
        endif
        d(i)   = v(i,i)
        v(i,i) = one
        if (l.ne.0) then
          do j=1,l
            v(i,j) = zero
            v(j,i) = zero
          enddo
        endif
      enddo

      do i=2,n
        wd(i-1) = wd(i)
      enddo

      wd(n) = zero
      b     = zero
      f     = zero
      do 212 l=1,n
        j = 0
        h = eps * ( abs(d(l)) + abs(wd(l)) )
        if (b.lt.h) b = h
        m = l - 1
  202   m = m + 1
        if (m.gt.n) go to 203
        if (abs(wd(m))-b) 203,203,202
  203   continue
        if (m.eq.l) go to 211
  204   continue
        if (j.eq.jstop) stop 'diagon jstop'
        j = j + 1
        p = (d(l+1)-d(l))/(two*wd(l))
        r = sqrt(p*p+one)
        h =  d(l)-wd(l)/(p+sign(r,p))
        do i=l,n
          d(i) = d(i) - h
        enddo
        f  = f+h
        p  = d(m)
        c  = one
        s  = zero
        m1 = m  - 1
        ml = m1 + l
        do ii=l,m1
          i = ml - ii
          g = c*wd(i)
          h = c*p
          if (abs(p).ge.abs(wd(i))) then
            c       = wd(i)/p
            r       = sqrt(c*c+one)
            wd(i+1) = s*p*r
            s       = c/r
            c       = one/r
          else
            c       = p/wd(i)
            r       = sqrt(c*c+one)
            wd(i+1) = s*wd(i)*r
            s       = one/r
            c       = c/r
          endif
          p      = c*d(i) - s*g
          d(i+1) = h + s*(c*g+s*d(i))
          do k=1,n
            h        = v(k,i+1)
            v(k,i+1) = s*v(k,i) + c*h
            v(k,i)   = c*v(k,i) - s*h
          enddo
        enddo
        wd(l) = s*p
        d(l)  = c*p
        if (abs(wd(l))-b) 211,211,204
  211   continue
        d(l) = d(l)+f
  212 continue

      n1 = n-1
      do i=1,n1
        k  = i
        p  = d(i)
        ii = i + 1
        do j=ii,n
          if (d(j).lt.p) then
            k = j
            p = d(j)
          endif
        enddo
        if (k.ne.i) then
          d(k) = d(i)
          d(i) = p
          do j=1,n
            p      = v(j,i)
            v(j,i) = v(j,k)
            v(j,k) = p
          enddo
        endif
      enddo

      do i=1,n
        a(i,i) = d(i)
      enddo

      return
      end subroutine diagon
end module diag
