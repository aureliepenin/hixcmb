!Module misc_tools2
!
!
!  Implicit None 
!  Public 
!Contains 

!-----------------------------------------------------------------------


! A USEFUL ROMBERG INTEGRATION ROUTINE
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        function rombint(f,a,b,tol)
!  Rombint returns the integral from a to b of using Romberg integration.
!  The method converges provided that f(x) is continuous in (a,b).
!  f must be DOUBLE PRECISION and must be declared external in the calling
!  routine.  tol indicates the desired relative accuracy in the integral.
!
        implicit none
        integer, parameter :: MAXITER=20,MAXJ=5
        dimension g(MAXJ+1)
        DOUBLE PRECISION f
        external f
        DOUBLE PRECISION :: rombint
        DOUBLE PRECISION, intent(in) :: a,b,tol
        integer :: nint, i, k, jmax, j
        DOUBLE PRECISION :: h, gmax, error, g, g0, g1, fourj
!
        h=0.5d0*(b-a)
        gmax=h*(f(a)+f(b))
        g(1)=gmax
        nint=1
        error=1.0d20
        i=0
10        i=i+1
          if (i.gt.MAXITER.or.(i.gt.5.and.abs(error).lt.tol)) &
            go to 40
!  Calculate next trapezoidal rule approximation to integral.

          g0=0.0
            do 20 k=1,nint
            g0=g0+f(a+(k+k-1)*h)
20        continue
          g0=0.5d0*g(1)+h*g0
          h=0.5d0*h
          nint=nint+nint
          jmax=min(i,MAXJ)
          fourj=1.0
            do 30 j=1,jmax
!  Use Richardson extrapolation.
            fourj=4.0*fourj
            g1=g0+(g0-g(j))/(fourj-1.0)
            g(j)=g0
            g0=g1
30        continue
          if (abs(g0).gt.tol) then
            error=1.0-gmax/g0
          else
            error=gmax
          end if
          gmax=g0
          g(jmax+1)=g0
        go to 10
40      rombint=g0
        if (i.gt.MAXITER.and.abs(error).gt.tol)  then
          write(*,*) 'Warning: Rombint failed to converge; '
          write (*,*)'integral, error, tol:', rombint,error, tol
        end if
        
        end function rombint


!-----------------------------------------------------------------------

   SUBROUTINE splint(xa, ya, y2a, n, x, y)
!   USE nrtype
!
! Given the arrays xa(1:n) and ya(1:n) of length n, which tabulate a function
! (with the xa(i) in order), and given the array y2a(1:n), which is the output
! from the subroutine spline, and given a value of x, this routine returns a
! cubic spline interpolated value y.
! (adopted from Numerical Recipes in FORTRAN 77)
!
   INTEGER, PARAMETER :: DP = KIND(1.0D0)
   INTEGER:: n
   REAL(DP):: x, y, xa(n), y2a(n), ya(n)
   INTEGER:: k, khi, klo
   REAL(DP):: a, b, h

     klo=1
     khi=n
1   if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if (xa(k).gt.x) then
           khi=k
        else
           klo=k
        endif
        goto 1
     endif

     h=xa(khi)-xa(klo)
     if (h.eq.0.) pause 'bad xa input in splint'

     a=(xa(khi)-x)/h
     b=(x-xa(klo))/h
     y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.

     return
   END SUBROUTINE splint

!-----------------------------------------------------------------------

   SUBROUTINE spline(x, y, n, yp1, ypn, y2)
!   use nrtype
!
! Given arrays x(1:n) and y(1:n) containing a tabulated function, i.e.
! y(i)=f(x(i)), with x(1)<x(2)<...<x(n), and given values yp1 and ypn for
! the first derivative of the interpolating function at points 1 and n,
! respectively, this routine returns an array y2(1:n) of length n which
! contains the second derivatives of the interpolating function at the
! tabulated points x(i).  If yp1 and/or ypn are equal to 1.e30 or larger,
! the routine is signaled to set the corresponding boundary condition for a
! natural spline with zero second derivative on that boundary.
! Parameter: nmax is the largest anticipiated value of n
! (adopted from Numerical Recipes in FORTRAN 77)
!
   INTEGER, PARAMETER :: DP=KIND(1.0D0)
   INTEGER:: n
   INTEGER, PARAMETER:: nmax=20000
   REAL(DP):: yp1, ypn, x(n), y(n), y2(n)
   INTEGER:: i, k
   REAL(DP):: p, qn, sig, un, u(nmax)

     if (yp1.gt..99e30) then
        y2(1)=0.
        u(1)=0.
     else
        y2(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
     endif

     do i=2, n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/&
             & (x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
     enddo

     if (ypn.gt..99e30) then
        qn=0.
        un=0.
     else
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
     endif

     y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)

     do k=n-1, 1, -1
        y2(k)=y2(k)*y2(k+1)+u(k)
     enddo

     return
   END SUBROUTINE spline







!End Module misc_tools2
