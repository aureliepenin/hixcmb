Subroutine cisi(x,ci,si)
  Use nrtype; Use nrutil, Only : nrerror
  Implicit None
  Real(DP), Intent(IN) :: x
  Real(DP), Intent(OUT) :: ci,si
  Integer(I4B), Parameter :: MAXIT=1000
  Real(DP), Parameter :: EPS=Epsilon(x),FPMIN=4.0_dp*Tiny(x),BIG=Huge(x)*EPS,TMIN=2.0
  Integer(I4B) :: i,k
  Real(DP) :: a,err,fact,sign,sum,sumc,sums,t,term
  Complex(DPC) :: h,b,c,d,del
  Logical(LGT) :: odd
  t=Abs(x)
  If (t == 0.0) Then
     si=0.0
     ci=-BIG
     Return
  End If
  If (t > TMIN) Then
     b=Cmplx(1.0_dp,t,kind=dpc)
     c=BIG
     d=1.0_dp/b
     h=d
     Do i=2,MAXIT
        a=-(i-1)**2
        b=b+2.0_dp
        d=1.0_dp/(a*d+b)
        c=b+a/c
        del=c*d
        h=h*del
        If (absc(del-1.0_dp) <= EPS) Exit
     End Do
     If (i > MAXIT) Call nrerror('continued fraction failed in cisi')
     h=Cmplx(Cos(t),-Sin(t),kind=dpc)*h
     ci=-Real(h)
     si=PIO2+Aimag(h)
  Else
     If (t < Sqrt(FPMIN)) Then
        sumc=0.0
        sums=t
     Else
        sum=0.0
        sums=0.0
        sumc=0.0
        sign=1.0
        fact=1.0
        odd=.True.
        Do k=1,MAXIT
           fact=fact*t/k
           term=fact/k
           sum=sum+sign*term
           err=term/Abs(sum)
           If (odd) Then
              sign=-sign
              sums=sum
              sum=sumc
           Else
              sumc=sum
              sum=sums
           End If
           If (err < EPS) Exit
           odd=.Not. odd
        End Do
        If (k > MAXIT) Call nrerror('MAXIT exceeded in cisi')
     End If
     si=sums
     ci=sumc+Log(t)+EULER
  End If
  If (x < 0.0) si=-si

Contains
  !BL
  Function absc(z)
    Implicit None
    Complex(DPC), Intent(IN) :: z
    Real(DP) :: absc
    absc=Abs(Real(z))+Abs(Aimag(z))
  End Function absc
End Subroutine cisi
