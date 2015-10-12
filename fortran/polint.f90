Subroutine polint(xa,ya,x,y,dy,n)
  Use nrtype; Use nrutil, Only : assert_eq,iminloc,nrerror
  Implicit None
  Real(DP), Dimension(1:n), Intent(IN) :: xa,ya
  Real(DP), Intent(IN) :: x
  Real(DP), Intent(OUT) :: y,dy
  Integer(I4B) :: m,n,ns
  Real(DP), Dimension(1:n) :: c,d,den,ho
  !n=assert_eq(Size(xa),Size(ya),'polint')
  c=ya
  d=ya
  ho=xa-x
  ns=iminloc(ABS(x-xa))
  y=ya(ns)
  ns=ns-1
  Do m=1,n-1
     den(1:n-m)=ho(1:n-m)-ho(1+m:n)
     If (Any(den(1:n-m) == 0.0)) &
          Call nrerror('polint: calculation failure')
     den(1:n-m)=(c(2:n-m+1)-d(1:n-m))/den(1:n-m)
     d(1:n-m)=ho(1+m:n)*den(1:n-m)
     c(1:n-m)=ho(1:n-m)*den(1:n-m)
     If (2*ns < n-m) Then
        dy=c(ns+1)
     Else
        dy=d(ns)
        ns=ns-1
     End If
     y=y+dy
  End Do
End Subroutine polint
