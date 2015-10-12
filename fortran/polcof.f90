Subroutine polcof(xa,ya,polc,n)
  Use nrtype; Use nrutil, Only : assert_eq,iminloc
  !Use nr, Only : polint
  Implicit None
  Real(DP), Dimension(1:n), Intent(IN)    :: xa,ya
  Real(DP), Dimension(1:n), Intent(INOUT) :: polc
  Integer(I4B) :: j,k,m,n
  Real(DP) :: dy
  Real(DP), Dimension(1:n) :: x,y
  !Write(*,*) ' n pol ',size(xa),size(ya)!n
  !n=assert_eq(Size(xa),Size(ya),'polcof')
  x=xa
  y=ya
  Do j=1,n
     m=n+1-j 
     Call polint(x(1:m),y(1:m),0.0_dp,polc(j),dy,m)
     k=iminloc(Abs(x(1:m)))
     Where (x(1:m) /= 0.0) y(1:m)=(y(1:m)-polc(j))/x(1:m)
     y(k:m-1)=y(k+1:m)
     x(k:m-1)=x(k+1:m)
  End Do
End Subroutine polcof
