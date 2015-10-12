Function splin2(x1a,x2a,n1,n2,ya,y2a,x1,x2)
  !USE nrtype
  ! USE nrutil, ONLY : assert_eq
  !USE nr, ONLY : spline,splint
  Implicit None
  Integer :: j,m,ndum,n1,n2
  Real*8, Dimension(1:n1)  , Intent(IN) :: x1a
  Real*8, Dimension(1:n2)  , Intent(IN) :: x2a
  Real*8, Dimension(1:n1,1:n2), Intent(IN) :: ya,y2a
  Real*8, Intent(IN) :: x1,x2
  Real*8  :: splin2
  Real*8, Dimension(1:n1) :: yytmp,y2tmp2
  m    = Size(x1a)!assert_eq(size(x1a),size(ya,1),size(y2a,1),'splin2: m')
  ndum = Size(x2a)!assert_eq(size(x2a),size(ya,2),size(y2a,2),'splin2: ndum')
  Do j=1,n1
     !yytmp(j)=splint(x2a,ya(j,:),y2a(j,:),x2)
     Call splint(x2a,ya(j,:),y2a(j,:),n2,x2,yytmp(j))
  End Do
  !call spline(x1a,yytmp,1.0e30_sp,1.0e30_sp,y2tmp2)
  Call spline(x1a,yytmp,n1,3.d30,3.d30,y2tmp2)
  !splin2=splint(x1a,yytmp,y2tmp2,x1)
  Call splint(x1a,yytmp,y2tmp2,n1,x1,splin2)
End Function splin2
