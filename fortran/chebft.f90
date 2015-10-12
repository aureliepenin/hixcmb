Function chebft(a,b,n,func)
  Use nrtype; Use nrutil, Only : arth,outerprod
  Implicit None
  Real(SP), Intent(IN) :: a,b
  Integer(I4B), Intent(IN) :: n
  Real(SP), Dimension(n) :: chebft
  Interface
     Function func(x)
       Use nrtype
       Implicit None
       Real(SP), Dimension(:), Intent(IN) :: x
       Real(SP), Dimension(Size(x)) :: func
     End Function func
  End Interface
  Real(DP) :: bma,bpa
  Real(DP), Dimension(n) :: theta
  bma=0.5_dp*(b-a)
  bpa=0.5_dp*(b+a)
  theta(:)=PI_D*arth(0.5_dp,1.0_dp,n)/n
  chebft(:)=Matmul(Cos(outerprod(arth(0.0_dp,1.0_dp,n),theta)), &
       func(Real(Cos(theta)*bma+bpa,sp)))*2.0_dp/n
  !Write(*,*) 'r0 ', bma,bpa
  !Write(*,*) 'r1 ', theta
  !Write(*,*) 'r2 ', Real(Cos(theta)*bma+bpa,sp)
End Function chebft
