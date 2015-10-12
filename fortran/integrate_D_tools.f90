!======================
Function integrate_D(a)
!======================

  ! a is the expansion factor

  Use nrtype
  Use nr
  Use ode_path

  Implicit None
  
  Real(DP),    Intent(IN) :: a
  Real(DP)                :: integrate_D

  Integer(I4B), Parameter :: Nvar=2
  Integer(I4B) :: i,nrhs
  Real(DP)     :: eps,h1,hmin,x1,x2,a0,ni
  Real(DP), Dimension(Nvar) :: ystart
  Common nrhs

  Interface
     Subroutine derivs(x,y,dydx)
       Use nrtype
       Implicit None
       Real(DP), Intent(IN) :: x
       Real(DP), Dimension(:), Intent(IN) :: y
       Real(DP), Dimension(:), Intent(OUT) :: dydx
     End Subroutine derivs
  End Interface
  ! 
  !===============================================================================
  ! Solve using RK of 4th order with adaptive step-size (\dot = \partial_t)
  ! \ddot D + 2H\dot D - 3/2 (H0^2\Om0/a^3) D f(k) = 0
  ! We make the change of variable u = ln a so that
  ! H^2 D'' + [2H^2 - 3/2 H0^2Om0/a^3 ]D' = 3/2 (H0^2Om0/a^3) D f(k)
  ! where H^2 = H0^2(Om0/a^2+Oml)
  ! Assume 
  !     y1 = D'
  !     y2 = D
  ! so that 
  !     \dot y1 = D'' = [3/2 H0^2Om0/a^3 f(k) y2 + y1[3/2 H0^2Om0/a^3 - 2H^2]]/H^2
  !     \dot y2 = D'  = y1
  !===============================================================================
  !
  !
  nrhs = 0
  !
  ! Integration limits
  ! ------------------
  !
  a0   = 1.e-6! (z=1.d4)
  x1   = log(a0)    ! a = 0
  x2   = log(a)     ! alog(1.)    ! a = 1
  !
  ! Set-up initial conditions
  ! -------------------------
  !
  ! We here use the fact that D \propto a when a->0 since H(a)^2->H0^2Om0/a^3
  ! Hence D'= D and we choose D = a^ni initially with
  ! n = 1/4 [ \sqrt{1+24f(k)} -1 ]
  ! The normalisation is arbitrary and will be absorbed eventually in sigma8 or As
  !
  ! Note that when f is introduced, we want to make sure that f->1 when a -> 0 otherwise the power-law
  ! solution does not hold...
  ! 
  ni = 1.
  !Write(*,*) 'F0 ',f_func(k,a0),ni
  !
  ystart(1) = ni*a0**ni ! D' = ni*D
  ystart(2) = a0**ni    ! D  = a^ni
  !
  ! Integration parameters
  ! ----------------------
  !  
  eps   = 1.0e-6_dp
  h1    = 0.01_dp
  hmin  = 0.0
  dxsav = (x2-x1)/20.0_dp
  save_steps=.True.
  !
  Call odeint(ystart,x1,x2,eps,h1,hmin,derivs,rkqs)
  !Write(*,'(/1x,a,t30,i3)') 'Successful steps          :',nok
  !Write(*,'(1x,a,t30,i3)')  'Bad steps                 :',nbad
  !Write(*,'(1x,a,t30,i3)')  'Function evaluations      :',nrhs
  !Write(*,'(1x,a,t30,i3)')  'Stored intermediate values:',kount
  !Write(*,'(/1x,t9,a,t20,a,t33,a)') 'X','Integral'
  !Do i=1,kount
  !   Write(*,'(1x,f10.4,2x,2f14.6)') xp(i),yp(2,i),ni
  !End Do
  ! The ouput results in in yp(2,kount)
  integrate_D = yp(2,kount)
  !
  deallocate(xp,yp) ! Important if not keep allocating xp and yp in odeint...
  !
End Function integrate_D
!=======================

!==========================
Subroutine derivs(x,y,dydx)
!==========================

  Use nrtype
  Use cosmo_tools

  Implicit None

  Real(DP), Intent(IN) :: x
  Real(DP), Dimension(:), Intent(IN) :: y
  Real(DP), Dimension(:), Intent(OUT) :: dydx
  Real(DP) :: H2_func,a,H2,cst
  External :: H2_func

  Integer(I4B) :: nrhs
  COMMON nrhs
  
  nrhs = nrhs+1
  
  a    = exp(x)
  H2   = H2_func(a)
  cst  = 3./2.*H02*Om0/a**3

  ! Define derivatives
  !dydx(1) = (cst*fk*y(2)+(cst-2*H2)*y(1))/H2 
  dydx(1) = (cst*y(2)+(cst-2*H2)*y(1))/H2 
  dydx(2) =  y(1)

  !Write(*,*) 'T ',a,cst,fk,(cst-2*H2),y(1),y(2)

End Subroutine derivs   
!====================

!==================
Function H2_func(a)
!==================

  Use nrtype
  Use cosmo_tools

  Implicit None
  
  Real(DP)  :: H2_func
  Real(DP), Intent(IN) :: a

  !!H2_func = H02*(Om0/a**3+Oml)
  !H2_func = H02*(Om0/a**3+Omk/a**2+Oml*Exp(3.*((1.+w0+wa)*Log(1./a)-wa*(1-a))))
  H2_func = H02*(Om0/a**3+Omr/a**4+Omk/a**2+Oml*Exp(3.*((1.+w0+wa)*Log(1./a)-wa*(1-a))))
  !Write(*,*) 'H2 ',Oml,Omk,Om0
  !Write(*,*) 'H2 ',Om0,Omr,Omk,Oml,w0,wa,Om0+Omr+Omk+Oml,sqrt(H02)

End Function H2_func
!===================




