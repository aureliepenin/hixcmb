!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	Real*8 Function spherbesselJD(L,x)
!!$c  Calculates the spherical bessel function j_l(x) 
!!$c  and optionally its derivative for real x and integer l>=0. 
!!$c  Asymptotic approximations 8.11.5, 8.12.5, and 8.42.7 from 
!!$c  G.N.Watson, A Treatise on the Theory of Bessel Functions,
!!$c  2nd Edition (Cambridge University Press, 1944).
!!$c  Higher terms in expansion for x near l given by
!!$c  Airey in Phil. Mag. 31, 520 (1916).
!!$
!!$c  This approximation is accurate to near 0.1% at the boundaries
!!$c  between the asymptotic regions; well away from the boundaries
!!$c  the accuracy is better than 10^{-5}. The derivative accuracy
!!$c  is somewhat worse than the function accuracy but still better
!!$c  than 1%.
!!$
!!$c  Point *jlp initially to a negative value to forego calculating
!!$c  the derivative; point it to a positive value to do the derivative
!!$c  also (Note: give it a definite value before the calculation
!!$c  so it's not pointing at junk.) The derivative calculation requires 
!!$
!!$c  only arithmetic operations, plus evaluation of one sin() for the
!!$c  x>>l region.  
!!$
!!$
!!$c  Original code by Arthur Kosowsky   akosowsky@cfa.harvard.edu
!!$c  This fortran version only computes j_l(x)

	Implicit Double Precision(a-h,o-z)
	Integer L
	Real*8 nu, nu2,ax,ax2,beta,beta2,beta4,beta6
  	Real*8 sx,sx2,cx,sum1,sum2,sum3,sum4,sum5,deriv1
 	Real*8 cotb,cot3b,cot6b,secb,sec2b,sec4b,sec6b
 	Real*8 trigarg,trigcos,expterm,prefactor,llimit,ulimit,fl
	Real*8 x,jl

	PI = 3.1415926536
	ROOTPI = 1.772453851
	GAMMA1 = 2.6789385347             !/* Gamma function of 1/3 */
	GAMMA2 = 1.3541179394             !/* Gamma function of 2/3 */

  	ax = Abs(Dble(x))
 	fl = l

	beta = fl**0.325
 	llimit=1.31*beta	!/* limits of asymptotic regions; fitted */
  	ulimit=1.48*beta
	
	nu= fl + 0.5        

	nu2=nu*nu
	
  	If (l .Lt. 0) Then
	   Print*, 'Bessel function index < 0\n'
	   Stop
	End If
	
!c       /************* Use closed form for l<6 **********/
	
 	If (l .Lt. 6) Then                

	   sx=Sin(ax)
	   cx=Cos(ax)
	   ax2=ax*ax
	   
	   If (l .Eq. 0) Then
	      If(ax .Gt. 0.001) Then
		 jl=Real(sx/ax)
	      Else 
		 jl=Real(1.0d0-ax2/6.0d0)
	      End If		!   /* small x trap */
	   Endif
 

	 If(l .Eq. 1) Then 
	    If(ax .Gt. 0.001) Then
	       jl=Real((sx/ax -cx)/ax)
	    Else 
	       jl=Real(ax/3.0d0)
	    End If
	 Endif

	 If(l .Eq. 2) Then
	    If(ax .Gt. 0.001) Then
	       jl=Real((-3.0d0*cx/ax-sx*(1.0d0-3.0d0/ax2))/ax)
	    Else 
	       jl=Real(ax2/15.0d0)
	    End If
	 Endif
	
	If(l .Eq. 3) Then
	   If(ax .Gt. 0.001) Then
	      jl=Real((cx*(1.0d0-15.0d0/ax2)-sx*(6.0d0-15.0d0/ax2)/ax)/ax)
	   Else 
	      jl=Real(ax*ax2/105.0d0)
	   Endif
	Endif
	
	If(l .Eq. 4) Then
	   If(ax .Gt. 0.001) Then
	      jl=Real((sx*(1.0d0-45.0d0/(ax*ax)+105.0d0/(ax*ax*ax*ax)) +cx*(10.0d0-105.0d0/(ax*ax))/ax)/ax)
	   Else 
	      jl=Real(ax2*ax2/945.0d0)
	   End If
	Endif
	
	If(l .Eq. 5) Then 
	   
	   If(ax .Gt. 0.001) Then
	      jl=Real((sx*(15.0d0-420.0d0/(ax*ax)+945.0d0 &
	&	   /(ax*ax*ax*ax))/ax -cx*(1.0-105.0d0/(ax*ax)+945.0d0 &
	&	   /(ax*ax*ax*ax)))/ax)
	   Else 
	      jl=Real(ax2*ax2*ax/10395.0d0)
	   Endif
	Endif


!c       /********************** x=0 **********************/
	  
  	Else If (ax .Lt. 1.d-30) Then
	   jl=0.0
	   
!c       /*************** Region 1: x << l ****************/
	   
	Else If (ax .Le. fl+0.5-llimit) Then
 

!c       beta=acosh(nu/ax)
	If (nu/ax .Lt. 1.d0) Print*, 'trouble with acosh'
	beta = dlog(nu/ax + dsqrt((nu/ax)**2 - 1.d0) ) 
			!(4.6.21)
	cotb=nu/Sqrt(nu*nu-ax*ax) ! /* cotb=coth(beta) */
    	cot3b=cotb*cotb*cotb
    	cot6b=cot3b*cot3b
        secb=ax/nu
    	sec2b=secb*secb
    	sec4b=sec2b*sec2b
    	sec6b=sec4b*sec2b
    	sum1=2.0+3.0*sec2b
    	expterm=sum1*cot3b/(24.0*nu)
    	sum2=4.0+sec2b
    	expterm = expterm - sum2*sec2b*cot6b/(16.0*nu2)
    	sum3=16.0-1512.0*sec2b-3654.0*sec4b-375.0*sec6b
   	expterm = expterm - sum3*cot3b*cot6b/(5760.0*nu*nu2)
    	sum4=32.0+288.0*sec2b+232.0*sec4b+13.0*sec6b
   	expterm = expterm - sum4*sec2b*cot6b*cot6b/(128.0*nu2*nu2)
    	expterm=Exp(-nu*beta+nu/cotb-expterm)
    	prefactor=Sqrt(cotb/secb)/(2.0*nu)
    	jl=Real(prefactor*expterm)

!c       /**************** Region 2: x >> l ****************/
                          

  	Else If (ax .Ge. fl+0.5+ulimit) Then         


    	beta=Acos(nu/ax)
    	cotb=nu/Sqrt(ax*ax-nu*nu)      !/* cotb=cot(beta) */
    	cot3b=cotb*cotb*cotb
    	cot6b=cot3b*cot3b
    	secb=ax/nu
    	sec2b=secb*secb
    	sec4b=sec2b*sec2b
    	sec6b=sec4b*sec2b
    	trigarg=nu/cotb - nu*beta - PI/4.0
    	sum1=2.0+3.0*sec2b
    	trigarg = trigarg - sum1*cot3b/(24.0*nu)
    	sum3=16.0-1512.0*sec2b-3654.0*sec4b-375.0*sec6b
    	trigarg = trigarg - sum3*cot3b*cot6b/(5760.0*nu*nu2)
    	trigcos=Cos(trigarg)
    	sum2=4.0+sec2b
    	expterm=sum2*sec2b*cot6b/(16.0*nu2)
    	sum4=32.0+288.0*sec2b+232.0*sec4b+13.0*sec6b
    	expterm = expterm - sum4*sec2b*cot6b*cot6b/(128.0*nu2*nu2)
    	expterm=Exp(-expterm)
    	prefactor=Sqrt(cotb/secb)/nu
    	jl=Real(prefactor*expterm*trigcos)

!c       /***************** Region 3: x near l ****************/
	
  	Else                       

    	beta=ax-nu        

    	beta2=beta*beta
    	beta4=beta2*beta2
    	beta6=beta2*beta4
    	sx=6.0/ax
    	sx2=sx*sx
    	cx=Sqrt(sx)                   

    	secb=sx**0.333333333    

    	sec2b=secb*secb

    	deriv1=GAMMA1*secb
    	deriv1= deriv1+ beta*GAMMA2*sec2b
    	sum1=(beta2/6.0-1.0/15.0)*beta
    	deriv1 = deriv1 - sum1*sx*secb*GAMMA1/3.0
    	sum2=beta4/24.0-beta2/24.0+1.0/280.0
    	deriv1 = deriv1 - 2.0*sum2*sx*sec2b*GAMMA2/3.0
    	sum3=beta6/720.0-7.0*beta4/1440.0+beta2/288.0-1.0/3600.0
    	deriv1 = deriv1 + 4.0*sum3*sx2*secb*GAMMA1/9.0
    	sum4=(beta6/5040.0-beta4/900.0+19.0*beta2/12600.0-13.0/31500.0)*beta
    	deriv1 = deriv1 + 10.0*sum4*sx2*sec2b*GAMMA2/9.0
    	sum5=(beta4*beta4/362880.0-beta6/30240.0+71.0*beta4/604800.0 &
     &               -121.0*beta2/907200.0 + 7939.0/232848000.0)*beta
    	deriv1 = deriv1 - 28.0*sum5*sx2*sx*secb*GAMMA1/27.0    

    	jl=Real(deriv1*cx/(12.0*ROOTPI))

	End If

        spherbesselJD=jl
	Return
	End

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 


!!$C////////////////////////////////////////
!!$C
!!$C	Function dbessel
!!$C
!!$C	Calculates the derivative of spherical
!!$C	bessel functions using 
!!$C
!!$C	j_l'(x) = [1/(2l + 1)] * [l*j_(l-1)(x) - (l+1)*j_(l+1)(x)] 
!!$C
!!$C////////////////////////////////////////
 Real*8 Function dbessel(l, x)
   Implicit None

   !C		// argument types
   Integer L
   Real*8 x

   !C		// External functions
   Real*8 spherbesselJD
   External spherbesselJD

   !C		// Local variables
   Real*8 Result

   Result = l * spherbesselJD(l-1, x)
   Result = Result - (l+1) * spherbesselJD(l+1, x)
   Result = Result / (2*l + 1)

   dbessel = Result
   Return

 End Function dbessel


