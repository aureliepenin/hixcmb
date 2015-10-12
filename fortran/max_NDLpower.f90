! By Xiaomin Wang & Max Tegmark, August 2003
! Computes P(k) rapidly using the fitting formulas from
! astro-ph/9811262:
! Title: An analytic approximation of MDM power spectra in four dimensional   parameter space
! Authors:  B. Novosyadlyj (L'viv Univ.),  R. Durrer (Geneva Univ.),  V.N. Lukash (AstroSpaceCenter, Moscow)
! Also makes the separability approximation of eq (C3) in astro-ph/0008167.

!       program test ! Comment out when using subroutines with another program
!       call NDLdemo
!       end
	
subroutine NDLdemo
  implicit none
  integer i, kbands
  parameter(kbands=1001)
  real*8 p(8), kvec(kbands), power(kbands), t,w
  logical verdict,depert
  p(1)    = 0.000 ! Ok
  p(2)    = 0.700 ! Ol
  p(3)    = 0.120 ! od
  p(4)    = 0.024 ! ob
  p(5)    = 0.000 ! fn
  p(6)    = 1.000 ! ns
  p(7)    = 2.2*10e-10 !0.650 ! As
  p(8)    = 0.000 ! al
  do i = 1,kbands
     ! lgk ranges from -5 to 5
     t = (i-1.)/(kbands-1.)			! Range = [0,1]
     t = 10.*t - 5.				! Range = [-5,5]
     kvec(i) = 10**t
  end do
  call compute_NDL_power(p,kbands,kvec,power,verdict)
  open(3,file='max_NDL_power_testout.dat')
  do i=1,kbands
     write(3,'(1f25.15,1f65.25)') kvec(i), power(i)
  end do
  close(3)
  return
end subroutine NDLdemo

subroutine compute_NDL_power(p,bands,kvec,power,verdict)
  ! Uses the Novosyadlyj et al. fitting formulas to compute P(k) rapidly.
  implicit none
  integer bands, i
  real*8 p(8), kvec(bands), power(bands)
  logical verdict
  
  !open(2,file='qaz_NDL.dat')
  do i=1,bands
     ! Notice input and output k are different
     call NDL_computeP(p,kvec(i),power(i),verdict)  	  
     ! To prevent overflow:
     ! if (power(i).gt.99999.) power(i) = 99999.
     ! print *, kvec(i),power(i)
     !write(2,'(2e18.8)') kvec(i), power(i)
  end do
  
  !close(2)
  return
end subroutine compute_NDL_power

subroutine NDL_computeP(p,k,power,verdict)
  ! Uses the Novosyadlyj et al. fitting formulas to compute P(k) rapidly.
  ! k is in h/Mpc, *not* 1/Mpc.
  ! 1   Ok = Omega_k
  ! 2   Ol = Omega_Lambda
  ! 3   od = Omega_dark h^2
  ! 4   ob = Omega_baryon h^2
  ! 5   fn = Omega_nu/Omega_dark, Omega_dark = Omega_cdm + Omega_nu
  ! 6   ns = scalar spectral index
  ! 7   As = primordial power spectrum amplitute
  ! 8   al = dn/dlnk    
  implicit none
  integer i
  real*8 T0
  parameter(T0=2.725d6) ! CMB temperature in uK
  real*8 k, kk, p(8), lastp(8), h, T_cmb, omhh, f_baryon
  real*8  Ok, Ol, Om, od, ob, fn, ns, As, al, Tnorm, cmbfastnorm
  real*8 transfer, tk, power, T, pi
  real*8 chub,omgb,omgh,annu,omgc,omega,zsp,zr
  real*8 hreal
  common /b1/ chub,omgb,omgh,annu,omgc,omega,zsp,zr
  
  logical initialized, died,verdict	
  save lastp, died	
  data lastp/666, 666, 666, 666, 666, 666, 666, 666/

  !p = (/0.000000000000000E+000,0.758732616901398,9.986562292960506E-002, &
  !     1.944022519335900E-002,0.000000000000000E+000,0.909507989883423, &
  !     2.205733900909667E-009,0.000000000000000E+000 /)
  !WritE(*,*) 'max ',p
  
  pi = 4.d0*atan(1.d0)
  ! Input Cosmological Parameters
  Ok	= 0.    ! Make Novos power work, then shifts k later
  Ol	= 0.    ! Make Novos power work, then shifts k later
  od	= p(3)
  ob	= p(4)
  fn	= p(5)
  ns	= p(6)
  As	= p(7)
  al	= p(8)
  if (fn.lt.1.d-3) fn=0   ! NDL can't handle tiny positive fn-values, giving NAN
  ! Atmospheric neutrino oscillations tell us that fn>0.01 or so.
  ! But there's an if-statement for the case fn=0!!!
  Om	= 1.-Ol-Ok
  omhh	= od + ob
  f_baryon= ob/omhh
  h	= sqrt(omhh/(1.-Ok-Ol))
  hreal	= sqrt(omhh/(1.-p(1)-p(2)))
  T_cmb	= T0/1d6	! CMB temperature in K
  ! parameters' names in NDL
  chub	= h
  zsp	= 0. ! redshift for which calculation of MDM transfer function will be done
  zr	= 1000.
  annu	= 1. ! number density of massive neutrino flavors
  omega	= Om
  omgb	= ob/h**2.	! Omega_b
  omgh	= fn*od/h**2.	! Omega_n
  omgc	= Om - omgb 	! Omega_c
  
  initialized = .true.

  !MARIE
  verdict=.true.

  do i=1,5 ! Need to recompute transfer function if any of 1st five parameters changed
     if (p(i).ne.lastp(i)) initialized = .false.
  end do
  if (.not.initialized) then
     do i=1,8
        lastp(i) = p(i)
     end do
     ! print *,'NDL Oc, Ob, On, h =',omgc,omgb,omgh,h
     call appcoef(omhh,ob)
     verdict = .true.
     ! Test if we're outside the range where the fits are valid:
     if (fn.gt.0.5) then   ! fit suitable only for Omega_nu<0.5
        write(*,*)"NDL death warning: Neutrino fraction greater than 0.5."
        verdict = .false.
     else if (omgh+omgb.gt.0.73) then ! fit suitable only for Omega_nu+Omega_b<0.73
        write(*,*)"NDL death warning: Neutrino + baryon density  greater than 0.75."
        verdict = .false.
     end if
  endif
  if (.not.verdict) then ! Avoid printing the same error messages thousands of times for the same model
!     WritE(*,*) 'Here verdict',verdict
     power = 0
     return
  end if
  
  ! NDL work in 1/Mpc, not h/Mpc
  kk = hreal*k		! kk has units of 1/Mpc, k has unit of h/Mpc
  call pksub(kk,transfer)
  if (.not.((transfer.ge.0).or.(transfer.le.0))) then
     write(*,*)"NDL death error: T(k) = NaN"
     write(*,*)p
     write(*,*)k,kk,transfer
     transfer = 666e15
     verdict = .false.
     return
     !	  stop
  end if
  
  Ok = p(1)
  Ol = p(2)
  Om = 1 - Ol - Ok
  call ComputeT0(hreal,Om,Ol,Tnorm) ! Could be moved to initialization part, since independent of k
  T = Tnorm*transfer
  cmbfastnorm = 800*(pi*5./3.)**2/T0**2 ! Converts to our old normalization As=1
  ! The next three lines are the same as in cmbfast_compute_P.f: (writing kk in place of k)
  power = 2*pi**2*(kk/0.05)**(ns-1+al*log(kk/0.05)/2)*kk*T**2 ! Now power has units of Mpc**3
  power = hreal**3 * power			       ! Now power has units of (Mpc/h)**3, normalized to A=1
  !Stop
!!$c rachel - normalize on output
!!$c	power = cmbfastnorm * As * power	
  power = power *As  ! As = my As *clnorm
  !WritE(*,*) 'power ',k,hreal,As,power,Tnorm,transfer
  
  if (.not.initialized) then
     !write(*,*)'h,Om,Ol: ',h,Om,Ol,norm
     !write(*,*)'k,T,P: ',k,transfer,power
  end if
  return
end subroutine NDL_computeP

subroutine ComputeT0(h,Om,Ol,Tnorm) ! Should improve to depend on w too
  ! Computes Tnorm, the CMBfast transfer function T(k) as k=0.
  implicit none
  real*8 h, Om, Ol, Tnorm, Tstar, growth
  Tstar = 3594770 ! This is empirical cmbfast T(0) for flat model with Om=1, h=1  
  growth = 2.5*Om/(Om**(4./7.)-Ol+(1+Om/2)*(1+Ol/70)) ! The Carroll et all approximation
  Tnorm = Tstar*growth/(Om*h**2)
  return
end subroutine ComputeT0

!	subroutine old_incorrect_ComputeNorm(h,Om,Ol,norm)
!	! Computes the normalization constant norm such that
!	! the normalization of 
!	!   P(k) = norm * (k/(h*0.05))**n * T(k)**2          [with k in h/Mpc, *not* h/Mpc?]  
!	! corresponds to the C_l-normalization of CMBfast.	
!	! Assumes that CMBfast does NOT
!	! call COBEnormalize, so that it normalizes to get psi_{initiall}=1.
!	implicit real*8(a-h,o-z)
!	real*8 pi, T0, Om, Ol, norm, den1, d1, d2norm, tfnorm, h
!	pi = 4.*atan(1.)
        !	T0 = 2.726e6 ! CMB teperature in uK
!	! CMBFAST has power=5.*(k/0.05)**ns
!	! Fitting formula for the growth rate
!	! THIS IS THE CARROLL ET ALL APPROXIMATION:
!	den1=Om**(4./7.)-Ol+(1.+0.5*Om)*(1.+Ol/70.)
!	d1=2.5/den1	! WHY ISN'T THIS d1=2.5*Om/den1?
!	tfnorm=d1*5.99e6/h**2
!	! Why the h^2 factor?
!	d2norm=2.*pi/T0**2
!	norm = 5.*d2norm*tfnorm**2*h**3*(2.*pi)**3
!	return
!	end

!!$********************************************************************************
!!$**************** Here is the tranfer function T_MDM from Novosyadlyj et. 1998************
!!$****************     astro-ph/9811262   **************
!!$********************************************************************************
!!$**      subroutine transfer_mdm
!!$**      implicit double precision (a-h,o-z) 
!!$**      common /b1/chub,omgb,omgh,annu,omgc,omega,zsp,zr 
!!$**      common /b2/alc,bc,sh,bnode,alb,betb,akj,aeq,bb,aks,bps,akeq,t0,
!!$**     #a1,a2,a3,a4,a5
!!$**
!!$**c Input data are following:
!!$**
!!$**c redshift 'zsp' for which calculation of MDM transfer function will be done
!!$**c 0\le zsp \le 30
!!$**!      zsp=0.0
!!$**c Hubble parameter 'chub'
!!$**c 30<chub<70
!!$**! It seems that chub can be greater than 70, but no less than 30 in nonflat
!!$**! model.
!!$**!      chub=0.5
!!$**c baryon density parameter 'omgb'
!!$**c omegab \le 0.3
!!$**!      omgb=0.2      
!!$**c hot dark matter density parameter 'omgh'
!!$**c omgh \le 0.5 but omgh+omegab \le 0.73
!!$**! It seems that omgh can not be exactly 0
!!$**!      omgh=0.2
!!$**c number density of massive neutrino flavors 'annu'
!!$**c annu=1, 2 or 3 for omgh \ne 0
!!$**      ! annu=1.0
!!$**      zsp=0.0
!!$**      chub=0.5
!!$**      omgb=0.06
!!$**      ! omgh=0.000001
!!$**      omgh=0.2
!!$**      annu=1.
!!$**      omega=1.0
!!$**          
!!$**c code works for matter dominated flat model only!, omega=1.
!!$**!      omega=1.
!!$**      omgc=omega-omgb 
!!$**      omh=omega*chub*chub 
!!$**      ombh=omgb*chub*chub 
!!$**      zr=1000. 
!!$**  
!!$**      call appcoef(omh,ombh) 
!!$**    
!!$**      open(2,file='tk.dat',status='unknown')    
!!$**      akmin=1.0D-04
!!$**      akmax=25. 
!!$**      ! nk=1000.
!!$**      nk=10000.
!!$**      ank=nk
!!$**      dkl=(dlog10(akmax)-dlog10(akmin))/ank 
!!$**      do 11 i1=1,nk+1  
!!$**      ak=10.**(dlog10(akmin)+dkl*(i1-1))
!!$**      akh=ak/chub 
!!$**      call pksub(ak,tk)  
!!$**      ! write (2,10) akh,tk    ! akh in h Mpc^(-1)
!!$**      write (2,10) ak,akh,tk
!!$**  10  format (6e12.4) 
!!$**  11  continue 
!!$**      close(2) 
!!$**      return
!!$**      end 

subroutine pksub(ak,tk) 
  implicit double precision (a-h,o-z) 
  common /b1/chub,omgb,omgh,annu,omgc,omega,zsp,zr 
  common /b2/alc,bc,sh,bnode,alb,betb,akj,aeq,bb,aks,bps,akeq,t0,a1,a2,a3,a4,a5
  
  qq=ak*(t0/2.7)**2/chub/chub/omega 
  cc=14.2/alc+386./(1.+69.9*qq**1.08) 
  cc1=14.2+386./(1.+69.9*qq**1.08) 
  t0kab=dlog(2.71828+1.8*bc*qq)/(dlog(2.71828+1.8*bc*qq)+cc*qq*qq) 
  t0k1b=dlog(2.71828+1.8*bc*qq)/(dlog(2.71828+1.8*bc*qq)+cc1*qq*qq) 
  t0k11=dlog(2.71828+1.8*qq)/(dlog(2.71828+1.8*qq)+cc1*qq*qq) 
  ff=1./(1.+(ak*sh/5.4)**4) 
  tkc=ff*t0k1b+(1.-ff)*t0kab 
  stk=sh/(1.+(bnode/ak/sh)**3)**0.333333 
  xxx=ak*stk  
  bes0=dsin(xxx)/xxx 
  akxs=(ak/aks)**1.4 
  albe=0. 
  if (akxs.gt.30.) go to 6550 
  albe=alb/(1.+(betb/ak/sh)**3)/dexp(akxs) 
6550 tkb=(t0k11/(1.+(ak*sh/5.2)**2)+albe)*bes0 
  tk1=omgb/omega*tkb+omgc/omega*tkc 
  tk = tk1 ! Xiaomin: surely this line was missing before?
  if (omgh.eq.0.0D0) return 
6545 du=(a1*ak/akj+a2*(ak/akj)**1.5+ a3*(ak/akj)**2+a4*(ak/akj)**2.5+a5*(ak/akj)**3.0) 
  
  amn=(1.-omgh*(1.+2.*(omgb-0.06)))**(1./bps)*aeq*(1.+zsp)* (bb*akj)**3.0 
  if (.not.((amn.ge.0).or.(amn.le.0))) then
     !	write(*,*)"NDL death error: amn = NaN"
     !	write(*,*)k, omghm omgb, amn, bps, aeq, zsp
     stop
  end if
  ch=1.+amn*du 
  zn=1.+(bb*ak)**3.0 
  ! ch can go slightly negative for small values of fn, giving NaN, in which 
  ! case we simply make the approximation that fn=0 and that there's no power suppression at all.
  if (ch.lt.0) return
  tkps=(ch/zn)**bps 
  if (.not.((tkps.ge.0).or.(tkps.le.0))) then
     !	write(*,*)"NDL death error: tkps = NaN"
     write(*,*)k, omghm, omgb, amn, bps, aeq, zsp
     stop
  end if
  tk=tk1*tkps 
  return 
end subroutine pksub

subroutine appcoef(omh,ombh) 
  implicit double precision (a-h,o-z)       
  common /b1/chub,omgb,omgh,annu,omgc,omega,zsp,zr 
  common /b2/alc,bc,sh,bnode,alb,betb,akj,aeq,bb,aks,bps,akeq,t0,a1,a2,a3,a4,a5 

  tg=0.9967 
  t0=2.735*tg 
  aeq=3.5D-05*(t0/2.7)**4/omh 
  zeq=1./aeq-1 
  akeq=dsqrt(2.*omh*zeq)/3000.
  rnr=2.7/chub/chub*tg*tg*sqrt(annu) 
  bb=rnr*(1.+0.1435)/(omgh+0.1435) 
  bps=(5.-dsqrt(25.-24.*omgh))/4. 
  akj=74.5*omgh*chub**3/dsqrt(1.+zsp)/annu 

  b1=0.313/omh**0.419*(1.+0.607*omh**0.674) 
  b2=0.238*omh**0.223 
  zd=1291.*omh**0.251/(1.+0.659*omh**0.828)*(1.+b1*ombh**b2) 
  req=31.5*ombh/(t0/2.7)**4*(1000./zeq) 
  rbd=31.5*ombh/(t0/2.7)**4*(1000./zd) 
  sh=2./3./akeq*dsqrt(6./req)*dlog((dsqrt(1.+rbd)+dsqrt(rbd+req))/(1.+dsqrt(req))) 
  
  aks=1.6*ombh**0.52*omh**0.73*(1.+(10.4*omh)**(-0.95)) 
  a1=(46.9*omh)**0.670*(1.+(32.1*omh)**(-0.532)) 
  a2=(12.0*omh)**0.424*(1.+(45.0*omh)**(-0.582)) 
  alc=1./a1**(omgb/omega)/a2**((omgb/omega)**3) 
  b1=0.944/(1.+(458.*omh)**(-0.708)) 
  b2=1./(0.395*omh)**0.0266 
  bc=1./(1.+b1*((omgc/omega)**b2-1.)) 
  yy=(1.+zeq)/(1.+zd) 
  gg=yy*(-6.*dsqrt(1.+yy)+(2.+3.*yy)*dlog((dsqrt(1.+yy)+1.)/(dsqrt(1.+yy)-1.))) 
  alb=2.07*akeq*sh/(1.+rbd)**0.75*gg 
  betb=0.5+omgb/omega+(3.-2.*omgb/omega)*dsqrt((17.2*omh)**2+1.) 
  bnode=8.41*omh**0.435 
  
  if (omgh.eq.0.) go to 8787 
  
  cnu1=1. 
  cnu2=1. 
  cnu3=1. 
  cnu4=1. 
  cnu5=1. 
  if (annu.eq.1.0) go to 139 
  if (annu.eq.3.0) go to 138 
  cnu1=(2.67-4.986*omgh/(1.+zsp)**0.066)* (chub/0.5)**(0.49+0.0137*zsp) 
  cnu2=(3.39-8.656*omgh/(1.+zsp)**0.19)* (chub/0.5)**(0.84+0.04*zsp) 
  cnu3=(2.91*(1.+zsp)**0.052-5.904*omgh/(1.+zsp)**0.061)* (chub/0.5)**(1.17+0.023*zsp) 
  cnu4=(3.04*(1.+zsp)**0.18-6.404*omgh*(1.+zsp)**0.15)* (chub/0.5)**(1.65+0.025*zsp) 
  cnu5=(2.61*(1.+zsp)**0.0135-4.148*omgh/(1.+zsp)**0.039)* (chub/0.5)**(0.565-0.03*zsp) 
  go to 139 
138 cnu1=(5.-12.06*omgh/(1.+zsp)**0.046)* (chub/0.5)**(0.41+0.017*zsp) 
  cnu2=(6.87*(1.+zsp)**0.044-20.51*omgh/(1.+zsp)**0.083)* (chub/0.5)**(0.71+0.041*zsp) 
  cnu3=(5.92*(1.+zsp)**0.13-14.70*omgh*(1.+zsp)**0.075)* (chub/0.5)**(1.03+0.02*zsp) 
  cnu4=(6.82*(1.+zsp)**0.29-17.22*omgh*(1.+zsp)**0.27)* (chub/0.5)**(1.41+0.008*zsp) 
  cnu5=(5.47*(1.+zsp)**0.023-11.4*omgh)* (chub/0.5)**(0.4-0.021*zsp) 
  
139 bo=(omgb-0.06) 
  ho=omgh*(1.+2.*bo) 
  corh=(1.0+0.15*(chub-0.5)*(omgh/0.2)**0.73* (omgb/0.06)**0.58)**(1./bps) 
  bps=(5.-dsqrt(25.-24.*omgh))/4. 
  bet1=0.0184629-0.72025*omgh+4.77036*omgh**2-10.9186*omgh**3 +11.62891*omgh**4-1.98952*omgh**5 
  
  be1=0.0647623-0.966035*omgh+5.52736*omgh**2-12.0666*omgh**3 +11.24344*omgh**4-0.488735*omgh**5 
  
  yz2=11.**be1 
  yz3=21.**bet1 
  az=1.155-0.21*yz2+0.055*yz3 
  bz=0.22*yz2-0.06*yz3-0.16 
  cz=0.005*(yz3+1.)-0.01*yz2 
  czsp=az+bz*(1.+zsp)+cz*(1.+zsp)**2 
  a1=(14.33879-110.752*ho+455.7737*ho**2-1053.96*ho**3+1271.372 *ho**4-627.584*ho**5)*czsp/corh*cnu1 
  
  bet2=-1.34913+24.69214*omgh-216.684*omgh**2+892.1813*omgh**3 -1696.36*omgh**4+1230.299*omgh**5 
  
  be2=-1.6563+34.34481*omgh-297.198*omgh**2+1196.067*omgh**3 -2243.6*omgh**4+1612.738*omgh**5 
  yz2=11.**be2 
  yz3=21.**bet2 
  az=1.155-0.21*yz2+0.055*yz3 
  bz=0.22*yz2-0.06*yz3-0.16 
  cz=0.005*(yz3+1.)-0.01*yz2 
  czsp=az+bz*(1.+zsp)+cz*(1.+zsp)**2 
  
  a2=(-24.2273+158.3797*ho-509.232*ho**2+912.4319*ho**3-811.954 *ho**4+280.4024*ho**5)*czsp/corh*cnu2 
  
  
  bet3=0.1657301-8.43583*omgh+47.33781*omgh**2-104.57*omgh**3 +102.2927*omgh**4-25.6059*omgh**5 
  
  be3=0.0937063-6.1168*omgh+33.33297*omgh**2-68.411*omgh**3 +57.25015*omgh**4-1.17787*omgh**5 
  yz2=11.**be3 
  yz3=21.**bet3 
  az=1.155-0.21*yz2+0.055*yz3 
  bz=0.22*yz2-0.06*yz3-0.16 
  cz=0.005*(yz3+1.)-0.01*yz2 
  czsp=az+bz*(1.+zsp)+cz*(1.+zsp)**2 
  a3=(31.14532-220.676*ho+868.4743*ho**2-1971.05*ho**3+2336.014* ho**4-1140.03*ho**5)*czsp/corh*cnu3 
  
  bet4=0.6073224-17.5771*omgh+84.96842*omgh**2-156.886*omgh**3 +106.8634*omgh**4 
  
  be4=0.4335517-14.5545*omgh+80.9298*omgh**2-196.306*omgh**3 +232.8642*omgh**4-102.861*omgh**5 
  yz2=11.**be4 
  yz3=21.**bet4 
  az=1.155-0.21*yz2+0.055*yz3 
  bz=0.22*yz2-0.06*yz3-0.16 
  cz=0.005*(yz3+1.)-0.01*yz2 
  czsp=az+bz*(1.+zsp)+cz*(1.+zsp)**2 
  a4=(-12.1073+88.2023*ho-389.51*ho**2+986.4919*ho**3-1290.4*ho**4 +687.0816*ho**5)*czsp/corh*cnu4 
  
  corh=1.0 
  bet5=0.0193762-1.88788*omgh+7.51947*omgh**2-11.4426*omgh**3 +7.33912*omgh**4 
  
  be5=0.0281617-2.1184*omgh+8.578*omgh**2-13.2839*omgh**3 +8.77538*omgh**4 
  yz2=11.**be5 
  yz3=21.**bet5 
  az=1.155-0.21*yz2+0.055*yz3 
  bz=0.22*yz2-0.06*yz3-0.16 
  cz=0.005*(yz3+1.)-0.01*yz2 
  czsp=az+bz*(1.+zsp)+cz*(1.+zsp)**2 
  a5=(5.89136-53.3049*ho+266.297*ho**2-701.69*ho**3+933.0745*ho**4 -495.85*ho**5)*czsp/corh*cnu5 
8787 continue 
  return 
end subroutine appcoef


