Module eh99_tools
  
  Real*8 :: omhh,f_nu,f_baryon,N_nu,y_d,alpha_nu,beta_c,sound_horizon,theta_cmb,omega,omegal,z_equality
  Real*8 :: scale,tilt,zeh99,ipower

Contains

  ! =========================
  Subroutine TFset_parameters
  ! =========================
      
      Implicit Real*8 (a-h,k,o-z)
      
      ! === Auxiliary variable
      obhh = omhh*f_baryon
      
      ! ===  Main variables
      z_equality = 2.50e4*omhh*theta_cmb**(-4.) - 1.D0
      k_equality = 0.0746*omhh*theta_cmb**(-2.)
      
      z_drag = 0.313*omhh**(-0.419)*(1.+0.607*omhh**(0.674))
      z_drag = 1e0 + z_drag*obhh**(0.238*omhh**(0.223))
      z_drag = 1291e0 * omhh**(0.251)/(1e0 + 0.659*omhh**(0.828)) * z_drag
      
      y_d = (1.+z_equality)/(1.+z_drag)
      

      R_drag = 31.5*obhh*theta_cmb**(-4.)*1000e0/(1e0 + z_drag)
      R_equality = 31.5*obhh*theta_cmb**(-4.)*1000e0/(1e0 + z_equality)
      
      sound_horizon = 2./3./k_equality*Sqrt(6./R_equality)* &
     &     Log(( Sqrt(1.+R_drag)+Sqrt(R_drag+R_equality) )  &
     &     /(1.+Sqrt(R_equality)))
      
      p_c  = -(5.-Sqrt(1.+24*(1.-f_nu-f_baryon)))/4.
      p_cb = -(5.-Sqrt(1.+24*(1.-f_nu)))/4.
      f_c  = 1.-f_nu-f_baryon
      f_cb = 1.-f_nu
      f_nub= f_nu+f_baryon
      
      alpha_nu = (f_c/f_cb)* (2.*(p_c+p_cb)+5.)/(4.*p_cb+5)
      alpha_nu = alpha_nu*(1.-0.553*f_nub+0.126*f_nub**3)
      alpha_nu = alpha_nu/(1.-0.193*Sqrt(f_nu)+0.169*f_nu)
      alpha_nu = alpha_nu*(1.+y_d)**(p_c-p_cb)
      alpha_nu = alpha_nu*(1.+ (p_cb-p_c)/2.*(1.+1./(4.*p_c+3.)/(4.*p_cb+7.))/(1.+y_d))
      beta_c   = 1./(1.-0.949*f_nub)

      !Write(*,*) 'cs ',sound_horizon
      Return
      
    End Subroutine TFset_parameters
    ! =============================

    ! ==========================
    Real*8 Function TF_master(k)
    ! ==========================
      
      Implicit Real*8 (a-h,k,o-z)
        
      q = k*theta_cmb**2/omhh
      gamma_eff=(Sqrt(alpha_nu) + (1.-Sqrt(alpha_nu))/(1.+(0.43*k*sound_horizon)**4))
      
      q_eff = q/gamma_eff
      TF_master= dlog(dexp(1.d0)+1.84*beta_c*Sqrt(alpha_nu)*q_eff)
      TF_master = TF_master/(TF_master + q_eff**2*(14.4 + 325./(1.+60.5*q_eff**1.11)))
      
      q_nu = 3.92*q*Sqrt(N_nu/f_nu)
      TF_master = TF_master*(1.+(1.2*f_nu**(0.64)*N_nu**(0.3+0.6*f_nu))/(q_nu**(-1.6)+q_nu**(0.8)))
      
      Return 
      
    End Function TF_master
    ! ====================
      
!!$    ! ======================================
!!$    Subroutine Growth(z,k,DD_cb,DD_cbnu,DD0)
!!$    ! ======================================
!!$      
!!$      Implicit Real*8 (a-h,k,o-z)
!!$      
!!$      q = k*theta_cmb**2/omhh
!!$      
!!$      y_fs = 17.2*f_nu*(1.+0.488*f_nu**(-7./6.))*(N_nu*q/f_nu)**2
!!$      
!!$      oz=omega*(1.+z)**3/(omegal+(1.-omegal-omega)*(1.+z)**2+omega*(1.+z)**3)
!!$      olz=omegal/(omegal+(1.-omegal-omega)*(1.+z)**2+omega*(1.+z)**3)
!!$      
!!$      D = (1.+z_equality)/(1.+z)*5.*oz/2.*(oz**(4./7.)-olz+(1.+oz/2.)*(1.+olz/70.))**(-1.)
!!$      
!!$      DD0= D/((1.+z_equality)*5.*omega/2.*(omega**(4./7.)-omegal+(1.+omega/2.)*(1.+omegal/70.))**(-1.))
!!$      
!!$      p_cb = -(5.-Sqrt(1.+24*(1.-f_nu)))/4.
!!$      
!!$      DD_cb = (1.+(D/(1.+y_fs))**(0.7))**(-p_cb/0.7)*D**(p_cb)
!!$      
!!$      DD_cbnu = ((1.-f_nu)**(-0.7/p_cb)+(D/(1.+y_fs))**(0.7))**(-p_cb/0.7)*D**(p_cb)
!!$      
!!$      Return
!!$
!!$    End Subroutine Growth
!!$    ! ===================
 
    !! =========================
    Real*8 Function sigmatop(kl)
    !! =========================

      ! Redefinition consistent with P_dd_ln

    Implicit None
    
    Real*8 :: kl,x,k,DoD0
    Real*8 :: hub_tmp,nt,dh,Tm,P_dd,c_over_H0

    k = Exp(kl)
    x = scale*k

    Tm = TF_master(k)
    !P_dd      = k**(3.+tilt)*Tm**2.*2.*3.14159**2/k**3
    !sigmatop  = P_dd*k**2/(2.*3.14159**2)*(3.*(x*dcos(x) - dsin(x))/x**3)**2*k
    sigmatop  = Tm**2*k**(3.+tilt)*(3.*(x*dcos(x) - dsin(x))/x**3)**2
    
!!$    sigmatop = k**(3.+tilt)*(DD0*DD_cbnu*TF_master(k))**2*(3.*(x*dcos(x) - dsin(x))/x**3)**2
!!$ Else
!!$    sigmatop = k**(3.+tilt)*(DD0*DD_cb*TF_master(k))**2*(3.*(x*dcos(x) - dsin(x))/x**3)**2
!!$ End If
    
    !Write(*,*) 's ',k,sigmatop

    Return
    
  End Function sigmatop
  !! ==================

!!$    ! ==========================
!!$    Real*8 Function sigmatop(kl)
!!$    ! ==========================
!!$      
!!$      Implicit Real*8 (a-h,k,o-z)
!!$      
!!$      !External Growth    
!!$      
!!$      k = Exp(kl)
!!$      x = scale*k
!!$      
!!$      Call Growth(zeh99,k,DD_cb,DD_cbnu,DD0)
!!$      
!!$      If (ipower.Eq.0) Then
!!$         sigmatop = k**(3.+tilt)*(DD0*DD_cbnu*TF_master(k))**2*(3.*(x*dcos(x) - dsin(x))/x**3)**2
!!$      Else
!!$         sigmatop = k**(3.+tilt)*(DD0*DD_cb*TF_master(k))**2*(3.*(x*dcos(x) - dsin(x))/x**3)**2
!!$      End If
!!$      
!!$      Return
!!$    End Function sigmatop
!!$    ! ===================
     
  End Module Eh99_tools
  !====================
