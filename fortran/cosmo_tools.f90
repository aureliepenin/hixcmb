Module cosmo_tools

  Use nrtype
  Use eh99_tools


  Implicit None

  Real(DP)       :: H0, H02, Oml, Omk, Om0, hub, N_eff
  Real(DP)       :: wg, T_cmb, wr, Omr, anorm, dr2, ldr2
  Real(DP)       :: w0, wa, wm, wb, ns, sig8, omegab, omeganu
  Real(DP)       :: mucc, ff, cc, ncc, cs, ncs, Mcut

  !! HOD model parameters
  Real(DP):: alpha_0, gamma_0, beta_0, T0_0, delta_0, log10_Mmin_0, log10_Meff_0
  Real(DP):: alpha, gamma, beta, T0, delta, log10_Mmin, log10_Meff

  !!Constants
  Real(DP), Parameter :: invhub  = 2997.9   ! c/H_0 in h^(-1) Mpc
  Real(DP), Parameter :: z_star  = 1090.51  ! z_\star decoupling
  Real(DP), Parameter :: th_star = 0.0104004! 0.5959* ! z_\star for WMAP5 (in radian)
  Real(DP), Parameter :: c_sol   = 299792.458 ! km/s

  !! Parameter vector
!  Integer , Parameter         :: npar = 13

  !! Conftau
  Integer, Parameter :: nz_tau = 1000
  Real(DP)  :: z_tau_min = 1.d-5
  Real(DP)  :: z_tau_max = 1.d4
  Real(DP)  :: conf_lnz_arr(1:nz_tau),conf_tau_arr(1:nz_tau),conf_tau2_arr(1:nz_tau)

  !! DoD0 related
  Integer, Parameter :: nd0 = 1000
  Real(DP) :: z_min_D0 = 1.d-8
  Real(DP)  :: z_max_D0 = 1.d4
  Real(DP)  :: dod0_lnz_arr(1:nd0),dod0_D_arr(1:nd0),dod0_D2_arr(1:nd0)

  !! Power spec related
  Real(DP)  :: Mstar



  Public 

Contains

  !==================================================
  Subroutine set_default_cosmology(iparam0)
!  Subroutine set_default_cosmology(iparam0,param_vec)
    !==================================================
    !! iparam : 0:reset all cosmology to default except the one selected by iparam
    !!           which has to be defined outside independantly
    !! values of the parameters cf http://lambda.gsfc.nasa.gov/product/map/dr4/params/lcdm_sz_lens_wmap7.cfm
    !! Normalized at sigma_8


    Integer                        , Optional :: iparam0
 !!   Real(SP)    , Dimension(1:npar), Optional :: param_vec

    Real(DP)  :: sigma_8,rombint,tol_i,tmp
    Integer   :: iparam

    External :: rombint

   If (.Not. Present(iparam0)) Then 
     iparam = 0
  Else
      iparam = iparam0
  End If

   !! Define cosmological parameters
    T_cmb = 2.725
    N_nu  = 0.0
    N_eff = 3.04
    wg    = 2.469d-5 ! \omega_\gamma h^2 for Tcmb = 2.725K as in Komatsu et al. arXiv:0803.0547 Eq. 7 

    !! Planck
    If (iparam .Ne. 1) w0   = -1.0
    If (iparam .Ne. 2) wa   =  0. 
    If (iparam .Ne. 3) wm   =  0.1430!0.1334
    If (iparam .Ne. 4) wb   =  0.022068!0.0258
    If (iparam .Ne. 5) ns   =  0.9624!0.963
    If (iparam .Ne. 6) sig8 =  0.8344!0.801 !



    hub = 0.6711

    Oml  = 0.6825!0.7259! 0.734             ! Omega_Lambda
    Omk  = 0.                ! Omega_curvature
!    Om0  = 1.- Oml - Omk     ! Omega_matter
    Om0 = wm/hub**2
    omega = Om0
    omegal= Oml
!    hub   =  Sqrt(wm/omega)

    omegab  = wb/hub**2      ! Omega_baryons
    omeganu = 0.0            ! Omega_neutrinos
    H0  = hub*100.
    H02 = H0**2
    Omr = wg/hub**2*(1.+0.2271*N_eff) 

    mucc = 0.9
    ff   = 0.1
    cc   = 7.5
    ncc  =-0.3
    cs   = 9.
    ncs  =-0.13
    Mcut = 1.d13

    tilt    = ns

    !! ======================
    !! EH99 stuff, Don't edit
    !! ======================
    !! Translate Parameters into forms GLOBALVARIABLES form
    f_nu     = omeganu/omega
    f_baryon = omegab/omega
    If (f_nu.Eq.0) f_nu = 1E-10
    If (f_baryon.Eq.0) f_baryon = 1E-10
    If (T_cmb.Le.0) T_cmb = 2.726
    If (N_nu.Lt.0)  N_nu=0.
    If (hub.Gt.10)   Write(6,*) 'WARNING: H_0/100km/s/Mpc needed'
    theta_cmb = T_cmb/2.7
    omhh      = omega*hub**2
    !Write(*,*) 'omhh ',omhh, 'wm ',wm,omega,hub,tilt,ns
    !Write(*,*) 'obhh ',omhh*f_baryon, 'wb ',wb
    Call TFset_parameters

    !! Sigma_8: Mass fluctuations in a top hat of 8./h Mpc at
    !! redshift z=0; total mass: CDM+baryon+neutrino
    scale  = 8./hub
    ipower = 0
    tol_i  = 1.d-6
    anorm  = 1.94E-5*omega**(-0.785-0.05*Log(omega))*Exp(-0.95*(tilt-1.)-0.169*(tilt-1.)**2)
    anorm  = anorm**2*(2997.9/hub)**(3.+tilt)

    !!  !! Sigma_8: Mass fluctuations in a top hat of 8./h Mpc at
    !!  !! redshift z=0; total mass: CDM+baryon+neutrino
    !!  scale  = 8./hub
    !!  ipower = 0
    !!  tol_i  = 1.d-6
    !!  anorm  = 1.94E-5*omega**(-0.785-0.05*Log(omega))*Exp(-0.95*(tilt-1.)-0.169*(tilt-1.)**2)
    !!  anorm  = anorm**2*(2997.9/hub)**(3.+tilt)
    !! ==
    !! Use this line if you normalize at \Delta_R^2
    dr2    = Exp(ldr2)
    !!    If (.Not. use_sig8_norm) anorm  = dr2*(0.796/4.396835056303888E-008)**2
    !! This value is chosen so that sigma8 and dr2 matches WMAP5 wmap only values
    !! ==
    !Write(*,*) 'Computing sig8'
    sigma_8= &
         &    Sqrt(anorm*( rombint(sigmatop,dlog(0.001/scale),dlog(0.1d0/scale),tol_i)+ &
         &                 rombint(sigmatop,dlog(0.1d0/scale),dlog(1d0/scale),tol_i)  + &
         &                 rombint(sigmatop,dlog(1d0/scale),dlog(10.d0/scale),tol_i)  + &
         &                 rombint(sigmatop,dlog(10.d0/scale),dlog(100.d0/scale),tol_i)))
    anorm = (sig8/sigma_8)**2*anorm
    !!
    !!    If (.Not. use_sig8_norm) sig8 = sigma_8 !<= Only if As normalized
    !! ==
    !! Use this line if you want to normalize at sig8
    !    If (use_sig8_norm) Then
    !       anorm = (sig8/sigma_8)**2*anorm
    !    End If
    !! ==
    !Write(*,*) 'sigma_8,anorm ',sigma_8,sig8/sigma_8,anorm,anorm**2*1.e-10,dr2*1.e10,anorm**2*1.e-10/(dr2*1.e10),anorm**2/dr2
    !Stop
    !! ==

  End Subroutine set_default_cosmology
  !=====================================

  !! ========================
  Real*8 Function conftime(z)
    !! ========================

    Implicit None
    Real*8 :: z,lnz,tau

    If (z .Lt. z_tau_min) Then
       lnz=Log(z_tau_min)
    Else
       lnz = Log(z)
    Endif
    Call splint(conf_lnz_arr,conf_tau_arr,conf_tau2_arr,nz_tau,lnz,tau)
    conftime = tau

    Return
  End Function conftime
  !! ==================      

  !! ===========================
  Subroutine initiate_conftime()
    !! ===========================

    Implicit None

    Real*8  :: z,tau
    Integer :: iz

    Do iz = 1,nz_tau
       conf_lnz_arr(iz) = Log(z_tau_min*0.1) + (iz-1.)*(Log(z_tau_max)-Log(z_tau_min*0.1))/(nz_tau-1.)
       z = Exp(conf_lnz_arr(iz))
       conf_tau_arr(iz) = conftime_integrate(z)
    End Do

    Call spline(conf_lnz_arr,conf_tau_arr,nz_tau,3.d30,3.d30,conf_tau2_arr)   

  End Subroutine initiate_conftime
  !! =============================

  !! =========================
  Real*8 Function dconftime(z)
    !! =========================

    Implicit None

    Real(DP) :: z

    dconftime  = dconftime_integrate(z)

    Return
  End Function dconftime
  !     =====================================      

  !! =====================================
  Real*8 Function conftime_integrate(z)
    !! =====================================
    ! Flat and non-Flat Universe
    ! See eg astro-ph/9709054 for our conventions

    Implicit None

    Real(DP) :: z,zmin,lzmin,zmax,lzmax,rombint

    External  rombint      

    If (z .Le. z_tau_min) Then
       zmin = z_tau_min
    Else
       zmin  = z
    Endif
    zmax  = 1.d50
    lzmin = Log(zmin)
    lzmax = Log(zmax)

    conftime_integrate = invhub*H0*rombint(conftime_integrate_int_ln,lzmin,lzmax,1.d-5)

    Return

  End Function conftime_integrate
  ! ==============================
  !!
  ! ==============================================
  Real*8 Function conftime_integrate_int_ln(lnz)
    ! ==============================================
    !     ! Flat Universe only

    Implicit None

    !       passed variables:
    Real*8   :: z,a,lnz,H2_func
    External :: H2_func

    z = Exp(lnz)
    a = 1./(1.+z)
    conftime_integrate_int_ln = z/Sqrt(H2_func(a))

    Return
  End Function conftime_integrate_int_ln
  !======================================
  !!
  ! =====================================
  Real*8 Function dconftime_integrate(z)
    ! =====================================
    ! Flat and non-Flat Universe
    ! See eg astro-ph/9709054 for our conventions

    Implicit None

    Real*8   :: z,a,H2_func
    External :: H2_func

    a     = 1./(1.+z)
    dconftime_integrate = invhub*H0/Sqrt(H2_func(a))

    Return
  End Function dconftime_integrate
  ! ===============================

  !!  !! =======================
  Subroutine initiate_dod0()
    !! =======================

    Implicit None

    Real*8  :: D0,z
    Integer :: iz

    D0 = Dlingrow(0.d0)
    Do iz = 1,nD0
       dod0_lnz_arr(iz) = Log(z_min_D0) + (iz-1.)*(Log(z_max_D0)-Log(z_min_D0))/(nD0-1.)
       z                = Exp(dod0_lnz_arr(iz))
       dod0_D_arr(iz)   = Dlingrow(z)/D0
       !Write(*,*) dod0_lnz_arr(iz),dod0_D_arr(iz)
    End Do
    Call spline(dod0_lnz_arr,dod0_D_arr,nD0,3.d30,3.d30,dod0_D2_arr)   
    !Write(*,*) 'D over D0 array computed'

  End Subroutine initiate_dod0
  !! =========================

  !! ==================
  Function fast_DoD0(z)
    !! ==================

    Implicit None

    Real*8 :: z,lnz,fast_DoD0,ratio

    If (z .Gt. z_min_D0) Then
       lnz = Log(z)
       Call splint(dod0_lnz_arr,dod0_D_arr,dod0_D2_arr,nD0,lnz,ratio)
       fast_DoD0 = ratio
    Else
       fast_DoD0 = 1.d0
    Endif

  End Function fast_DoD0
  !! ===================



  !! ========================
  Real*8 Function Dlingrow(z)
    !! ========================

    Implicit None

    Real*8   :: a,z,integrate_D
    External :: integrate_D

    a = 1./(1.+z)
    Dlingrow = integrate_D(a)

  End Function Dlingrow
  !! ==================

  !! =========================
  Real*8 Function rho_bar_m(z)
    !! =========================
    ! Compute the mean matter density rho bar as a function of z
    ! eg Kitayama and Suto 96 (~2.4 & A6) and Peacock (p. 664)
    ! units [ h^2 M_sun Mpc^-3 ]
    ! z ok

    Implicit None

    Real*8 :: z

    rho_bar_m = 2.775d11*(1.+z)**3*omhh 
    ! 2.775d11 [h^2 M_sun/Mpc^3]
    Return

  End Function rho_bar_m
  !! ===================


  !! =======================
  Real*8 Function omega_f(z)
    !! =======================
    !     Compute the evolution if Omega
    !     Kitayama,Suto 1996  (Eq. A7)
    !     Eisenstein & Hu   (Eq. A5)

    Implicit None

    Real*8 :: z,a,num,denom

    !num    = omega*(1.0+z)**3
    !denom  = omega*(1.0+z)**3 + (1.0-omega-omegal)*(1.0+z)**2 + omegal
    !Write(*,*) 'omega ',omega,omegal
    a      = 1./(1.+z)
    num    = omega
    denom  = a+omega*(1.0-a)+omegal*(a**3-a)

    omega_f = num/denom

  End Function omega_f
  !! =================

  !==================
  Function Hubble(z)
    !==================
    
    Implicit None
    
    Real(DP)  :: Hubble
    Real(DP)  :: z, a
    
    a = 1.d0/(1. + z)
    
    Hubble = sqrt(H02*(Om0/a**3+Omr/a**4+Omk/a**2+Oml*Exp(3.*((1.+w0+wa)*Log(1./a)-wa*(1-a)))))
    
  End Function Hubble
!===================








End Module cosmo_tools
