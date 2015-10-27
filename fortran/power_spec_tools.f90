Module power_spec_tools

  Use nrtype
  Use eh99_tools
  Use cosmo_tools

  Implicit None

  !! Sigma interp
  Integer, Parameter :: n_sig = 500
  Real(DP)  :: M_min_sig = 1.d0!10.d0
  Real(DP)  :: M_max_sig = 1.d22
  Real(DP)  :: sig_arr(n_sig),sig2_arr(n_sig),M_arr(n_sig),sigg_arr(n_sig),sigg2_arr(n_sig)

  !! Pk array 
  Integer, Parameter :: na_arr = 10!8!20!5
  Integer, Parameter :: nk_arr = 100!80!100!50
  Real*8  :: kmin_arr  = 1.d-4
  Real*8  :: kmax_arr  = 1.d+4
!  Real*8  :: zmax_arr  = 7.!1100.d0
  Real*8  :: ku_buf,Mu_buf
  Real*8  :: lp3d_lna_arr(na_arr),lp3d_lnk_arr(nk_arr),lp3d_arr2(na_arr,nk_arr),lp3d_arr(na_arr,nk_arr)
  Real*8  :: lp3d_arr_1g(na_arr,nk_arr),lp3d_arr_2g(na_arr,nk_arr)
  Real*8  :: f_ng_bar_arr(na_arr), f_ng_bar_arr2(na_arr)


  !! Integrands sigma ln_gau
  Real(DP) :: R_sig,z_sig

 !! Global halo
  Real*8  :: M_min,M_max,z_Mstar,M_min0,M_max0,lM_min0,lM_max0
  Real*8  :: tol  = 1.d-3
  Real*8  :: dM   = 1.d3

!! 1h and 2h term
  Integer :: mass_function
!!  Real(DP), parameter :: A_ST=0.322, aa_ST=0.75, p_ST=0.3
  Real(DP), parameter :: A_ST=0.622, aa_ST=0.75, p_ST=0.33
  Real*8  :: k_buf_2h,z_buf_2h
  Real*8  :: k_buf_1h,z_buf_1h
  Real*8  :: k_2h_norm,P_2h_norm,z_2h_norm

  Public 
Contains

  !! =========================================
  Subroutine fast_pk_nl(k,pnl,pnl_1g,pnl_2g,z)
  !! =========================================

    Implicit None

    Real*8 pnl,lpnl,k,z,lnk,lna, lpnl_1g, lpnl_2g,pnl_1g,pnl_2g
    Real*8 splin2,pk
    External splin2
    lnk  =  Log(k)
    lna  = -Log(1.+z)
    lpnl =  splin2(lp3d_lna_arr,lp3d_lnk_arr,na_arr,nk_arr,lp3d_arr,lp3d_arr2,lna,lnk)
    lpnl_1g =  splin2(lp3d_lna_arr,lp3d_lnk_arr,na_arr,nk_arr,lp3d_arr_1g,lp3d_arr2,lna,lnk)
    lpnl_2g =  splin2(lp3d_lna_arr,lp3d_lnk_arr,na_arr,nk_arr,lp3d_arr_2g,lp3d_arr2,lna,lnk)
    pnl     =  Exp(lpnl)
    pnl_1g  =  Exp(lpnl_1g)
    pnl_2g  =  Exp(lpnl_2g)
!    write(*,*) 'ds fast_pk_nl', lpnl_1g, pnl_1g

  End Subroutine fast_pk_nl
  !! =====================

  !! =============================
  Real*8 Function P_dd_ln(k,z_arg)
  !! =============================
    !    Compute the linear power spectrum using the TF from Eisenstein & Hu 99
    !    Norm from Bunn & White 97
    !    Flat cosmology Omega_l = 1.0 - Omega_m
    !    k in input is in [h Mpc-1]
    !    Pk in [Mpc/h]^3

    Implicit None
    
    Real*8   :: z_arg,k
    Real*8   :: DoD0,Tm

    Tm      = TF_master(k*hub)
    P_dd_ln = anorm*Tm**2.*(k*hub)**ns*2.*pi**2*hub**3 !(EH99 A1)
    !! Agrees as is with CAMB with the hub**3 factor (OD 12/09/10)...
    !! So does the NL correction with halofit
    !! This works for flat or non-flat cosmology, but not DoD0

    If (z_arg .Ne. 0.0) Then
       DoD0    = fast_DoD0(z_arg)
       P_dd_ln = P_dd_ln*DoD0**2
    Endif

    Return

  End Function P_dd_ln
  !! =================

  !! ==========================
  Real(DP) Function bias(Mv,z)
  !! ==========================
    !     Compute the bias parameter 
    !     Mo & White 1997, ST99 & HK02

    !     Compute as a function of M_180 but input is Mv

    Implicit None

    Real(DP) :: Mv,z,sig,a,p,dc,M180,nu,nu2,M200,y,Delta!,mv2mh

    Real(DP) :: B,C,aa,bb,cc

    !external delta_c,splint,mv2mh,fast_sigma_ln

!!$    !! HK02
!!$    M180 = mv2mh(Mv,180.d0,z)
!!$    sig  = fast_sigma_ln(M180,z)
!!$    a   = 0.75
!!$    p   = 0.3
!!$    dc  = delta_c(z)
!!$    nu  = dc/sig
!!$    nu2 = nu**2.
!!$    bias = 1. + (a*nu2 - 1.)/dc + 2.*p/(dc*(1.+(a*nu2)**p))
!!$
    !! === Tinker et al. 1001.3162, Eq 6, Table 2
    !! ==========================================
    !! No observed redshift dependence till z~2
    Delta = 200
    y     = log10(Delta)

    !! --------- Table 2
    A  = 1.0+0.24*y*exp(-(4/y)**4)
    aa = 0.44*y-0.88
    B  = 0.183
    bb = 1.5
    C  = 0.019+0.107*y+0.19*exp(-(4/y)**4)
    cc = 2.4
    !! --------- 

    M200 = Mv ! IF we use the Tinker dN/dM, then the mass definition is consistent and we are using M200 basically
    sig  = fast_sigma_ln(M200,z)
    dc   = 1.686
    nu   = dc/sig
    bias = 1. - A*nu**aa/(nu**aa+dc**aa) + B*nu**bb + C*nu**cc
    !! ====

    Return

  End Function bias
  !! ==============


  !! ==============================
  Real*8 Function bias1_ST(Mv,zarg)
  !! ==============================
!! Order 1 bias as in Manera et al 2011

    Implicit None

    Real(DP)  :: Mv, zarg, M180,sig, epsilon1, E1
    Real(DP)  :: d_c, aa, p, nu

    d_c = delta_c(zarg)
    M180 = mv2mh(Mv   ,180.d0,zarg) 
    sig  = fast_sigma_ln(M180 ,zarg)
    nu = (d_c/sig )**2

    epsilon1 = (aa_ST*nu-1)/d_c
    E1 = 2*p_ST/d_c/(1+(aa_ST*nu)**p_ST)
    bias1_ST = 1 + epsilon1 + E1

    Return

  End Function Bias1_ST
  !! =======================



  !! ==============================
  Real(DP) Function bias2_ST(Mv,zarg)
  !! ==============================
!! Order 2 bias as in Manera et al 2011

    Implicit None

    Real(DP)  :: Mv, zarg, M180, sig, epsilon1, E1
    Real(DP)  :: d_c, nu, epsilon2, E2, a2

    d_c = delta_c(zarg)
    M180 = mv2mh(Mv   ,180.d0,zarg) 
    sig  = fast_sigma_ln(M180 ,zarg)
    nu = (d_c/sig)**2

    epsilon1 = (aa_ST*nu-1)/d_c
!    epsilon2 = (aa_ST*nu)/d_c*(((aa_ST*nu)**2-6*aa_ST*nu+3.)/d_c)
    epsilon2 = aa_ST*nu*(aa_ST*nu-3)/d_c**2 
    E1 = 2*p_ST/d_c/(1+(aa_ST*nu)**p_ST)
!    E2 = E1*((1+2*p_ST)/d_c+2*epsilon1)
    E2 = E1*(2*p_ST+2*aa_ST*nu-1)/d_c
    a2 = -17./21.
    bias2_ST = 2*(1+a2)*(epsilon1+E1)+epsilon2+E2
!    write(*,*)aa_ST, p_ST, bias2_ST, epsilon1, E1, epsilon2, E2
    Return

  End Function Bias2_ST
  !! =====================



  !! =======================
  Real*8 Function conc(Mv,z)
  !! =======================
    !     Compute the concentration parameter adequate for the NFW profile
    !     From Komatsu & Seljak 02 (Eq. 10) or Hu & Kratsov 02 or Cooray, Hu & Miralda-Escude 01
    !     Flat cosmology (fiducial)

    !     Definition valid for Mv

    Implicit None

    Real*8   Mv,z

    If (Mstar .Eq. 0.0 .Or. z_Mstar .Ne. z) Then 
       Mstar   = Mstar_finder(z)
       z_Mstar = z
    Endif

    If (Mv .Gt. Mcut) Then
       conc = 9.0/(1.0+z)*(Mv/Mstar)**(-0.13)    ! Dolney 04 (6) & HK C14
    Else
       conc = cs/(1.0+z)*(Mv/Mstar)**(ncs)    ! Dolney 04 (6) & HK C14
    End If

!!$!     -------------------------------------------------------------------
!!$     Formula of Cooray, Hu & Miralda-Escude 00 (Eq. 8)
!!$      a = 10.3*(1.+z)**(-0.3)
!!$      b = 0.24*(1.+z)**(-0.3)!!$      Mstar = Mstar_finder(z)
!!$      conc  = a * (Mv/Mstar)**(-b)
!!$      WRITE(*,*) 'conc0,conc1,mstar',conc,9.0/(1.0+z)*(Mv/Mstar_finder(0.d0))**(-0.13),mstar
!!$!     ------------------------------------------------------------------

  End Function conc
  !! ===============
  
 !! ==============================
  Real*8 Function u_rho_nfw(k,Mv,z)
  !! ==============================
    ! Compute analyticaly the Fourier transform of a NFW density profile
    ! Mv input is virial mass
    ! Use eg formula (81) from Sheth & Cooray KH02

    Implicit None

    Real*8 :: k,Mv,z,mh
    Real*8 :: dv,c,rho_bar,r_v3,r_v,r_s,x,f,rho_s
    Real*8 :: cikrs,sikrs,cickrs,sickrs,krs,ckrs,krv,ckrv,cikrv,cickrv,sikrv,sickrv

    dv       = delta_v(z)
    c        = conc(Mv,z)
    Mh       = mv2mh(Mv,180.d0,z) 
    rho_bar  = rho_bar_m(z) ! Fix rho_ba at z = 0

    ! Choose Mh or Mv here...
    r_v3     = 3.*Mv/(4.*3.14159*rho_bar*dv) ! Follow definition from HK (C5)
    !r_v3     = 3.*Mh/(4.*3.14159*rho_bar*dv)
    r_v      = (r_v3)**(1./3.)
    r_s      = r_v/c

    krs      = k*r_s
    !krv      = k*r_v
    ckrs     = (1+c)*krs  
    !ckrv     = (1+c)*krv  

    Call cisi(krs ,cikrs ,sikrs )
    Call cisi(ckrs,cickrs,sickrs)
    !Call cisi(krv/c ,cikrv ,sikrv)
    !Call cisi(ckrv/c,cickrv,sickrv)
    x        = 1./c
    f        = (x**3)*(Log(1.+1./x)-1./(1.+x))
    rho_s    = Mv/(4.*3.14159*r_v3*f)              ! HK (C5)  (Checked that both are equivalent)
    !rho_s    = Mh/(4.*3.14159*r_v3*f)              ! HK (C5)  (Checked that both are equivalent)

    u_rho_nfw = 4.*3.14159*rho_s*r_s**3/Mv
    !u_rho_nfw = 4.*3.14159*rho_s*r_s**3/Mh
    u_rho_nfw = u_rho_nfw*(Sin(krs)*(sickrs-sikrs)+Cos(krs)*(cickrs-cikrs)-Sin(c*krs)/ckrs)
    ! SC (81)
    !u_rho_nfw = -sin(krv)/ckrv+cos(krv/c)*(cickrv-cikrv)

    Return

  End Function u_rho_nfw
  !! ===================

  !! =======================
  Real*8 Function delta_v(z) 
  !! =======================
    !     Gives a definition of delta_v necessary to compute M_vir
    !     Inspired from BN 98 quoted by KH (C6) 

    Implicit None

    Real*8 :: z,x

    x = omega_f(z) - 1.0

    delta_v = 18.*3.14159**2+82.*x-39.*x**2.
    delta_v = delta_v/(1.+x)

    Return
  End Function Delta_v
  !! ==================

  !! ===========================
  Real*8 Function mv2mh(Mv,dh,z) 
  !! ===========================
    !     Convert Mv to Mh defined as Mh = 4pi/3 r_h^3 dh rho_m 
    !     For our purpose, dh = 180
    !     Inspired from KH (C12)

    Implicit Real*8 (a-h,k,o-z)

    Real*8 Mv,g,mh

!!$    dv = delta_v(z)
!!$    x  = conc(Mv,z)
!!$    y  = (dh/dv)**(1./3.)
!!$
!!$    g  = -3.*(x+y)
!!$    g  =  g * (x-x*y+(1.+x)*(x+y)*Log(1.+x) - (1.+x)*(x+y)*Log(1.+x/y))
!!$    g  =  g / ( 3.*(1.+x)*((x+y)**2)*Log(1.+x) - x * (x + 4.*(x**2.) + 6.*x*y + 3.*y**2) )
!!$
!!$    ratio = 1. + g
!!$
!!$    mv2mh = ratio * Mv
!!$
!!$    Write(*,*) 'mv2mh 1',mv2mh,Mv,ratio,g,x,y,dv!Mv,dh

    !! == Test an alternative way using HK C8-10
    !     Test ok at low z... but different at high z. To be investigated more
    dv = delta_v(z)
    c  = conc(Mv,z)
    x  = 1./c
    f  = (x**3)*(log(1.+1./x)-1./(1.+x))
    fh = dh/dv*f
    
    a1 =  0.5116
    a2 = -0.4283
    a3 = -3.13*10.**(-3)
    a4 = -3.52*10.**(-5)
    
    p  = a2 + a3*log(fh) + a4*(log(fh))**2.
    xf = ( a1*(fh)**(2.*p) + (3./4.)**2. )**(-1./2.) + 2.*fh   
    mh = Mv * dh/dv * (1./(c*xf))**3.
    mv2mh = mh
    !! ==
    !Write(*,*) 'mv2mh 1',mv2mh,Mv,ratio,g,x,y,dv!Mv,dh

    Return

  End Function mv2mh
  !! ===============


 !     ============================================================================
  Real*8 Function Mstar_finder(z)
    !     ============================================================================
    !     Compute the linear mass scale
    !     by simple dichotomy if needed or use the formual of Komatsu & Seljak 02
    !     Compute as a function of M_180
    !     Should output M_vir

    Implicit Real*8 (a-h,k,o-z)

    Real(DP) ::  M,Mo,Mt,nu

    !external delta_c,fast_sigma_ln

    tol_loc  = 0.01
    zref     = z

    !     First guess
    d_c  = delta_c(zref)
    M    = 1.0d6
    sig  = fast_sigma_ln(M,zref)

    nu   = d_c/sig

    !     In case of a bad first guess
    Do While (nu .Gt. 1.0)
       Mo  = M
       M   = Mo/10.0
       sig = fast_sigma_ln(M,zref)
       !sig = sigma_ln(M,zref)
       nu  = d_c/sig
    End Do

    !      Rough search
    Do While (nu .Lt. 1.0)
       Mo  = M
       M   = Mo*5.0
       sig = fast_sigma_ln(M,zref) 
       !sig = sigma_ln(M,zref) 
       nu  = d_c/sig
    End Do

    !      Proper dichotomy
    Do While ((M-Mo)/M .Gt. tol_loc)
       Mt  = (M + Mo)/2.
       sig = fast_sigma_ln(Mt,zref)
       nu  = d_c/sig
       If  (nu .Lt. 1) Then
          Mo = Mt
       Else
          M = Mt
       Endif
    End Do

    !write(*,*) 'Mstar_finder, nu :', Mt, nu

    Mstar_finder = Mt

!!$    WRITE(*,*) ''
!!$    WRITE(*,*)'==================================================='
!!$    WRITE(*,*) 'Mstar computed at z = ', zref
!!$    WRITE(*,*) 'Mstar                 ', Mstar_finder
!!$    WRITE(*,*)'==================================================='
!!$    WRITE(*,*) ''

  End Function Mstar_finder
  !      ============================================================================

  !! ===============================
  Real*8 Function fast_sigma_ln(M,z)
  !!================================
    !    Compute sigma_m at redshift M and z using precompiled table at z = 0.0
    !    Top hat smoothing
    !    Flat Universe
    !    M in [h-1 Ms]

    Implicit None

    Real(DP)   :: M,D,D0,M0,z,rho0,rho,sig

    External splint

    ! Convert Mass scales from 0 (where sig table is calculated) to z
    ! Mass corresponding to the same smoothing scale R
    !rho0 =  rho_bar_m(0.d0)
    !rho  =  rho_bar_m(z)
    !M0   = (rho0/rho) * M 
    M0 = M ! OD No need for scaling for Tinker 12/11/10
    M0 = M/hub ! Ad hoc rescaling here. Allow all the masses to be Msun/h later. 12/11/10
    !Write(*,*) 'rho ',z,M,M0,rho0,rho

    ! Test range of interpolation
    If (M0 .Le. M_min_sig .Or. M0 .Ge. M_max_sig) Then
       Write(*,*) ''
       Write(*,*) 'Sigma interpolation out of range '
       Write(*,*) 'Stop here'
       Write(*,*) 'M_min_sig :  ',M_min_sig
       Write(*,*) 'M_max_sig :  ',M_max_sig
       Write(*,*) 'M         :  ',M0
       Write(*,*) ''
       Stop
    Endif

    ! Interpolate
    Call splint(M_arr,sig_arr,sig2_arr,n_sig,M0,sig)
    fast_sigma_ln = sig

    ! Evolve linearly
    If (z .Ne. 0.0) Then
       fast_sigma_ln = sig*fast_DoD0(z)
    Endif

    Return

  End Function fast_sigma_ln
  !! =======================

  !! =========================
  Subroutine initiate_sigmas()
  !! =========================

    Implicit None

    Real(DP)  :: M,sig,sig_g,M_lnstep,lMin,lMax
    Integer   :: i

    M_lnstep  = Log10(M_max_sig/M_min_sig)/(n_sig-1.)
    lMin      = Log10(M_min_sig)
    lMax      = Log10(M_max_sig)
    Do i = 1,n_sig
       M           = 10.**(lMin+(i-1.)*M_lnstep)
       sig         = sigma_ln(M,0.d0)
       sig_g       = sigma_ln_gauss(M,0.d0)
       M_arr(i)    = M
       sig_arr(i)  = sig
       sigg_arr(i) = sig_g 
    End Do
    Call spline(M_arr,sig_arr ,n_sig,3.d30,3.d30,sig2_arr )
    Call spline(M_arr,sigg_arr,n_sig,3.d30,3.d30,sigg2_arr)
    !Write(*,*) 'Sigma array computed'
61  Format(1X,3E13.5)

  End Subroutine initiate_sigmas
  !! ===========================
  !! ==========================
  Real*8 Function sigma_ln(M,z)
  !! ==========================
    !    Compute sigma_m at redshift M and z, follow EH99 (A7)
    !    Top hat smoothing
    !    Flat Universe
    !    M in input is in [Ms]

    Implicit None

    Real*8 :: z,rho_bar,M,R,R3,int,rombint,kmin,kmax,sigma_ln2

    External rombint

    z_sig   = z
    rho_bar = rho_bar_m(z_sig)
    R3      = 3.0*M/(4.0*3.14159*rho_bar)  
    R_sig   = R3**(1./3.)

    kmin = 1.0d-5    ![Mpc-1]
    kmax = 1.0d5

    int  = rombint(pwk2,Log(kmin),Log(kmax),tol)

    sigma_ln2 = int!/(2.*3.14159**2.)       ! cf eg EH (A7)
    sigma_ln  = Sqrt(sigma_ln2)

    Return

  End Function sigma_ln
  !! ==================

  !! ======================
  Real*8 Function pwk2(lnk)
  !! ======================
    !    From EH99 (A7) 

    Implicit None

    Real*8 :: wr,kr,N_nu,lnk,k,dln,pln,Tm

    k    = Exp(lnk)

    Tm   = TF_master(k)
    !pln  = P_dd_ln(k,z_sig)*(2.*pi)**3 ! Cf comment in P_dd_ln
    dln  = anorm*k**(3+ns)*Tm**2     
    !Call NDL_computeP(par_ndl,k,pln,verdict)
    kr   = k*R_sig  
    wr   = 3.*(Sin(kr) - kr*Cos(kr))/(kr**3)
    pwk2 = dln*wr**2
    ! Extra k due to the log integration
    Return

  End Function pwk2
  !! ==============

  !! ================================
  Real*8 Function sigma_ln_gauss(M,z)
  !! ================================
    !    Compute sigma_m at redshift M and z, follow EH99 (A7)
    !    Gaussian smoothing
    !    Flat Universe
    !    M in input is in [Ms]

    Implicit None

    Real*8 :: M,z,rho_bar,R,R3,int,kmin,kmax,sigma_ln2,rombint

    External rombint

    z_sig   = z
    rho_bar = rho_bar_m(z_sig)
    R3      = 3.0*M/(4.0*3.14159*rho_bar)   
    R_sig   = R3**(1./3.)

    kmin = 1.0d-5     ![Mpc-1]
    kmax = 1.0d5

    int  = rombint(pwk2_gauss,Log(kmin),Log(kmax),tol)

    sigma_ln2 = int!/(2.*3.14159**2.)       ! cf eg EH (A7)
    sigma_ln_gauss  = Sqrt(sigma_ln2)

    Return

  End Function sigma_ln_gauss
  !! ========================

  !! ============================
  Real*8 Function pwk2_gauss(lnk)
  !! ============================
    !    From Smith et al. (56) 

    Implicit None

    Real*8 :: k,delta2,wr,kr,pln,lnk,Tm,dln

    k          = Exp(lnk)
    !pln  = P_dd_ln(k,z_sig) ! Cf comment in P_dd_ln
    Tm         = TF_master(k)
    dln        = anorm*k**(3+ns)*Tm**2
    !delta2     = k**3*pln
    kr         = k*R_sig  
    wr         = Exp(-kr**2.)
    pwk2_gauss = dln*wr
    ! extra k due to the log integration

    Return

  End Function pwk2_gauss
  !! ====================

  !      ============================================================================
  Real*8 Function delta_c(z)
    !      ============================================================================
    !      Compute the critical density
    !      Kitayama,Suto 1996 (Eq. A6)
    !      Flat Universe

    Implicit Real*8 (a-h,k,o-z)

    Real(DP) :: number,om_f
    !,omega_f
    !external omega_f

    om_f    = omega_f(z)
    number  = (3./20.)*(12*3.14159)**(2./3.)

    delta_c = number*(1.0 + 0.0123*Log10(om_f))

  End Function delta_c
  !      ============================================================================

 
  !! =============================
  Function gg(om_m_loc,om_v_loc)        
  !! =============================
  !! growth factor for linear fluctuations 
    Implicit None
    
    Real*8 gg,om_m_loc,om_v_loc
    
    gg=2.5*om_m_loc/(om_m_loc**(4./7.)-om_v_loc+(1d0+om_m_loc/2.)*(1.+om_v_loc/70.))
    
    Return
  End Function gg
  !! ============

  !! ============================
  Real*8 Function dndm_f(Mv,zarg)
  !! ============================
    !     Compute dndM for mass M at z, using Jenkins et al. 01
    !     or Sheth & Tormen 99
    !     Top hat smoothing (e.g. HK02)
    !     Flat Universe
    !     Jenkins def. valid for M180 and here Mv in input

    Implicit None
    
    Real*8 :: Mv,zarg,M180,nu,M180b,d_c,nub,nup,A,M,b,c,aa,dmlnsigdMv,al 
    Real*8 :: f,pdn_dm,sigb,sig,dnudMv,dlns_dM,sig2,p,dM180dMv,rho_bar
    Real*8 :: dnu_dm,dc,siga,Ma,Mb


    Select Case(mass_function)
       
       Case(1,0)
          !! ================================================================
          !!     Jenkins part
          rho_bar  = rho_bar_m(zarg)  
          
          M180     = mv2mh(Mv   ,180.d0,zarg)
          M180b    = mv2mh(Mv+dM,180.d0,zarg)     
          
          !Write(*,*) 'dndm ',Mv,M180,zarg
          
          dM180dMv = (M180b-M180)/dM

          ! Compute dlns_dM by straightforward finite differencing
          sig  = fast_sigma_ln(M180 ,zarg)
          sig2 = fast_sigma_ln(M180b,zarg)

          
          dlns_dM = Log(sig/sig2)/dM
          
          dndm_f = 0.301*Exp(-(Abs(-1.*Log(sig)+0.64))**3.82)
          dndm_f = dndm_f * (rho_bar/M180) * dlns_dM * dM180dMv
          !! ===============================================================

       Case(2)
          !! ===============================
          !! Press & Schechter
          !! ===============================
          !!   This part is the classical PS
          d_c = delta_c(0.d0)          
          dnu_dm = dc*(1./sig2 - 1./sig)/dM
          nu     = dc/sig
          dndm_f = rho_bar/M * dnu_dm * Sqrt(2./3.14)*Exp(-nu**2/2.)
          !!  ===============================================================

       Case(3)
          !! ================================
          !! Sheth & Tormen as in Sheth and tormen 1999
          !! Watch out. Diff definition of nu
          !! ================================
          d_c  = delta_c(zarg) 
          M180 = mv2mh(Mv   ,180.d0,zarg) 
          M180b= mv2mh(Mv+dM,180.d0,zarg)
          dM180dMv = (M180b-M180)/dM
          sig  = fast_sigma_ln(M180 ,zarg)
          sigb = fast_sigma_ln(M180b,zarg)
          nu   = (d_c/sig )**2
          nub  = (d_c/sigb)**2
          !A    = 0.116!1.0
          !a    = 0.67 !0.707
          !p    = 0.33 !0.3
!as in the original paper (sheth)
!!---------------------------
!          A_ST    = 0.322
!          aa_ST   = 0.75!0.707
!          p_ST    = 0.3
!!---------------------------
          nup  = aa_ST*nu
          f    = A_ST/nu*(1.+nup**(-p_ST))*Sqrt(nup)*Exp(-nup/2.)/sqrt(2*3.14159)
          rho_bar = rho_bar_m(0.d0)
          dnudMv  = (nub-nu)/dM
          dndm_f  = rho_bar/M180*f*dnudMv*dM180dMv

          !!Input M is Msun/h. Output dndM is h^4/Ms/Mpc^3 
          !! =====================================

       Case(4)
          !! =================================
          !! Tinker et al. 0803.2706 fitting formula
          !! Eq 3, \Delta = 200
          !! And parameters given in Eq. 5,6,7 and Tab 2
          !! =================================
          Ma   = 0.9*Mv
          Mb   = 1.1*Mv 
          sig  = fast_sigma_ln(Mv,zarg)
          siga = fast_sigma_ln(Ma,zarg)
          sigb = fast_sigma_ln(Mb,zarg)
          A    = 0.186*(1+zarg)**(-0.14) !A
          aa   = 1.47*(1+zarg)**(-0.06)  !a
          al   = 10.**(-(0.75/log10(200./75.))**1.2)
          b    = 2.57*(1+zarg)**(-al)
          c    = 1.19
          f    = A*((sig/b)**(-aa)+1.)*exp(-c/sig**2)

          rho_bar    =  rho_bar_m(0.d0)
          dmlnsigdMv = -(log(sigb)-log(siga))/(Mb-Ma)

          dndm_f     =  f*rho_bar/Mv*dmlnsigdMv
          dndm_f     =  dndm_f/hub**2
          If (dndm_f .LT. 1.d-90) dndm_f = 0.d0
!          write(*,*) zarg, log10(Mv), f, rho_bar, dmlnsigdMv, log10(dndm_f)
!stop
!          Write(*,50) M,dndm_f,A,sig/b,aa,c/sig**2,exp(-c/sig**2),(sig/b)**a
!50        Format(1X,20E17.9)
          !! OD : This has been thoroughly tested to agree with Tinker, including z scaling now 12/11/10
          !! Input M is Msun/h. Output dndM is h^4/Ms/Mpc^3 
          !! ===============================================

       Case(5)
          !! ================================
          !! Tinker et al. 0803.2706 fitting formula
          !! Eq 3, \Delta = 500
          !! And parameters given in Eq. 5,6,7 and Tab 2
          !! ================================
          !! Computed with average between 400 - 600 fit would work better
          !! A      0.215205
          !! a       1.57375
          !! b       1.94432
          !! c       1.39670
          !d_c  = delta_c(z) 
          Ma   = 0.9*Mv
          Mb   = 1.1*Mv        
          sig  = fast_sigma_ln(Mv,zarg)
          siga = fast_sigma_ln(Ma,zarg)
          sigb = fast_sigma_ln(Mb,zarg)
          A    = 0.215205*(1+zarg)**(-0.14) !A
          aa   = 1.57375*(1+zarg)**(-0.06) !a
          al   = 10.**(-(0.75/log10(500./75.))**1.2)
          b    = 1.94432*(1+zarg)**(-al)
          c    = 1.39670
          f    = A*((sig/b)**(-aa)+1.)*exp(-c/sig**2)
          rho_bar    = rho_bar_m(0.d0)
          dmlnsigdMv = -(sigb-siga)/sig/(Mb-Ma)
          dndm_f     = f*rho_bar/Mv*dmlnsigdMv
          dndm_f     = dndm_f/hub**2
          If (dndm_f .LT. 1.d-90) dndm_f = 0.d0
          !! ===============================

       Case default
          Write(*,*) 'Wrong mass function'
          Stop

       End Select
       Return

  End Function dndm_f
  !! ================


  !! =====================================================
  Function dndm_sub_f(msub,Mhalo,zarg)
!! Subhalos mass functions from Tinker & Wetzel 2010 Eq. 10 
!! Mhalo : M200
!! Checked!

    Real(DP)   :: msub, Mhalo, zarg, ratio_mass, inexp, dndm_sub_f

    ratio_mass = msub/Mhalo

    inexp = -9.9 * ratio_mass**2.5
    dndm_sub_f = 0.13 * ratio_mass**(-0.7) * exp(inexp) / msub
!!    dndm_sub_f = ratio_mass
    Return


  END Function dndm_sub_f
  !! =====================================================

End Module power_spec_tools
