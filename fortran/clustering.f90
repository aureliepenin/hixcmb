Module clustering

  Use nrtype
  Use eh99_tools
  Use cosmo_tools
  Use power_spec_tools
  Use initiate_all
  Use HI_stuffs
  Use clustering_limber

  Implicit None

  Real(DP) :: kperp_for_int, zmean_for_int, zmin_for_int, zmax_for_int
  Real(DP) :: eta_for_int, delta_eta_for_int, kpara_for_int
  Real(DP) :: bias_for_fisher, f_for_fisher, zmean_for_fisher
  Real(DP) :: sigma_NL = 7 
!!!ATTENTION AU H!!!!!!!!!!!!!!!!!!

  Public 
  
Contains



 !======================================================
  Subroutine Compute_cl_hi_kap_for_fisher()

    Integer  :: ibin, il
    Real(DP) :: zmin, zmax

    Do ibin = 1, nbins_freq
       zmin = 1.0!redshift_from_freq(freq_min + dfreq * ibin)
       zmax = 1.005!redshift_from_freq(freq_min + dfreq * (ibin-1))
       write(*,*) "zmin, zmax = ", freq_min + dfreq * ibin, zmin, zmax 
       Allocate(freq_arr(ibin)%cl_hi_kap_arr(1:nl_arr))
       Do il = 1, nl_arr
          freq_arr(ibin)%cl_hi_kap_fish_arr(il) = cl_hi_kap_fish_func(l_arr(il),zmin,zmax)
!          write(*,*) l_arr(il), freq_arr(ibin)%cl_hi_kap_fish_arr(il)
       End Do
    End Do

  End Subroutine Compute_cl_hi_kap_for_fisher
  !======================================================
 
 !======================================================
  Subroutine Compute_cl_hi_kap()

    Integer  :: ibin, il
    Real(DP) :: zmin, zmax

    Do ibin = 1, nbins_freq
       zmin = 1.0!redshift_from_freq(freq_min + dfreq * ibin)
       zmax = 1.005!redshift_from_freq(freq_min + dfreq * (ibin-1))
       write(*,*) "zmin, zmax = ", freq_min + dfreq * ibin, zmin, zmax 
       Allocate(freq_arr(ibin)%cl_hi_kap_arr(1:nl_arr))
       Do il = 1, nl_arr
          freq_arr(ibin)%cl_hi_kap_arr(il) = cl_hi_kap_func(l_arr(il),zmin,zmax)
!          write(*,*) l_arr(il), freq_arr(ibin)%cl_hi_kap_arr(il)
       End Do
    End Do


  End Subroutine Compute_cl_hi_kap
  !======================================================
 
  !======================================================
  Function cl_hi_kap_func(ell,zmin,zmax)
    
    Real(DP) :: cl_hi_kap_func, ell, zmin, zmax, rombint, tol_i
    Real(DP) :: eta, eta0, Temp
    Real(DP) :: kpara_min, kpara_max, log_kpara_min, log_kpara_max
    external :: rombint

    tol_i = 1.d-4
    eta0          = conftime(0.d0)
    eta_for_int   = eta0 - conftime(zmin)

    kperp_for_int = ell / eta_for_int
    zmin_for_int  = zmin 
    zmax_for_int  = zmax 


   kpara_min = 2.0 * pi / eta_for_int * (1.d0 + zmin_for_int) 
!    kpara_max = 0.21 ! cf Bull+2015
    kpara_max = 0.15 * 1.d0 ! cf Bull+2015 h^-1 taken into account

! to check
!    kpara_min =  1.e-5!2.0 * pi / eta_for_int * (1.d0 + zmean_for_int)
!    kpara_max = pi/2./delta_eta_for_int 
!    write(*,*) "klim = ", kpara_min, kpara_max

    log_kpara_min = log10(kpara_min) * 1.d0
    log_kpara_max = log10(kpara_max) * 1.d0

    cl_hi_kap_func = rombint(cl_hi_kap_int,log_kpara_min,log_kpara_max,tol_i) 

    Temp = T_HI_mean(zmin)
    cl_hi_kap_func = cl_hi_kap_func * Temp / eta_for_int 
!    write(*,*) "in cl_cross", Temp, eta_for_int
    
  End Function cl_hi_kap_func
  !======================================================


  !======================================================
  Function cl_hi_kap_rsd_func(ell,zmin,zmax)
    
    Real(DP) :: cl_hi_kap_rsd_func, ell, zmin, zmax, rombint, tol_i
    Real(DP) :: eta, eta0, Temp
    Real(DP) :: kpara_min, kpara_max, log_kpara_min, log_kpara_max
    external :: rombint

    tol_i = 1.d-4
    eta0          = conftime(0.d0)
    eta_for_int   = eta0 - conftime(zmin)

    kperp_for_int = ell / eta_for_int
    zmin_for_int  = zmin 
    zmax_for_int  = zmax 

    kpara_min = 2.0 * pi / eta_for_int * (1.d0 + zmin_for_int) 
!    kpara_max = 0.21 ! cf Bull+2015
    kpara_max = 0.15 * 1.d0 ! cf Bull+2015 h^-1 taken into account
    log_kpara_min = log10(kpara_min) * 1.d0
    log_kpara_max = log10(kpara_max) * 1.d0

    cl_hi_kap_rsd_func = rombint(cl_hi_kap_rsd_int,log_kpara_min,log_kpara_max,tol_i) 

    Temp = T_HI_mean(zmin)
    cl_hi_kap_rsd_func = cl_hi_kap_rsd_func * Temp / eta_for_int 
!    write(*,*) "in cl_cross", Temp, eta_for_int
    
  End Function cl_hi_kap_rsd_func
  !======================================================

  !======================================================
  Function cl_hi_kap_bias_func(ell,zmin,zmax)
    
    Real(DP) :: cl_hi_kap_bias_func, ell, zmin, zmax, rombint, tol_i
    Real(DP) :: eta, eta0, Temp
    Real(DP) :: kpara_min, kpara_max, log_kpara_min, log_kpara_max
    external :: rombint

    tol_i = 1.d-4
    eta0          = conftime(0.d0)
    eta_for_int   = eta0 - conftime(zmin)

    kperp_for_int = ell / eta_for_int
    zmin_for_int  = zmin 
    zmax_for_int  = zmax 

    kpara_min = 2.0 * pi / eta_for_int * (1.d0 + zmin_for_int) 
!    kpara_max = 0.21 ! cf Bull+2015
    kpara_max = 0.15 * 1.d0 ! cf Bull+2015 h^-1 taken into account
    log_kpara_min = log10(kpara_min) * 1.d0
    log_kpara_max = log10(kpara_max) * 1.d0

    cl_hi_kap_bias_func = rombint(cl_hi_kap_bias_int,log_kpara_min,log_kpara_max,tol_i) 

    Temp = T_HI_mean(zmin)
    cl_hi_kap_bias_func = cl_hi_kap_bias_func * Temp / eta_for_int 
!    write(*,*) "in cl_cross", Temp, eta_for_int
    
  End Function cl_hi_kap_bias_func
  !======================================================

  !======================================================
  Function cl_hi_kap_int(logkpara)

    Real(DP) :: logkpara, cl_hi_kap_int, tol_i
    Real(DP) :: rombint
    external :: rombint

    tol_i = 1.d-4 
    kpara_for_int = 10.d0**logkpara
    cl_hi_kap_int  = rombint(cross_over_r_int,zmin_for_int,zmax_for_int,tol_i)
    cl_hi_kap_int  = cl_hi_kap_int * kpara_for_int * log(10.)

!! extra kara because of the log integral

  End Function cl_hi_kap_int
  !======================================================

  !======================================================
  Function cl_hi_kap_bias_int(logkpara)

    Real(DP) :: logkpara, cl_hi_kap_bias_int, tol_i
    Real(DP) :: rombint
    external :: rombint

    tol_i = 1.d-4 
    kpara_for_int = 10.d0**logkpara
    cl_hi_kap_bias_int  = rombint(cross_bias_over_r_int,zmin_for_int,zmax_for_int,tol_i)
    cl_hi_kap_bias_int  = cl_hi_kap_bias_int * kpara_for_int * log(10.)

!! extra kara because of the log integral

  End Function cl_hi_kap_bias_int
  !======================================================

 !======================================================
  Function cl_hi_kap_rsd_int(logkpara)

    Real(DP) :: logkpara, cl_hi_kap_rsd_int, tol_i
    Real(DP) :: rombint
    external :: rombint

    tol_i = 1.d-4 
    kpara_for_int = 10.d0**logkpara
    cl_hi_kap_rsd_int  = rombint(cross_rsd_over_r_int,zmin_for_int,zmax_for_int,tol_i)
    cl_hi_kap_rsd_int  = cl_hi_kap_rsd_int * kpara_for_int * log(10.)

!! extra kara because of the log integral

  End Function cl_hi_kap_rsd_int
  !======================================================

  !======================================================
  Function cl_hi_kap_fish_func(ell,zmin,zmax)
    
    Real(DP) :: cl_hi_kap_fish_func, ell, zmin, zmax, rombint, tol_i
    Real(DP) :: eta, eta0, Temp
    Real(DP) :: kpara_min, kpara_max, log_kpara_min, log_kpara_max
    external :: rombint


    tol_i = 1.d-4
    eta0          = conftime(0.d0)
    eta_for_int   = eta0 - conftime(zmin)

    kperp_for_int = ell / eta_for_int
    zmin_for_int  = zmin 
    zmax_for_int  = zmax 

    kpara_min = 2.0 * pi / eta_for_int * (1.d0 + zmin_for_int) 
!    kpara_max = 0.21 ! cf Bull+2015
    kpara_max = 0.15 * 1.d0 ! cf Bull+2015 h^-1 taken into account
    log_kpara_min = log10(kpara_min) * 1.d0
    log_kpara_max = log10(kpara_max) * 1.d0

    cl_hi_kap_fish_func = rombint(cl_cross_for_fisher_int,log_kpara_min,log_kpara_max,tol_i) 

    Temp = T_HI_mean(zmin)
    cl_hi_kap_fish_func = cl_hi_kap_fish_func * Temp / eta_for_int 
!    write(*,*) "in cl_cross", Temp, eta_for_int
    
  End Function cl_hi_kap_fish_func
  !======================================================

  !======================================================
  Function cl_cross_for_fisher_int(logkpara)

    Real(DP) :: logkpara, cl_cross_for_fisher_int, tol_i
    Real(DP) :: rombint
    external :: rombint

    tol_i = 1.d-4 
    kpara_for_int = 10.d0**logkpara

    cl_cross_for_fisher_int  = rombint(cross_over_r_for_fisher_int,zmin_for_int,zmax_for_int,tol_i)
    cl_cross_for_fisher_int  = cl_cross_for_fisher_int * kpara_for_int * log(10.)

!! extra kara because of the log integral

  End Function cl_cross_for_fisher_int
  !======================================================

  !======================================================
  Function cross_over_r_for_fisher_int(z)

    Real(DP)  :: a, eta0, eta, z_source, eta_star, Omegam, frac, kappa_kernel
    Real(DP)  :: knorm, Plin, delta_eta, in_cos, bias, mu, f, bias_term
    Real(DP)  :: cross_over_r_for_fisher_int, z

    a        = 1. / (1. + z)
    eta0     = conftime(0.d0)
    eta      = eta0 - conftime(z)
    z_source = z_cmb
    eta_star = eta0 - conftime(z_source)
    Omegam   = Om0
    frac     = (eta_star - eta) / (eta_star * eta)
    kappa_kernel = eta / a * frac * 3.d0 / 2.d0 / invhub**2  * Omegam

    knorm    = sqrt(kpara_for_int**2 + kperp_for_int**2)
    Plin     = P_lin_cross(knorm,zmin_for_int,z)

    delta_eta = abs(eta_for_int - eta)
    in_cos    = kpara_for_int * delta_eta

    mu    = kpara_for_int / knorm

    bias_term = bias_for_fisher +  &
                f_for_fisher * kpara_for_int**2/knorm**2

    cross_over_r_for_fisher_int = kappa_kernel * cos(in_cos) * bias_term * Plin
!    write(*,*) kappa_kernel, cos(in_cos), bias_term, cross_over_r_int

  End Function Cross_over_r_for_fisher_int
  !======================================================

 !======================================================
  Subroutine Compute_cl_hi()
    
    Integer  :: il, ibin 
    Real(DP) :: zmin, zmax
    
    Do ibin = 1, nbins_freq
       zmin = 1.0!redshift_from_freq(freq_min + dfreq * ibin)
       zmax = 1.005!redshift_from_freq(freq_min + dfreq * (ibin-1))
       write(*,*) "zmin, zmax = ", freq_min + dfreq * ibin, zmin, zmax 
       Allocate(freq_arr(ibin)%cl_hi_arr(1:nl_arr))
       Do il = 1, nl_arr
          freq_arr(ibin)%cl_hi_arr(il) = cl_hi_func(l_arr(il),zmin,zmax)
!          write(*,*) l_arr(il), freq_arr(ibin)%cl_hi_arr(il)
       End Do
    End Do

  End Subroutine Compute_cl_hi
    !======================================================

  !======================================================
  Function cl_hi_func(ell,zmin,zmax)
!! Appendix A in Shaw+2014 for eta = eta_mean 
    
    Real(DP) :: cl_hi_func, ell, zmin, zmax, rombint, tol_i
    Real(DP) :: zmean, f, bOmega_HI, eta1, eta2, eta0, Tmin, Tmax, eta_mean_for_int
    Real(DP) :: kpara_min, kpara_max, log_kpara_min, log_kpara_max
    external :: rombint

    tol_i = 1.d-4
    zmean_for_int = (zmin+zmax)/2.

!    bOmega_HI = (Bias_HI(zmin)*Omega_HI(zmin) + Bias_HI(zmax)*Omega_HI(zmax))/2.
    eta0      = conftime(0.d0)
    eta1      = eta0 - conftime(zmin)
    eta2      = eta0 - conftime(zmax)
    eta_mean_for_int = (eta1 + eta2)/2.
    delta_eta_for_int = eta2 - eta1

    kperp_for_int = ell/eta_mean_for_int 
    zmin_for_int  = zmin 
    zmax_for_int  = zmax 

    kpara_min = 2.0 * pi / eta_mean_for_int * (1.d0 + zmean_for_int)
    kpara_max = 0.15*1.d0 ! cf Bull+2015

! to check
!    kpara_min =  1.d-5!2.0 * pi / eta_mean_for_int * (1.d0 + zmean_for_int)
!    kpara_max = pi/2./delta_eta_for_int
!    write(*,*) "kmin, kmax for hi only = ", kpara_min, kpara_max

    log_kpara_min = log10(kpara_min)
    log_kpara_max = log10(kpara_max)

    cl_hi_func = rombint(cl_hi_int,log_kpara_min,log_kpara_max,tol_i)

    Tmin = T_HI_mean(zmin)
    Tmax = T_HI_mean(zmax)

    cl_hi_func = cl_hi_func * Tmin * Tmax / eta1 / eta2 / pi 
!   write(*,*) "in cl_hi end ", kperp_for_int

  End Function cl_hi_func
  !======================================================


  !======================================================
  Function cl_hi_int(logkpara)
!! eta = eta mean 
!! Cf appendix A Shaw+2014

    Real(DP) :: cl_hi_int, bias, kpara, knorm, mu
    Real(DP) :: logkpara, f, eta_mean_for_int 
    Real(DP) :: in_cos, Plin

!!! ATTENTION AU ZMEAN
    bias  = Bias_HI(zmean_for_int) 
    kpara = 10.d0**logkpara
    knorm = sqrt(kpara**2 + kperp_for_int**2)
    mu    = kpara / knorm
    f     = growth_factor(zmean_for_int)

!    write(*,*) kpara, kperp_for_int, knorm
    Plin = P_lin_cross(knorm,zmin_for_int,zmax_for_int)
    in_cos = kpara * delta_eta_for_int
    cl_hi_int = (bias + f * mu**2)**2 * cos(in_cos) * Plin * log(10.) * kpara
!    cl_hi_int = (bias + f * mu**2)**2  * Plin * log(10.) * kpara

!! Extra kpara * log 10 because log integration

  End Function cl_hi_int
  !======================================================

  !======================================================
  Function cross_over_r_int(z)

    Real(DP)  :: a, eta0, eta, z_source, eta_star, Omegam, frac, kappa_kernel
    Real(DP)  :: knorm, Plin, delta_eta, in_cos, bias, mu, f, bias_term
    Real(DP)  :: cross_over_r_int, z

    a        = 1. / (1. + z)
    eta0     = conftime(0.d0)
    eta      = eta0 - conftime(z)
    z_source = z_cmb
    eta_star = eta0 - conftime(z_source)
    Omegam   = Om0
    frac     = (eta_star - eta) / (eta_star * eta)
    kappa_kernel = eta / a * frac * 3.d0 / 2.d0 / invhub**2  * Omegam

    knorm    = sqrt(kpara_for_int**2 + kperp_for_int**2)
    Plin     = P_lin_cross(knorm,zmin_for_int,z)

    delta_eta = abs(eta_for_int - eta)
    in_cos    = kpara_for_int * delta_eta

    bias  = Bias_HI(z) 
    mu    = kpara_for_int / knorm
    f     = growth_factor(z)
    bias_term = bias + f * kpara_for_int**2/knorm**2

!    cross_over_r_int = kappa_kernel * cos(in_cos) * bias_term * Plin
    cross_over_r_int = kappa_kernel *  bias_term * Plin
!    write(*,*) kappa_kernel, cos(in_cos), bias_term, cross_over_r_int

  End Function Cross_over_r_int
  !======================================================

 !======================================================
  Function cross_bias_over_r_int(z)

    Real(DP)  :: a, eta0, eta, z_source, eta_star, Omegam, frac, kappa_kernel
    Real(DP)  :: knorm, Plin, delta_eta, in_cos, bias, mu, f, bias_term
    Real(DP)  :: cross_bias_over_r_int, z

    a        = 1. / (1. + z)
    eta0     = conftime(0.d0)
    eta      = eta0 - conftime(z)
    z_source = z_cmb
    eta_star = eta0 - conftime(z_source)
    Omegam   = Om0
    frac     = (eta_star - eta) / (eta_star * eta)
    kappa_kernel = eta / a * frac * 3.d0 / 2.d0 / invhub**2  * Omegam

    knorm    = sqrt(kpara_for_int**2 + kperp_for_int**2)
    Plin     = P_lin_cross(knorm,zmin_for_int,z)

    delta_eta = abs(eta_for_int - eta)
    in_cos    = kpara_for_int * delta_eta

    bias  = Bias_HI(z) 
    bias_term = bias 

    cross_bias_over_r_int = kappa_kernel * cos(in_cos) * bias_term * Plin
!    write(*,*) kappa_kernel, cos(in_cos), bias_term, cross_over_r_int

  End Function Cross_bias_over_r_int
  !======================================================

 !======================================================
  Function cross_rsd_over_r_int(z)

    Real(DP)  :: a, eta0, eta, z_source, eta_star, Omegam, frac, kappa_kernel
    Real(DP)  :: knorm, Plin, delta_eta, in_cos, bias, mu, f, bias_term
    Real(DP)  :: cross_rsd_over_r_int, z

    a        = 1. / (1. + z)
    eta0     = conftime(0.d0)
    eta      = eta0 - conftime(z)
    z_source = z_cmb
    eta_star = eta0 - conftime(z_source)
    Omegam   = Om0
    frac     = (eta_star - eta) / (eta_star * eta)
    kappa_kernel = eta / a * frac * 3.d0 / 2.d0 / invhub**2  * Omegam

    knorm    = sqrt(kpara_for_int**2 + kperp_for_int**2)
    Plin     = P_lin_cross(knorm,zmin_for_int,z)

    delta_eta = abs(eta_for_int - eta)
    in_cos    = kpara_for_int * delta_eta

    mu    = kpara_for_int / knorm
    f     = growth_factor(z)
    bias_term =  f * kpara_for_int**2/knorm**2

    cross_rsd_over_r_int = kappa_kernel * cos(in_cos) * bias_term * Plin
!    write(*,*) kappa_kernel, cos(in_cos), bias_term, cross_over_r_int

  End Function Cross_rsd_over_r_int
  !======================================================

 !======================================================
  Subroutine Compute_cl_kap()

    Integer  :: ibin, il
    Real(DP) :: zmin, zmax

    Do ibin = 1, nbins_freq
       zmin = 1.0!redshift_from_freq(freq_min + dfreq * ibin)
       zmax = 1.005!redshift_from_freq(freq_min + dfreq * (ibin-1))
       write(*,*) "zmin, zmax = ", freq_min + dfreq * ibin, zmin, zmax 
       Allocate(freq_arr(ibin)%cl_kap_arr(1:nl_arr))
       Do il = 1, nl_arr
          freq_arr(ibin)%cl_kap_arr(il) = cl_kap_func(l_arr(il),zmin,zmax)
!          write(*,*) l_arr(il), freq_arr(ibin)%cl_kap_arr(il)
       End Do
    End Do


  End Subroutine Compute_cl_kap
  !======================================================


 !======================================================
  Function cl_kap_func(ell,zmin,zmax)
    
    Real(DP) :: cl_kap_func, ell, zmin, zmax, rombint, tol_i
    Real(DP) :: eta, eta0, Temp
    Real(DP) :: kpara_min, kpara_max, log_kpara_min, log_kpara_max
    external :: rombint

    tol_i         = 1.d-4
    eta0          = conftime(0.d0)
    eta_for_int   = eta0 - conftime(zmin)

    kperp_for_int = ell / eta_for_int
    zmin_for_int  = zmin 
    zmax_for_int  = zmax 

    kpara_min = 2.0 * pi / eta_for_int * (1.d0 + zmin_for_int) 
!    kpara_max = 0.21 ! cf Bull+2015
    kpara_max = 0.15 * 1.d0 ! cf Bull+2015 h^-1 taken into account
    log_kpara_min = log10(kpara_min) * 1.d0
    log_kpara_max = log10(kpara_max) * 1.d0
!!wrong 
!!!    cl_kap_func = rombint(cl_hi_kap_int,log_kpara_min,log_kpara_max,tol_i) 

    cl_kap_func = cl_kap_func / eta_for_int **2 /pi
!    write(*,*) "in cl_cross", Temp, eta_for_int
    
  End Function cl_kap_func
  !======================================================


 !======================================================
  Function cl_kap_int(logkpara)

    Real(DP) :: logkpara, cl_kap_int, tol_i
    Real(DP) :: rombint
    external :: rombint

    tol_i = 1.d-4 
    kpara_for_int  = 10.d0**logkpara
    cl_kap_int  = rombint(kappa_over_r_int,zmin_for_int,zmax_for_int,tol_i)
    cl_kap_int  = cl_kap_int * kpara_for_int * log(10.)

!! extra kara because of the log integral

  End Function cl_kap_int
  !======================================================


  !======================================================
  Function kappa_over_r_int(z)

    Real(DP)  :: a, eta0, eta, z_source, eta_star, Omegam, frac, kappa_kernel
    Real(DP)  :: knorm, Plin, delta_eta, in_cos, bias, mu, f, bias_term
    Real(DP)  :: kappa_over_r_int, z

    a        = 1. / (1. + z)
    eta0     = conftime(0.d0)
    eta      = eta0 - conftime(z)
    z_source = z_cmb
    eta_star = eta0 - conftime(z_source)
    Omegam   = Om0
    frac     = (eta_star - eta) / (eta_star * eta)
    kappa_kernel = eta / a * frac * 3.d0 / 2.d0 / invhub**2  * Omegam

    knorm    = sqrt(kpara_for_int**2 + kperp_for_int**2)
    Plin     = P_lin_cross(knorm,zmin_for_int,z)

    delta_eta = abs(eta_for_int - eta)
    in_cos    = kpara_for_int * delta_eta

    kappa_over_r_int = kappa_kernel**2 * cos(in_cos) * Plin
!    write(*,*) kappa_kernel, cos(in_cos), bias_term, cross_over_r_int

  End Function Kappa_over_r_int
  !======================================================

  !======================================================
  Function growth_factor(z)
!!    Following Bull+2015 page 3

    Real(DP) :: growth_factor, z
    Real(DP) :: H_of_z, Omega_m, gamma 
    
    H_of_z = Hubble(z)
    Omega_m = Om0 * (1.0 + z)**3 * H02 / H_of_z**2
    gamma = 0.55
    growth_factor = Omega_m ** gamma 
!    write(*,*) "in growth factor", z, H_of_z, Om0, Omega_m, gamma, growth_factor

  End Function growth_factor
  !======================================================
















End Module clustering
