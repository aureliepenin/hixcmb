Module clustering_limber

  Use nrtype
  Use eh99_tools
  Use cosmo_tools
  Use power_spec_tools
  Use initiate_all
  Use HI_stuffs

  Implicit None

  Real(DP) :: l_for_integrand
  Real(DP) ::  z_cmb   = 1089.


  Public 
  
Contains

  !======================================================
  Subroutine Compute_cl_hi_kappa_for_fisher(zmin,zmax)

    Integer  :: il
    Real(DP) :: rombint, tol_i, zmin, zmax, tmp, bOm
    external :: rombint

    tol_i = 1.d-4
    Allocate(cl_hi_kappa_for_fisher(1:nl_arr))
    bOm = (Bias_HI(zmin)*Omega_HI(zmin) + Bias_HI(zmax)*Omega_HI(zmax))/2.
!    write(*,*) "bOm = ", bOm
    cl_hi_kappa_for_fisher = 0.
    Do il = 1, nl_arr 
       l_for_integrand = l_arr(il)
       tmp = rombint(cl_hi_kappa_for_fisher_int,zmin,zmax,tol_i)
       cl_hi_kappa_for_fisher(il) = tmp * bOm
!       write(*,*) il, l_for_integrand, cl_hi_kappa(il)
    End Do

  End Subroutine Compute_cl_hi_kappa_for_fisher
  !======================================================

 
  !======================================================
  Function cl_hi_kappa_for_fisher_int(z)
    !! integration on z for the moment
    !! detadz       =>  Mpc h^{-1}
    !! Kernel_kappa => (Mpc h^-1)^-2 
    !! Plin         => (Mpc h^-1)^3 
    !! Kernel_HI    => (Mpc h^-1)^-2  mK
    !! Cl           => mK



    Real(DP) :: cl_hi_kappa_for_fisher_int, z, Plin, k_buf, detadz, eta, eta0
    Real(DP) :: w_kappa, w_hi

    w_kappa = kernel_kappa(z,l_for_integrand)
    w_hi    = kernel_hi_for_fisher(z)
    eta0    = conftime(0.d0)
    eta     = eta0 - conftime(z)
    k_buf   = (l_for_integrand + 0.5) / eta 
    Plin    = P_dd_ln(k_buf,z)
    detadz  = dconftime(z)

    cl_hi_kappa_for_fisher_int = detadz * w_hi * w_kappa * Plin

  End Function cl_hi_kappa_for_fisher_int
  !======================================================


  !======================================================
  Subroutine Compute_cl_hi_kappa(zmin,zmax)

    Integer  :: il
    Real(DP) :: rombint, tol_i, zmin, zmax
    external :: rombint

    tol_i = 1.d-4
    Allocate(cl_hi_kappa(1:nl_arr))
    cl_hi_kappa = 0.
    Do il = 1, nl_arr 
       l_for_integrand = l_arr(il)
       cl_hi_kappa(il) = rombint(cl_hi_kappa_int,zmin,zmax,tol_i)
       write(*,*) "old = ", l_arr(il), cl_hi_kappa(il), nl_Arr
    End Do

  End Subroutine Compute_cl_hi_kappa
  !======================================================

 
  !======================================================
  Function cl_hi_kappa_int(z)
    !! integration on z for the moment
    !! detadz       =>  Mpc h^{-1}
    !! Kernel_kappa => (Mpc h^-1)^-2 
    !! Plin         => (Mpc h^-1)^3 
    !! Kernel_HI    => (Mpc h^-1)^-2  mK
    !! Cl           => mK


    Real(DP) :: cl_hi_kappa_int, z, Plin, k_buf, detadz, eta, eta0
    Real(DP) :: w_kappa, w_hi
    
    w_kappa = kernel_kappa(z,l_for_integrand)
    w_hi    = kernel_hi(z)
    eta0    = conftime(0.d0)
    eta     = eta0 - conftime(z)
    k_buf   = (l_for_integrand + 0.5) / eta 
    Plin    = P_dd_ln(k_buf,z)
    detadz  = dconftime(z)

    cl_hi_kappa_int = detadz * w_hi * w_kappa * Plin

  End Function cl_hi_kappa_int
  !======================================================

 !======================================================
  Subroutine Compute_cl_kappa_kappa(zmin,zmax)
  !! cl_kk        => nothing
    
    Integer  :: il
    Real(DP) :: rombint, tol_i, zmin, zmax
    external :: rombint

    tol_i = 1.d-4
    Allocate(cl_kappa_kappa(1:nl_arr))
    cl_kappa_kappa = 0.
    Do il = 1, nl_arr 
       l_for_integrand = l_arr(il)
       cl_kappa_kappa(il) = rombint(cl_kappa_kappa_int,zmin,zmax,tol_i)
!      write(*,*) "il = ", il, l_arr(il), cl_kappa_kappa(il)
    End Do

  End Subroutine Compute_cl_kappa_kappa
  !======================================================

  !======================================================
  Function cl_kappa_kappa_int(z)
    !! detadz       =>  Mpc h^{-1}
    !! Kernel_kappa => (Mpc h^-1)^-2 
    !! Plin         => (Mpc h^-1)^3 
    !! cl_kk        => nothing

!integration on z for the moment

    Real(DP) :: cl_kappa_kappa_int, z, Plin, k_buf, detadz
    Real(DP) :: eta0, eta, w_kappa
    
    w_kappa = kernel_kappa(z,l_for_integrand)
    eta0    = conftime(0.d0)
    eta     = eta0 - conftime(z)
    k_buf   = (l_for_integrand + 0.5) / eta 
    Plin    = P_dd_ln(k_buf,z)
    detadz  = dconftime(z)

    cl_kappa_kappa_int = detadz * w_kappa**2 * Plin

  End Function cl_kappa_kappa_int
  !======================================================


  !======================================================
  Function Kernel_kappa(z,ell)
!!!CHECK LES C et autres trucs 

!! eta     => Mpc h^{-1}
!! invhub  => Mpc h^-1
!! Kernel_kappa => (Mpc h^-1)^-2 * Mpc h^-1 * (Mpc h^-1)^-1
!! Kernel_kappa => (Mpc h^-1)^-2 
    

    Real(DP) :: Kernel_kappa, z, Omegam, a, ell, z_source
    Real(DP) :: detadz, eta0, eta, eta_star, frac

    Omegam   = Om0
    a        = 1. / (1. + z)
    eta0     = conftime(0.d0)
    eta      = eta0 - conftime(z)
    z_source = z_cmb
    eta_star = eta0 - conftime(z_source)

    frac = (eta_star - eta) / (eta_star * eta)
    Kernel_kappa = 3.0 / 2. / invhub**2 * Omegam * (eta / a) * frac
!    write(*,*) 'eta_star = ', eta_star


  End Function Kernel_kappa
  !======================================================

  !======================================================
  Function Kernel_HI(z)
!! eta  => Mpc h^{-1}
!! T    => mK (Mpc h^-1)^-1
!! Kernel_HI => (Mpc h^-1)^-2  mK

    Real(DP) :: z, Kernel_HI, a, eta0, eta, temp, bias

    a      = 1. / (1. + z)
    eta0   = conftime(0.d0)
    eta    = eta0 - conftime(z)

    temp    = T_HI_mean(z)
    bias    = bias_HI(z)

    Kernel_HI = (a/eta) * temp * bias

  End Function Kernel_HI
  !======================================================

 !======================================================
  Function Kernel_HI_for_fisher(z)
!! eta  => Mpc h^{-1}
!! T    => mK (Mpc h^-1)^-1
!! Kernel_HI => (Mpc h^-1)^-2  mK

    Real(DP) :: z, Kernel_HI_for_fisher, a, eta0, eta, temp, num, den

    a      = 1. / (1. + z)
    eta0   = conftime(0.d0)
    eta    = eta0 - conftime(z)

!!temp = temp / OmegaHI
    num = (1.0 + z)**2
    den = hubble(z)/H0
    temp    = 180.d0 * num / den * hub 

    Kernel_HI_for_fisher = (a/eta) * temp

  End Function Kernel_HI_for_fisher
  !======================================================


  !======================================================
  Subroutine Compute_cl_hi_hi(zmin,zmax)

    Integer  :: il
    Real(DP) :: zmin, zmax, rombint, tol_i, tmp
    External :: rombint

    Allocate(cl_hi_hi(1:nl_arr))
    cl_hi_hi = 0.
    tol_i = 1.d-4
    Do il = 1, nl_arr 
       l_for_integrand = l_arr(il)
       cl_hi_hi(il) = rombint(cl_hi_hi_int,zmin,zmax,tol_i)
       write(*,*) "old = ", l_arr(il), cl_hi_hi(il)
    End Do


  End Subroutine Compute_cl_hi_hi
    !======================================================

  !======================================================
  Function cl_hi_hi_int(z)
!! Kernel_HI => (Mpc h^-1)^-2  mK
!! detadz    => Mpc h^-1
!! Cl HIHI   => 

    Real(DP)  :: cl_hi_hi_int,  w_HI, geo, Plin, z, a, detadz
    Real(DP)  :: k_buf, eta, eta0

    a      = 1. / (1. + z)
    detadz = dconftime(z)
    eta0   = conftime(0.d0)
    eta    = eta0 - conftime(z)
    w_HI   = kernel_HI(z)
    k_buf  = (l_for_integrand + 0.5) / eta 
    Plin   = P_dd_ln(k_buf,z)

    cl_hi_hi_int = detadz * w_HI**2 * Plin 

  End Function Cl_hi_hi_int
  !======================================================

  !======================================================
  Function P_HI(z,k)
!! Input k in (Mpc h^-1)^-1
!! Plin in (Mpc h^-1)^3
!! Temp in mK here (divided by hub)
!! P_HI in mK^2 (Mpc h^-1)^3

    Real(DP)  :: z, k, Plin, P_HI, temp, bias

    temp   = T_HI_mean(z) / hub
    bias   = bias_HI(z)
    Plin   = P_dd_ln(k,z)

    P_HI = bias**2 * temp**2 * Plin 

  End Function P_HI
  !======================================================


  !======================================================
  Function geom(z)
!! in h^-1 pour le moment
    Real(DP)  :: a, z, geom, detadz, eta0, eta

    a      = 1. / (1. + z)
    detadz = dconftime(z)
    eta0   = conftime(0.d0)
    eta    = eta0 - conftime(z)
    geom   =  detadz * (a / eta)**2 

  End Function geom
  !======================================================









End Module clustering_limber
