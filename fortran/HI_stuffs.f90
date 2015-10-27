Module HI_stuffs


  Use nrtype
  Use eh99_tools
  Use cosmo_tools
  Use power_spec_tools
  Use initiate_all


  Real(DP)  :: z_buf 

  Public 
  
Contains
  !======================================================
  Function T_HI_mean(z)
!! From Pourtsidou+2015 Eq. 12 
!! and Battye+2103 Eq. 9 (derivation)
!! in mK (Mpc h^-1)^-1, yes because there is a hub....
!! in this formula h has units, km/s/Mpc
!! otherwise, unit issue in Cl

    Real(DP) :: z, T_HI_mean
    Real(DP) :: num, den 

    num = (1.0 + z)**2
    den = hubble(z)/H0

    T_HI_mean = 180.d0 * Omega_HI(z) * num / den * hub 

  End Function T_HI_mean
  !======================================================
    
  !======================================================
  Function Omega_HI(z)
!! From Pourtsidou+2015 Eq 13
!! From Crighton+2015 
!! (1+z)^-3 because should be in proper units
!! normalisation Omega_HI0 = 6.5d-4 from Bull+2015

    Real(DP)  :: Omega_HI, z, rho_crit0, norma
    write(*,*) "Attention a la normalisation de Omega_HI"
!!    Omega_HI = 4.0d-4 * (1. + z)**0.6
    norma = 37.157521791021466
    rho_crit0 = rho_bar_m(z)
    Omega_HI = (1. + z)**(-3.) * rho_HI(z) / rho_crit0 * norma

  End Function Omega_HI
  !======================================================

  !======================================================
  Function bias_HI(z)
!! From Pourtsidou+2015 Eq 13
!! From Crighton+2015 

    Real(DP)  :: bias_HI, z_buf, rombint, tol_i, rho, z
    Real(DP)  :: logMmin, logMmax
    external :: rombint

    tol_i   = 1.d-4
    z_buf   = z
    rho     = rho_HI(z)
    logMmin = log10(Mmin_for_hosting_HI(z))
    logMmax = log10(Mmax_for_hosting_HI(z))
    bias_HI = rombint(bias_HI_int, logMmin, logMmax,tol_i)
    bias_HI = bias_HI / rho
!    write(*,*) "z, bias HI, rhoHI = ", z, bias_HI, rho


  End Function Bias_HI
  !======================================================

  !======================================================
  Function bias_HI_int(logM)
!! in proper units
    Real(DP) :: M, dndm, bias_HI_int, logM

    M         = 10.d0**logM
    dndm      = dndm_f(M,z_buf)! * (1. + z_buf)**(-3.)
    bias_halo = bias(M,z_buf)
    MHI       = MHI_of_M(M)

    bias_HI_int = log(10.) * M * dndm * MHI * bias_halo
!! extra M for log integration

  End Function bias_HI_int
  !======================================================

  !======================================================
  Function rho_HI(z)

    Real(DP) :: rombint, rho_HI, z, tol_i, logMmin, logMmax
    external :: rombint

    tol_i = 1.d-4 
    z_buf = z
    logMmin = log10(Mmin_for_hosting_HI(z))
    logMmax = log10(Mmax_for_hosting_HI(z))
    rho_HI = rombint(rho_HI_int, logMmin, logMmax, tol_i)

  End Function rho_HI
  !======================================================

  !======================================================
  Function rho_HI_int(logM)
!! in proper units

    Real(DP)  :: logM, M, MHI, rho_HI_int

    M         = 10.d0**logM
    dndm      = dndm_f(M,z_buf)! * (1.+ z)**(-3.)
    MHI       = MHI_of_M(M)

    rho_HI_int  = log(10.) * M * dndm * MHI 
!! extra M for log integration

  End Function rho_HI_int
  !======================================================

  !======================================================
  Function Mmin_for_hosting_HI(z)
!! Cf Appendix B Bull+2015
!! Bagla+2010

    Real(DP) ::  z, Mmin_for_hosting_HI, vcirc

    vcirc = 30. !!km/s
    Mmin_for_hosting_HI = M_from_vcirc(vcirc,z)

  End Function Mmin_for_hosting_HI
  !======================================================

 !======================================================
  Function Mmax_for_hosting_HI(z)
!! Cf Appendix B Bull+2015
!! Bagla+2010

    Real(DP) ::  z, Mmax_for_hosting_HI, vcirc

    vcirc = 200. !!km/s
    Mmax_for_hosting_HI = M_from_vcirc(vcirc,z)

  End Function Mmax_for_hosting_HI
  !======================================================

 !======================================================
  Function M_from_vcirc(vcirc,z)
!! Cf Appendix B Bull+2015
!! Bagla+2010
!! v in km s^-1

    Real(DP) :: M_from_vcirc, vcirc, z

    M_from_vcirc = 1.d10 * (vcirc / 30. / sqrt(1. + z))**3.

  End Function M_from_vcirc
  !======================================================
  Function MHI_of_M(M)
!! See Appendix B of Bull+2015
!! normalisation with the constraint at z=0.8 of Switzer+2013
!! the normalisation to match 6.2d-4
!! M in Msol h^-1 

    Real(DP) :: M, alpha, norm, MHI_of_M

    norm = 2.9558233886481609
    alpha = 0.6
    MHI_of_M = M**alpha * norm
!    write(*,*) "in func MHI", M, MHI_of_M, alpha, norm

  End Function MHI_of_M
  !======================================================







End Module HI_stuffs
