Module HI_stuffs


  Use nrtype
  Use eh99_tools
  Use cosmo_tools
  Use power_spec_tools
  Use initiate_all


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

    Real(DP)  :: Omega_HI, z

    Omega_HI = 4.0d-4 * (1. + z)**0.6

  End Function Omega_HI
  !======================================================

  !======================================================
  Function bias_HI(z)
!! From Pourtsidou+2015 Eq 13
!! From Crighton+2015 

    Real(DP)  :: bias_HI, z

    Bias_HI = 1.0
!    write(*,*) "===================================="
!    write(*,*) "BE CAREFUL WRONG VALUE OF BIAS HI"
!    write(*,*) "===================================="

  End Function Bias_HI
  !======================================================



End Module HI_stuffs
