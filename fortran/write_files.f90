Module write_files

  Use nrtype
  Use eh99_tools
  Use cosmo_tools
  Use power_spec_tools
  Use initiate_all
  Use HI_stuffs
  Use clustering_limber
  Use clustering

  Implicit None


  Public 
  
Contains

  !======================================================
  Subroutine write_cl_for_fisher_versus_cosmo()

    Integer, parameter :: np_fisher = 4
    Integer, parameter :: nfac      = 3
    Integer :: ipar, ifac
    Real(DP), Dimension(1:nfac)      :: fac_vec


!! you can check if 10% is fine for the derivative 
    fac_vec = [0.9, 1., 1.1]

    Do ipar = 1, np_fisher
       Do ifac = 1, nfac
          Call set_default_cosmology()
          if (ipar .eq. 1)  wm   = wm   * fac_vec(ifac)
          if (ipar .eq. 2)  wb   = wb   * fac_vec(ifac)
          if (ipar .eq. 3)  ns   = ns   * fac_vec(ifac)
          if (ipar .eq. 4)  sig8 = sig8 * fac_vec(ifac)
          Call initiate_conftime()
          Call initiate_dod0()
          Call initiate_sigmas()
          Call write_cl_for_fisher()
       End Do
    End Do
    
End Subroutine write_cl_for_fisher_versus_cosmo
  !======================================================


  !======================================================
  Subroutine write_cl_for_fisher()

    Real(DP) :: dfreq, freq_min, zmin, zmax
    Integer  :: ibin, nbin

    nbin     = 4
    dfreq    = 100 ! in MHz
    freq_min = 400

    Do ibin = 1, nbin
       zmin = redshift_from_freq(freq_min + dfreq * ibin)
       zmax = redshift_from_freq(freq_min + dfreq * (ibin-1))
       write(*,*) "zmin, zmax = ", zmin, zmax 
       Call compute_cl_hi_kappa_for_fisher(zmin,zmax)
       Call write_cls_for_cosmo(zmin,zmax)
       Deallocate(cl_hi_kappa_for_fisher)
    End Do

  End Subroutine write_cl_for_fisher
  !======================================================

  !======================================================
  Subroutine write_Plin()

    Integer   :: ikpara, ikperp, nk 
    Real(DP)  :: logkmin, logkmax, dlogk, z, kpara, kperp, knorm, Plin
    Real(DP)  :: delta_eta_for_int, in_cos
    Character(Len=120) :: file_Plin

    nk = 100 
    logkmin = -5.
    logkmax =  2. 
    dlogk = (logkmax - logkmin)/nk
    z = 1.0

    Write(file_Plin,"(2A)") Trim(dir_out), "Plin_kk.dat"
    write(*,*) file_Plin
    Open(unit=15,file=file_Plin)


    delta_eta_for_int = 40.
    Do ikpara = 1, nk
       kpara = 10.d0**(logkmin + (ikpara - 1.)*dlogk)
       in_cos = kpara * delta_eta_for_int
       Do ikperp = 1, nk
          kperp = 10.d0**(logkmin + (ikperp - 1.)*dlogk)
          knorm = sqrt(kpara**2 + kperp**2)
          Plin = P_dd_ln(knorm,z) * cos(in_cos)
          write(15,*) kpara, knorm, Plin
       End Do
    End Do


  End Subroutine write_Plin
  !======================================================



  !======================================================
  Subroutine write_HI_stuffs()

    Integer   :: iz 
    Real(DP)  :: z
    Character(Len=120) :: file_hi

    Write(file_HI,"(2A)") Trim(dir_out), 'HI_inputs.dat'
    write(*,*) Trim(file_HI)
    Open(unit=15,file=file_hi)
    Do iz = 1, nz_arr
       z = z_arr(iz)
       write(15,*) z , T_HI_mean(z), Omega_HI(z), Bias_HI(z), rho_HI(z)
    End Do
    Close(15)

  End Subroutine write_HI_stuffs
  !======================================================

  !======================================================
  Subroutine write_cls(zmin,zmax)

    Integer   :: il 
    Real(DP)  :: zmin, zmax, bOm
    Character(Len=120) :: file_cl

    bOm = (Bias_HI(zmin)*Omega_HI(zmin) + Bias_HI(zmax)*Omega_HI(zmax))/2.

    Write(file_cl,"(2A,F4.2,A,F4.2,A)") Trim(dir_out), 'cls_zmin=', zmin, &
         '_zmax=', zmax, '.dat'
    write(*,*) Trim(file_cl)
    Open(unit=10,file=file_cl)
    write(10,*) bOm
    Do il = 1, nl_arr 
       write(10,*) l_arr(il), cl_hi_kappa_for_fisher(il) 
    End Do
    Close(10)

  End Subroutine write_cls
  !======================================================


  !======================================================
  Subroutine write_cls_for_cosmo(zmin,zmax)

    Integer   :: il 
    Real(DP)  :: zmin, zmax, bOm
    Character(Len=120) :: file_cl

    bOm = (Bias_HI(zmin)*Omega_HI(zmin) + Bias_HI(zmax)*Omega_HI(zmax))/2.

    Write(file_cl,"(2A,F4.2,A,F4.2,A,F4.2,A,F4.2,A,F4.2,A,F4.2,A)") &
         Trim(dir_out), 'cls_zmin=', zmin, &
         '_zmax=', zmax, &
         '_wm=', wm, '_wb=', wb, '_ns=', ns, '_sig8=', sig8, '.dat'
    write(*,*) Trim(file_cl)
    Open(unit=10,file=file_cl)
    write(10,*) bOm
    Do il = 1, nl_arr 
       write(10,*) l_arr(il), cl_hi_kappa_for_fisher(il) 
    End Do
    Close(10)

  End Subroutine write_cls_for_cosmo
  !======================================================


  !======================================================
  Subroutine write_cl_kappakappa()

    Integer   :: il 
    Real(DP)  :: zmin, zmax
    Character(Len=120) :: file_cl

    zmin = 0.07d0
    zmax = 9.d0
    Call Compute_cl_kappa_kappa(zmin,zmax)
 
    Write(file_cl,"(2A,F4.2,A,F4.2,A)") Trim(dir_out), 'cl_kappakappa_zmin=', zmin, &
         '_zmax=', zmax, '.dat'
    write(*,*) Trim(file_cl)
    Open(unit=10,file=file_cl)

    Do il = 1, nl_arr 
       write(10,*) l_arr(il),  cl_kappa_kappa(il) 
    End Do
    Close(10)
    Deallocate(cl_kappa_kappa)
  End Subroutine write_cl_kappakappa
  !======================================================

  !======================================================
  Subroutine write_PHI()
!! output P_HI in mK^2 (Mpc h^-1)^3 and k in (Mpc h^-1)^-1

    Integer   :: ik, nk
    Real(DP)  :: z, k, log_kmin, log_kmax, dlogk
    Character(Len=120) :: file_PHI

    log_kmin = -3.0
    log_kmax =  3.0
    nk = 100
    dlogk = (log_kmax - log_kmin) / nk
    z = 0.1    

    Write(file_PHI,"(2A,F3.1,A)") Trim(dir_out), 'P_HI_z=', z, 'PourtsidouFig1.dat'
    write(*,*) Trim(file_PHI)
    Open(unit=10,file=file_PHI)

    Do ik = 1, nk
       k  = 10.d0**(log_kmin + dlogk * (ik - 1))
       write(10,*), k, P_HI(z,k)
    End Do
    Close(10)

  End Subroutine write_PHI
  !======================================================

  
  !======================================================
  Subroutine write_over_r_integrand()

    Integer   :: iz, nz, ik, nk
    Real(DP)  :: eta0, eta, ell, dz, func, z, dk
    Real(DP)  :: z_source, eta_star, Omegam, frac, kappa_kernel, a
    Real(DP)  :: delta_eta, in_cos, bias, f, growth_factor, knorm, bias_term
    Character(Len=120) :: file


    zmin_for_int  = 1.d0
    eta0          = conftime(0.d0)
    eta_for_int   = eta0 - conftime(zmin_for_int)
    kpara_for_int = 0.20 
    ell           = 1000.d0
    dz = 0.00001 
    nz = 50
    dk = 0.05
    nk = 3

    write(file,"(2A)") Trim(dir_out), "cross_integrand_over_r.dat"
    write(*,*) Trim(file)
    Open(unit=10,file=file)
    Do ik = 1, nk
       kpara_for_int = 5.d-3 + dk * (ik - 1)
       Do iz = 1, nz 
          z    = zmin_for_int + dz * (iz - 1)
          eta  = eta0 - conftime(z)
          kperp_for_int = ell / eta
          eta_for_int = eta 
          func = cross_over_r_int(z)

          a        = 1. / (1. + z)
          z_source = z_cmb
          eta_star = eta0 - conftime(z_source)
          Omegam   = Om0
          frac     = (eta_star - eta) / (eta_star * eta)
          kappa_kernel = eta / a * frac * 3.d0 / 2.d0 / invhub**2  * Omegam

          delta_eta = abs(eta_for_int - eta)
          in_cos    = kpara_for_int * delta_eta

          knorm    = sqrt(kpara_for_int**2 + kperp_for_int**2)
          bias  = Bias_HI(z) 
          f     = 1.!growth_factor(z)
!          bias_term =  f * kpara_for_int**2/knorm**2
          bias_term = bias + f * kpara_for_int**2/knorm**2

          write(10,*) kpara_for_int, z, func, in_cos, kappa_kernel, bias_term
          write(*,*) kperp_for_int, ell, knorm
       End Do
    End Do
    Close(10)

  End Subroutine write_over_r_integrand
  !======================================================


  !======================================================
  !Subroutine write_full_integrand()



  !End Subroutine write_full_integrand
  !======================================================



End Module write_files
