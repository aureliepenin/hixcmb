Program Main 

  Use nrtype
  Use eh99_tools
  Use cosmo_tools
  Use power_spec_tools
  Use initiate_all
  Use HI_stuffs
  Use clustering
  Use write_files

  Implicit None 
  Real(DP) :: zmin, zmax, bias_main, Omega_main, rho_main, zz, pp1, pp2, kk 
  Real(DP) :: t1, t2, freq, redshift

  write(*,*) "starting"
  Call cpu_time(t1)
  write(*,*) "cosmo"
  Call initiate_cosmo()
  write(*,*) "initiate_parameters"
  Call initiate_parameters()
  write(*,*) "initiate_arrays"
  Call initiate_arrays()
  write(*,*) "write HI stuffs"
  Call  write_cl_for_fisher()
  Call write_cl_for_fisher_versus_cosmo()
stop
!  Call write_HI_stuffs()
!stop
  zz = 0.8
  kk = 0.01
  pp1 = P_dd_ln(kk, zz)
  pp2 = growth_factor(zz)
  write(*,*) zz, pp2 * sig8
!Call write_Plin()
stop
!  write(*,*) 6.2d-4/(bias_main*omega_main)
!  write(*,*) rho_bar_m(0.d0) * 4.86d-4 / rho_HI(0.d0)

  Call write_PHI()
!  Call  write_cl_kappakappa()
!  Call  write_cl_for_fisher()

!  write(*,*) "compute_cl HIHI"
!  zmin = 1.
!  zmax = 1.005

  zmin = 1.
  zmax = 2.

  write(*,*) "before limber one"
  Call write_over_r_integrand()
!stop
 Call compute_cl_hi_hi(zmin,zmax)
!  write(*,*) "before new one!"
 Call compute_cl_hi(zmin,zmax)
!  write(*,*) "before limber one"
  Call compute_cl_hi_kappa(zmin,zmax)
!  write(*,*) "before new one!"
  Call compute_cl_cross(zmin, zmax)

!  Call Compute_cl_kappa_kappa(zmin,zmax)
  Call write_cls(zmin,zmax)
  Call cpu_time(t2)
  write(*,*), "time = ", t2 - t1 
  write(*,*) "That's all folks! "

End Program Main
