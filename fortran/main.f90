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
  Real(DP) :: t1, t2, redshift

  Call cpu_time(t1)
  Call initiate_cosmo()
  Call initiate_parameters()
  Call initiate_arrays()
  write(*,*) "Initialisation done"

!  Call Compute_cl_hi_kap()
!  Call Compute_cl_kap()
!  Call compute_cl_hi()
!  Call write_all_cls()

  zmin = 1.
  zmax = 1.005
  
  Call  Compute_cl_hi_limber(zmin,zmax)
  Call  write_cls_check(zmin,zmax)

  Call cpu_time(t2)
  write(*,*), "time = ", t2 - t1 
  write(*,*) "That's all folks! "

End Program Main
