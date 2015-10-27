Program Main 

  Use nrtype
  Use eh99_tools
  Use cosmo_tools
  Use power_spec_tools
  Use initiate_all
  Use HI_stuffs
  Use clustering
  Use write_files

  Real(DP) :: zmin, zmax, zz, plouf1, plouf2


  write(*,*) "starting"
  Call cpu_time(t1)
  write(*,*) "cosmo"
  Call initiate_cosmo()
  write(*,*) "initiate_parameters"
  Call initiate_parameters()
  write(*,*) "initiate_arrays"
  Call initiate_arrays()
  write(*,*) "write HI stuffs"
  zz = 0.8
  write(*,*) 'at z = 0.8 = '
  plouf1 =  bias_HI(zz)
  plouf2  = Omega_HI(zz)
  write(*,*) plouf2, plouf1, plouf2*plouf1, 0.62d-4/(plouf2*plouf1)
stop
  Call write_HI_stuffs()
  Call write_PHI()
  zmin = 1.
  zmax = 1.3
  write(*,*) "compute_cl HIHI"
  Call compute_cl_hi_hi(zmin,zmax)
  write(*,*) "compute_cl HI Kappa"
  Call compute_cl_hi_kappa(zmin,zmax)
  write(*,*) "compute_cl Kappa Kappa"
  Call Compute_cl_kappa_kappa(zmin,zmax)
  Call write_cls()

  Call cpu_time(t2)
  write(*,*), "time = ", t2 - t1 
  write(*,*) "That's all folks! "


End Program Main
