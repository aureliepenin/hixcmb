Module Initiate_all


  Use nrtype
  Use eh99_tools
  Use cosmo_tools
  Use power_spec_tools

  Integer, parameter :: nl_arr = 100
  Integer, parameter :: nz_arr = 100

  Real(DP) :: lmin_arr, lmax_arr, zmin_arr, zmax_arr
  Real(DP), Dimension(:)  , Allocatable  :: z_arr, l_arr, cl_arr, cl_hi_kappa, cl_hi_hi, cl_kappa_kappa, cl_hi_kappa_for_fisher
  Real(DP), Dimension(:)  , Allocatable  :: cl_cross_arr, cl_hi_arr
  Real(DP), Dimension(:,:), Allocatable  :: k_arr, Plin_arr

  Character(Len=50) :: dir_out='../Output/'


  Public 
  
Contains
  !======================================================
  Subroutine initiate_cosmo()
    
    Call set_default_cosmology()
    Call initiate_conftime()
    Call initiate_dod0()
    Call initiate_sigmas()

 End Subroutine Initiate_cosmo
  !======================================================


 !======================================================
 Subroutine initiate_parameters()
   
   lmin_arr = 30.d0
   lmax_arr = 200.d0

   zmin_arr = 0.07
   zmax_arr = 5.
   
!! Use Tinker mass function (Cf power_spec_tools.f90)
   mass_function = 4

 End Subroutine initiate_parameters
  !======================================================


 !======================================================
 Subroutine Initiate_arrays()

   Integer   :: iz, il, ik
   Real(DP)  :: log_zmin, log_zmax, dlogz, dz_arr
   Real(DP)  :: log_lmin, log_lmax, dlogl, tmp

   !! redshift array
   Allocate(z_arr(1:nz_arr))
   log_zmin = log10(zmin_arr)
   log_zmax = log10(zmax_arr)
   dlogz = (log_zmax - log_zmin)/nz_arr

    Do iz = 1, nz_arr
       dz_arr = log_zmin + dlogz*(iz - 1.)
       z_arr(iz) = 10**dz_arr
!       write(*,*), iz, z_arr(iz)
    End do
    


  !! multipole arrays
    Allocate(l_arr(1:nl_arr))
    log_lmin = log10(lmin_arr)
    log_lmax = log10(lmax_arr)

    dlogl = (log_lmax - log_lmin)/nl_arr
    Do il = 1, nl_arr
       tmp = log_lmin + (il - 1)*dlogl
       l_arr(il) = 10**tmp
    End Do

   

!!    !! k arrays
!!    !! The output k is in (Mpc/h)^(-1) => divided by h
!!    !! Plin array also 
!!    !! The output Plin is in (Mpc/h)^3 => multiplied by h**3
!!    !! Also geometry term for the Cl
!!
!!    Allocate(k_perp_arr(1:nz_arr,1:nl_arr))
!!    Allocate(Plin_arr(1:nz_arr,1:nl_arr))
!!
!!    !! geom = dr/dz * (a/r)**2 => in (Mpc/h)**(-1) from the cosmo routines => / hub
!!
!!    Do iz = 1, nz_arr 
!!       a = 1. / ( 1 + z_arr(iz))
!!       detadz = dconftime(z_arr(iz))
!!       eta0   = conftime(0.d0)
!!       eta    = eta0 - conftime(z_arr(iz))
!!
!!       Do il = 1, nl_arr
!!          k_arr(iz,il)    = (l_arr(il) + 0.5 ) / eta / hub
!!          Plin_arr(iz,il) = P_dd_ln(k_arr(iz,il)*hub,z_arr(iz)) * hub**3
!!       End Do
!!    End Do
!!

 End Subroutine Initiate_arrays

 !======================================================













End Module Initiate_all
