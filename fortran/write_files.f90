Module write_files

  Use nrtype
  Use eh99_tools
  Use cosmo_tools
  Use power_spec_tools
  Use initiate_all
  Use HI_stuffs
  Use clustering

  Implicit None



  Public 
  
Contains

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
       write(15,*) z , T_HI_mean(z), Omega_HI(z), Bias_HI(z)
    End Do
    Close(15)

  End Subroutine write_HI_stuffs
  !======================================================

  !======================================================
  Subroutine write_cls()

    Integer   :: il 
!    Real(DP)  :: 
    Character(Len=120) :: file_cl

    Write(file_cl,"(2A)") Trim(dir_out), 'cls.dat'
    write(*,*) Trim(file_cl)
    Open(unit=10,file=file_cl)

    Do il = 1, nl_arr 
       write(10,*) l_arr(il), cl_hi_kappa(il), cl_hi_hi(il), cl_kappa_kappa(il) 
    End Do
    Close(10)

  End Subroutine write_cls
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

    Write(file_PHI,"(2A,F3.1,A)") Trim(dir_out), 'P_HI_z=', z, '.dat'
    write(*,*) Trim(file_PHI)
    Open(unit=10,file=file_PHI)

    Do ik = 1, nk
       k  = 10.d0**(log_kmin + dlogk * (ik - 1))
       write(10,*), k, P_HI(z,k)
    End Do
    Close(10)

  End Subroutine write_PHI
  !======================================================

  
End Module write_files
