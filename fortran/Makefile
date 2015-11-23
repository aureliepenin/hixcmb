F90C    = gfortran
OPTIM   = -O3
DEBUG   = -g -fbacktrace -fbounds-check
FFLAGS  = ${OPTIM}

IFLAG   = -I/. 
LIB     = 

F90FLAGS = $(FFLAGS) $(IFLAG)

FC       = $(F90C) $(F90FLAGS)
OBJ_FI_CIB = nrtype.o nr.o nrutil.o cisi.o rkck.o rkqs.o gammln.o gcf.o gser.o gammp.o erf.o odeint.o chebpc.o chebev.o \
	chebft.o chder.o polint.o polcof.o gaussj.o covsrt.o  bessel_AK.o misc_tools2.o splin2.o \
	 eh99_tools.o cosmo_tools.o  integrate_D_tools.o  max_NDLpower.o power_spec_tools.o \
	initiate_all.o HI_stuffs.o clustering_limber.o clustering.o write_files.o
#compute_dj_dM.o clustering.o write_files.o other_routines.o

cross:$(OBJ_FI_CIB) main.o
	$(FC) -o $@ $(OBJ_FI_CIB) main.o $(LIB)

.f.o:
	$(F90C) $(F90FLAGS) -c $<

%.o: %.f90
	$(F90C) $(F90FLAGS) -c $*.f90 -o $*.o

clean:
	rm -f *.o *.mod *.d *.pc *.obj core
