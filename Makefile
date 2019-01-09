F90=mpif90
#F90=/usr/local/bin/gfortran
#F90=/opt/local/bin/gfortran-mp-7
FFLAGS=-g -cpp -fbacktrace -fcheck=all

PRODUCT=mcrad.x 

SRCS := mod_mympi.f90 mod_optic.f90 mod_ioset.f90 mod_utest.f90
SRCS += mcr_main.f90
OBJS=$(SRCS:.f90=.o)

$(PRODUCT): $(OBJS)
	$(F90) $(FFLAGS) $^ -o $@
%.o: %.f90
	$(F90) $(FFLAGS) -c $< -o $@

clean: 
	rm -rf *.o $(PRODUCT) *.mod *.dSYM 
clear:
	rm -f *.o *.mod
