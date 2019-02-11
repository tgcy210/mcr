F90=mpif90
#F90=/usr/local/bin/gfortran
#F90=/opt/local/bin/gfortran-mp-7
FFLAGS := -g -cpp -fbacktrace -fcheck=all

PRODUCT=mcrad.x 

SRCS := m01_mympi.f90 m02_optic.f90 m03_ioset.f90 m11_utest.f90
SRCS += mcr_main.f90
OBJS=$(SRCS:.f90=.o)

all: $(PRODUCT)

mpi: FFLAGS += -D_USE_MPI
mpi: clean
mpi: $(PRODUCT)

$(PRODUCT): $(OBJS)
	$(F90) $(FFLAGS) $^ -o $@
%.o: %.f90
	$(F90) $(FFLAGS) -c $< -o $@

clean: 
	rm -rf *.o $(PRODUCT) *.mod *.dSYM 
clear:
	rm -f *.o *.mod
