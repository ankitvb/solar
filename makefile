## Makefile
## -----------------------------------------------------
## MPI Version of the Spherical Acoustic Sun Simulator.
## Copyright 2006, Shravan Hanasoge
##
## Hansen Experimental Physics Laboratory
## 455 Via Palou way, Stanford
## CA 94305, USA
## Email: shravan@stanford.edu
## -----------------------------------------------------
##

# add 
OBJS=   cube.o init.o derivative.o interpolation.o metric.o output.o transform.o \
        init_mpi.o grad.o filter.o timestep.o auxiliary.o dirs.o 
OBJ_IC= ic_acoustic_wave.o #ic_gaussian_pulse.o
OBJ_RHS= compute_rhs_acoustic.o #compute_rhs_gaussian.o
OBJ_BC= bc_acoustic_wave.o sponge.o

FC= mpixlf90
CC= mpicc
FFLAGS=  #-C -g -traceback
LIBS=	-L/home/ankitb/local/fftw3/lib -lfftw3 -lm -L/home/ankitb/local/cfitsio/lib -lcfitsio 
INCLUDE= -I/home/ankitb/local/fftw3/include fftw3.f 

COMMAND= glass

$(COMMAND): $(OBJS) $(OBJ_IC) $(OBJ_RHS) $(OBJ_BC)
	$(FC) -o $(COMMAND) $(OBJS) $(OBJ_IC) $(OBJ_RHS) $(OBJ_BC) $(LIBS) 

###	rm core.*
%.o : %.f90
	$(FC) $(FFLAGS) -c $< 
%.o : %.c
	$(CC) -c $<

clean:
	rm *.o *.mod $(COMMAND)

# add to cube.o
cube.o:	header init.o derivative.o interpolation.o metric.o output.o transform.o grad.o filter.o auxiliary.o
grad.o:	header init.o derivative.o	interpolation.o	transform.o	output.o
derivative.o: header init.o
interpolation.o: header init.o
transform.o: header init.o
compute_rhs_acoustic.o:	header grad.o init.o 
compute_rhs_gaussian.o:  header grad.o init.o
timestep.o: header init.o filter.o
ic_acoustic_wave.o: header init.o 
ic_gaussian_pulse.o: header init.o
bc_acoustic_wave.o: header init.o
metric.o: header init.o
filter.o: header interpolation.o output.o init.o
init_mpi.o: header 
init.o: header
auxiliary.o: header init.o
sponge.o: header init.o
output.o: init.o
dirs.o: dirs.h
