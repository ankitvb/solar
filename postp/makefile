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
OBJS=	cube.o init.o bspline90_22.o interpolation.o output.o \
	metric.o transform.o 	

##OBJS=	sphturb.o bspline90_22.o

FC= 		ifort
FFLAGS= -O3 -r8 -zero
LIBS=	-L/data/seismo1/hanasoge/cfitsio -lcfitsio

COMMAND=	tester


$(COMMAND): $(OBJS) 
	$(FC) -I $(INCLUDE) $(FFLAGS) -o $(COMMAND) $(OBJS) $(LIBS) 

%.o : %.f90
	$(FC) $(FFLAGS) -c $< 

clean:
	rm *.o *.mod $(COMMAND)

##sphturb.o:	bspline90_22.o
cube.o:	header init.o interpolation.o metric.o output.o transform.o bspline90_22.o 
interpolation.o: header	output.o init.o
transform.o: header init.o
metric.o: header init.o
output.o: header init.o
init.o: header
