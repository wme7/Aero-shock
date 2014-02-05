# Specify search path for subroutines that perform 1D ppmlr hydrodynamics
VPATH = ../PPMLR 
#
# System-dependent parameters 
#
# F90           FORTRAN 90 Compiler
# LDR           Program to load the objects into an executable, typically the same as F90
# LDFLAGS       Flags to pass to the compiler during linking
# LIBS          A list of libraries to link with the executable, normally only netCDF
#
#  the following is set for gfortran compiler  http://gcc.gnu.org/wiki/GFortran
#
F90 = gfortran 
FFLAGS = -c 
LDR=	gfortran              
LDRFLAGS=              
LIBS=  
#
#
# List of objects to build multidimensional VH-1 hydro code:

VHOBJS = vh1mods.o vhone.o \
         ppmlr.o forces.o flatten.o evolve.o remap.o \
         states.o boundary.o volume.o riemann.o parabola.o
         

# Lines from here on down should not need to be changed.  They are the
# actual rules which make uses to build the executable
#
.SUFFIXES: .f90 .o

.f90.o:
	$(F90) $(FFLAGS) $< 

vhzero:	$(VHOBJS)
	$(LDR) $(LDRFLAGS) -o vh1-starter $(VHOBJS) $(LIBS); mv vh1-starter ../../
#
#
clean:
	rm -f *.o *.l

clobber:	clean
	rm -f ../../vh1-starter
	rm -f *.mod


# Dependencies for the object files
vhone.o:    vhone.f90    global.mod sweeps.mod sweepsize.mod 

evolve.o:   evolve.f90   global.mod sweeps.mod sweepsize.mod
flatten.o:  flatten.f90  global.mod sweeps.mod sweepsize.mod
forces.o:   forces.f90   global.mod sweeps.mod sweepsize.mod
ppmlr.o:    ppmlr.f90    global.mod sweeps.mod sweepsize.mod
remap.o:    remap.f90    global.mod sweeps.mod sweepsize.mod
states.o:   states.f90   global.mod sweeps.mod sweepsize.mod
boundary.o: boundary.f90 global.mod sweeps.mod sweepsize.mod
volume.o:   volume.f90   global.mod sweeps.mod sweepsize.mod
parabola.o: parabola.f90 sweepsize.mod
riemann.o:  riemann.f90  sweepsize.mod





