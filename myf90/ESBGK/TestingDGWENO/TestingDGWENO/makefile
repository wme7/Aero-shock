# Make File to Compile my Fortran Programs
# By Manuel Diaz, NTU, 2013.07.10
#
# Defining Variables
objects = main.o dispmodule.o mathmodule.o IDmodule.o DGmodule.o quadraturemodule.o tecplotmodule.o initmodule.o WENOmodule.o
#f90comp = ifort
f90comp = gfortran
# Start the makefile:
runme: $(objects)
	$(f90comp) -o runme $(objects)
dispmodule.o: dispmodule.f90
	$(f90comp) -c dispmodule.f90
mathmodule.o: mathmodule.f90
	$(f90comp) -c mathmodule.f90
IDmodule.o: IDmodule.f90
	$(f90comp) -c IDmodule.f90
DGmodule.o: DGmodule.f90
	$(f90comp) -c DGmodule.f90
quadraturemodule.o: quadraturemodule.f90
	$(f90comp) -c quadraturemodule.f90
tecplotmodule.o: tecplotmodule.f90
	$(f90comp) -c tecplotmodule.f90
initmodule.o: initmodule.f90
	$(f90comp) -c initmodule.f90
WENOmodule.o: WENOmodule.f90
	$(f90comp) -c WENOmodule.f90

main.o: dispmodule.o mathmodule.o IDmodule.o DGmodule.o quadraturemodule.o tecplotmodule.o initmodule.o WENOmodule.o main.f90
	$(f90comp) -c main.f90
# Cleaning everything
clean:
	rm *.o
	rm *.mod
	rm runme
	#rm $(objects)
	#rm some.mod runme
# End of the makefile
