# Make File to Compile my Fortran Programs
# By Manuel Diaz, NTU, 2013.07.10
#
# Defining Variables
objects = main.o dispmodule.o mathmodule.o 
#f90comp = ifort
f90comp = gfortran
# Start the makefile:
runme: $(objects)
	$(f90comp) -o runme $(objects)
dispmodule.o: dispmodule.f90
	$(f90comp) -c dispmodule.f90
mathmodule.o: mathmodule.f90
	$(f90comp) -c mathmodule.f90
main.o: dispmodule.o mathmodule.o main.f90
	$(f90comp) -c main.f90
# Cleaning everything
clean:
	rm *.o
	rm *.mod
	rm runme
	#rm $(objects)
	#rm some.mod runme
# End of the makefile
