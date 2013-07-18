all : DG1D_kinetic

PROG = DG1D_kinetic

DIR    = .
SRCDIR = ./src
BIN    = ./bin
FBIN   = ./bin


#MOD_SFX= d
MOD_SFX=mod
.SUFFIXES:
.SUFFIXES: .f90 .o .$(MOD_SFX)

F90=ifort
#F90=gfortran
#F90=pgf95

FC = ifort

LINK	= $(F90) -O3
COMP90	= $(F90) -O3
COMP77	= $(F90) -O3 #-132
FFLAGS  =  -m64 #-heap-arrays  -fp-stack-check -g -check all
FFLAGS2 =   -m64 #-heap-arrays  -fp-stack-check -g -check all
#FFLAGS	= -unroll8 -autodouble -W1 -mtune=core2 -g -heap-arrays
#FFLAGS2 = -unroll8 -autodouble -W1 -mtune=core2 -g -heap-arrays
#LIBDIR  = -L/usr/local/Matlab6p1/extern/lib/glnx86
#LIB     = -leng -lmat -lmx -lut
#LIBDIR2 = -L/opt/intel/compiler70/ia32/lib
#LIB2    = -lintrins

# ABSOFT fortran
#COMP90  = f90 -O2
#COMP77  = f77 -O2 -W
#LINK    = f90 -O2  
#FFLAGS  =  -cpu:p7  
#FFLAGS2 =  -cpu:p7
#LIBDIR  = -L/usr/local/Matlab6p1/extern/lib/glnx86
#LIB     = -leng -lmat -lmx -lut 
#-llibmat.so -llibmx.so -llibut.so  


##LIBDIR  = -L/usr/local/matlab6p1/extern/lib/glnx86/
#LIB     =  libeng.so libmat.so libmx.so libut.so  

#LINK    = f77 -fast -xO4 -xsafe=mem -xarch=v8plusa -xchip=ultra
#COMP    = f77 -e -fast -xO4 -xsafe=mem -xarch=v8plusa -xchip=ultra

OBJS = 	$(FBIN)/splibsrc1.o\
	$(FBIN)/universal.o\
	$(FBIN)/Legendre.o\
	$(FBIN)/State_Var.o\
	$(FBIN)/RK_Var.o\
	$(FBIN)/RK_Pack.o\
	$(FBIN)/Kinetic_Pack.o\
	$(FBIN)/Output_Pack.o\
	$(FBIN)/LinearSolvers.o\
	$(FBIN)/NSOLI.o\
	$(FBIN)/DG1D_kinetic.o
	
#	$(FBIN)/Filter_Var.o\
#	$(FBIN)/Filter_Pack.o\
#	$(FBIN)/Time_Pack.o\
#	$(FBIN)/OutPut_Tec.o\
#	$(FBIN)/Visual_Graph.o\
#	$(FBIN)/MG_Var.o\

                        
# Rules to build main programs.

$(PROG):  $(OBJS)  #nitsol_lib lapack_lib blas_lib
	echo $(BIN)
	@echo ""
	@echo "Linking object files $^ into $@"
	$(LINK) $(OBJS) -o $(PROG) $(LIBDIR) $(LIB) #$(NITSOL) $(LAPACK) $(BLAS) 
#	mv $(PROG) ../.

#$(LIBDIR) $(LIB)

#Source file
$(FBIN)/%.o: $(SRCDIR)/%.f90
	@echo ""
	@echo "Compiling source f90 file $< into $@"
	$(COMP90) $(FFLAGS) -c $<
	mv *.o $(FBIN)
	@echo "Done"


$(FBIN)/%.o: $(SRCDIR)/%.for
	@echo ""
	@echo "Compiling source for file $< into $@"
	$(COMP77) $(FFLAGS2) -c $<
	mv *.o $(FBIN)
	@echo "Done"

$(FBIN)/%.o: $(SRCDIR)/%.f
	@echo ""
	@echo "Compiling source f file $< into $@"
	$(COMP77) $(FFLAGS2) -c $<
	mv *.o $(FBIN)
	@echo "Done"



#Module 
#%.$(MOD_SFX):
#	@echo ""
#	@if [ ! -s $@ ]; then \
#		echo "Error, module file deleted?"; exit 1; fi
#	@echo "Module file $@ considered updated through dependence on: $^"
#	touch $@

#Dependencies for creating objectfile (USE Statement)
#$(FBIN)/DG1D_scalar.o:\
#	$(FBIN)/UNIVERSAL_CONST.$(MOD_SFX)\
#	$(FBIN)/LEGENDRE.$(MOD_SFX)\
#	$(FBIN)/STATE_VAR.$(MOD_SFX)\	
#	$(FBIN)/RK_VAR.$(MOD_SFX)






#Location of modules in source file
$(FBIN)/UNIVERSAL_CONST.$(MOD_SFX) :$(FBIN)/universal.o
$(FBIN)/LEGENDRE.$(MOD_SFX)        :$(FBIN)/Legendre.o
$(FBIN)/RK_VAR.$(MOD_SFX)          :$(FBIN)/RK_Var.o
$(FBIN)/STATE_VAR.$(MOD_SFX)       :$(FBIN)/State_Var.o

#FILTER_VAR.$(MOD_SFX):$(FBIN)/Filter_Var.o
#MG_VAR.$(MOD_SFX):$(FBIN)/MG_Var.o 

clean:
	rm -f $(FBIN)/*.o
	rm -f $(PROG) *.mod


