SHELL = /bin/bash
SRCS := $(wildcard *.f90)
CC = gfortran

#Compiler flags
ifeq ($(mode),release)
   FLAGS = -O3 -fdefault-real-8
else
   mode = debug
   FLAGS = -g -O0 -fdefault-real-8
endif

#LOPTION = -static
ifeq ($(arch),32)
   EXE = dg-weno-bgk
   arch = 32
   LIB = 
   LIB_PATH= 
else
   EXE = dg-weno-bgk-64
   LIB = 
   LIB_PATH=
endif

MODULE = data_module.mod
OBJS = dg_weno_bgk_1d.o

$(EXE): ${OBJS}
	$(CC) $(OBJS) com.o $(LOPTION) $(OPTION) $(FLAGS) $(LIB) $(LIB_PATH) -o $(EXE)
#	mv $(EXE) $(HOME)/bin/$(EXE)


$(MODULE):com.f90
	$(CC) -c $(FLAGS) com.f90

$(OBJS):dg_weno_bgk_1d.f90 $(MODULE)
	$(CC) -c $(FLAGS) dg_weno_bgk_1d.f90

%.o: %.cpp *.h
	$(CC) $(INC) $(CFLAGS) $(OPTION) -c $<

.PHONY: clean
clean:
	rm -fv $(EXE) 
	rm -fv ${OBJS} *.o *.dat *.plt *.mod

information:
ifneq ($(mode),release)
ifneq ($(mode),debug)
	@echo "Invalid build mode." 
	@echo "Please use 'make mode=release' or 'make mode=debug'"
	@exit 1
endif
endif
	@echo "Building on "$(mode)" mode"
	@echo ".........................."
