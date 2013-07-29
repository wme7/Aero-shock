FC=gfortran
FC_FLAGS=-O3
FC_LIBS=-lstdc++ 
# -lcudart -lcublas
#FC_LD_LIB=/usr/local/cuda-5.0/lib64
#FC_INCL=/usr/local/cuda-5.0/include

CC=g++
CC_FLAGS=-O3
CC_LIBS=
CC_LD_LIB=
CC_INCL=


#CU_CC=nvcc
#CU_FLAGS=-arch compute_20 
#CU_LIBS=-lcublas


TARGET=stuDG_F90
SOURCE=stuDG.f90


$(TARGET): $(SOURCE) lgl_interp.o mesh_tools_2d_dg.o problem_setup_2d_dg.o \
lin_alg.o dg_solvers.o dg_c_solvers.o 
	$(FC) -o $(TARGET) $(SOURCE) $(FC_FLAGS) lgl_interp.o \
	mesh_tools_2d_dg.o problem_setup_2d_dg.o lin_alg.o \
	dg_solvers.o dg_c_solvers.o -L$(FC_LD_LIB) -I$(FC_INCL) $(FC_LIBS)

lgl_interp.o: lgl_interp.f90
	$(FC) -c lgl_interp.f90 $(FC_FLAGS)

mesh_tools_2d_dg.o: mesh_tools_2d_dg.f90
	$(FC) -c mesh_tools_2d_dg.f90 $(FC_FLAGS)

problem_setup_2d_dg.o: problem_setup_2d_dg.f90
	$(FC) -c problem_setup_2d_dg.f90 $(FC_FLAGS)

lin_alg.o: lin_alg.f90
	$(FC) -c lin_alg.f90 $(FC_FLAGS)

dg_solvers.o: dg_solvers.f90
	$(FC) -c dg_solvers.f90 $(FC_FLAGS)

dg_c_solvers.o: dg_c_solvers.cpp
	$(CC) -c dg_c_solvers.cpp $(CC_FLAGS)

#dg_cuda_solvers.o: dg_cuda_solvers.cu 
#	$(CU_CC) -c dg_cuda_solvers.cu $(CU_LIBS) $(CU_FLAGS)


clean:
	rm *.o *.mod $(TARGET)
