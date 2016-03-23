CUDA_HOME = /usr/local/cuda
CC = g++ -O3 -c
NVCC = nvcc -DUNIX -O3 -c --use_fast_math -arch sm_20
LINK=g++ -fPIC
INCLUDE := -I$(CUDA_HOME)/include
LIB := -L$(CUDA_HOME)/lib64 -lcudart


DSMCSG_FUNCTIONS_C = dsmcsg_main.o dsmcsg_class.o dsmcsg_init.o dsmcsg_output.o dsmcsg_readfile.o
DSMCSG_FUNCTIONS_CUDA = dsmcsg_cudafunction.o scanLargeArray_kernel.o

all : dsmcsg dsmcsg_postprocessor

dsmcsg : ${DSMCSG_FUNCTIONS_C} ${DSMCSG_FUNCTIONS_CUDA}
	${LINK} -o dsmcsg ${DSMCSG_FUNCTIONS_C} ${DSMCSG_FUNCTIONS_CUDA} ${LIB}

dsmcsg_postprocessor : dsmcsg_postprocessor.o
	${LINK} -o dsmcsg_postprocessor dsmcsg_postprocessor.o

dsmcsg_main.o : dsmcsg_main.cpp
	${CC} -c dsmcsg_main.cpp ${INCLUDE}

dsmcsg_class.o : dsmcsg_class.cpp
	${CC} -c dsmcsg_class.cpp

dsmcsg_init.o : dsmcsg_init.cpp
	${CC} -c dsmcsg_init.cpp

dsmcsg_cudafunction.o : dsmcsg_cudafunction.cu
	${NVCC} -c dsmcsg_cudafunction.cu ${INCLUDE}

scanLargeArray_kernel.o : scanLargeArray_kernel.cu
	${NVCC} -c scanLargeArray_kernel.cu ${INCLUDE}

dsmcsg_output.o : dsmcsg_output.cpp
	${CC} -c dsmcsg_output.cpp

dsmcsg_readfile.o : dsmcsg_readfile.cpp
	${CC} -c dsmcsg_readfile.cpp

dsmcsg_postprocessor.o : dsmcsg_postprocessor.cpp
	${CC} -c dsmcsg_postprocessor.cpp

clean : 
	rm -f dsmcsg dsmcsg_postprocessor *.o

clean_dat :
	rm -f *.dat
