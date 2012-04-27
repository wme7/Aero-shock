#include "stdio.h"

__global__ void add_arrays_gpu( float *in1, float *in2, float *out, int Ntot)
{
	int idx=blockIdx.x*blockDim.x+threadIdx.x;
	if ( idx )
		out[idx]=in1[idx]+in2[idx];
}

int main()
{
	/* pointers to host memory */
	float *a, *b, *c;
	/* pointers to device memory */
	float *a_d, *b_d, *c_d;
	int N=18;
	int i;
	
	/* Allocate arrays a, b and c on host*/
	a = (float*) malloc(N*sizeof(float));
	b = (float*) malloc(N*sizeof(float));
	c = (float*) malloc(N*sizeof(float));
	
	/* Allocate arrays a_d, b_d and c_d on device*/
	cudaMalloc ((void **) &a_d, sizeof(float)*N);
	cudaMalloc ((void **) &b_d, sizeof(float)*N);
	cudaMalloc ((void **) &c_d, sizeof(float)*N);
	
	/* Initialize arrays a and b */
	for (i=0; i<N; i++)
	{
		a[i] = (float) i;
		b[i] =(float) i;
	}
	
	
	/* Copy data from host memory to device memory */
	cudaMemcpy(a_d, a, sizeof(float)*N, cudaMemcpyHostToDevice);
	cudaMemcpy(b_d, b, sizeof(float)*N, cudaMemcpyHostToDevice);
	
	/* Compute the execution configuration */
	int block_size=8;
	dim3 dimBlock(block_size);
	dim3 dimGrid ( (N/dimBlock.x) + (!(N%dimBlock.x)?0:1) );
	
	/* Add arrays a and b, store result in c */
	add_arrays_gpu<<<dimGrid,dimBlock>>>(a_d, b_d, c_d, N);
	
	/* Copy data from deveice memory to host memory */
	cudaMemcpy(c, c_d, sizeof(float)*N, cudaMemcpyDeviceToHost);
	
	/* Print c */
	for (i=0; i<N; i++)
		printf(" c[%d]=%f\n",i,c[i]);
	
	/* Free the memory */
	free(a); free(b); free(c);
	cudaFree(a_d); cudaFree(b_d);cudaFree(c_d);
	
}

