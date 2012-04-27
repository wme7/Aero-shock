//------------------------------------------------------------ 
// VectorAddition.cu
//------------------------------------------------------------

// Includes
#include <stdio.h>

// Variables
float* h_A;
float* h_B;
float* h_C;
float* d_A;
float* d_B;
float* d_C;
bool noprompt = false;

// Functions
void Cleanup(void);
void RandomInit(float*, int);

// Device Code
__global__ void VecAdd(const float* A, const float* B, float* C, int N)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if (i < N) {
		C[i] = A[i] + B[i];
	}
}

// Host Code
int main(int nNumberofArgs, char** psArgs)
{
	printf("Vector addition\n");
	int N = 50000;
	size_t size = N*sizeof(float);

	// Allocate input vectors h_A and h_B in host memory
	h_A = (float*)malloc(size);
	if (h_A == 0) Cleanup();
	h_B = (float*)malloc(size);
	if (h_B == 0) Cleanup();
	h_C = (float*)malloc(size);
	if (h_C == 0) Cleanup();

	// Initializa input vectors
	RandomInit(h_A, N);
	RandomInit(h_B, N);

	// Allocate vectors in device memory
	cudaMalloc((void**)&d_A, size);
	cudaMalloc((void**)&d_B, size);
	cudaMalloc((void**)&d_C, size);

	// Copy vectors from host memory to device memory
	cudaMemcpy(d_A, h_A, size, cudaMemcpyHostToDevice);
	cudaMemcpy(d_B, h_B, size, cudaMemcpyHostToDevice);

	// Invoke kernel
	int threadsPerBlock = 256;
	int blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;
	VecAdd<<<blocksPerGrid, threadsPerBlock>>>(d_A, d_B, d_C, N);

	//Copy result from Device to host memory
	// Remember: h_C contains the result in host memory
	cudaMemcpy(h_C, d_C, size, cudaMemcpyDeviceToHost);

	// Verify Result
	int i;
	for (i = 0; i < N; i++) {
		float sum = h_A[i] + h_B[i];
		if (fabs(h_C[i] - sum) > 1e-5)
		break;
	}
	printf("%s \n", (i == N) ? "PASSED" : "FAILED");
	Cleanup();
}

void Cleanup (void)
{
	// Free Device Memory
	if (d_A)
		cudaFree(d_A);
	if (d_B)
		cudaFree(d_B);
	if (d_C) 
		cudaFree(d_C);

	// Free host memory
	if (h_A)
		free(h_A);
	if (h_B)
		free(h_B);
	if (h_C)
		free(h_C);

	cudaThreadExit();
	
	if (!noprompt) {
	printf("\nPress ENTER to exit ... \n");
	fflush( stdout);
	fflush( stderr);
	getchar();
	}
	exit(0);
}

// Allocates an array with random float entries.
void RandomInit(float* data, int n)
{
	for (int i = 0; i < n; ++i)
	data[i] = rand() / (float)RAND_MAX;
}
