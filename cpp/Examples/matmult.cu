// CUDA version of Matrix multiplication by Manuel Diaz (c)

#include <iostream>
#include <fstream>
#include "../common/book.h"
using namespace std;

__global__ void MatrixMultiplicationKernel(float* d_A, float* d_B, float* d_C, int n) 
{
	int x = blockIdx.x * blockDim.x + threadIdx.x;
	int y = blockIdx.y + blockDim.y + threadIdx.y;

	int temp = 0;
	d_C[y*n+x]=0;
	for (int k=0; k<n; k++)
		d_C[y*n+x]=d_A[y*n+k]*d_B[k*n+x];
}

int main(void) 
{
	int i, j, n, m, p, q;

	// load matrix A
	cout << "...loading file A.mat\n";
	cout << "Enter the size n of matrix A: ";
	cin >> n;
	p=n*n;
	float A[p];
		
	// load matrix B
	cout << "...loading file B.mat\n";
	cout << "Enter the size m of matrix B: ";
	cin >> m;
	q=m*m;
	float B[q];
	if (m == n)
	{
		cout << "Computing A*B = \n";
		ifstream Amat ("A.mat");
		if (Amat.is_open())
		{
			while (Amat.good())
			{
				for (i=0; i<n; i++)
				{
				Amat >> A[i];
				}
			}
			Amat.close();
			cout << "Data loaded in vector B \n";
		}
		else cout << "Unable to open file A.mat";
		ifstream Bmat ("B.mat");		
		if (Bmat.is_open())
		{
			while (Bmat.good())
			{
				for (j=0; j<m; j++)
				{
				Bmat >> B[j];
				}
			}
			Bmat.close();
			cout << "Data loaded in vector A \n";
		}
		else cout << "Unable to open file B.mat";
		
		// If matrices A and B we start the result matrix C
		float C[p];

		// Load Matrices on the Device
		float* d_A;
		float* d_B;
		float* d_C;
		size_t matsize = n*n*sizeof(float);
		cudaMalloc((void**)&d_A,matsize);
		cudaMemcpy(d_A, A, matsize, cudaMemcpyHostToDevice);
		cudaMalloc((void**)&d_B,matsize);
		cudaMemcpy(d_B, B, matsize, cudaMemcpyHostToDevice);
		cudaMalloc((void**)&d_C,matsize);

		// Invoque cuda kernel
		int tile_width = 16; // Because I want 16x16 = 256 threads per block 
		dim3 dimGrid(n/tile_width,n/tile_width,1); // 2D grid
		dim3 dimBlock(tile_width,tile_width,1); // 2D blocks
		MatrixMultiplicationKernel<<<dimGrid,dimBlock>>>(d_A,d_B,d_C,n);

		// Copy result to device
		cudaMemcpy(C, d_C, matsize, cudaMemcpyDeviceToHost);

		// Free Device Memory
		cudaFree(d_A); 	cudaFree(d_B); 	cudaFree(d_C);
		
		// write to file matrix C
		ofstream myfile ("C2.mat");
		if (myfile.is_open())
		{
			for (j=0; j<n; j++)
			{
				for (i=0; i<n; i++)
				{
				myfile << C[j*n+i] << " ";
				}
			myfile << endl;
			}
			myfile.close();
		}
		else cout << "Unable to open file";
	}
	else cout << "Matrix A and B have different sizes"; 
return 0;
}
