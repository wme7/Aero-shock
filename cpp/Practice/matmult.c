/* c source code for matrix multiplication */

# include <stdio.h>
# include <stdlib.h>
# include <mkl.h>

/* To be compiled as: $ icc -mkl matmult.c */

/* Consider adjusting LOOP_COUNT based on the performance of your computer to 
make sure that a total run time is at least 1 second */
/* LOOP_COUNT as a constant macro value */ 
# define LOOP_COUNT 1
/* define minimun and maximun macro functions */
# define min(x,y) ( ((x)<(y)) ? (x):(y) )
# define max(x,y) ( ((x)>(y)) ? (x):(y) )

int main()
{
	double *A, *B, *C;
	int m, n, k, i, j, l, r;
	double alpha, beta;
	double sum;
	double s_initial, s_elapsed;

	printf(" compute real matrix C = alpha*A*B+beta*C using \n"
		"Intel(R) MKL dgemm. Here A, B and C are matrices \n"
		"and alpha beta are double precision scalars \n\n");

	m=1500; k=2000; n=1000;
	printf(" Initializing data for matrix multiplication C = A*B for \n"
		"matrix A(%ix%i) and matrix B(%ix%i) \n\n",m,k,k,n);

	alpha=1.0; beta=0.0;

	printf(" Allocating memory for matrices aligned on 64-byte boundary \n"
		"for better performance \n\n");
	A = (double *)mkl_malloc(m*k*sizeof(double),64);
	B = (double *)mkl_malloc(k*n*sizeof(double),64);
	C = (double *)mkl_malloc(m*n*sizeof(double),64);
	if (A==NULL||B==NULL||C==NULL){
	printf("\n ERROR: Can't allocate memory for matrices. Aborting ... \n\n");
	mkl_free(A);
	mkl_free(B);
	mkl_free(C);
	return 1;
	}

	printf(" Initializing matrix data \n\n");
	for (i=0;i<(m*k);i++){
		A[i] = (double)(i+1);
	}

	for (i=0;i<(k*n);i++){
		B[i] = (double)(-(i+1));
	}

	for (i=0;i<(m*n);i++){
		C[i] = 0.0;
	}

	printf(" Computing matrix product using tripe nested loop \n\n");
	s_initial = dsecnd();
	for (r=0;r<LOOP_COUNT;r++){
		for (i=0;i<m;i++){
			for (j=0;j<n;j++){
				sum = 0.0;
				for (l=0;l<k;l++)
					sum += A[k*i+l] * B[n*l+j];
				C[n*i+j] = sum;
			}
		}
	}
	s_elapsed = (dsecnd() - s_initial)/LOOP_COUNT;

	printf(" Matrix A: \n");
	for (i=0;i<min(m,5);i++){
		for (j=0;j<min(k,5);j++){
			printf("%12.0f", A[j+i*k]);
		}
		printf("\n");
	}

	printf(" Matrix B: \n");
	for (i=0;i<min(k,5);i++){
		for (j=0;j<min(n,5);j++){
			printf("%12.0f", B[j+i*n]);
		}
		printf("\n");
	}

	printf(" Matrix C: \n");
	for (i=0;i<min(m,5);i++){
		for (j=0;j<min(n,5);j++){
			printf("%12.5G", C[j+i*n]);
		}
		printf("\n");
	}

	printf("Matrix multiplication using triple nested loop completed == \n"
		"== at %.5f milliseconds == \n\n",s_elapsed);

	/* end of program */
	printf("\n Deallocatting memory \n\n");
	mkl_free(A);
	mkl_free(B);
	mkl_free(C);

	printf("Example completed. \n\n");
	return 0;

}