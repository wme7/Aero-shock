#include <stdio.h>
#include <mkl.h>

int main(void) {
  double *a, *b, *c;
  int n, i;
  double alpha, beta;
  MKL_INT64 AllocatedBytes;
  int N_AllocatedBuffers;
 
  alpha = 1.1; beta = -1.2;
  n = 1000;
  mkl_peak_mem_usage(MKL_PEAK_MEM_ENABLE);
  a = (double*)mkl_malloc(n*n*sizeof(double),64);
  b = (double*)mkl_malloc(n*n*sizeof(double),64);
  c = (double*)mkl_calloc(n*n,sizeof(double),64);
  for (i=0;i<(n*n);i++) {
     a[i] = (double)(i+1);
     b[i] = (double)(-i-1);
  }
 
  dgemm("N","N",&n,&n,&n,&alpha,a,&n,b,&n,&beta,c,&n);
  AllocatedBytes = mkl_mem_stat(&N_AllocatedBuffers);
  printf("\nDGEMM uses %d bytes in %d buffers",AllocatedBytes,N_AllocatedBuffers);
 
  mkl_free_buffers();
  mkl_free(a);
  mkl_free(b);
  mkl_free(c);
 
  AllocatedBytes = mkl_mem_stat(&N_AllocatedBuffers);
  if (AllocatedBytes > 0) {
      printf("\nMKL memory leak!");
      printf("\nAfter mkl_free_buffers there are %d bytes in %d buffers",
         AllocatedBytes,N_AllocatedBuffers);
  }
  printf("\nPeak memory allocated by Intel MKL memory allocator %d bytes. Start to count new memory peak",
         mkl_peak_mem_usage(MKL_PEAK_MEM_RESET));
  a = (double*)mkl_malloc(n*n*sizeof(double),64);
  a = (double*)mkl_realloc(a,2*n*n*sizeof(double));
  mkl_free(a);
  printf("\nPeak memory allocated by Intel MKL memory allocator after reset of peak memory counter %d bytes\n",
         mkl_peak_mem_usage(MKL_PEAK_MEM));
 
  return 0;
}