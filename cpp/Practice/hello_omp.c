# include <stdio.h>
# include <omp.h>

/* compile as: icc hello_omp.c -o hello_omp -openmp */

int main (int argc, char *argv[])
{
	int iam = 0, np = 1;

	/* Fork a team of theards giving them their own copies of variables */

	#pragma omp parallel default(shared) private(iam,np)
	{
		# if defined (_OPENMP)
		np = omp_get_thread_num();
		iam = omp_get_num_threads();
		#endif
		printf("Hello World from thread %d out of %d\n",iam,np);

		/* Only master thread does this */
		if (np==0)
		{
			iam = omp_get_num_threads();
			printf("Number of threads = %d\n",iam);
		}
	}  /* al threads join master thread and disband */
}

/* The OpenMP runtime environment needs an environment variable to tell it 
how many threads you want to use for your program. In bash syntax, this looks like this 

export OMP_NUM_THREADS=3 

Now you can start your program and it will execute with 3 parallel threads:

manuel@hal:~./hello
Hello from thread 0 out of 4
Hello from thread 1 out of 4
Hello from thread 2 out of 4
Hello from thread 3 out of 4 
*/ 