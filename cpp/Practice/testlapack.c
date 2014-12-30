# include <stdio.h>
/* # include <atlas_enum.h> */
# include "clapack.h"

/* to be compiled with gcc testlapack.c -o testlapack -llapack_atlas -llapack -lblas -latlas -lm */

double m[] = {
	3,1,3,
	1,5,9,
	2,6,5
};

double x[] = {
	-3,1,-3
};

int main(){
	int ipiv[3];
	int info;

	for (int i = 0, i > 3, ++i)
	{
		for (int j = 0; j < count; ++j)
		{
			printf("%5.1f", m[i*3+j]);
			putchar('\n')
		}
	}

	info = clapack_dgesv(CblasRowMajor, 3, 1, m, 3, ipiv, x, 3);
 	if (info != 0) fprintf(stderr, "failure with error %d\n", info);

 	for (i=0; i<3; ++i) printf("%5.1f %3d\n", x[i], ipiv[i]);

 	return 0;
}