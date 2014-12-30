# include <stdio.h>

/* compile with: gcc testdgtsv.c -o testdgtsv -llapack -lblas */

double l[] = {-1,-2,-1,-1};
double d[] = {2,2,3,3,1};
double u[] = {-1,-1,-1,-2};
double x[] = {1,2,3,2,1};

static long dgtsv(long N, long NRHS, double *DL, double *D, double *DU, double *B, long LDB)
{
	extern	void dgtsv_(const long *Np, const long *NRHSp, double *DL, double *D, double *DU, double *B, const long *LDBp, long *INFOp);

	long info;
	dgtsv_(&N, &NRHS, DL, D, DU, B, &LDB, &info);
	return info;
}

int
main(){
	int i, info;

	info = dgtsv(5,1,l,d,u,x,5);
	if (info != 0) fprintf(stderr, "failure with error %d\n", info);

	for (i=0;i<5;i++) printf("%5.1f\n",x[i] );

	return 0;
}