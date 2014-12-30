#include <stdio.h>
#include <math.h>

/* root.c: computes an approximate square root of two */

main(){
	float X, Y, epsilon;
	printf("Enter epsilon: "); scanf("%f", &epsilon);

	X = 1;
	do {
		Y = 2/X;
		X = (X+Y)/2;
	} while ( fabs(X-Y) > epsilon );

	printf("Approximate square root of 2 = %f\n", X);
}