#include <stdio.h>
#include "foo.h" /* Include the header here, to obtain the function declaration */ 
#include "data.h" /* Include the vector */

/* compile as: gcc -o my_app main.c foo.c*/

int main(void){
	int y = foo(3); /* Use the function here */
	printf("%d\n", y);
	double x = matrix[0][0];
	printf("%2.4f\n", x);
	return	0;
}