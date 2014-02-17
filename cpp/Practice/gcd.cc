#include "gcd.h" 	// our header file with procedure 'gcd(a,b)'
#include <iostream> // input and output streams are need.
using namespace std;

/**
 * Find the Greates Common Divisor (GCD) using C/C++
 * Manuel Diaz, NTU, 2012.09.07
 * 
 * Refs. 
 * 1. E. Scheinerman; C++ for mathematiticians, 2006.
 * 2. P. Geteruer; Writing Matlab C/MEX code, 2010.
 * 3. J. Sanders, E. Kandrot; CUDA by Example, 2009.
 */

long gcd(long a, long b) {
	
	// if a and b both zero, print an error and return 0
	if ((a==0) && (b==0)) {
		cerr << "Warning: gcd called with both arguments equal to zero." << endl;
		return 0;
	}
	
	// make sure a and b are both nonnegative
	if (a<0) {
		a = -a;
	}
	if (b<0) {
		b = -b;
	}
	
	// if a is zero, the answer is b
	if (a==0) { 
		return b;
	}

	//otherise, we check all possibilities from 1 to a
	long d; // will hold the answer

	for (long t=1; t<=a; t++) {
		if ((a%t==0) && (b%t==0)) {
			d = t;
		}
	}

	return d;
}