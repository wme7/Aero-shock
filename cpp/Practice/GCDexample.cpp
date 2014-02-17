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

 /**
  * The following formulation is correct but highly inefficient.
  * Problem in Chapter 3 of reference [1]. 
  * Is only used to introduce some of the fundamental concepts
  * in c/c++. MD 2012
  */

int main() {
	cout << "Greates Common Divisor program " << endl; 

	long a,b;

	cout << "Enter the first number --> ";
	cin >> a;
	cout << "Enter the second number --> ";
	cin >> b;

	cout << "The gcd of " << a << " and " << b << " is " << gcd(a,b) << endl;

	/**
	 * End of Chapters 3.
	 */

	return 0;
}