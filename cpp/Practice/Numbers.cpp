#include <iostream> // input and output streams
#include <climits>	// max & min size of interger types
#include <cfloat>   // max & min size of real types
using namespace std;

/*
 * This is a simple example to help myself become familiar with C/C++ language
 * Manuel Diaz, NTU, 2012.09.07
 * 
 * Refs. 
 * 1. E. Scheinerman; C++ for mathematiticians, 2006.
 * 2. P. Geteruer; Writing Matlab C/MEX code, 2010.
 * 3. J. Sanders, E. Kandrot; CUDA by Example, 2009.
 */

int main() {
	cout << "Hello C++, this is very cool!";
	cout << endl; // this is to finish the print of my line

	int x; 	int y;

	// set:
	x = 3; 	y = 4;
	
	/* 
	 * show me x and y:
	 */
	cout << "x is : " << x << "." << endl;
	cout << "y is : " << y << "." << endl;

	/*
	 * Let me see the result of this sum
	 */
	cout << "then x+y is " << x+y << endl;

	/*
	 * Integer tets
	 */
	int million = 1000000; 	int trillion = million*million;
	cout << "1 million is " << 1000000 << endl;
	cout << "According to this computer, " << million << " squared is " << trillion << "." << endl;
	cout << "which cannot be right, or wallstreet will be stealing my money!" << endl;

	/*
	 * Lets report the size of the different C/C++ variables
	 */
	 // integers type
	cout << "the size of short 	is " << sizeof(short) << "bytes" << endl;
	cout << "the size of int  	is " << sizeof(int)	<< "bytes" << endl;
	cout << "the size of long 	is " << sizeof(long) << "bytes" << endl;
	cout << "the size of long long	is " << sizeof(long long) << "bytes" << endl;
	
	// characters type
	cout << "the size of char 	is " << sizeof(char) << "bytes" << endl;
	cout << "the size of bool 	is " << sizeof(bool) << "bytes" << endl;

	// floating point types
	cout << "the size of float 	is " << sizeof(float) << "bytes" << endl;
	cout << "the size of double	is " << sizeof(double) << "bytes" << endl;

	// long double might not exist in all computers
	cout << "the size of long double	is " << sizeof(long double)  << "bytes" << endl;
	cout << endl;

	/*
	 * The range of values for each type of variable is,
	 */

	// Integer types
	cout << "The maximun size of short is " << SHRT_MAX << endl;
	cout << "The minimum size of short is " << SHRT_MIN << endl;

	cout << "The minimum size of int is " << INT_MIN << endl;
	cout << "The maximun size of int is " << INT_MAX << endl;

	cout << "The maximun size of long is " << LONG_MAX << endl;
	cout << "The minimum size of long is " << LONG_MIN << endl;

	cout << "The maximun size of long long is " << LLONG_MAX << endl;
	cout << "The minimum size of long long is " << LLONG_MIN << endl;	
	cout << endl;

	// Float types
	cout << "The maximun value of a float is " << FLT_MAX << endl;
	cout << "The minimum epsilon size of float is " << FLT_EPSILON << endl;
	cout << "The minimum value of float is " << FLT_MIN << endl;

	cout << "The maximun value of a double is " << DBL_MAX << endl;
	cout << "The minimum epsilon size of double is " << DBL_EPSILON << endl;
	cout << "The minimum value of double is " << DBL_MIN << endl;

	// Long double might not be defined on some systems
	cout << "The maximun value of a long double is " << LDBL_MAX << endl;
	cout << "The minimum epsilon size of long double is " << LDBL_EPSILON << endl;
	cout << "The minimum value of long double is " << LDBL_MIN << endl;
	cout << endl;

	/*
	 * About standar operations
	 */

	int numerator = 13;
	int denominator = 5;
	double quotient;

	quotient = numerator/denominator;

	cout << "a: " << numerator << ", b: " <<  denominator << ", the quotient a/b is: " << quotient;
	cout << endl;
	cout << "This cannot be!, wrong selection of variable type!, using double" << endl;

	quotient = double(numerator)/double(denominator);

	cout << "a: " << numerator << ", b: " <<  denominator << ", the quotient a/b ins: " << quotient;
	cout << endl;
	cout << "now, yes!" << endl; cout << endl;

	/*
	 * The modulus '%' operator between integers
	 */

	cout << "The expression 237%100 gives: " << 237%100 << endl;
	cout << endl;

	/*
	 * C++'s combination of arithmetic operations
	 */

	int a = x;

	// Recall
	cout << "x is : " << x << "." << endl;
	cout << "y is : " << y << "." << endl;
	//operation
	a += y; cout << " x += y gives: " << a << endl; a=x;
	a -= y; cout << " x -= y gives: " << a << endl; a=x;
	a *= y; cout << " x *= y gives: " << a << endl; a=x;
	a /= y; cout << " x /= y gives: " << a << endl; a=x;
	a %= y; cout << " x %= y gives: " << a << endl; //a=x;
	cout << "finally x is " << a << endl;

	/*
	 * End of Chapters 1 & 2.
	 */

	return 0;
}