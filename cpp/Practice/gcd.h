#ifndef GCD_H
#define GCD_H

/**
 * Find the Greates Common Divisor (GCD) using C/C++
 * Manuel Diaz, NTU, 2012.09.07
 * 
 * Refs. 
 * 1. E. Scheinerman; C++ for mathematiticians, 2006.
 * 2. P. Geteruer; Writing Matlab C/MEX code, 2010.
 * 3. J. Sanders, E. Kandrot; CUDA by Example, 2009.
 */

/*
 * Calculate the greastest common divisor of two integers
 * Note: gcd(0,0) will return 0 and print an error message.
 * @param a the first integer
 * @param b the second integer
 * @return the greastest common divisor of a and b
 */

 long gcd(long a, long b);

 #endif