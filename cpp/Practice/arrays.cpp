#include <iostream>
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

int input(istream& in=cin)
{
	int x;
	in >> x;
	return x;
}

int main()
{
	int board[3][3]; //creates a 3*3 matrix or a 2d array.

	cout << "Enter the values for a 3x3 matrix: " << endl;

	for(int i=0; i<3; i++)    //This loops on the rows.
	{
		for(int j=0; j<3; j++) //This loops on the columns
		{
			board[i][j] = input(); //you can also connect to the file
			//and place the name of your ifstream in the input after opening the file will
			//let you read from the file.
		}
	}
	cout << endl;

	cout << "printing matrix array: " << endl;

	for(int i=0; i<3; i++)    //This loops on the rows.
	{
		for(int j=0; j<3; j++) //This loops on the columns
		{
			cout << board[i][j]  << "  ";
		}
		cout << endl;
	}
	cout << endl;

	/*
	 * This a little bit aside from the original problem in arrays, but I doit for the sake of pure curiosity
	 * What will happen if I do: vector++?
	 */

	int vector[3][3];

	vector[0][0] = 1; 	vector[0][1] = 2; 	vector[0][2] = 3;
	vector[1][0] = 4; 	vector[1][1] = 5; 	vector[1][2] = 6;
	vector[2][0] = 7; 	vector[2][1] = 8; 	vector[2][2] = 9;

	// Method 1, this method fails! it modifies the pointer instead of adding 1 to each element!.

	/*cout << vector++ << enld;*/

	// Method 2, the correct procedure:

	for (int i=0; i<3; i++) {
		for (int j=0; j<3; j++) {
			vector[i][j]++;
		}
	}

	for(int i=0; i<3; i++) {
		for(int j=0; j<3; j++) {
			cout << vector[i][j]  << "  ";
		}
		cout << endl;
	}
	cout << endl;

}