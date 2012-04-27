// Random Square Matrix Generator by Manuel Diaz (C)

/* we want to produce numbers in a specific range, 
 * rather than between 0 and RAND_MAX, we can use 
 * the modulo operator. It's not the best way to 
 * generate a range but it's the simplest. If we 
 * use rand()%n we generate a number from 0 to n-1. 
 * By adding an offset to the result we can produce 
 * a range that is not zero based.
 */

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
using namespace std;

int main() {
	char filename[6];
	int i, j, n, nrand;
	cout << "Square matrix Generator \n";
	cout << "Enter the name of the matrix file: ";
	cin >> filename;
	cout << "\n";
	cout << "Enter the desired size n: " ;
	cin >> n;
	cout << "\n";
	srand((unsigned)time(0)); //srand is called just one time.
	//To see this PC RAND_MAX number erase the comment lines:
	//cout <<"my RAND_MAX value is:" << RAND_MAX << endl;
	ofstream myfile;
	myfile.open (filename);
	for(j=0; j<n; j++)
	 {
	 for(i=0; i<n; i++)
	  {
	  int nrand = (rand()%10)+1;
	  myfile << nrand <<" "; 
	  }
	 myfile << endl;
	 }
myfile.close();
return 0;
}

