// Read matrix file in to memory algorithm by Manuel Diaz (C)

#include <iostream>
#include <fstream>
using namespace std;

int main() {
	double x;
	double sum = 0;
	ifstream myfile ("data2.dat");
	if (myfile.is_open())
	{
	  while ( myfile.good() )
	  {
 	    while (myfile >> x)
	    {
	    sum += x;
	    }
	  }
	cout <<"the sum is: " << sum << "\n";
	myfile.close();
	}
	else 
	cout << "Unable to open file";

return 0;
}
