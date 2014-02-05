// Read integers from file and compute its sum

#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;

int main() {
	int sum = 0;
	int x;
	ifstream inFile;

	inFile.open("test.txt");
	if (!inFile) {
		cout << "Unable to open file";
		exit (1); // terminate with error
	}
	
	while (inFile >> x) {
		sum = sum + x;
	}

	inFile.close();
	cout << "sum = " << sum << endl;
	return 0;
}
