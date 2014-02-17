// write string to a file

#include <iostream>
#include <fstream>
using namespace std;

int main() {
	ofstream myFile;
	myFile.open ("example.txt");
	myFile << "Writing this to a file.\n";
	myFile.close();
	return 0;
}
