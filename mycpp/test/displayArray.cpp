// diplay the members of an arra of length nsize

#include <iostream>
using namespace std;

void displayArray( int intArray[], int nSize) {
	cout << "The value of the array is: \n";
	
	int* pArray = intArray;
	for( int n=0; n<nSize; n++)
	{
		cout << n << ": " << intArray[n] << "\n";
	}
	cout << "\n";
}

int main() {
	int array[] = {4, 3, 2, 1};
	displayArray(array,4);

return 0;
}
