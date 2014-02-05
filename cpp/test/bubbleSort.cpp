// This is the Classic Bubble Sort Algotirhm by Manuel Diaz (C)

#include <iostream>
#include <fstream>
using namespace std;

void bubbleSort(float array[], int array_size) 
{
	cout << "Starting the bubble sort! \n";
	bool done = true;
	float temp;
	while (done)
	{ 
		done = false;
		for (int i=0; i<array_size; i++)
		{
			if ( array[i]> array[i+1] )
			{
			temp = array[i];
			array[i] = array[i+1];
			array[i+1] = temp;
			done = true;
			}
		}
	}
}

int main () 
{
	float mylist[] = { 10.2, 9.3, 8.4, 7.6, 6.1, 5.5, 4.4, 3.8, 2.9, 1.3 };
	int nelements = sizeof(mylist) / sizeof(mylist[1]);
	bubbleSort(mylist, nelements);
	cout << "Printing last element: \n";
	cout << mylist[10] << endl;

return 0;
}
