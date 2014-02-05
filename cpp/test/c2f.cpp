//
// Example Program to convert tempererature from Celcious to Fahrenheit degree units
//

// Includes
#include <cstdio>
#include <cstdlib>
#include <iostream>
using namespace std;

// Enter the temperature in Celsius
int main(int nNumberofArg, char* pszArgs[])
{
	int celsius;
	cout << "Enter the temperature in Celsius:";
	cin >> celsius;

	// Compute fahrenheit values
	int factor;
	factor = 212 - 32;
	int fahrenheit;
	fahrenheit = factor * celsius / 100 + 32;

	// Output the results
	cout << "Fahrenheit value is:";
	cout << fahrenheit << endl;

	// Wait for the user to see the result under windows
	//system("PAUSE");
	//return 0;
}
