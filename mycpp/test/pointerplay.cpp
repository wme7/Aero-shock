// File for Playing with pointers by Manuel Diaz (C)

#include <iostream>
using namespace std;

void double_it(int *p) {
	*p = *p * 2;
}

int main() {
	int home = 10;
	int* phome = &home;
	*phome = 15;

	cout << "my home number is: " << *phome << endl;
	double_it(phome);
	cout << "the double of my home number is: " << *phome << endl;
	cout << "my home number is: " << home << endl;

return 0;
}
