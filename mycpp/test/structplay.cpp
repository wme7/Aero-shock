// Learning to use struct in c++ by Manuel Diaz (C)

#include <iostream>
#include <string>
#include <sstream>
using namespace std;

struct product {
	int weigth;
	float price;
	string name;
	} prod1, prod2, prod3;

void printprod ( product fruit );

int main() 
{
	string temp;
	cout << "2 products are in memory, please enter one more: \n";	

	// self entered
	prod1.name = "banana";
	prod1.price = 1.50;
	prod1.weigth = 1;

	prod2.name = "apples";
	prod2.price = 3.25;
	prod2.weigth = 2;
	
	// ask user to enter one:
	cout << "Enter product name: \n";
	getline (cin,prod3.name);
	cout << "Enter product price (UD$): \n";
	getline (cin,temp); //careful!!! getline get you a string, not numbers!
	stringstream(temp) >> prod3.price;

	cout << "Enter product weigth (kg): \n";
/*	getline (cin,temp);
	stringstream(temp) >> prod3.weigth;
*/
	cin >> prod3.weigth;	// Easier!

	// Print those 3 products
	cout << "list of things to buy: \n";
	cout << "Name | " << "USD | " << "kg" << endl;
	printprod ( prod1 );
	printprod ( prod2 );
	printprod ( prod3 );

return 0;
}

void printprod ( product fruit ) {
	cout << fruit.name <<" "<< fruit.price <<" "<< fruit.weigth <<endl;
}
