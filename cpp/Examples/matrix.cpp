// Simple matrix multiplication

#include <iostream>
using namespace std;

int main() {
	int a[3][3], b[3][3], c[3][3];
	int i, j, k;
	
	cout << "Enter matrix A: \n";
	for(i=0; i<3; i++)
	 for(j=0; j<3; j++)
	  {
	  cout <<"A("<< i <<","<< j <<"): ";
	  cin >> a[i][j];
	  }

	cout << "Enter matrix B: \n";
	for(i=0; i<3; i++)
	 for(j=0; j<3; j++)
	  {
	  cout <<"B("<< i <<","<< j <<"): ";
	  cin >> b[i][j];
	  }

	for(i=0; i<3; i++)
	 for(j=0; j<3; j++)
	   for(k=0; k<3; k++)
	   c[i][j] += a[i][k] * b[k][j];
	  
	cout << "The resultant matrix is: \n";
	for(i=0; i<3; i++)
	 {
	 for(j=0; j<3; j++)
	 cout << c[i][j] <<" ";
	 cout << endl;
	 }
	return 0;
}
