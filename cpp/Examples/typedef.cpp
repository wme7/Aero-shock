// test file

#include <iostream>
using namespace std;

#define MATRIX_SIZE 4

typedef struct {
	int width;
	int height;
	int stride;
	float* elements;
	} Matrix;

int main()
{
	Matrix A;
	A.width = MATRIX_SIZE;
	A.height = MATRIX_SIZE;
	A.stride = A.width;
	for (int i=0; i<A.width*A.height ;i++)
	{	
		A.elements[i]=1; //random();
	}

	cout << A.width << endl;
	cout << A.height << endl;
	cout << A.stride << endl;
	for (int j=0; j<A.width*A.height ;j++)
	{	
		cout << A.elements[j] <<" "; //random();
	}
	cout << endl;

return 0;
}

