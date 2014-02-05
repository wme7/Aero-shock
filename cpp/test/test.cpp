// test file

#include <iostream>
using namespace std;

void double_it(int *p) {

	*p = *p * 2;

}

int main() {
	int home[] = {9, 8, 7, 6, 5, 4, 3, 2, 1};
	int *phome = &home;
	cout << *phome << endl;
	double_it(phome);
	cout << *phome << endl;

return 0;
}
