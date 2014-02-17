// Median algorithm by Manuel Diaz (C)

#include <iostream>
using namespace std;

/*
 * Algorithm from N. Wirth's book, implementation by N. Devillard.
 * This code in public domain.
 */

typedef float elem_type ;

#define ELEM_SWAP(a,b) { register elem_type t=(a);(a)=(b);(b)=t; }

elem_type kth_smallest(elem_type a[], int n, int k)
{
    // register i,j,l,m ;
	int i, j, l, m;
    register elem_type x ;

    l=0 ; m=n-1 ;
    while (l<m) {
        x=a[k] ;
        i=l ;
        j=m ;
        do {
            while (a[i]<x) i++ ;
            while (x<a[j]) j-- ;
            if (i<=j) {
                ELEM_SWAP(a[i],a[j]) ;
                i++ ; j-- ;
            }
        } while (i<=j) ;
        if (j<k) l=i ;
        if (k<i) m=j ;
    }
    return a[k] ;
}

#define median(a,n) kth_smallest(a,n,(((n)&1)?((n)/2):(((n)/2)-1)))

int main() {
	elem_type array[] = { 5, 4, 6, 7, 8, 9, 6, 11, 7};
	int elements = 9;
	cout <<" data loaded " << endl;
	cout << median(array,elements) << endl;

return 0;
}
