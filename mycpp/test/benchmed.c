
/*----------------------------------------------------------------------------
   Comparison of qsort-based and optimized median search methods
   Nicolas Devillard <ndevilla@free.fr> August 1998
   This code in public domain.
 ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

/*
 * A note about this benchmarking method:
 * the reported times indicate the actual time spent on CPU, they are
 * quite dependent on CPU load. However, the test is quite fast and it
 * is reasonable to assume a constant machine load throughout the test.
 */

#include <time.h>

/* Number of pixels in the array */
#define BIG_NUM (1024*1024)
/* #define BIG_NUM (1001) */

/* Generated numbers are in [0...MAX_ARRAY_VALUE-1] */
#define MAX_ARRAY_VALUE     1024    


/* Change this to whatever suits your needs */
typedef float pixelvalue ;


/* Macro to determine an integer's oddity */
#define odd(x) ((x)&1)

/*------------------- Functions: ANSI C prototypes -------------------------*/

int compare(const void *, const void*) ;

void pixel_qsort(pixelvalue *, int) ;

pixelvalue median_AHU(pixelvalue *, int) ;
pixelvalue select_k(int, pixelvalue *, int);

pixelvalue kth_smallest(pixelvalue *, int, int);
#define median_WIRTH(a,n) kth_smallest(a,n,(((n)&1)?((n)/2):(((n)/2)-1)))


pixelvalue torben(pixelvalue a[], int n) ;

pixelvalue quick_select(pixelvalue a[], int n) ;

void bench(int, size_t);


/*-------------------------- main() ----------------------------------------*/

int main(int argc, char * argv[]) 
{
    int     i ;
    int     from, to, step ;
    int     count ;

    if (argc<2) {
        printf("usage:\n") ;
        printf("%s <n>\n", argv[0]) ;
        printf("\tif n=1 the output is verbose for one attempt\n") ;
        printf("\tif n>1 the output reads:\n") ;
        printf("\t# of elements | method1 | method2 | ...\n") ;
        printf("\n") ;
        printf("%s <from> <to> <step>\n", argv[0]) ;
        printf("\twill loop over the number of elements in input\n") ;
        printf("\n") ;
        return 0 ;
    }

    if (argc==2) {
        count = atoi(argv[1]) ;
        if (count==1) {
            bench(1, BIG_NUM) ;
        } else {
            for (i=0 ; i<atoi(argv[1]) ; i++) {
                bench(0, BIG_NUM) ;
            }
        }
    } else if (argc==4) {
        from = atoi(argv[1]) ;
        to   = atoi(argv[2]) ;
        step = atoi(argv[3]) ;
        for (count=from ; count<=to ; count+=step) {
            bench(0, count) ;
        }
    }

    return 0 ;
}


#define N_METHODS   5

void bench(int detailed_printout, size_t array_size)
{
    int         i ;
    int         mednum ;
    pixelvalue  med[N_METHODS] ;

    clock_t     chrono ;
    double       elapsed ;

    pixelvalue  *   array_init,
                *   array ;

    /*
     * initialize random generator with PID
     * This is the only Unix-ish thing, replace this by any
     * initialization scheme you wish.
     */
    srand48(getpid()) ;

    if (detailed_printout) {
        printf("generating numbers...\n") ;
        fflush(stdout) ;
    }
    if (array_size<1) array_size = BIG_NUM ;

    if (detailed_printout) {
        printf("array size: %ld\n", (long)array_size) ;
    } else {
        printf("%ld\t", (long)array_size) ;
    }
    
    array_init = malloc(array_size * sizeof(pixelvalue)) ;
    array      = malloc(array_size * sizeof(pixelvalue)) ;

    if (array_init==NULL || array==NULL) {
        printf("memory allocation failure: aborting\n") ;
        return ;
    }

    for (i=0 ; i<array_size; i++) {
        array_init[i] = (pixelvalue)(lrand48() % MAX_ARRAY_VALUE) ;
    }
    mednum = 0 ;

    /* benchmark the quick select sort */
    memcpy(array, array_init, array_size * sizeof(pixelvalue)) ;
    if (detailed_printout) {
        printf("quick select    :\t") ;
        fflush(stdout) ;
    }
    chrono = clock() ;
    med[mednum] = quick_select(array, array_size) ;
    elapsed = (double)(clock() - chrono) / (double)CLOCKS_PER_SEC ;
    if (detailed_printout) {
        printf("%5.3f sec\t", elapsed) ;
        fflush(stdout) ;
        printf("med %g\n", (double)med[mednum]) ;
        fflush(stdout) ;
    } else {
        printf("%5.3f\t", elapsed) ;
        fflush(stdout) ;
    }
    mednum++ ;

    /* benchmark WIRTH */
    memcpy(array, array_init, array_size * sizeof(pixelvalue)) ;
    if (detailed_printout) {
        printf("WIRTH median    :\t") ;
        fflush(stdout) ;
    }
    chrono = clock() ;
    med[mednum] = median_WIRTH(array, array_size) ;
    elapsed = (double)(clock() - chrono) / (double)CLOCKS_PER_SEC ;
    if (detailed_printout) {
        printf("%5.3f sec\t", elapsed) ;
        fflush(stdout) ;
        printf("med %g\n", (double)med[mednum]) ;
        fflush(stdout) ;
    } else {
        printf("%5.3f\t", elapsed) ;
        fflush(stdout) ;
    }
    mednum++ ;

    /* benchmark the AHU sort */
    memcpy(array, array_init, array_size * sizeof(pixelvalue)) ;
    if (detailed_printout) {
        printf("AHU median      :\t") ;
        fflush(stdout) ;
    }
    chrono = clock() ;
    med[mednum] = median_AHU(array, array_size) ;
    elapsed = (double)(clock() - chrono) / (double)CLOCKS_PER_SEC ;
    if (detailed_printout) {
        printf("%5.3f sec\t", elapsed) ;
        fflush(stdout) ;
        printf("med %g\n", (double)med[mednum]) ;
        fflush(stdout) ;
    } else {
        printf("%5.3f\t", elapsed) ;
        fflush(stdout) ;
    }
    mednum++ ;

    /* benchmark torben's method */
    memcpy(array, array_init, array_size * sizeof(pixelvalue)) ;
    if (detailed_printout) {
        printf("torben          :\t") ;
        fflush(stdout) ;
    }
    chrono = clock() ;
    med[mednum] = torben(array, array_size) ;
    elapsed = (double)(clock() - chrono) / (double)CLOCKS_PER_SEC ;
    if (detailed_printout) {
        printf("%5.3f sec\t", elapsed) ;
        fflush(stdout) ;
        printf("med %g\n", (double)med[mednum]) ;
        fflush(stdout) ;
    } else {
        printf("%5.3f\t", elapsed) ;
        fflush(stdout) ;
    }
    mednum++ ;

    /* benchmark the eclipse fast pixel sort */
    memcpy(array, array_init, array_size * sizeof(pixelvalue)) ;
    if (detailed_printout) {
        printf("fast pixel sort :\t") ;
        fflush(stdout) ;
    }
    chrono = clock() ;
    pixel_qsort(array, array_size) ;
    elapsed = (double)(clock() - chrono) / (double)CLOCKS_PER_SEC ;

    if (odd(array_size)) {
        med[mednum] = array[array_size/2] ;
    } else {
        med[mednum] = array[(array_size/2) -1] ;
    }

    if (detailed_printout) {
        printf("%5.3f sec\t", elapsed) ;
        fflush(stdout) ;
        printf("med %g\n", (double)med[mednum]) ;
        fflush(stdout) ;
    } else {
        printf("%5.3f\t", elapsed) ;
        fflush(stdout) ;
    }
    mednum++ ;

    free(array) ;
    free(array_init) ;

    for (i=1 ; i<N_METHODS ; i++) {
        if (med[i-1]!=med[i]) {
            printf("diverging median values!\n") ;
            fflush(stdout) ;
        }
    }
    printf("\n") ;
    fflush(stdout) ;
    return ;
}

/*
 * This function only useful to the qsort() routine
 */

int compare(const void *f1, const void *f2)
{ return ( *(pixelvalue*)f1 > *(pixelvalue*)f2) ? 1 : -1 ; } 




/*----------------------------------------------------------------------------
   Function :   pixel_qsort()
   In       :   pixel array, size of the array
   Out      :   void
   Job      :   sort out the array of pixels
   Notice   :   optimized implementation, unreadable.
 ---------------------------------------------------------------------------*/


#define PIX_SWAP(a,b) { pixelvalue temp=(a);(a)=(b);(b)=temp; }
#define PIX_STACK_SIZE 50

void pixel_qsort(pixelvalue *pix_arr, int npix)
{
    int         i,
                ir,
                j,
                k,
                l;
    int *       i_stack ;
    int         j_stack ;
    pixelvalue  a ;

    ir = npix ;
    l = 1 ;
    j_stack = 0 ;
    i_stack = malloc(PIX_STACK_SIZE * sizeof(pixelvalue)) ;
    for (;;) {
        if (ir-l < 7) {
            for (j=l+1 ; j<=ir ; j++) {
                a = pix_arr[j-1];
                for (i=j-1 ; i>=1 ; i--) {
                    if (pix_arr[i-1] <= a) break;
                    pix_arr[i] = pix_arr[i-1];
                }
                pix_arr[i] = a;
            }
            if (j_stack == 0) break;
            ir = i_stack[j_stack-- -1];
            l  = i_stack[j_stack-- -1];
        } else {
            k = (l+ir) >> 1;
            PIX_SWAP(pix_arr[k-1], pix_arr[l])
            if (pix_arr[l] > pix_arr[ir-1]) {
                PIX_SWAP(pix_arr[l], pix_arr[ir-1])
            }
            if (pix_arr[l-1] > pix_arr[ir-1]) {
                PIX_SWAP(pix_arr[l-1], pix_arr[ir-1])
            }
            if (pix_arr[l] > pix_arr[l-1]) {
                PIX_SWAP(pix_arr[l], pix_arr[l-1])
            }
            i = l+1;
            j = ir;
            a = pix_arr[l-1];
            for (;;) {
                do i++; while (pix_arr[i-1] < a);
                do j--; while (pix_arr[j-1] > a);
                if (j < i) break;
                PIX_SWAP(pix_arr[i-1], pix_arr[j-1]);
            }
            pix_arr[l-1] = pix_arr[j-1];
            pix_arr[j-1] = a;
            j_stack += 2;
            if (j_stack > PIX_STACK_SIZE) {
                printf("stack too small in pixel_qsort: aborting");
                exit(-2001) ;
            }
            if (ir-i+1 >= j-l) {
                i_stack[j_stack-1] = ir;
                i_stack[j_stack-2] = i;
                ir = j-1;
            } else {
                i_stack[j_stack-1] = j-1;
                i_stack[j_stack-2] = l;
                l = i;
            }
        }
    }
    free(i_stack) ;
}
#undef PIX_STACK_SIZE
#undef PIX_SWAP



/*---------------------------------------------------------------------------
   Function :   select_k()
   In       :   list of pixelvalues, # of values, kth smallest element
                which value is searched for (median for k=n/2)
   Out      :   pixelvalue
   Job      :   find out the kth smallest value of the list 
   Notice   :   recursively called by median_AHU()
 ---------------------------------------------------------------------------*/


pixelvalue select_k(int k, pixelvalue * list, int n)
{
    int             n1 = 0,
                    n2 = 0,
                    n3 = 0 ;
    pixelvalue  *   S ;
    int             i, j ;
    pixelvalue      p ;

    if (n==1) return list[0] ;
    p = list[(n>>1)] ;
    for (i=0 ; i<n ; i++) {
        if (list[i]<p) {
            n1++ ;
        } else if (list[i]==p) {
            n2++ ;
        } else {
            n3++ ;
        }
    }
    if (n1>=k) {
        S = malloc(n1*sizeof(pixelvalue)) ;
        j = 0 ;
        for (i=0 ; i<n ; i++) {
            if (list[i]<p) S[j++] = list[i] ;
        }
        p = select_k(k, S, n1) ;
        free(S) ;
    } else {
        if ((n1+n2)<k) {
            S = malloc(n3 * sizeof(pixelvalue)) ;
            j = 0 ;
            for (i=0 ; i<n ; i++) {
                if (list[i]>p) S[j++] = list[i] ;
            }
            p = select_k(k-n1-n2, S, n3) ;
            free(S) ;
        }
    }
    return p ;
}


pixelvalue median_AHU(pixelvalue * list, int n)
{
    if (odd(n)) {
        return select_k((n/2)+1, list, n) ;
    } else {
        return select_k((n/2), list, n) ;
    }
}

/*---------------------------------------------------------------------------
   Function :   kth_smallest()
   In       :   array of elements, # of elements in the array, rank k
   Out      :   one element
   Job      :   find the kth smallest element in the array
   Notice   :   use the median_WIRTH() macro to get the median. 

                Reference:

                  Author: Wirth, Niklaus 
                   Title: Algorithms + data structures = programs 
               Publisher: Englewood Cliffs: Prentice-Hall, 1976 
    Physical description: 366 p. 
                  Series: Prentice-Hall Series in Automatic Computation 

 ---------------------------------------------------------------------------*/


#define PIX_SWAP(a,b) { register pixelvalue t=(a);(a)=(b);(b)=t; }

pixelvalue kth_smallest(pixelvalue a[], int n, int k)
{
    register i,j,l,m ;
    register pixelvalue x ;

    l=0 ; m=n-1 ;
    while (l<m) {
        x=a[k] ;
        i=l ;
        j=m ;
        do {
            while (a[i]<x) i++ ;
            while (x<a[j]) j-- ;
            if (i<=j) {
                PIX_SWAP(a[i],a[j]) ;
                i++ ; j-- ;
            }
        } while (i<=j) ;
        if (j<k) l=i ;
        if (k<i) m=j ;
    }
    return a[k] ;
}

#undef PIX_SWAP

pixelvalue torben(pixelvalue m[], int n)
{
    int         i, less, greater, equal, half;
    pixelvalue  min, max, guess, maxltguess, mingtguess;

    half = (n+1)/2 ;
    min = max = m[0] ;
    for (i=1 ; i<n ; i++) {
        if (m[i]<min) min=m[i];
        if (m[i]>max) max=m[i];
    }

    while (1) {
        guess = (min+max)/2;
        less = 0; greater = 0; equal = 0;
        maxltguess = min ;
        mingtguess = max ;
        for (i=0; i<n; i++) {
            if (m[i]<guess) {
                less++;
                if (m[i]>maxltguess) maxltguess = m[i] ;
            } else if (m[i]>guess) {
                greater++;
                if (m[i]<mingtguess) mingtguess = m[i] ;
            } else equal++;
        }
        if (less <= half && greater <= half) break ;
        else if (less>greater) max = maxltguess ;
        else min = mingtguess; 
    }
    if (less >= half) return maxltguess;
    else if (less+equal >= half) return guess;
    else return mingtguess;
}   

#define PIX_SWAP(a,b) { pixelvalue temp=(a);(a)=(b);(b)=temp; }

pixelvalue quick_select(pixelvalue a[], int n) 
{
    int low, high ;
    int median;
    int middle, ll, hh;

    low = 0 ; high = n-1 ; median = (low + high) / 2;
    for (;;) {
        if (high <= low) /* One element only */
            return a[median] ;

        if (high == low + 1) {  /* Two elements only */
            if (a[low] > a[high])
                PIX_SWAP(a[low], a[high]) ;
            return a[median] ;
        }

    /* Find median of low, middle and high items; swap into position low */
        middle = (low + high) / 2;
        if (a[middle] > a[high])    PIX_SWAP(a[middle], a[high]) ;
        if (a[low] > a[high])       PIX_SWAP(a[low], a[high]) ;
        if (a[middle] > a[low])     PIX_SWAP(a[middle], a[low]) ;

    /* Swap low item (now in position middle) into position (low+1) */
        PIX_SWAP(a[middle], a[low+1]) ;

    /* Nibble from each end towards middle, swapping items when stuck */
        ll = low + 1;
        hh = high;
        for (;;) {
            do ll++; while (a[low] > a[ll]) ;
            do hh--; while (a[hh]  > a[low]) ;

            if (hh < ll)
            break;

            PIX_SWAP(a[ll], a[hh]) ;
        }
        
        /* Swap middle item (in position low) back into correct position */
        PIX_SWAP(a[low], a[hh]) ;
        
        /* Re-set active partition */
        if (hh <= median)
            low = ll;
        if (hh >= median)
            high = hh - 1;
    }
}   
#undef PIX_SWAP





