
/*----------------------------------------------------------------------------
   Comparison of qsort-based and optimized median search
   Nicolas Devillard <ndevilla@free.fr> August 1997
   This code in public domain.
 ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*
 * A note about this benchmarking method:
 * the reported times indicate the actual time spent on CPU, they are
 * quite dependent on CPU load. However, the test is quite fast and it
 * is reasonable to assume a constant machine load throughout the test.
 */

#include <time.h>

/* Number of sorts to perform */
#define BIG_NUM 100000

/* Change this to whatever suits your needs */
typedef double pixelvalue ;

/* Global variable: contains the arrays to be sorted */
pixelvalue   array[BIG_NUM][9] ;

/*------------------- Functions: ANSI C prototypes -------------------------*/

int compare(const void *, const void*) ;
pixelvalue median_sort9(pixelvalue *) ;
pixelvalue opt_med9(pixelvalue *) ;

/*-------------------------- main() ----------------------------------------*/

int main(int argc, char *argv[])
{
    pixelvalue  med;
    int         i, j ;

    clock_t     chrono ;
    double       elapsed ;

    /* initialize random generator with PID */
    srand48(getpid()) ;
    printf("generating numbers...\n") ;
    for (i=0 ; i<BIG_NUM; i++)
        for (j=0 ; j<9 ; j++)
            array[i][j] = (pixelvalue)(lrand48() % 1024) ;

    /* benchmark the qsort median */
    printf("qsort median:\n") ;
    chrono = clock() ;
    for (i=0 ; i<BIG_NUM; i++) {
        med = median_sort9(array[i]) ;
    }
    elapsed = (double)(clock() - chrono) / (double)CLOCKS_PER_SEC ;
    printf("elapsed time: %g seconds\n", elapsed) ;

    printf("generating numbers...\n") ;
    for (i=0 ; i<BIG_NUM; i++)
        for (j=0 ; j<9 ; j++)
            array[i][j] = lrand48() % 1024 ;

    /* benchmark the fast median */
    chrono = clock() ;
    printf("fast median:\n") ;
    for (i=0 ; i<BIG_NUM; i++) {
        med = opt_med9(array[i]) ;
    }
    elapsed = (double)(clock() - chrono) / (double)CLOCKS_PER_SEC ;
    printf("elapsed time: %g seconds\n", elapsed) ;

    return 0 ;
}

/*
 * This function only useful to the qsort() routine
 */

int compare(const void *f1, const void *f2)
{ return ( *(pixelvalue*)f1 > *(pixelvalue*)f2) ? 1 : -1 ; } 

/*
 * Standard median search: by use of the qsort() function
 */

pixelvalue median_sort9(pixelvalue *array)
{
    qsort(array, 9, sizeof(pixelvalue), compare) ;
    return array[4] ;
}

/*
 * Optimized median search on 9 values
 */

#define PIX_SORT(a,b) { if ((a)>(b)) PIX_SWAP((a),(b)); }
#define PIX_SWAP(a,b) { pixelvalue temp=(a);(a)=(b);(b)=temp; }

pixelvalue opt_med9(pixelvalue * p)
{
    PIX_SORT(p[1], p[2]) ; PIX_SORT(p[4], p[5]) ; PIX_SORT(p[7], p[8]) ; 
    PIX_SORT(p[0], p[1]) ; PIX_SORT(p[3], p[4]) ; PIX_SORT(p[6], p[7]) ; 
    PIX_SORT(p[1], p[2]) ; PIX_SORT(p[4], p[5]) ; PIX_SORT(p[7], p[8]) ; 
    PIX_SORT(p[0], p[3]) ; PIX_SORT(p[5], p[8]) ; PIX_SORT(p[4], p[7]) ; 
    PIX_SORT(p[3], p[6]) ; PIX_SORT(p[1], p[4]) ; PIX_SORT(p[2], p[5]) ; 
    PIX_SORT(p[4], p[7]) ; PIX_SORT(p[4], p[2]) ; PIX_SORT(p[6], p[4]) ; 
    PIX_SORT(p[4], p[2]) ; return(p[4]) ;
}

#undef PIX_SWAP
#undef PIX_SORT


/*------------------------------ end of file -------------------------------*/
