#include "dsmcsg_class.h"

#if !defined(__DSMCSG_OUTPUT_H)
#define __DSMCSG_OUTPUT_H


// Output simulation results (Result-XXXXX.dat -- Note: XXXXX is timestep number).
void OutputResult( DSMC_CELL *h_Cell , DSMC_OUTPUT *h_Result , int CellNum , 	int TimestepNo ) ;


// Output simulation time.
void OutputTime( DSMC_TIME Time ) ;


#endif
