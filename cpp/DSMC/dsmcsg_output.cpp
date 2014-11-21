//===================================================================================================
//
// These functions are used to output DSMC simulation result (macroscopic properties) and 
// computational time.
//
//===================================================================================================

#include <iostream>
#include <fstream>
#include <iomanip>
#include "dsmcsg_output.h"


using namespace std ;


// Output simulation results (Result-XXXXX.dat -- Note: XXXXX is timestep number).
void OutputResult( DSMC_CELL *h_Cell , DSMC_OUTPUT *h_Result , int CellNum , 	int TimestepNo ){

	ofstream		OutFile ;
	string			FileName , buffer1 ;
	char			buffer2[8] ;


	sprintf( buffer2 , "%d" , TimestepNo ) ;
	buffer1		=	buffer2 ;
	FileName	=	"Result-"+buffer1+".dat" ;

	OutFile.open( FileName.c_str() , ios::out | ios::trunc ) ;


	if ( !OutFile.fail() ){
		OutFile << setw(20) << "CellNo" << setw(20) << "XCenter" << setw(20) << "YCenter" << setw(20) << "ZCenter" << setw(20) << "Density"
			<< setw(20) << "NumDensity" << setw(20) << "U-Vel" << setw(20) << "V-Vel" << setw(20) << "W-Vel" << setw(20) << "Temp"  
			<< setw(20) << "AveParticlesPerCell" << '\n' ;

		for ( int i=0 ; i<CellNum ; i++ ){
			OutFile << setw(20) << (i+1) << setw(20) << h_Cell[i].XCenter << setw(20) << h_Cell[i].YCenter << setw(20) << "0" << setw(20) << h_Result[i].Density
				<< setw(20) << h_Result[i].NumDensity << setw(20) << h_Result[i].XVel << setw(20) << h_Result[i].YVel << setw(20) << h_Result[i].ZVel 
				<< setw(20) << h_Result[i].Temp << setw(20) << h_Result[i].AveParticleNum << '\n' ;
		}
	}else{
		cout << "Fail to open " << FileName << '\n' ;
	}

	OutFile.clear() ;
	OutFile.close() ;
}

//------------------------------------------------------------------------------------------------------------

// Output simulation time.
void OutputTime( DSMC_TIME Time ){
	ofstream		OutFile ;

	OutFile.open( "SimTime.dat" , ios::out | ios::trunc ) ;

	if ( !OutFile.fail() ){
		OutFile << setw(20) << "DSMC Simulation Time: (Unit: second)\n" ;
		OutFile << setw(20) << "Total Time:"	<< setw(20) << Time.Total 	<< '\n' ;
		OutFile << setw(20) << "Initial:"	<< setw(20) << Time.Init 	<< '\n' ;
		OutFile << setw(20) << "Move:"		<< setw(20) << Time.Move 	<< '\n' ;
		OutFile << setw(20) << "Index:"		<< setw(20) << Time.Index 	<< '\n' ;
		OutFile << setw(20) << "Collision:"	<< setw(20) << Time.Collision 	<< '\n' ;
		OutFile << setw(20) << "Sample:"	<< setw(20) << Time.Sample 	<< '\n' ;
	}else{
		cout << "Fail to open SimTime.dat\n" ;
	}

	OutFile.clear() ;
	OutFile.close() ;
}