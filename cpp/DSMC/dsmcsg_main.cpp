//===================================================================================================
//
// This program is a two-dimensional direct simulation Monte Carlo (DSMC) code based on Cartesian
// structured grids on single graphics processing unit (GPU).
//
//===================================================================================================

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cuda_runtime.h>
#include <cuda.h>
#include "parameter.h"
#include "dsmcsg_init.h"
#include "dsmcsg_cudafunction.h"
#include "dsmcsg_output.h"
#include "dsmcsg_class.h"
#include "dsmcsg_readfile.h"

using namespace std ;

int main( int argc , char* argv[] ){
	
	// set gpu device
	cudaSetDevice(0) ;

	ofstream	Outputfile ;
	DSMC_TIME	Time ;
	DSMC_DOMAIN	h_Domain , *d_Domain ;
	DSMC_NODE	*h_Node , *d_Node ;
	DSMC_CELL	*h_Cell , *d_Cell ;
	DSMC_BLOCK	*h_Block , *d_Block ;
	DSMC_OUTPUT	*h_Result ;
	int		TimestepNo ; 

	// Index data
	int		*h_IndexParticle , *h_IndexCell1 , *h_IndexCell2 ;

	// Particle data
	float		*h_ParticleXCoord , *h_ParticleYCoord ;
	float		*h_ParticleXVel , *h_ParticleYVel , *h_ParticleZVel ;
	int		*h_ParticleCellNo , *h_ParticleLastCollide ;
	int		*h_ParticleIn ;

	// Inlet Face data
	int		*h_InletCellNo , *h_InletEdgeNo ;
	float		*h_InletNum , *h_InletNumR , *h_InletXVel , *h_InletYVel , *h_InletZVel , *h_InletNumDen , *h_InletTemp , *h_InletArea ;

	// Sample data
	int		*h_SampleParticleNum ;
	float		*h_SampleXVel , *h_SampleYVel , *h_SampleZVel , *h_SampleXVelSq , *h_SampleYVelSq , *h_SampleZVelSq ;

	// Random Number
	RANDOM		*h_Rand , *d_Rand ;

	//----------------------------------------------------------------------------------------------------------------

	// Index data
	int		*d_IndexParticle , *d_IndexCell1 , *d_IndexCell2 ;

	// Particle data
	float		*d_ParticleXCoord , *d_ParticleYCoord ;
	float		*d_ParticleXVel , *d_ParticleYVel , *d_ParticleZVel ;
	int		*d_ParticleCellNo , *d_ParticleLastCollide ;
	int		*d_ParticleIn , *d_ParticleWrite ;

	// Inlet face data
	int		*d_InletCellNo , *d_InletEdgeNo ;
	float		*d_InletNum , *d_InletNumR , *d_InletXVel , *d_InletYVel , *d_InletZVel , *d_InletNumDen , *d_InletTemp , *d_InletArea ;

	// Sample data
	int		*d_SampleParticleNum ;
	float		*d_SampleXVel , *d_SampleYVel , *d_SampleZVel , *d_SampleXVelSq , *d_SampleYVelSq , *d_SampleZVelSq ;

	//----------------------------------------------------------------------------------------------------------------


	// Start DSMC Simulation
	Time.Start() ;


	// read simulation conditions (input.txt)
	ReadInput( &h_Domain ) ;


	// allocate memory on host
	h_Node		=	new DSMC_NODE[h_Domain.NodeNum] ;
	h_Cell		=	new DSMC_CELL[h_Domain.CellNum] ;
	h_Block		=	new DSMC_BLOCK[h_Domain.BlockNum] ;
	h_Result	=	new DSMC_OUTPUT[h_Domain.CellNum] ;

	h_IndexParticle	=	new int[h_Domain.MaxParticleNum] ;
	h_IndexCell1	=	new int[h_Domain.CellNum] ;
	h_IndexCell2	=	new int[h_Domain.CellNum] ;

	h_ParticleXCoord	=	new float[h_Domain.MaxParticleNum] ;
	h_ParticleYCoord	=	new float[h_Domain.MaxParticleNum] ;
	h_ParticleXVel		=	new float[h_Domain.MaxParticleNum] ;
	h_ParticleYVel		=	new float[h_Domain.MaxParticleNum] ;
	h_ParticleZVel		=	new float[h_Domain.MaxParticleNum] ;
	h_ParticleCellNo	=	new int[h_Domain.MaxParticleNum] ;
	h_ParticleLastCollide	=	new int[h_Domain.MaxParticleNum] ;
	h_ParticleIn		=	new int[h_Domain.MaxParticleNum] ;

	h_InletCellNo		=	new int[h_Domain.InletFaceNum] ;
	h_InletEdgeNo		=	new int[h_Domain.InletFaceNum] ;
	h_InletNum		=	new float[h_Domain.InletFaceNum] ;
	h_InletNumR		=	new float[h_Domain.InletFaceNum] ;
	h_InletXVel		=	new float[h_Domain.InletFaceNum] ;
	h_InletYVel		=	new float[h_Domain.InletFaceNum] ;
	h_InletZVel		=	new float[h_Domain.InletFaceNum] ;
	h_InletNumDen		=	new float[h_Domain.InletFaceNum] ;
	h_InletTemp		=	new float[h_Domain.InletFaceNum] ;
	h_InletArea		=	new float[h_Domain.InletFaceNum] ;

	h_SampleParticleNum	=	new int[h_Domain.CellNum] ;
	h_SampleXVel		=	new float[h_Domain.CellNum] ;
	h_SampleYVel		=	new float[h_Domain.CellNum] ;
	h_SampleZVel		=	new float[h_Domain.CellNum] ;
	h_SampleXVelSq		=	new float[h_Domain.CellNum] ;
	h_SampleYVelSq		=	new float[h_Domain.CellNum] ;
	h_SampleZVelSq		=	new float[h_Domain.CellNum] ;

	h_Rand			=	new RANDOM[8192*256] ;
	
	//----------------------------------------------------------------------------------------------------------------
	
	// allocate memory on device
	cudaMalloc((void**)&d_Domain		,	sizeof(DSMC_DOMAIN)) ;
	cudaMalloc((void**)&d_Node		,	sizeof(DSMC_NODE)*h_Domain.NodeNum) ;
	cudaMalloc((void**)&d_Cell		,	sizeof(DSMC_CELL)*h_Domain.CellNum) ;
	cudaMalloc((void**)&d_Block		,	sizeof(DSMC_BLOCK)*h_Domain.BlockNum) ;

	cudaMalloc((void**)&d_IndexParticle	,	sizeof(int)*h_Domain.MaxParticleNum) ;
	cudaMalloc((void**)&d_IndexCell1	,	sizeof(int)*h_Domain.CellNum) ;
	cudaMalloc((void**)&d_IndexCell2	,	sizeof(int)*h_Domain.CellNum) ;

	cudaMalloc((void**)&d_ParticleXCoord	,	sizeof(float)*h_Domain.MaxParticleNum) ;
	cudaMalloc((void**)&d_ParticleYCoord	,	sizeof(float)*h_Domain.MaxParticleNum) ;
	cudaMalloc((void**)&d_ParticleXVel	,	sizeof(float)*h_Domain.MaxParticleNum) ;
	cudaMalloc((void**)&d_ParticleYVel	,	sizeof(float)*h_Domain.MaxParticleNum) ;
	cudaMalloc((void**)&d_ParticleZVel	,	sizeof(float)*h_Domain.MaxParticleNum) ;
	cudaMalloc((void**)&d_ParticleCellNo	,	sizeof(int)*h_Domain.MaxParticleNum) ;
	cudaMalloc((void**)&d_ParticleLastCollide,	sizeof(int)*h_Domain.MaxParticleNum) ;
	cudaMalloc((void**)&d_ParticleIn	,	sizeof(int)*h_Domain.MaxParticleNum) ;
	cudaMalloc((void**)&d_ParticleWrite	,	sizeof(int)*h_Domain.MaxParticleNum) ;

	cudaMalloc((void**)&d_InletCellNo	,	sizeof(int)*h_Domain.InletFaceNum) ;
	cudaMalloc((void**)&d_InletEdgeNo 	,	sizeof(int)*h_Domain.InletFaceNum) ;
	cudaMalloc((void**)&d_InletNum 		,	sizeof(float)*h_Domain.InletFaceNum) ;
	cudaMalloc((void**)&d_InletNumR 	,	sizeof(float)*h_Domain.InletFaceNum) ;
	cudaMalloc((void**)&d_InletXVel 	,	sizeof(float)*h_Domain.InletFaceNum) ;
	cudaMalloc((void**)&d_InletYVel 	,	sizeof(float)*h_Domain.InletFaceNum) ;
	cudaMalloc((void**)&d_InletZVel 	,	sizeof(float)*h_Domain.InletFaceNum) ;
	cudaMalloc((void**)&d_InletNumDen 	,	sizeof(float)*h_Domain.InletFaceNum) ;
	cudaMalloc((void**)&d_InletTemp 	,	sizeof(float)*h_Domain.InletFaceNum) ;
	cudaMalloc((void**)&d_InletArea 	,	sizeof(float)*h_Domain.InletFaceNum) ;

	cudaMalloc((void**)&d_SampleParticleNum	,	sizeof(int)*h_Domain.CellNum) ;
	cudaMalloc((void**)&d_SampleXVel	,	sizeof(float)*h_Domain.CellNum) ;
	cudaMalloc((void**)&d_SampleYVel	,	sizeof(float)*h_Domain.CellNum) ;
	cudaMalloc((void**)&d_SampleZVel	,	sizeof(float)*h_Domain.CellNum) ;
	cudaMalloc((void**)&d_SampleXVelSq	,	sizeof(float)*h_Domain.CellNum) ;
	cudaMalloc((void**)&d_SampleYVelSq	,	sizeof(float)*h_Domain.CellNum) ;
	cudaMalloc((void**)&d_SampleZVelSq	,	sizeof(float)*h_Domain.CellNum) ;

	cudaMalloc((void**)&d_Rand		,	sizeof(RANDOM)*8192*256) ;


	cudaMemcpy(d_Rand , h_Rand , sizeof(RANDOM)*8192*256 , cudaMemcpyHostToDevice) ;

	cudaMemset(d_ParticleWrite , 0 , sizeof(int)*h_Domain.MaxParticleNum) ;


	// set initial value
	InitValue( h_IndexParticle , h_IndexCell1 , h_IndexCell2 , h_ParticleXCoord , h_ParticleYCoord , h_ParticleXVel , h_ParticleYVel , 
		   h_ParticleZVel , h_ParticleCellNo , h_ParticleLastCollide , h_ParticleIn , 
		   h_SampleParticleNum , h_SampleXVel , h_SampleYVel , h_SampleZVel , h_SampleXVelSq , h_SampleYVelSq , h_SampleZVelSq , 
		   h_InletCellNo , h_InletEdgeNo , h_InletNum , h_InletNumR , h_InletXVel , h_InletYVel , 
		   h_InletZVel , h_InletNumDen , h_InletTemp , h_InletArea , &h_Domain ) ;


	// set Boundary condition
	BoundaryCondition( &h_Domain , h_Block ) ;


	// create Mesh for simulation
	Time.Time() ;
	InitGrid( h_InletCellNo , h_InletEdgeNo , h_InletXVel , h_InletYVel , h_InletZVel , h_InletNumDen , h_InletTemp , h_InletArea , &h_Domain , h_Node ,
		  h_Cell , h_Block ) ;

	// setup initial particles
	InitFlow( h_ParticleXCoord , h_ParticleYCoord , h_ParticleXVel , h_ParticleYVel , h_ParticleZVel , h_ParticleCellNo , 
		  h_ParticleLastCollide , h_ParticleIn , &h_Domain , h_Cell ) ;


	// setup inlet face
	if ( h_Domain.InletFaceNum > 0 ){
		SetupInletFace( h_InletCellNo , h_InletEdgeNo , h_InletNum , h_InletNumR , h_InletXVel , h_InletYVel , h_InletZVel , h_InletNumDen , h_InletTemp , 
				h_InletArea , &h_Domain ) ;
	}


	// setup internal soild cell
	if ( h_Domain.InternalBlockNum > 0 ){
		SetInternalBlock( &h_Domain , h_Cell , h_Block ) ;
	}


	// output cell data.
	PrintCellInfo( h_Cell , h_Domain.CellNum ) ;
	PrintNeighborInfo( h_Cell , h_Domain.CellNum ) ;

	Time.Time(&Time.Init) ;


	// memory copied from host to device after initialization
	MemoryCopy( 0 , d_IndexParticle , h_IndexParticle , d_IndexCell1 , h_IndexCell1 , d_IndexCell2 , h_IndexCell2 , d_ParticleXCoord , 
			h_ParticleXCoord , d_ParticleYCoord , h_ParticleYCoord , d_ParticleXVel , h_ParticleXVel , d_ParticleYVel , 
			h_ParticleYVel , d_ParticleZVel , h_ParticleZVel , d_ParticleCellNo , h_ParticleCellNo , d_ParticleLastCollide , 
			h_ParticleLastCollide , d_ParticleIn , h_ParticleIn , d_SampleParticleNum , h_SampleParticleNum , d_SampleXVel , 
			h_SampleXVel , d_SampleYVel, h_SampleYVel , d_SampleZVel , h_SampleZVel , d_SampleXVelSq , h_SampleXVelSq , 
			d_SampleYVelSq , h_SampleYVelSq , d_SampleZVelSq , h_SampleZVelSq , d_InletCellNo , h_InletCellNo , d_InletEdgeNo , 
			h_InletEdgeNo , d_InletNum , h_InletNum , d_InletNumR , h_InletNumR , d_InletXVel , h_InletXVel , d_InletYVel , 
			h_InletYVel , d_InletZVel , h_InletZVel , d_InletNumDen , h_InletNumDen , d_InletTemp , h_InletTemp , d_InletArea , 
			h_InletArea , d_Domain , &h_Domain , d_Node , h_Node , d_Cell , h_Cell , d_Block , h_Block ) ;

	Outputfile.open( "particles.dat" , ios::out | ios::trunc ) ;

	
	for ( TimestepNo=1 ; TimestepNo<=h_Domain.TimestepNum ; TimestepNo++ ){
		cout << "TimestepNo:" << setw(20) << TimestepNo << " , ParticleNum:" << setw(20) << h_Domain.ParticleNum  << '\n' ;
		Outputfile << setw(20) << TimestepNo << setw(20) << h_Domain.ParticleNum << '\n' ;

	
		// particle movement
		Time.Time() ;
		ParticleMove( d_ParticleXCoord , d_ParticleYCoord , d_ParticleXVel , d_ParticleYVel , d_ParticleZVel , d_ParticleCellNo , 
			      d_ParticleLastCollide , d_ParticleIn , d_ParticleWrite , d_InletCellNo , d_InletEdgeNo , d_InletXVel , 
			      d_InletYVel , d_InletTemp , d_InletNum , d_InletNumR , &h_Domain , d_Domain , d_Node , d_Cell , d_Block , 
			      d_Rand ) ;
		Time.Time(&Time.Move) ;


		// index
		Time.Time() ;
		Index( d_IndexParticle , d_IndexCell1 , d_IndexCell2 , d_ParticleCellNo , h_Domain.ParticleNum , h_Domain.CellNum , &Time ) ;
		Time.Time(&Time.Index) ;


		// collision
		Time.Time() ;
		Collision( d_IndexParticle , d_IndexCell1 , d_IndexCell2 , d_ParticleLastCollide , d_ParticleXVel , d_ParticleYVel , 
			   d_ParticleZVel , d_SampleParticleNum , d_Domain , d_Cell , d_Rand ) ;
		Time.Time(&Time.Collision) ;
		

		if ( TimestepNo>h_Domain.SampleNo && (TimestepNo%2)==0 ){
			// sample particle properties into cell
			Time.Time() ;
			Sample( d_IndexParticle , d_IndexCell1 , d_IndexCell2 , d_SampleParticleNum , d_SampleXVel , d_SampleYVel , 
				d_SampleZVel , d_SampleXVelSq , d_SampleYVelSq , d_SampleZVelSq , d_ParticleXVel , d_ParticleYVel , 
				d_ParticleZVel , d_Cell , d_Domain ) ;
			cudaMemcpy( &h_Domain ,d_Domain , sizeof(DSMC_DOMAIN) , cudaMemcpyDeviceToHost ) ;
			Time.Time(&Time.Sample) ;


			if ( TimestepNo%8000 == 0 && TimestepNo != h_Domain.TimestepNum ){
				// momory copied from device to host for sampling data
				MemoryCopy( d_SampleParticleNum	, h_SampleParticleNum , d_SampleXVel , h_SampleXVel , d_SampleYVel , h_SampleYVel , 
					    d_SampleZVel , h_SampleZVel , d_SampleXVelSq , h_SampleXVelSq , d_SampleYVelSq , h_SampleYVelSq , 
					    d_SampleZVelSq , h_SampleZVelSq , d_Domain , &h_Domain ) ;

				// calculate macroscopic properties (simulation results)
				CalculateResult( h_SampleParticleNum , h_SampleXVel , h_SampleYVel , h_SampleZVel , h_SampleXVelSq , 
						 h_SampleYVelSq , h_SampleZVelSq , h_Result , &h_Domain , h_Cell ) ;
	
				// output macroscopic properties (simulation results)
				OutputResult( h_Cell , h_Result , h_Domain.CellNum , TimestepNo ) ;
			}
		}


		// ouput executed time (unit: sec.)
		cout << "Move:       " << setw(15) << Time.Move << '\n' ;
		cout << "Index:      " << setw(15) << Time.Index << '\n' ;
		cout << "Collision:  " << setw(15) << Time.Collision << '\n' ;
		cout << "Sample:     " << setw(15) << Time.Sample << '\n' ;
	}
	Outputfile.clear() ;
	Outputfile.close() ;

	// momory copied from device to host for sampling data
	MemoryCopy( d_SampleParticleNum	, h_SampleParticleNum , d_SampleXVel , h_SampleXVel , d_SampleYVel , h_SampleYVel , d_SampleZVel , 
		    h_SampleZVel , d_SampleXVelSq , h_SampleXVelSq , d_SampleYVelSq , h_SampleYVelSq , d_SampleZVelSq , h_SampleZVelSq , 
		    d_Domain , &h_Domain ) ;
	
	// calculate macroscopic properties (simulation results)
	CalculateResult( h_SampleParticleNum , h_SampleXVel , h_SampleYVel , h_SampleZVel , h_SampleXVelSq , h_SampleYVelSq , 
			 h_SampleZVelSq , h_Result , &h_Domain , h_Cell ) ;

	// output macroscopic properties (simulation results)
	OutputResult( h_Cell , h_Result , h_Domain.CellNum , (TimestepNo-1) ) ;


	// delete allocated memory on device
	cudaFree(d_Domain) ;
	cudaFree(d_Node) ;
	cudaFree(d_Cell) ;
	cudaFree(d_Block) ;
	
	cudaFree(d_IndexParticle) ;
	cudaFree(d_IndexCell1) ;
	cudaFree(d_IndexCell2) ;

	cudaFree(d_ParticleXCoord) ;
	cudaFree(d_ParticleYCoord) ;
	cudaFree(d_ParticleXVel) ;
	cudaFree(d_ParticleYVel) ;
	cudaFree(d_ParticleZVel) ;
	cudaFree(d_ParticleCellNo) ;
	cudaFree(d_ParticleLastCollide) ;
	cudaFree(d_ParticleIn) ;
	cudaFree(d_ParticleWrite) ; 

	cudaFree(d_InletCellNo) ;
	cudaFree(d_InletEdgeNo) ;
	cudaFree(d_InletNum) ;
	cudaFree(d_InletNumR) ;
	cudaFree(d_InletXVel) ;
	cudaFree(d_InletYVel) ;
	cudaFree(d_InletZVel) ;
	cudaFree(d_InletNumDen) ;
	cudaFree(d_InletTemp) ;
	cudaFree(d_InletArea) ;

	cudaFree(d_SampleParticleNum) ;
	cudaFree(d_SampleXVel) ;
	cudaFree(d_SampleYVel) ;
	cudaFree(d_SampleZVel) ;
	cudaFree(d_SampleXVelSq) ;
	cudaFree(d_SampleYVelSq) ;
	cudaFree(d_SampleZVelSq) ;
	
	cudaFree(d_Rand) ;


	// delete allocated meomry on host
	delete [] h_Node ;
	delete [] h_Cell ;
	delete [] h_Block ;
	delete [] h_Result ;

	delete [] h_IndexParticle ;
	delete [] h_IndexCell1 ;
	delete [] h_IndexCell2 ;

	delete [] h_ParticleXCoord ;
	delete [] h_ParticleYCoord ;
	delete [] h_ParticleXVel ;
	delete [] h_ParticleYVel ;
	delete [] h_ParticleZVel ;
	delete [] h_ParticleCellNo ;
	delete [] h_ParticleLastCollide ;
	delete [] h_ParticleIn ;

	delete [] h_InletCellNo ;
	delete [] h_InletEdgeNo ;
	delete [] h_InletNum ;
	delete [] h_InletNumR ;
	delete [] h_InletXVel ;
	delete [] h_InletYVel ;
	delete [] h_InletZVel ;
	delete [] h_InletNumDen ;
	delete [] h_InletTemp ;
	delete [] h_InletArea ;

	delete [] h_SampleParticleNum ;
	delete [] h_SampleXVel ;
	delete [] h_SampleYVel ;
	delete [] h_SampleZVel ;
	delete [] h_SampleXVelSq ;
	delete [] h_SampleYVelSq ;
	delete [] h_SampleZVelSq ;

	delete [] h_Rand ;


	// End DSMC Simulation
	Time.End() ;

	// output simulation time to file (SimTime.dat).
	OutputTime( Time ) ;

	cout << "==============================DSMC Simulation End==============================\n" ;
	cout << "Unit: second\n" ;
	cout << "Total Time: " << setw(15) << Time.Total << '\n' ;
	cout << "Init:       " << setw(15) << Time.Init << '\n' ;
	cout << "Move:       " << setw(15) << Time.Move << '\n' ;
	cout << "Index:      " << setw(15) << Time.Index << '\n' ;
	cout << "Collision:  " << setw(15) << Time.Collision << '\n' ;
	cout << "Sample:     " << setw(15) << Time.Sample << '\n' ;

	return	0 ;
}
