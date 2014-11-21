#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include "dsmcsg_class.h"

using namespace std ;


DSMC_DOMAIN::DSMC_DOMAIN(){
	TimestepNum 	=	0 ;
	SampleNo	=	0 ;
	ParticleNum	=	0 ;
	MaxParticleNum	=	0 ;

	NodeNum		=	0 ;
	CellNum 	=	0 ;
	XNodeNum 	=	0 ;
	YNodeNum 	=	0 ;
	XCellNum	=	0 ;
	YCellNum	=	0 ;
	BlockNum	=	0 ;
	InternalBlockNum=	0 ;
	InletFaceNum	=	0 ;

	SampleNum	=	0 ;
	XCellSize	=	0. ;
	YCellSize	=	0. ;
	CellVolume	=	0. ;
	Timestep	=	0. ;
	ParticleWeight	=	0. ;
	WeightRatio	=	0. ;
	XVel		=	0. ;
	YVel		=	0. ;
	ZVel		=	0. ;
	Temp		=	0. ;
	Density		=	0. ;
	NumDensity	=	0. ;
	XL		=	0. ;
	XH		=	0. ;
	YL		=	0. ;
	YH		=	0. ;
}

//------------------------------------------------------------------------------------------------------------

void DSMC_DOMAIN::Dump(){
	cout << "TimestepNum	" << setw(15) << TimestepNum 	<< '\n' ;
	cout << "SampleNo	" << setw(15) << SampleNo	<< '\n' ;
	cout << "ParticleNum	" << setw(15) << ParticleNum	<< '\n' ;
	cout << "MaxParticleNum	" << setw(15) << MaxParticleNum	<< '\n' ;
	cout << "NodeNum	" << setw(15) << NodeNum	<< '\n' ;
	cout << "CellNum	" << setw(15) << CellNum 	<< '\n' ;
	cout << "XNodeNum	" << setw(15) << XNodeNum 	<< '\n' ;
	cout << "YNodeNum	" << setw(15) << YNodeNum 	<< '\n' ;
	cout << "XCellNum	" << setw(15) << XCellNum	<< '\n' ;
	cout << "YCellNum	" << setw(15) << YCellNum	<< '\n' ;
	cout << "BlockNum	" << setw(15) << BlockNum	<< '\n' ;
	cout << "InletFaceNum	" << setw(15) << InletFaceNum	<< '\n' ;
	cout << "BCBlockNum	" << setw(15) << BCBlockNum	<< '\n' ;
	cout << "InternalBlockNum" << setw(15) << InternalBlockNum<< '\n' ;
	cout << "SampleNum	" << setw(15) << SampleNum	<< '\n' ;
	cout << "XCellSize	" << setw(15) << XCellSize	<< '\n' ;
	cout << "YCellSize	" << setw(15) << YCellSize	<< '\n' ;
	cout << "CellVolume	" << setw(15) << CellVolume	<< '\n' ;
	cout << "Timestep	" << setw(15) << Timestep	<< '\n' ;
	cout << "ParticleWeight	" << setw(15) << ParticleWeight	<< '\n' ;
	cout << "WeightRatio	" << setw(15) << WeightRatio	<< '\n' ;
	cout << "XVel		" << setw(15) << XVel		<< '\n' ;
	cout << "YVel		" << setw(15) << YVel		<< '\n' ;
	cout << "ZVel		" << setw(15) << ZVel		<< '\n' ;
	cout << "Temp		" << setw(15) << Temp		<< '\n' ;
	cout << "Density	" << setw(15) << Density	<< '\n' ;
	cout << "NumDensity	" << setw(15) << NumDensity	<< '\n' ;
	cout << "XL		" << setw(15) << XL		<< '\n' ;
	cout << "XH		" << setw(15) << XH		<< '\n' ;
	cout << "YL		" << setw(15) << YL		<< '\n' ;
	cout << "YH		" << setw(15) << YH		<< '\n' ;
}

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

DSMC_NODE::DSMC_NODE(){
	XCoord	=	0. ;
	YCoord	=	0. ;
}

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

DSMC_CELL::DSMC_CELL(){
	Type	=	0 ;
	XCenter		=	0. ;
	YCenter		=	0. ;
	MaxCrossSectionSpeed	=	0. ;
	RemainderCollisionPair	=	0. ;

	for ( int i=0 ; i<4 ; i++ ){
		Node[i]		=	0 ;
		Neighbor[i]	=	-999 ;
		EdgeFA[i]	=	0. ;
		EdgeFB[i]	=	0. ;
		EdgeFC[i]	=	0. ;
	}
}

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

DSMC_BLOCK::DSMC_BLOCK(){
	Type		=	0 ;
	XVel		=	0. ;
	YVel		=	0. ;
	NumDensity	=	0. ;
	Temp		=	0. ;
	XL		=	0. ;
	XH		=	0. ;
	YL		=	0. ;
	YH		=	0. ;
	NodeNum		=	0 ;
}

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

DSMC_OUTPUT::DSMC_OUTPUT(){
	NumDensity	=	0. ;
	Density		=	0. ;
	XVel		=	0. ;
	YVel		=	0. ;
	ZVel		=	0. ;
	Temp		=	0. ;
	AveParticleNum	=	0. ;
}

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

DSMC_TIME::DSMC_TIME(){
	Total	=	0. ;
	Init	=	0. ;
	Move	=	0. ;
	Index	=	0. ;
	Collision=	0. ;
	Sample	=	0. ;
	time	=	0. ;
	tstart	=	0. ;
}

//------------------------------------------------------------------------------------------------------------

void DSMC_TIME::Start(){
	gettimeofday(&start,0) ;
	tstart = (long double)(start.tv_sec+start.tv_usec*1e-6) ;
}

//------------------------------------------------------------------------------------------------------------

void DSMC_TIME::End(){
	gettimeofday(&end,0) ;
	Total = (long double)(end.tv_sec+end.tv_usec*1e-6) - tstart ; 
}

//------------------------------------------------------------------------------------------------------------

void DSMC_TIME::Time(){
	struct timeval	tv ;
	gettimeofday(&tv,0) ;
	time = (long double)(tv.tv_sec+tv.tv_usec*1e-6) ;
}

//------------------------------------------------------------------------------------------------------------

void DSMC_TIME::Time( long double *t ){
	struct timeval	tv ;
	gettimeofday(&tv,0) ;
	(*t) += (long double)(tv.tv_sec+tv.tv_usec*1e-6) - time ;  
}

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

RANDOM::RANDOM(){
	inext	=	0 ;
	inextp	=	0 ;
	iff	=	0 ;

	for ( int i=0 ; i<56 ; i++ )
		ma[i]   =       0 ;
}

//------------------------------------------------------------------------------------------------------------

float randn(){
	float		random ;

	random		=	(float)rand()/(float)RAND_MAX ;

	return random*(1. - 1.e-10) + 1.e-10 ;
}

//------------------------------------------------------------------------------------------------------------

void RandVel( float *u , float *v , float most_probable_speed ){
	float		a ;
	float		b ;

	(*u)	=	0. ;
	(*v)	=	0. ;
	
	a	=	sqrt(-log( randn() )) ;
	b	=	6.283185308 * randn() ;

	(*u)	=	a*sin(b)*most_probable_speed ;
	(*v)	=	a*cos(b)*most_probable_speed ;
}