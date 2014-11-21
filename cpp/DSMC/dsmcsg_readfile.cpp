//===================================================================================================
//
// This function is uesd to read simulation conditions and initial conditions from the input file 
// "input.txt".
//
//===================================================================================================


#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include "dsmcsg_readfile.h"


using namespace std ;


// Read initial condition and simulation condition (input.txt).
void ReadInput( DSMC_DOMAIN *Domain ){
	int			i ;
	ifstream		Input ;
	string			get_line ;
	string			word ;
	string			number ;
	int			wordstart ;
	int			wordend ;

	Input.open( "input.txt" , ios::in ) ;

	if ( !Input ){
		cout << "Failed to open file input.txt" << endl ;
		exit(1) ;
	}

	while ( getline(Input, get_line) ){
		if ( get_line[0] != '#' ){
			wordend = get_line.find(' ') ;
			word = get_line.substr(0,wordend) ;


			if ( word == "SAMPLING_TIME"){
				wordstart = get_line.find('=') ;
				number = get_line.substr(wordstart+2) ;

				Domain->SampleNo = atoi( number.c_str() ) ;

			}else if ( word == "NUMBER_OF_TIMESTEP"){
				wordstart = get_line.find('=') ;
				number = get_line.substr(wordstart+2) ;

				Domain->TimestepNum = atoi( number.c_str() ) ;

			}else if ( word == "PARTICLE_PER_CELL"){
				wordstart = get_line.find('=') ;
				number = get_line.substr(wordstart+2) ;

				Domain->WeightRatio = atof( number.c_str() ) ;

			}else if ( word == "MAX_PARTICLE_NUMBER"){
				wordstart = get_line.find('=') ;
				number = get_line.substr(wordstart+2) ;

				Domain->MaxParticleNum = atoi( number.c_str() ) ;

			}else if ( word == "X_LOWER"){
				wordstart = get_line.find('=') ;
				number = get_line.substr(wordstart+2) ;

				Domain->XL = atof( number.c_str() ) ;

			}else if ( word == "X_HIGHER" ){
				wordstart = get_line.find('=') ;
				number = get_line.substr(wordstart+2) ;

				Domain->XH = atof( number.c_str() ) ;

			}else if ( word == "Y_LOWER" ){
				wordstart = get_line.find('=') ;
				number = get_line.substr(wordstart+2) ;

				Domain->YL = atof( number.c_str() ) ;

			}else if ( word == "Y_HIGHER" ){
				wordstart = get_line.find('=') ;
				number = get_line.substr(wordstart+2) ;

				Domain->YH = atof( number.c_str() ) ;

			}else if ( word == "X_CELL_NUMBER" ){
				wordstart = get_line.find('=') ;
				number = get_line.substr(wordstart+2) ;

				Domain->XCellNum = atoi( number.c_str() ) ;

			}else if ( word == "Y_CELL_NUMBER" ){
				wordstart = get_line.find('=') ;
				number = get_line.substr(wordstart+2) ;

				Domain->YCellNum = atoi( number.c_str() ) ;

			}else if ( word == "X_VELOCITY" ){
				wordstart = get_line.find('=') ;
				number = get_line.substr(wordstart+2) ;

				Domain->XVel = atof( number.c_str() ) ;

			}else if ( word == "Y_VELOCITY" ){
				wordstart = get_line.find('=') ;
				number = get_line.substr(wordstart+2) ;

				Domain->YVel = atof( number.c_str() ) ;

			}else if ( word == "Z_VELOCITY" ){
				wordstart = get_line.find('=') ;
				number = get_line.substr(wordstart+2) ;

				Domain->ZVel = atof( number.c_str() ) ;

			}else if ( word == "NUMBER_DENSITY" ){
				wordstart = get_line.find('=') ;
				number = get_line.substr(wordstart+2) ;

				Domain->NumDensity = atof( number.c_str() ) ;

			}else if ( word == "TEMPERATURE" ){
				wordstart = get_line.find('=') ;
				number = get_line.substr(wordstart+2) ;

				Domain->Temp = atof( number.c_str() ) ;

			}else if ( word == "BOUNDARY_BLOCK_NUMBER" ){
				wordstart = get_line.find('=') ;
				number = get_line.substr(wordstart+2) ;

				Domain->BCBlockNum = atoi( number.c_str() ) ;

			}else if ( word == "INLET_CELL_FACE_NUMBER" ){
				wordstart = get_line.find('=') ;
				number = get_line.substr(wordstart+2) ;

				Domain->InletFaceNum = atoi( number.c_str() ) ;

			}else if ( word == "INTERNAL_BLOCK_NUMBER" ){
				wordstart = get_line.find('=') ;
				number = get_line.substr(wordstart+2) ;

				Domain->InternalBlockNum = atoi( number.c_str() ) ;		
			}
		}
	}


	// To close the input.txt file.
	Input.clear() ;
	Input.close() ;


	Domain->XNodeNum	=	Domain->XCellNum + 1 ;
	Domain->YNodeNum	=	Domain->YCellNum + 1 ;
	Domain->NodeNum		=	Domain->XNodeNum*Domain->YNodeNum ;
	Domain->CellNum		=	Domain->XCellNum*Domain->YCellNum ;

	Domain->XCellSize	=	(Domain->XH-Domain->XL )/Domain->XCellNum ;
	Domain->YCellSize	=	(Domain->YH-Domain->YL )/Domain->YCellNum ;
	Domain->CellVolume	=	Domain->XCellSize*Domain->YCellSize ;

	Domain->BlockNum	=	Domain->BCBlockNum + Domain->InternalBlockNum ;
}
