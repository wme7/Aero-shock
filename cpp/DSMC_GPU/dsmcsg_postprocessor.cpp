//===================================================================================================
//
// This program is used to convert simulation results (Result.dat) to Tecplot format (Result-tec.dat) 
// for visualization.
//
//===================================================================================================

#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>


using namespace std ;

class NODE{
	public:
		int		Id ;
		double		Coord[3] ;

		NODE() ;
} ;

NODE::NODE(){
	Id	=	0 ;
	for ( int i=0 ; i<3 ; i++ )
		Coord[i]	=	0. ;
}


class CELL{
	public:
		int		Id , Type ;
		NODE		*Node[4] ;
		double		Coord[3] , Density , NumDensity , Velocity[3] , Temp ;
		double		AveParticleNum ;

		CELL() ;
} ;


CELL::CELL(){
	Id	=	0 ;
	Type	=	0 ;
	for ( int i=0 ; i<4 ; i++ ) Node[i] = NULL ;
	for ( int i=0 ; i<3 ; i++ ) Coord[i] = 0. ;
	Density	=	0. ;
	NumDensity =	0. ;
	for ( int i=0 ; i<3 ;i ++ ) Velocity[i] = 0. ;
	Temp	=	0. ;
}


int FileNum( string Filename ) ;
void ReadGrid( NODE *Node , CELL *Cell , int NodeNum , int CellNum ) ;
void ReadResult( CELL *Cell , int CellNum ) ;
void OutputResultTec( NODE *Node , CELL *Cell , int NodeNum , int CellNum ) ;
int SortCell( CELL *Cell , int CellNum ) ;


using namespace std ;


int main(){
	int		NodeNum , CellNum ;
	NODE		*Node ;
	CELL		*Cell ;


	// Get number of nodes and cells.
	NodeNum = FileNum("node.dat") - 1 ;
	CellNum = FileNum("cell.dat") - 1 ;

	// Allocate the memory for cell and node.
	Node = new NODE[NodeNum] ;
	Cell = new CELL[CellNum] ;

	// Read mesh information, including the coordinate of nodes and cells (node.dat and cell.dat).
	ReadGrid( Node , Cell , NodeNum , CellNum ) ;

	// Read simulation results, including density, number density, velocities and temperature (Result.dat).
	ReadResult( Cell , CellNum ) ;
	
	// Remove solid cells.
	CellNum = SortCell( Cell , CellNum ) ;
	
	// Output the results for Tecplot format (Result-tec.dat).
	OutputResultTec( Node , Cell , NodeNum , CellNum ) ;


	return	0 ;
}

//------------------------------------------------------------------------------------------------------------

void ReadGrid( NODE *Node , CELL *Cell , int NodeNum , int CellNum ){
	int		NodeNo ;
	string		buffer ;
	ifstream	Input ;

	Input.open( "node.dat" , ios::in ) ;

	if ( !Input ){
 		cout << "Fail to open file incoord.txt \n" ;
		exit(1) ;
	}


	getline( Input , buffer ) ;
	for ( int i=0 ; i<NodeNum ; i++ )
		Input >> Node[i].Id >> Node[i].Coord[0] >> Node[i].Coord[1] >> Node[i].Coord[2] ;

	Input.clear() ;
	Input.close() ;
	

	Input.open( "cell.dat" , ios::in ) ;

	if ( !Input ){
 		cout << "Fail to open file node.dat \n" ;
		exit(1) ;
	}

	getline( Input , buffer ) ;
	for ( int i=0 ; i<CellNum ; i++ ){
		Input >> Cell[i].Id >> Cell[i].Type ;
		for ( int j=0 ; j<4 ; j++ ){
			Input >> NodeNo ;
			Cell[i].Node[j]	= &Node[NodeNo-1] ;
		}
	}

	Input.clear() ;
	Input.close() ;
}

//------------------------------------------------------------------------------------------------------------

void ReadResult( CELL *Cell , int CellNum ){
	ifstream	Input ;
	string		buffer ;


	Input.open( "Result.dat" , ios::in ) ;


	if ( !Input ){
 		cout << "Fail to open Result.dat\n" ;
		exit(1) ;
	}
	
	getline( Input , buffer ) ;
	for ( int i=0 ; i<CellNum ; i++ ){
		Input >> Cell[i].Id >> Cell[i].Coord[0] >> Cell[i].Coord[1] >> Cell[i].Coord[2] >> Cell[i].Density >> Cell[i].NumDensity 
		      >> Cell[i].Velocity[0] >> Cell[i].Velocity[1] >> Cell[i].Velocity[2] >> Cell[i].Temp >> Cell[i].AveParticleNum ;
	}

	Input.clear() ;
	Input.close() ;
}

//------------------------------------------------------------------------------------------------------------

void OutputResultTec( NODE *Node , CELL *Cell , int NodeNum , int CellNum ){
	ofstream	Output ;
	int		OutNum=0 ;

	Output.open( "Result-tec.dat" , ios::out | ios::trunc ) ;

	if ( !Output.fail() ){
		Output << "TITLE = \"Finite volume dataset\"\n" ;
		Output << "VARIABLES = \"X\", \"Y\", \"Density\", \"NumDensity\", \"U-Vel\", \"V-Vel\", \"W-Vel\", \"Temp\", \"#P\"\n" ;
		Output << "ZONE N=" << NodeNum << ", E=" << CellNum << ", DATAPACKING=BLOCK, ZONETYPE=FEQUADRILATERAL\n" ;
		Output << "VARLOCATION=([1-2]=NODAL, [3-9]=CELLCENTERED)\n" ;

		for ( int i=0 ; i<2 ; i++ ){
			for ( int j=0 ; j<NodeNum ; j++ ){
				OutNum++ ;
				Output << setw(18) << Node[j].Coord[i] ;
				if ( OutNum%6 == 0 ) Output << '\n' ;
			}
		}

		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << Cell[i].Density ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}


		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << Cell[i].NumDensity ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}

		for ( int i=0 ; i<3 ; i++ ){
			for ( int j=0 ; j<CellNum ; j++ ){
				OutNum++ ;
				Output << setw(18) << Cell[j].Velocity[i] ;
				if ( OutNum%6 == 0 ) Output << '\n' ;
			}
		}

		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << Cell[i].Temp ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}


		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << Cell[i].AveParticleNum ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}


		Output << "\n\n" ;
		for ( int i=0 ; i<CellNum ; i++ ){
			for ( int j=0 ; j<4 ; j++ )
				Output << setw(18) << Cell[i].Node[j]->Id ;
			Output << '\n' ;
		}
	}else{
		cout << "Fail to open Result-tec.dat \n" ;
		exit(1) ;
	}

	Output.clear() ;
	Output.close() ;
}

//------------------------------------------------------------------------------------------------------------

int SortCell( CELL *Cell , int CellNum ){

	for ( int i=0 ; i<CellNum ; i++ ){
		if ( Cell[i].Type == 2 ){
			Cell[i]	=	Cell[CellNum-1] ;
			i-- ;
			CellNum-- ;
		}
	}
	return	CellNum ;
}

//------------------------------------------------------------------------------------------------------------

int FileNum( string Filename ){
	int		Num=0  ;
	string		buffer ;
	ifstream	Input ;

	Input.open( Filename.c_str() , ios::in ) ;
	
	if ( !Input ){
 		cout << "Fail to open file " << Filename << '\n';
		exit(1) ;
	}
	
	while ( getline( Input , buffer ) ){
		Num+=1 ;
	}

	Input.clear() ;
	Input.close() ;

	return	Num ;
}