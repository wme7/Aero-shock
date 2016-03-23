//===================================================================================================
//
// These functions are to initialize DSMC simulation, including setting grid data, and initial 
// particles.
//
//===================================================================================================

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include "dsmcsg_init.h"
#include "parameter.h"


using namespace std ;

// Setup boundary conditions and internal blocks inside simulation domain.
void BoundaryCondition( DSMC_DOMAIN *h_Domain , DSMC_BLOCK *h_Block ){
	// h_Block[].Type:  -3: inlet 
	//                  -4: vacuum
	//                 -21: fully specular wall
	//                 -22: fully diffusive wall

	
	// M-4 Supersonic Flow with attack angle of 30 degree over a Block.	
	// BC 0 (inlet boundary: -3)
	h_Block[0].Type		=	-3 ;
	h_Block[0].XVel		=	1117.14 ;
	h_Block[0].YVel		=	644.98 ;
	h_Block[0].NumDensity	=	3.24e20 ;
	h_Block[0].Temp		=	300. ;

	h_Block[0].XL		=	0. ;
	h_Block[0].XH		=	0.8 ;
	h_Block[0].YL		=	0. ;
	h_Block[0].YH		=	0. ;


	// BC 1 (inlet boundary: -3)
	h_Block[1].Type		=	-3 ;
	h_Block[1].XVel		=	1117.14 ;
	h_Block[1].YVel		=	644.98 ;
	h_Block[1].NumDensity	=	3.24e20 ;
	h_Block[1].Temp		=	300. ;

	h_Block[1].XL		=	0.8 ;
	h_Block[1].XH		=	0.8 ;
	h_Block[1].YL		=	0. ;
	h_Block[1].YH		=	1. ;
	

	// BC 2 (inlet boundary: -3)
	h_Block[2].Type		=	-3 ;
	h_Block[2].XVel		=	1117.14 ;
	h_Block[2].YVel		=	644.98 ;
	h_Block[2].NumDensity	=	3.24e20 ;
	h_Block[2].Temp		=	300. ;

	h_Block[2].XL		=	0. ;
	h_Block[2].XH		=	0.8 ;
	h_Block[2].YL		=	1. ;
	h_Block[2].YH		=	1. ;


	// BC 3 (inlet boundary: -3)
	h_Block[3].Type		=	-3 ;
	h_Block[3].XVel		=	1117.14 ;
	h_Block[3].YVel		=	644.98 ;
	h_Block[3].NumDensity	=	3.24e20 ;
	h_Block[3].Temp		=	300. ;

	h_Block[3].XL		=	0. ;
	h_Block[3].XH		=	0. ;
	h_Block[3].YL		=	0. ;
	h_Block[3].YH		=	1. ;


	// Internal Block (fully diffusive wall: -22)
	h_Block[4].Type		=	-22 ;
	h_Block[4].XVel		=	0. ;
	h_Block[4].YVel		=	0. ;
	h_Block[4].NumDensity	=	0. ;
	h_Block[4].Temp		=	1200. ;

	h_Block[4].XL		=	0.3 ;
	h_Block[4].XH		=	0.5 ;
	h_Block[4].YL		=	0.4 ;
	h_Block[4].YH		=	0.6 ;

	h_Block[4].NodeNum      =       4 ;
	h_Block[4].Node[0].XCoord=      0.3 ;
	h_Block[4].Node[0].YCoord=      0.4 ;
	h_Block[4].Node[1].XCoord=      0.5 ;
	h_Block[4].Node[1].YCoord=      0.4 ;
	h_Block[4].Node[2].XCoord=      0.5 ;
	h_Block[4].Node[2].YCoord=      0.6 ;
	h_Block[4].Node[3].XCoord=      0.3 ;
	h_Block[4].Node[3].YCoord=      0.6 ;
	
	
	
	/*
	// M-5 Supersonic Mach reflection problem.
	// BC 0 (inlet boundary: -3)
	h_Block[0].Type		=	-3 ;
	h_Block[0].XVel		=	1612.45 ;
	h_Block[0].YVel		=	0. ;
	h_Block[0].NumDensity	=	1.294e20 ;
	h_Block[0].Temp		=	300. ;

	h_Block[0].XL		=	0. ;
	h_Block[0].XH		=	0. ;
	h_Block[0].YL		=	0. ;
	h_Block[0].YH		=	2. ;


	// BC 1 (fully specular wall: -21)
	h_Block[1].Type		=	-21 ;
	h_Block[1].XVel		=	0. ;
	h_Block[1].YVel		=	0. ;
	h_Block[1].NumDensity	=	0. ;
	h_Block[1].Temp		=	0. ;

	h_Block[1].XL		=	0. ;
	h_Block[1].XH		=	3.2 ;
	h_Block[1].YL		=	0. ;
	h_Block[1].YH		=	0. ;
	

	// BC 2 (vacuum boundary: -4)
	h_Block[2].Type		=	-4 ;
	h_Block[2].XVel		=	0. ;
	h_Block[2].YVel		=	0. ;
	h_Block[2].NumDensity	=	0. ;
	h_Block[2].Temp		=	0. ;

	h_Block[2].XL		=	3.2 ;
	h_Block[2].XH		=	3.2 ;
	h_Block[2].YL		=	0. ;
	h_Block[2].YH		=	2. ;


	// BC 3 (fully diffusive wall: -22)
	h_Block[3].Type		=	-22 ;
	h_Block[3].XVel		=	0. ;
	h_Block[3].YVel		=	0. ;
	h_Block[3].NumDensity	=	0. ;
	h_Block[3].Temp		=	600. ;

	h_Block[3].XL		=	0.1 ;
	h_Block[3].XH		=	3.2 ;
	h_Block[3].YL		=	2. ;
	h_Block[3].YH		=	2. ;
	
	
	// BC 4 (inlet boundary: -3)
	h_Block[4].Type		=	-3 ;
	h_Block[4].XVel		=	1612.45 ;
	h_Block[4].YVel		=	0. ;
	h_Block[4].NumDensity	=	1.294e20 ;
	h_Block[4].Temp		=	300. ;

	h_Block[4].XL		=	0. ;
	h_Block[4].XH		=	0.1 ;
	h_Block[4].YL		=	2. ;
	h_Block[4].YH		=	2. ;


	// Internal Block (fully diffusive wall: -22)
	h_Block[5].Type		=	-22 ;
	h_Block[5].XVel		=	0. ;
	h_Block[5].YVel		=	0. ;
	h_Block[5].NumDensity	=	0. ;
	h_Block[5].Temp		=	600. ;

	h_Block[5].XL		=	0.1 ;
	h_Block[5].XH		=	3.2 ;
	h_Block[5].YL		=	1.5 ;
	h_Block[5].YH		=	2. ;

	h_Block[5].NodeNum      =       4 ;
	h_Block[5].Node[0].XCoord=      0.1 ;
	h_Block[5].Node[0].YCoord=      2. ;
	h_Block[5].Node[1].XCoord=      0.9660254 ;
	h_Block[5].Node[1].YCoord=      1.5 ;
	h_Block[5].Node[2].XCoord=      3.2 ;
	h_Block[5].Node[2].YCoord=      1.5 ;
	h_Block[5].Node[3].XCoord=      3.2 ;
	h_Block[5].Node[3].YCoord=      2. ;
	*/
	
	

	/*
	// M-22 Hypersonic flow over a ramp.
	// BC 0 (inlet boundary: -3)
	h_Block[0].Type		=	-3 ;
	h_Block[0].XVel		=	1831.87 ;
	h_Block[0].YVel		=	0. ;
	h_Block[0].NumDensity	=	1.294e20 ;
	h_Block[0].Temp		=	20. ;

	h_Block[0].XL		=	0. ;
	h_Block[0].XH		=	0. ;
	h_Block[0].YL		=	0. ;
	h_Block[0].YH		=	1. ;


	// BC 1 (inlet boundary: -3)
	h_Block[1].Type		=	-3 ;
	h_Block[1].XVel		=	1831.87 ;
	h_Block[1].YVel		=	0. ;
	h_Block[1].NumDensity	=	1.294e20 ;
	h_Block[1].Temp		=	20. ;

	h_Block[1].XL		=	0. ;
	h_Block[1].XH		=	0.2 ;
	h_Block[1].YL		=	0. ;
	h_Block[1].YH		=	0. ;
	

	// BC 2 (fully diffusive wall: -22)
	h_Block[2].Type		=	-22 ;
	h_Block[2].XVel		=	0. ;
	h_Block[2].YVel		=	0. ;
	h_Block[2].NumDensity	=	0. ;
	h_Block[2].Temp		=	500. ;

	h_Block[2].XL		=	0.2 ;
	h_Block[2].XH		=	2.066 ;
	h_Block[2].YL		=	0. ;
	h_Block[2].YH		=	0. ;


	// BC 3 (vacuum boundary: -4)
	h_Block[3].Type		=	-4 ;
	h_Block[3].XVel		=	0. ;
	h_Block[3].YVel		=	0. ;
	h_Block[3].NumDensity	=	0. ;
	h_Block[3].Temp		=	0. ;

	h_Block[3].XL		=	2.066 ;
	h_Block[3].XH		=	2.5 ;
	h_Block[3].YL		=	0. ;
	h_Block[3].YH		=	0. ;
	
	
	// BC 4 (vacuum boundary: -4)
	h_Block[4].Type		=	-4 ;
	h_Block[4].XVel		=	0. ;
	h_Block[4].YVel		=	0. ;
	h_Block[4].NumDensity	=	0. ;
	h_Block[4].Temp		=	0. ;

	h_Block[4].XL		=	2.5 ;
	h_Block[4].XH		=	2.5 ;
	h_Block[4].YL		=	0. ;
	h_Block[4].YH		=	1. ;


	// BC 5 (inlet boundary: -3)
	h_Block[5].Type		=	-3 ;
	h_Block[5].XVel		=	1831.87 ;
	h_Block[5].YVel		=	0. ;
	h_Block[5].NumDensity	=	1.294e20 ;
	h_Block[5].Temp		=	20. ;

	h_Block[5].XL		=	0 ;
	h_Block[5].XH		=	2.5 ;
	h_Block[5].YL		=	1. ;
	h_Block[5].YH		=	1. ;


	// Internal Block (fully diffusive wall: -22)
	h_Block[6].Type		=	-22 ;
	h_Block[6].XVel		=	0. ;
	h_Block[6].YVel		=	0. ;
	h_Block[6].NumDensity	=	1.294e20 ;
	h_Block[6].Temp		=	400. ;

	h_Block[6].XL		=	1.2 ;
	h_Block[6].XH		=	2.066 ;
	h_Block[6].YL		=	0. ;
	h_Block[6].YH		=	0.5 ;

	h_Block[6].NodeNum      =       3 ;
	h_Block[6].Node[0].XCoord=      1.2 ;
	h_Block[6].Node[0].YCoord=      0. ;
	h_Block[6].Node[1].XCoord=      2.066 ;
	h_Block[6].Node[1].YCoord=      0. ;
	h_Block[6].Node[2].XCoord=      2.066 ;
	h_Block[6].Node[2].YCoord=      0.5 ;
	*/


	for ( int i=0 ; i<h_Domain->BlockNum ; i++ ){
		h_Block[i].XL	-=	1.e-6 ;
		h_Block[i].XH	+=	1.e-6 ;
		h_Block[i].YL	-=	1.e-6 ;
		h_Block[i].YH	+=	1.e-6 ;
	}
}

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

// Setup solid cell inside domain.
void SetInternalBlock( DSMC_DOMAIN *h_Domain , DSMC_CELL *h_Cell , DSMC_BLOCK *h_Block ){
	int		Neighbor1 , Neighbor2 ;


	for ( int i=0 ; i<h_Domain->CellNum ; i++ ){
		for ( int j=h_Domain->BCBlockNum ; j<h_Domain->BlockNum ; j++ ){
			// Setup solid cell (Type: 2).
			if ( InBlock( &h_Cell[i] , &h_Block[j] ) ){
				h_Cell[i].Type	=	2 ;

				// Setup tarcking cell (Type: 1).
				for ( int m=0 ; m<4 ; m++ ){
					Neighbor1	=	h_Cell[i].Neighbor[m] ;
					if ( Neighbor1 < 0 ) continue ;

					if ( h_Cell[Neighbor1].Type != 2 && h_Cell[Neighbor1].Type >= 0 ) h_Cell[Neighbor1].Type = 1 ;
					

					for ( int n=0 ; n<4 ; n++ ){
						Neighbor2	=	h_Cell[Neighbor1].Neighbor[n] ;
						if ( Neighbor2 < 0 ) continue ;

						if ( h_Cell[Neighbor2].Type != 2 && h_Cell[Neighbor2].Type >= 0 ){
							h_Cell[Neighbor2].Type = 1 ;
						}else{
							h_Cell[Neighbor1].Neighbor[n] = -(j+1) ;
						}
					}
				}

				break ;
			}
		}
	}
}

//------------------------------------------------------------------------------------------------------------

bool InBlock( DSMC_CELL *_Cell , DSMC_BLOCK *_Block ){
	bool		inblock = true ;
	double		XCenter , YCenter ;
	double		vector_a[2] , vector_b[2] , cross ;

	XCenter	=	_Cell->XCenter ;
	YCenter	=	_Cell->YCenter ;

	if ( XCenter > _Block->XL && XCenter < _Block->XH && YCenter > _Block->YL && YCenter < _Block->YH ){
		if ( _Block->NodeNum > 0 ){
			for ( int i=0 ; i<_Block->NodeNum ; i++ ){
				vector_a[0]	=	_Block->Node[i].XCoord - XCenter ;
				vector_a[1]	=	_Block->Node[i].YCoord - YCenter ;
				if ( i==(_Block->NodeNum-1) ){
					vector_b[0]	=	_Block->Node[0].XCoord - XCenter ;
					vector_b[1]	=	_Block->Node[0].YCoord - YCenter ;
				}else{
					vector_b[0]	=	_Block->Node[i+1].XCoord - XCenter ;
					vector_b[1]	=	_Block->Node[i+1].YCoord - YCenter ;
				}

				cross	=	vector_a[0]*vector_b[1] - vector_a[1]*vector_b[0] ;
				if ( cross < 0. ){
					inblock = false ;
					break ;
				}
			}
		}
	}else{
		inblock	=	false ;
	}

	return	inblock ;
}

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

// Create mesh system.
void InitGrid(	int *h_InletCellNo ,
		int *h_InletEdgeNo ,
		float *h_InletXVel ,
		float *h_InletYVel ,
		float *h_InletZVel ,
		float *h_InletNumDen ,
		float *h_InletTemp ,
		float *h_InletArea ,
		DSMC_DOMAIN *h_Domain ,
		DSMC_NODE *h_Node ,
		DSMC_CELL *h_Cell , 
		DSMC_BLOCK *h_Block ){


	CreateGrid( h_Domain , h_Node , h_Cell ) ;
	

	SetNeighborCell( h_InletCellNo , h_InletEdgeNo , h_InletXVel , h_InletYVel , h_InletZVel , h_InletNumDen , h_InletTemp , h_InletArea , h_Domain , h_Node ,
			 h_Cell , h_Block ) ;
}

//------------------------------------------------------------------------------------------------------------

// Setup coordinate of nodes and cell data.
void CreateGrid( DSMC_DOMAIN *h_Domain , DSMC_NODE *h_Node , DSMC_CELL *h_Cell ){
	int		node_no , cell_no , node1 , node2 ;
	float		coordinate[2] , buffer[2] ;


	// Set coordinate of node in x/y-direction.
	for ( int j=0 ; j<h_Domain->YNodeNum ; j++ ){
		coordinate[1]	=	h_Domain->YL + h_Domain->YCellSize*j ;

		if ( j==h_Domain->YCellNum )
			coordinate[1] = h_Domain->YH ;

		for ( int i=0 ; i<h_Domain->XNodeNum ; i++  ){
			node_no	=	j*h_Domain->XNodeNum + i ;
			coordinate[0]   =       h_Domain->XL + h_Domain->XCellSize*i ;

			if ( i==h_Domain->XCellNum )
				coordinate[0] = h_Domain->XH ;

			h_Node[node_no].XCoord = coordinate[0] ;
			h_Node[node_no].YCoord = coordinate[1] ;
		}
	}


	//Set node no. of cell
	for ( int j=0 ; j<h_Domain->YCellNum ; j++ ){
		for ( int i=0 ; i<h_Domain->XCellNum ; i++ ){
			cell_no = j*h_Domain->XCellNum + i ;
			node_no = j*h_Domain->XNodeNum + i ;


			h_Cell[cell_no].Node[0]	=	node_no ;
			h_Cell[cell_no].Node[1]	=	node_no+1 ;
			h_Cell[cell_no].Node[2]	=	node_no+1+h_Domain->XNodeNum ;
			h_Cell[cell_no].Node[3]	=	node_no+h_Domain->XNodeNum ;

			if ( i<2 || i>=(h_Domain->XCellNum-2) || j<2 || j>=(h_Domain->YCellNum-2) )
				h_Cell[cell_no].Type	=	1 ;
			else
				h_Cell[cell_no].Type	=	0 ;


			for ( int k=0 ; k<4 ; k++ ){
				h_Cell[cell_no].XCenter	+=	h_Node[h_Cell[cell_no].Node[k]].XCoord ;
				h_Cell[cell_no].YCenter	+=	h_Node[h_Cell[cell_no].Node[k]].YCoord ;
			}

			h_Cell[cell_no].XCenter /= 4. ;
			h_Cell[cell_no].YCenter /= 4. ;


			// Set Edge Info.
			for ( int k=0 ; k<4 ; k++ ){
				if ( k<=2 ){
					node1	=	h_Cell[cell_no].Node[k] ;
					node2	=	h_Cell[cell_no].Node[k+1] ;
				}else{
					node1	=	h_Cell[cell_no].Node[k] ;
					node2	=	h_Cell[cell_no].Node[0] ;
				}

				buffer[0]	=	h_Node[node2].YCoord - h_Node[node1].YCoord ;
				buffer[1]	=	h_Node[node2].XCoord - h_Node[node1].XCoord ;

				h_Cell[cell_no].EdgeFA[k] = buffer[0] ;
				h_Cell[cell_no].EdgeFB[k] = -buffer[1] ;
				h_Cell[cell_no].EdgeFC[k] = buffer[1] * h_Node[node1].YCoord - buffer[0] * h_Node[node1].XCoord ;
			}
		}
	}



	// Output coordinate of nodes to node.dat file.
	PrintNodeInfo( h_Node , h_Domain->NodeNum ) ;
	//PrintCellInfo( h_Cell , h_Domain->CellNum ) ;
	//PrintEdgeInfo( h_EdgeFA , h_EdgeFB , h_EdgeFC , h_Domain->CellNum ) ;
}

//------------------------------------------------------------------------------------------------------------

// Setup neighbor cells of each cell.
void SetNeighborCell(	int *h_InletCellNo ,
			int *h_InletEdgeNo ,
			float *h_InletXVel ,
			float *h_InletYVel ,
			float *h_InletZVel ,
			float *h_InletNumDen ,
			float *h_InletTemp ,
			float *h_InletArea ,
			DSMC_DOMAIN *h_Domain ,
			DSMC_NODE *h_Node ,
			DSMC_CELL *h_Cell ,
			DSMC_BLOCK *h_Block ){

	int	*node_cell_num , *node_cell_relation ;
	int	max_node_cell_num , node_no , m ;
	int	cell_1 , cell_2 , edge_num , a_node[2] , b_node[2] ;
	int	a1 , node1 , node2 ;
	int	FaceNum ;	


        max_node_cell_num       =       4 ;
	FaceNum			=	0 ;

        node_cell_num           =       new int[h_Domain->NodeNum] ;
        node_cell_relation      =       new int[max_node_cell_num*h_Domain->NodeNum] ;

	

        for ( int i=0 ; i<h_Domain->NodeNum ; i++ )
                node_cell_num[i]        =       0 ;

        for ( int i=0 ; i<max_node_cell_num*h_Domain->NodeNum ; i++ )
                node_cell_relation[i]   =       -1 ;

        cell_1                  =       0 ;
        cell_2                  =       0 ;
        edge_num                =       4 ;
        m                       =       0 ;

        // Set relation between node and cell.  
	for ( int i=0 ; i<h_Domain->CellNum ; i++ ){
                for ( int j=0 ; j<4 ; j++ ){
                        node_no =       h_Cell[i].Node[j] ;
                        
			node_cell_num[node_no]  +=      1 ;
                        node_cell_relation[node_no*max_node_cell_num + (node_cell_num[node_no]-1)]      =       i ;
                }
        }

	
	// Set neighbor cells to each cell.
        for ( int i=0 ; i<h_Domain->NodeNum ; i++ ){
                for ( int j=0 ; j<node_cell_num[i] ; j++ ){
                        cell_1  =       node_cell_relation[j+max_node_cell_num*i] ;

                        for ( int k=0 ; k<edge_num ; k++ ){
				if ( h_Cell[cell_1].Neighbor[k] == -999 ){
                                        m       =       0 ;
					while ( m<node_cell_num[i] ){
                                                if ( j != m ){
                                                        cell_2  =       node_cell_relation[m+max_node_cell_num*i] ;

                                                        for ( int n=0 ; n<edge_num ; n++ ){
								a_node[0] = h_Cell[cell_1].Node[k] ;
								if ( k==3 )
									a_node[1] = h_Cell[cell_1].Node[0]  ;
								else
									a_node[1] = h_Cell[cell_1].Node[k+1]  ;
								

								b_node[0] = h_Cell[cell_2].Node[n] ;
								if ( n==3 )
									b_node[1] = h_Cell[cell_2].Node[0] ;
								else
									b_node[1] = h_Cell[cell_2].Node[n+1] ;
			

								if ( a_node[0] == b_node[1] && a_node[1] == b_node[0] ){
                                                                        h_Cell[cell_1].Neighbor[k]      =       cell_2 ;
                                                                        h_Cell[cell_2].Neighbor[n]      =       cell_1 ;

                                                                        m       =       node_cell_num[i] ;
                                                                        break ;
                                                                }
                                                        }

						}

                                                m       +=      1 ;
                                        }

                                }

                        }
                }
        }
	

	// Set boundary.
	for ( int i=0 ; i<h_Domain->CellNum ; i++ ){
                for ( int j=0 ; j<4 ; j++ ){
                        if (  h_Cell[i].Neighbor[j] == -999 ){
				node1	=	h_Cell[i].Node[j] ;
				if ( j==3 )
					node2	=	h_Cell[i].Node[0] ;
				else
					node2	=	h_Cell[i].Node[j+1] ;


				SetBC( &h_Cell[i].Neighbor[j] , h_Node[node1].XCoord , h_Node[node1].YCoord , h_Node[node2].XCoord , h_Node[node2].YCoord , 
					h_InletCellNo , h_InletEdgeNo , h_InletXVel , h_InletYVel , h_InletZVel , h_InletNumDen , h_InletTemp , h_InletArea , 
					&FaceNum , i , j , h_Domain , h_Block ) ;
                        }
                }
        }

	h_Domain->InletFaceNum	=	FaceNum ;

	delete [] node_cell_num ;
	delete [] node_cell_relation ;
}

//------------------------------------------------------------------------------------------------------------

// Setup boundary and inlet faces.
void SetBC(	int *_CellNeighbor ,
		float Node1XCoord ,
		float Node1YCoord , 
		float Node2XCoord , 
		float Node2YCoord ,
		int *h_InletCellNo ,
		int *h_InletEdgeNo ,
		float *h_InletXVel ,
		float *h_InletYVel ,
		float *h_InletZVel ,
		float *h_InletNumDen ,
		float *h_InletTemp ,
		float *h_InletArea ,
		int *FaceNum ,
		int CellNo ,
		int EdgeNo ,
		DSMC_DOMAIN *h_Domain ,
		DSMC_BLOCK *h_Block ){


	for ( int i=0 ; i<h_Domain->BlockNum ; i++ ){
		if ( Node1XCoord > h_Block[i].XL && Node1XCoord < h_Block[i].XH && Node1YCoord > h_Block[i].YL && Node1YCoord < h_Block[i].YH && 
		     Node2XCoord > h_Block[i].XL && Node2XCoord < h_Block[i].XH && Node2YCoord > h_Block[i].YL && Node2YCoord < h_Block[i].YH ){
			(*_CellNeighbor)	=	-(i+1) ;


			// Set Inlet Face.
			if ( h_Block[i].Type == -3 ){
				if ( (*FaceNum) >= h_Domain->InletFaceNum ){
					cout << "Inlet Face Number is ERROR. " << h_Domain->InletFaceNum << '\n' ;
					exit(0) ;
				}

				h_InletCellNo[*FaceNum]	=	CellNo ;
				h_InletEdgeNo[*FaceNum]	=	EdgeNo ;
				h_InletXVel[*FaceNum]	=	h_Block[i].XVel ;
				h_InletYVel[*FaceNum]	=	h_Block[i].YVel ;
				h_InletZVel[*FaceNum]	=	0. ;
				h_InletNumDen[*FaceNum]	=	h_Block[i].NumDensity ;
				h_InletTemp[*FaceNum]	=	h_Block[i].Temp ;
				h_InletArea[*FaceNum]	=	sqrt((Node1XCoord-Node2XCoord)*(Node1XCoord-Node2XCoord) + 
								     (Node1YCoord-Node2YCoord)*(Node1YCoord-Node2YCoord) ) ;

				(*FaceNum)++ ;
			}

			break ;
		}
	}
}

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

// Calculate number of enter particles in each timestep.
void SetupInletFace( int *h_InletCellNo , int *h_InletEdgeNo , float *h_InletNum , float *h_InletNumR , float *h_InletXVel ,float *h_InletYVel , float *h_InletZVel ,
		     float *h_InletNumDen , float *h_InletTemp , float *h_InletArea , DSMC_DOMAIN *h_Domain ){

	float		MostProbableSpeed , Mass , SpeedRatio , Flux ;

	Mass	=	MASS ;

	for ( int i=0 ; i<h_Domain->InletFaceNum ; i++ ){
		MostProbableSpeed	=	sqrt(2.*BOLTZ*h_InletTemp[i]/Mass) ;

		if ( h_InletEdgeNo[i] == 0 ) SpeedRatio = h_InletYVel[i]/MostProbableSpeed ;
		if ( h_InletEdgeNo[i] == 1 ) SpeedRatio = -h_InletXVel[i]/MostProbableSpeed ;
		if ( h_InletEdgeNo[i] == 2 ) SpeedRatio = -h_InletYVel[i]/MostProbableSpeed ;
		if ( h_InletEdgeNo[i] == 3 ) SpeedRatio = h_InletXVel[i]/MostProbableSpeed ;


		if ( fabs(SpeedRatio) < 10.1 )	Flux = (exp(-SpeedRatio*SpeedRatio) + sqrt(PI)*SpeedRatio*(1.+ERF(SpeedRatio))) /(2.*sqrt(PI)) ;
		if ( SpeedRatio > 10. )		Flux = SpeedRatio ;
		if ( SpeedRatio < -10. )	Flux = 0. ;

		h_InletNum[i]	=	h_InletNumDen[i]*Flux*MostProbableSpeed*h_Domain->Timestep*h_InletArea[i]/h_Domain->ParticleWeight ;
	}
}

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

// Setup initial particles.
void InitFlow(	float *h_ParticleXCoord ,
		float *h_ParticleYCoord ,
		float *h_ParticleXVel ,
		float *h_ParticleYVel ,
		float *h_ParticleZVel , 
		int *h_ParticleCellNo ,
		int *h_ParticleLastCollide ,
		int *h_ParticleIn ,
		DSMC_DOMAIN *h_Domain ,
		DSMC_CELL *h_Cell ){

	int		k = 0 , ParticleNum = 0 ;
	float		Temp , CrossSection , Mass , Vel[3] , MostProbableSpeed , Buffer ;
	

	Temp			=	h_Domain->Temp ;
	CrossSection		=	MIX_CROSS_SECTION ;
	Mass			=	MASS ;
	MostProbableSpeed	=	sqrt( 2*BOLTZ*Temp/Mass ) ;
	Vel[0]			=	h_Domain->XVel ;
	Vel[1]			=	h_Domain->YVel;
	Vel[2]			=	0. ;
	Buffer			=	0. ;

	// Calculate timestep and weighting.
	h_Domain->Timestep	=	h_Domain->XCellSize/( 3*(sqrt(pow(Vel[0],2)+pow(Vel[1],2)+pow(Vel[2],2)) + sqrt(2*BOLTZ*Temp/Mass)) );
	h_Domain->ParticleWeight=	h_Domain->NumDensity*h_Domain->CellVolume / h_Domain->WeightRatio ;


	for ( int i=0 ; i<h_Domain->CellNum ; i++ ){
		h_Cell[i].MaxCrossSectionSpeed	=	CrossSection * (300.*sqrt(Temp/300.));
		h_Cell[i].RemainderCollisionPair=	randn() ;


		ParticleNum	=	(h_Domain->NumDensity*h_Domain->CellVolume) / h_Domain->ParticleWeight ;
		
		for ( int j=0 ; j<ParticleNum ; j++ ){
			h_ParticleXCoord[k]	=	h_Cell[i].XCenter + (randn()-0.5) * h_Domain->XCellSize ;
			h_ParticleYCoord[k]	=	h_Cell[i].YCenter + (randn()-0.5) * h_Domain->YCellSize ;

			RandVel( &h_ParticleXVel[k] , &Buffer , MostProbableSpeed ) ;
			RandVel( &h_ParticleYVel[k] , &Buffer , MostProbableSpeed ) ;
			RandVel( &h_ParticleZVel[k] , &Buffer , MostProbableSpeed ) ;

			h_ParticleXVel[k]	+=	Vel[0] ;
			h_ParticleYVel[k]	+=	Vel[1] ;
			h_ParticleZVel[k]	+=	Vel[2] ;


			h_ParticleCellNo[k]	=	i ;
			h_ParticleLastCollide[k]=	-1 ;

			h_ParticleIn[k]		=	1 ;

			k++ ;
			if ( k>= h_Domain->MaxParticleNum ){
				cout << "Number of particles is more than maximum particle number." << setw(15) << k << setw(15) << h_Domain->MaxParticleNum << '\n' ;
				exit(1) ;
			}
		}
	}

	h_Domain->ParticleNum	=	k ;

	// Debug
	//PrintParticleInfo( h_ParticleXCoord , h_ParticleYCoord , h_ParticleXVel , h_ParticleYVel , h_ParticleZVel , h_ParticleCellNo , 
	//		   h_ParticleLastCollide , h_Domain->ParticleNum ) ;

}

//------------------------------------------------------------------------------------------------------------

void PrintNodeInfo( DSMC_NODE *Node , int NodeNum ){
	ofstream	outfile ;

	outfile.open( "node.dat" , ios::out | ios::trunc ) ;

	if ( !outfile.fail() ){
		outfile << setw(18) << "Node No." << setw(18) << "X Coord." << setw(18) << "Y Coord." << setw(18) << "Z Coord." << '\n' ;

		for ( int i=0 ; i<NodeNum ; i++ )
			outfile << setw(18) << (i+1) << setw(18) << Node[i].XCoord << setw(18) << Node[i].YCoord << setw(18) << "0" << '\n' ;
	}else{
		cout << "Fail to open \"node.dat\" file.\n" ;
	}
}

//------------------------------------------------------------------------------------------------------------

void PrintCellInfo( DSMC_CELL *Cell , int CellNum ){
	int		k = 0 ;
	ofstream	outfile ;

	outfile.open( "cell.dat" , ios::out | ios::trunc ) ;


	if ( !outfile.fail() ){
		outfile << setw(15) << "Cell No." << setw(8) << "Type" << setw(15) << "Node 1" << setw(15) << "Node 2" << setw(15) << "Node 3" 
			<< setw(15) << "Node 4" << '\n' ;

		for ( int i=0 ; i<CellNum ; i++ ){
			outfile << setw(15) << (i+1) << setw(8) << Cell[i].Type << setw(15) << (Cell[i].Node[k]+1) << setw(15) << (Cell[i].Node[k+1]+1) 
				<< setw(15) << (Cell[i].Node[k+2]+1) << setw(15) << (Cell[i].Node[k+3]+1) << '\n' ;
		}
	}else{
		cout << "Fail to open \"cell.dat\" file.\n" ;
	}
}

//------------------------------------------------------------------------------------------------------------

void PrintCellInfo( int *CellNode , float *XCenter , float *YCenter , int CellNum ){
	int		k = 0 ;
	ofstream	outfile ;

	outfile.open( "cell.dat" , ios::out | ios::trunc ) ;

	if ( !outfile.fail() ){
		outfile << setw(12) << "Cell No." << setw(15) << "XCenter" << setw(15) << "YCenter" << setw(12) << "Node 1" << setw(12) << "Node 2" << setw(12) << "Node 3" << setw(12) << "Node 4" << '\n' ;

		for ( int i=0 ; i<CellNum ; i++ ){
			k = i*4 ;
			outfile << setw(12) << (i+1) << setw(15) << XCenter[i] << setw(15) << YCenter[i] << setw(12) << (CellNode[k]+1) 
			        << setw(12) << (CellNode[k+1]+1) << setw(12) << (CellNode[k+2]+1) << setw(12) << (CellNode[k+3]+1) << '\n' ;
		}
	}else{
		cout << "Fail to open \"cell.dat\" file.\n" ;
	}
}

//------------------------------------------------------------------------------------------------------------

void PrintNeighborInfo( DSMC_CELL *Cell , int CellNum ){
	int		k = 0 ;
	ofstream	outfile ;

	outfile.open( "neighbor.dat" , ios::out | ios::trunc ) ;

	if ( !outfile.fail() ){
		outfile << setw(15) << "Cell No." << setw(8) << "Type" << setw(15) << "NBR 1" << setw(15) << "NBR 2" << setw(15) << "NBR 3" << setw(15) << "NBR 4" << '\n' ;

		for ( int i=0 ; i<CellNum ; i++ ){
			outfile << setw(15) << (i+1) << setw(8) << Cell[i].Type << setw(15) << (Cell[i].Neighbor[0]+1) << setw(15) << (Cell[i].Neighbor[1]+1) << setw(15) << (Cell[i].Neighbor[2]+1) 
			        << setw(15) << (Cell[i].Neighbor[3]+1) << '\n' ;
		}
	}else{
		cout << "Fail to open \"neighbor.dat\" file.\n" ;
	}
}

//------------------------------------------------------------------------------------------------------------

void PrintParticleInfo( float *h_ParticleXCoord , float *h_ParticleYCoord , float *h_ParticleXVel , float *h_ParticleYVel , float *h_ParticleZVel , int *h_ParticleCellNo ,
			int *h_ParticleLastCollide , int ParticleNum ){

	ofstream	outfile ;

	outfile.open( "particle_info.dat" , ios::out | ios::trunc ) ;

	if ( !outfile.fail() ){
		outfile << ParticleNum << '\n' ;
		outfile << setw(15) << "Particle No." << setw(15) << "Cell No." << setw(15) << "XCoord" << setw(15) << "YCoord" << setw(15) << "XVel" 
			<< setw(15) << "YVel" << setw(15) << "ZVel" << setw(15) << "LastCollide" << '\n' ;

		for ( int i=0 ; i<ParticleNum ; i++ ){
			outfile << setw(15) << i << setw(15) << h_ParticleCellNo[i] << setw(15) << h_ParticleXCoord[i] << setw(15) << h_ParticleYCoord[i] 
				<< setw(15) << h_ParticleXVel[i] << setw(15) << h_ParticleYVel[i] << setw(15) << h_ParticleZVel[i] 
				<< setw(15) << h_ParticleLastCollide[i] << '\n' ;
		}
	}else{
		cout << "Fail to open \"particle_info.dat\" file.\n" ;
	}
}

//------------------------------------------------------------------------------------------------------------

void PrintEdgeInfo( float *h_EdgeFA , float *h_EdgeFB , float *h_EdgeFC , int CellNum ){
	ofstream	outfile ;

	outfile.open( "edge_info.dat" , ios::out | ios::trunc ) ;

	if ( !outfile.fail() ){
		outfile << setw(12) << "Cell No." << setw(12) << "Edge No." << setw(15) << "Factor A" << setw(15) << "Factor B" << setw(15) << "Factor C" << '\n' ;

		for ( int i=0 ; i<CellNum ; i++ ){
			for ( int j=0 ; j<4 ; j++ ){
				outfile << setw(12) << i << setw(12) << j << setw(15) << h_EdgeFA[i*4+j] << setw(15) << h_EdgeFB[i*4+j] 
					<< setw(15) << h_EdgeFC[i*4+j] << '\n' ;
			}
		}
	}else{
		cout << "Fail to open \"particle_info.dat\" file.\n" ;
	}
}

//------------------------------------------------------------------------------------------------------------

void PrintInletInfo( int *h_InletCellNo , int *h_InletEdgeNo , float *h_InletNum , float *h_InletNumR , float *h_InletXVel , float *h_InletYVel , float *h_InletZVel ,
		     float *h_InletNumDen , float *h_InletTemp , float *h_InletArea , int InletFaceNum ){


	cout << setw(15) << "Inlet No" << setw(15) << "CellNo" << setw(15) << "EdgeNo" << setw(15) << "Num" << setw(15) << "NumR" << setw(15) << "XVel" 
	     << setw(15) << "YVel" << setw(15) << "ZVel" << setw(15) << "NumDen" << setw(15) << "Temp" << setw(15) << "Area\n" ;
	for ( int i=0 ; i<InletFaceNum ; i++ ){
		cout << setw(15) << i << setw(15) << h_InletCellNo[i] << setw(15) << h_InletEdgeNo[i] << setw(15) << h_InletNum[i] << setw(15) << h_InletNumR[i] 
		     << setw(15) << h_InletXVel[i] << setw(15) << h_InletYVel[i] << setw(15) << h_InletZVel[i] << setw(15) << h_InletNumDen[i] 
		     << setw(15) << h_InletTemp[i] << setw(15) << h_InletArea[i] << endl ;
	}
}

//------------------------------------------------------------------------------------------------------------

// Setup initial valuve before DSMC initilization.
void InitValue( int *h_IndexParticle , int *h_IndexCell1 , int *h_IndexCell2 , float *h_ParticleXCoord , float *h_ParticleYCoord , 
		float *h_ParticleXVel , float *h_ParticleYVel , float *h_ParticleZVel , int *h_ParticleCellNo , int *h_ParticleLastCollide , 
		int *h_ParticleIn , int *h_SampleParticleNum , float *h_SampleXVel , float *h_SampleYVel , float *h_SampleZVel , 
		float *h_SampleXVelSq , float *h_SampleYVelSq , float *h_SampleZVelSq , int *h_InletCellNo , int *h_InletEdgeNo , 
		float *h_InletNum , float *h_InletNumR , float *h_InletXVel , float *h_InletYVel , float *h_InletZVel , 
		float *h_InletNumDen , float *h_InletTemp , float *h_InletArea , DSMC_DOMAIN *h_Domain ){

	for ( int i=0 ; i<h_Domain->MaxParticleNum ; i++ )
		h_IndexParticle[i]	=	0 ;
	
	for ( int i=0 ; i<h_Domain->CellNum ; i++ ){
		h_IndexCell1[i]	=	0 ;
		h_IndexCell2[i]	=	0 ;
	}

	for ( int i=0 ; i<h_Domain->MaxParticleNum ; i++ ){
		h_ParticleXCoord[i]		=	0. ;
		h_ParticleYCoord[i]		=	0. ;
		h_ParticleXVel[i]		=	0. ;
		h_ParticleYVel[i]		=	0. ;
		h_ParticleZVel[i]		=	0. ;
		h_ParticleCellNo[i]		=	0 ;
		h_ParticleLastCollide[i]	=	-1 ;
		h_ParticleIn[i]			=	0 ;
	}

	
	for ( int i=0 ; i<h_Domain->InletFaceNum ; i++ ){
		h_InletCellNo[i]	=	-1 ;
		h_InletEdgeNo[i]	=	-1 ;
		h_InletNum[i]		=	0. ;
		h_InletNumR[i]		=	0. ;
		h_InletXVel[i]		=	0. ;
		h_InletYVel[i]		=	0. ;
		h_InletZVel[i]		=	0. ;
		h_InletNumDen[i]	=	0. ;
		h_InletTemp[i]		=	0. ;
		h_InletArea[i]		=	0. ;
	}


	for ( int i=0 ; i<h_Domain->CellNum ; i++ ){
		h_SampleParticleNum[i]	=	0 ;
		h_SampleXVel[i]		=	0. ;
		h_SampleYVel[i]		=	0. ;
		h_SampleZVel[i]		=	0. ;
		h_SampleXVelSq[i]	=	0. ;
		h_SampleYVelSq[i]	=	0. ;
		h_SampleZVelSq[i]	=	0. ;
	}
}

//------------------------------------------------------------------------------------------------------------

float ERF( float S ){
	float	error_x ;
	float	B , C , D , T ;

	B	=	fabs(S) ;
	if ( B > 4. ){
		D	=	1. ;
	}else{
		C	=	exp(-B*B) ;
		T	=	1./(1.+0.3275911*B) ;
		D	=	1. - (0.254829592*T - 0.284496736*T*T + 1.421413741*T*T*T - 1.453152027*T*T*T*T + 1.061405429*T*T*T*T*T )*C ;
	}

	if ( S < 0. ) D = -D ;
	error_x = D ;
	
	return	error_x ;
}