#include "dsmcsg_class.h"

#if !defined(__DSMCSG_INIT_H)
#define __DSMCSG_INIT_H


// Setup boundary conditions and internal blocks inside simulation domain.
void BoundaryCondition( DSMC_DOMAIN *h_Domain , DSMC_BLOCK *h_Block ) ;


// Setup solid cell inside domain.
void SetInternalBlock( DSMC_DOMAIN *h_Domain , DSMC_CELL *h_Cell , DSMC_BLOCK *h_Block ) ;


bool InBlock( DSMC_CELL *_Cell , DSMC_BLOCK *_Block ) ;


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
		DSMC_BLOCK *h_Block ) ;


// Setup coordinate of nodes and cell data.
void CreateGrid( DSMC_DOMAIN *h_Domain , DSMC_NODE *h_Node , DSMC_CELL *h_Cell ) ;


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
			DSMC_BLOCK *h_Block ) ;



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
		DSMC_BLOCK *h_Block ) ;


// Calculate number of enter particles in each timestep.
void SetupInletFace( int *h_InletCellNo , int *h_InletEdgeNo , float *h_InletNum , float *h_InletNumR , float *h_InletXVel ,float *h_InletYVel , float *h_InletZVel ,
		     float *h_InletNumDen , float *h_InletTemp , float *h_InletArea , DSMC_DOMAIN *h_Domain ) ;


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
		DSMC_CELL *h_Cell ) ;


void PrintNodeInfo( DSMC_NODE *Node , int NodeNum ) ;

void PrintCellInfo( DSMC_CELL *Cell , int CellNum ) ;

void PrintCellInfo( int *CellNode , float *XCenter , float *YCenter , int CellNum ) ;

void PrintNeighborInfo( DSMC_CELL *Cell , int CellNum ) ;

void PrintParticleInfo( float *h_ParticleXCoord , float *h_ParticleYCoord , float *h_ParticleXVel , float *h_ParticleYVel , float *h_ParticleZVel , int *h_ParticleCellNo ,
			int *h_ParticleLastCollide , int ParticleNum ) ;

void PrintEdgeInfo( float *h_EdgeFA , float *h_EdgeFB , float *h_EdgeFC , int CellNum ) ;

void PrintInletInfo( int *h_InletCellNo , int *h_InletEdgeNo , float *h_InletNum , float *h_InletNumR , float *h_InletXVel , float *h_InletYVel , float *h_InletZVel ,
		     float *h_InletNumDen , float *h_InletTemp , float *h_InletArea , int InletFaceNum ) ;


// Setup initial valuve before DSMC initilization.
void InitValue( int *h_IndexParticle , int *h_IndexCell1 , int *h_IndexCell2 , float *h_ParticleXCoord , float *h_ParticleYCoord , 
		float *h_ParticleXVel , float *h_ParticleYVel , float *h_ParticleZVel , int *h_ParticleCellNo , int *h_ParticleLastCollide , 
		int *h_ParticleIn , int *h_SampleParticleNum , float *h_SampleXVel , float *h_SampleYVel , float *h_SampleZVel , 
		float *h_SampleXVelSq , float *h_SampleYVelSq , float *h_SampleZVelSq , int *h_InletCellNo , int *h_InletEdgeNo , 
		float *h_InletNum , float *h_InletNumR , float *h_InletXVel , float *h_InletYVel , float *h_InletZVel , 
		float *h_InletNumDen , float *h_InletTemp , float *h_InletArea , DSMC_DOMAIN *h_Domain ) ;


float ERF( float S ) ;


#endif
