//===================================================================================================
//
// These functions are used to accelerated each component of DSMC simulation on single GPU, including 
// particle movement, index, collision, sampling, and transfer data between host and device. In addition, 
// some functions are used to generate uniform random number and calculate macroscopic properties.
//
//===================================================================================================


#include <cstdio>
#include <cmath>
#include <cuda.h>
#include <cuda_runtime.h>
#include "dsmcsg_cudafunction.h"
#include "parameter.h"


// Particle movement.
void ParticleMove(	float *d_ParticleXCoord , 
			float *d_ParticleYCoord ,
			float *d_ParticleXVel ,
			float *d_ParticleYVel ,
			float *d_ParticleZVel , 
			int *d_ParticleCellNo ,
			int *d_ParticleLastCollide ,
			int *d_ParticleIn ,
			int *d_ParticleWrite ,
			int *d_InletCellNo ,
			int *d_InletEdgeNo ,
			float *d_InletXVel ,
			float *d_InletYVel ,
			float *d_InletTemp ,
			float *d_InletNum ,
			float *d_InletNumR ,
			DSMC_DOMAIN *h_Domain ,
			DSMC_DOMAIN *d_Domain ,
			DSMC_NODE *d_Node ,
			DSMC_CELL *d_Cell ,
			DSMC_BLOCK *d_Block ,
			RANDOM *d_Rand ){

	float		*d_BufferParticleProperty ;
	int		ParticleNum , num1 , num2 ;


	// Update location of all particles.
	Kernel_ParticleMove<<< 8192 , 256 >>>( d_ParticleXCoord , d_ParticleYCoord , d_ParticleXVel , d_ParticleYVel , d_ParticleZVel , 
						d_ParticleCellNo , d_ParticleLastCollide , d_ParticleIn , d_Domain , d_Node , d_Cell , 
						d_Block , d_Rand )  ;
	cudaThreadSynchronize() ;
	

	// Enter New Particle and Move them.
	if ( h_Domain->InletFaceNum > 0 ){
		EnterParticle( d_ParticleXCoord , d_ParticleYCoord , d_ParticleXVel , d_ParticleYVel , d_ParticleZVel , d_ParticleCellNo , 
				d_ParticleLastCollide , d_ParticleIn , d_InletCellNo , d_InletEdgeNo , d_InletXVel , d_InletYVel , 
				d_InletTemp , d_InletNum , d_InletNumR , h_Domain , d_Domain , d_Node , d_Cell , d_Rand ) ;
	}
	cudaThreadSynchronize() ;


	cudaMalloc( (void**)&d_BufferParticleProperty , sizeof(float)*h_Domain->ParticleNum ) ;

	// Scan.
	preallocBlockSums( h_Domain->ParticleNum ) ;
	prescanArray( d_ParticleWrite , d_ParticleIn , h_Domain->ParticleNum ) ;
	cudaThreadSynchronize() ;
	deallocBlockSums() ;


	// Calculate new particle number.
	cudaMemcpy( &num1 , &d_ParticleWrite[h_Domain->ParticleNum-1] , sizeof(int) , cudaMemcpyDeviceToHost ) ;
	cudaMemcpy( &num2 , &d_ParticleIn[h_Domain->ParticleNum-1] , sizeof(int) , cudaMemcpyDeviceToHost ) ;
	ParticleNum	=	num1+num2 ;

	
	// Remove particle
	ParticleRemove<float>( d_ParticleIn , d_ParticleWrite , d_ParticleXCoord , d_BufferParticleProperty , h_Domain->ParticleNum , ParticleNum ) ;
	ParticleRemove<float>( d_ParticleIn , d_ParticleWrite , d_ParticleYCoord , d_BufferParticleProperty , h_Domain->ParticleNum , ParticleNum ) ;
	ParticleRemove<float>( d_ParticleIn , d_ParticleWrite , d_ParticleXVel , d_BufferParticleProperty , h_Domain->ParticleNum , ParticleNum ) ;
	ParticleRemove<float>( d_ParticleIn , d_ParticleWrite , d_ParticleYVel , d_BufferParticleProperty , h_Domain->ParticleNum , ParticleNum ) ;
	ParticleRemove<float>( d_ParticleIn , d_ParticleWrite , d_ParticleZVel , d_BufferParticleProperty , h_Domain->ParticleNum , ParticleNum ) ;
	ParticleRemove<int>( d_ParticleIn , d_ParticleWrite , d_ParticleCellNo , d_BufferParticleProperty , h_Domain->ParticleNum , ParticleNum ) ;
	ParticleRemove<int>( d_ParticleIn , d_ParticleWrite , d_ParticleLastCollide , d_BufferParticleProperty , h_Domain->ParticleNum , ParticleNum ) ;
	ParticleRemove<int>( d_ParticleIn , d_ParticleWrite , d_ParticleIn , d_BufferParticleProperty , h_Domain->ParticleNum , ParticleNum ) ;


	// Update particle number
	h_Domain->ParticleNum	=	ParticleNum ;
	cudaMemcpy( d_Domain , h_Domain , sizeof(DSMC_DOMAIN) , cudaMemcpyHostToDevice ) ;

	cudaFree( d_BufferParticleProperty ) ;
}

//------------------------------------------------------------------------------------------------------------

// Particle momvement (kernel function).
__global__ void Kernel_ParticleMove(	float *d_ParticleXCoord , 
					float *d_ParticleYCoord ,
					float *d_ParticleXVel ,
					float *d_ParticleYVel ,
					float *d_ParticleZVel , 
					int *d_ParticleCellNo ,
					int *d_ParticleLastCollide ,
					int *d_ParticleIn ,
					DSMC_DOMAIN *d_Domain ,
					DSMC_NODE *d_Node ,
					DSMC_CELL *d_Cell ,
					DSMC_BLOCK *d_Block ,
					RANDOM *d_Rand ){

	int	tidx = threadIdx.x , bidx = blockIdx.x , tdimx = blockDim.x , bdimx = gridDim.x , tdim = tdimx*bdimx , idx = bidx*tdimx + tidx ;

	int	CellNo , CellNum[2] , FromCellNo , EdgeNo=0 , TrackNum=0 ; 
	float	MoveTime , RMoveTime , BCLower[2] , BeforeXCoord , BeforeYCoord ;
	bool	Tracking = true ;


	CellNum[0]	=	d_Domain->XCellNum ;
	CellNum[1]	=	d_Domain->YCellNum ;
	BCLower[0]	=	d_Domain->XL ;
	BCLower[1]	=	d_Domain->YL ;

	for ( int i=idx ; i<d_Domain->ParticleNum ; i+=tdim ){
		FromCellNo	=	d_ParticleCellNo[i] ;
		RMoveTime	=	d_Domain->Timestep ;
		MoveTime	=	0. ;
		Tracking	=	true ;
		TrackNum	=	0 ;

		// Particle tracking.
		do {
			TrackNum++ ;
			
			if ( TrackNum == 5 ){
				d_ParticleXCoord[i]	=	d_Cell[d_ParticleCellNo[i]].XCenter ;
				d_ParticleYCoord[i]	=	d_Cell[d_ParticleCellNo[i]].YCenter ;
				break ;
			}

			BeforeXCoord	=	d_ParticleXCoord[i] ;
			BeforeYCoord	=	d_ParticleYCoord[i] ;
			CellNo		=	d_ParticleCellNo[i] ;
			RMoveTime	=	RMoveTime - MoveTime ;


			if ( d_Cell[CellNo].Type == 2 ){
				d_ParticleIn[i]	=	0 ;
				break ;
			}


			// Move i-th particle
			Move( &d_ParticleXCoord[i] , &d_ParticleYCoord[i] , &d_ParticleXVel[i] , &d_ParticleYVel[i] , RMoveTime ) ;
			

			if ( d_Cell[CellNo].Type == 0 ){
				// Update i-th cell no. of particle.
				d_ParticleCellNo[i]	=	GetCellNo( d_ParticleXCoord[i] , d_ParticleYCoord[i] , BCLower[0] , BCLower[1] , d_Domain->XCellSize , 
									   d_Domain->YCellSize , CellNum[0] , CellNum[1] ) ;
				break ;
			}

			// Check particle within the same cell.
			if ( InCell(d_ParticleXCoord[i] , d_ParticleYCoord[i] , d_Node , &d_Cell[CellNo]) ){
				break ;
			}

			// Get edge no. and moving time.
			GetOutEdge( BeforeXCoord , BeforeYCoord , d_ParticleXCoord[i] , d_ParticleYCoord[i] , FromCellNo , &EdgeNo , &MoveTime , 
				    RMoveTime , &d_Cell[CellNo] ) ;


			if ( EdgeNo == -1 ){
				d_ParticleXCoord[i]	=	d_Cell[CellNo].XCenter ;
				d_ParticleYCoord[i]	=	d_Cell[CellNo].YCenter ;
				break ;
			}

			FromCellNo	=	CellNo ;
			CellNo		=	d_Cell[FromCellNo].Neighbor[EdgeNo] ;

			if ( CellNo >= 0 ){
				d_ParticleXCoord[i]	=	BeforeXCoord ;
				d_ParticleYCoord[i]	=	BeforeYCoord ;

				Move( &d_ParticleXCoord[i] , &d_ParticleYCoord[i] , &d_ParticleXVel[i] , &d_ParticleYVel[i] , MoveTime ) ;

				d_ParticleCellNo[i]	=	CellNo ;
			}else{
				d_ParticleXCoord[i]	=	BeforeXCoord ;
				d_ParticleYCoord[i]	=	BeforeYCoord ;

				Move( &d_ParticleXCoord[i] , &d_ParticleYCoord[i] , &d_ParticleXVel[i] , &d_ParticleYVel[i] , MoveTime ) ;

				Tracking = ParticleCollideWall( d_ParticleXVel , d_ParticleYVel , d_ParticleZVel , d_ParticleLastCollide , 
								i , d_Block , EdgeNo , CellNo , &d_Rand[idx] , idx ) ;
			
	
				FromCellNo	=	CellNo ;

				// Removing particle
				if ( !Tracking ){
					d_ParticleIn[i]	=	0 ;
				}
			}
		}while ( Tracking ) ;
	}
}

//------------------------------------------------------------------------------------------------------------

// Update location of particle over a time.
__device__ void Move( float *XCoord , float *YCoord , float *XVel , float *YVel , float time ){
	(*XCoord) += (*XVel) * time ;
	(*YCoord) += (*YVel) * time ;
}

//------------------------------------------------------------------------------------------------------------

// Calculte cell no. of particle.
__device__ int GetCellNo( float XCoord , float YCoord , float XLower , float YLower , float XCellSize , float YCellSize , int XCellNum , int YCellNum ){
	int	XCellNo , YCellNo ;

	XCellNo	=	(XCoord-XLower)/XCellSize ;
	YCellNo	=	(YCoord-YLower)/YCellSize ;

	if ( XCellNo == XCellNum ) XCellNo -= 1 ;
	if ( YCellNo == YCellNum ) YCellNo -= 1 ;

	return	(YCellNo*XCellNum)+XCellNo ;
}

//------------------------------------------------------------------------------------------------------------

__device__ bool InCell( float XCoord , float YCoord , DSMC_NODE *d_Node , DSMC_CELL *_Cell ){
	bool		incell = true ;

	if ( XCoord < d_Node[_Cell->Node[0]].XCoord || XCoord > d_Node[_Cell->Node[1]].XCoord || 
	     YCoord < d_Node[_Cell->Node[0]].YCoord || YCoord > d_Node[_Cell->Node[3]].YCoord ){
		incell = false ;
	}

	return	incell ;
}

//------------------------------------------------------------------------------------------------------------

// Calculate time and edge number for particle collide edge of cell.
__device__ void GetOutEdge(	float BeforeXCoord , 
				float BeforeYCoord , 
				float XCoord , 
				float YCoord ,
				int FromCellNo ,
				int *EdgeNo , 
				float *MoveTime , 
				float RMoveTime , 
				DSMC_CELL *_Cell ){

	int		NBRCellNo ;
	float 		Distance1 , Distance2 , Time=0. ;

	*MoveTime	=	1.e5 ;
	*EdgeNo		=	-1 ;
	
	for ( int i=0 ; i<4 ; i++ ){
		NBRCellNo	=	_Cell->Neighbor[i] ;

		if ( NBRCellNo != FromCellNo ){
			Distance1	=	BeforeXCoord*_Cell->EdgeFA[i] + BeforeYCoord*_Cell->EdgeFB[i] + _Cell->EdgeFC[i] ;
			Distance2	=	XCoord*_Cell->EdgeFA[i] + YCoord*_Cell->EdgeFB[i] + _Cell->EdgeFC[i] ;

			if ( Distance1*Distance2 <= 0. ){
				Time	=	RMoveTime*fabsf( Distance1/(Distance1-Distance2) ) ;
				
				if ( Time < (*MoveTime) ){
					(*MoveTime)	=	Time ;
					(*EdgeNo)	=	i ;
				}
			}
		}
	}
}

//------------------------------------------------------------------------------------------------------------

// Calculate particle properties after reflection.
__device__ bool ParticleCollideWall(	float *d_ParticleXVel ,
					float *d_ParticleYVel ,
					float *d_ParticleZVel ,
					int *d_ParticleLastCollide ,
					int ParticleNo ,
					DSMC_BLOCK *d_Block ,
					int EdgeNo ,
					int BlockNo ,
					RANDOM* _Rand ,
					int idx ){

	bool		Tracking = true ;
	float		Vel[3] , WallVel[2] , MostProbableSpeed=0. , Temp , Mass ; 


	BlockNo		=	-1*(BlockNo+1) ;
	Vel[0]		=	d_ParticleXVel[ParticleNo] ;
	Vel[1]		=	d_ParticleYVel[ParticleNo] ;
	Vel[2]		=	0. ;
	WallVel[0]	=	d_Block[BlockNo].XVel ;
	WallVel[1]	=	d_Block[BlockNo].YVel ;
	Temp		=	d_Block[BlockNo].Temp ;
	Mass		=	MASS ;

	// particle collide inlet boundary.
	if ( d_Block[BlockNo].Type == -3 || d_Block[BlockNo].Type == -4 ){
		Tracking = false ;
		
	// particle collide fully-specular wall.
	}else if ( d_Block[BlockNo].Type == -21 ){
		if ( EdgeNo == 0 || EdgeNo == 2 )
			Vel[1]	=	-Vel[1] ;
		else if ( EdgeNo == 1 || EdgeNo == 3 )
			Vel[0]	=	-Vel[0] ;

		d_ParticleXVel[ParticleNo] 	=	Vel[0] ;
		d_ParticleYVel[ParticleNo]	=	Vel[1] ;
		
		d_ParticleLastCollide[ParticleNo]	=	-1 ;

	// particle collide fully-diffusive wall.
	}else if ( d_Block[BlockNo].Type == -22 ){
		MostProbableSpeed	=	sqrtf( 2.*BOLTZ*Temp/Mass ) ;

		if ( EdgeNo == 0 ){
			RandVel( &Vel[0] , &Vel[2] , MostProbableSpeed , _Rand , idx) ;
			
			Vel[0]	+=	WallVel[0] ;
			Vel[1]	=	sqrtf(-logf(Randn(_Rand,idx)))*MostProbableSpeed + WallVel[1] ;
		}else if ( EdgeNo == 1 ){
			RandVel( &Vel[1] , &Vel[2] , MostProbableSpeed , _Rand , idx) ;
			
			Vel[0]	=	-sqrtf(-logf(Randn(_Rand,idx)))*MostProbableSpeed + WallVel[0] ;
			Vel[1]	+=	WallVel[1] ;
		}else if ( EdgeNo == 2 ){
			RandVel( &Vel[0] , &Vel[2] , MostProbableSpeed , _Rand , idx) ;
			
			Vel[0]	+=	WallVel[0] ;
			Vel[1]	=	-sqrtf(-logf(Randn(_Rand,idx)))*MostProbableSpeed + WallVel[1] ;
		}else if ( EdgeNo == 3 ){
			RandVel( &Vel[1] , &Vel[2] , MostProbableSpeed , _Rand , idx) ;
			
			Vel[0]	=	sqrtf(-logf(Randn(_Rand,idx)))*MostProbableSpeed + WallVel[0] ;
			Vel[1]	+=	WallVel[1] ;
		}

		d_ParticleXVel[ParticleNo] = Vel[0] ;
		d_ParticleYVel[ParticleNo] = Vel[1] ;
		d_ParticleZVel[ParticleNo] = Vel[2] ;
	
		d_ParticleLastCollide[ParticleNo]	=	-1 ;
	}

	return	Tracking ;
}

//------------------------------------------------------------------------------------------------------------

// Remove particle.
template <class TYPE>
void ParticleRemove( int *d_ParticleIn , int *d_ParticleWrite , TYPE *d_ParticleProperty , float *d_BufferParticleProperty , int ParticleNum , int RParticleNum ){
	Kernel_ParticleRemove1<TYPE><<< 8192 , 256 >>>( d_ParticleIn , d_ParticleWrite , d_ParticleProperty , d_BufferParticleProperty , ParticleNum ) ;
	cudaThreadSynchronize() ;
	Kernel_ParticleRemove2<TYPE><<< 8192 , 256 >>>( d_ParticleProperty , d_BufferParticleProperty , RParticleNum ) ;
	cudaThreadSynchronize() ;
}

//------------------------------------------------------------------------------------------------------------

// Move particle to buffer.
template <class TYPE>
__global__ void Kernel_ParticleRemove1( int *d_ParticleIn , int *d_ParticleWrite , TYPE *d_ParticleProperty , float *d_BufferParticleProperty , int ParticleNum ){
	int	tidx = threadIdx.x , bidx = blockIdx.x , tdimx = blockDim.x , bdimx = gridDim.x , tdim = tdimx*bdimx , idx = bidx*tdimx + tidx ;

	for ( int i=idx ; i<ParticleNum ; i+=tdim )
		if ( d_ParticleIn[i] == 1 ) d_BufferParticleProperty[d_ParticleWrite[i]]	=	d_ParticleProperty[i] ;
}

//------------------------------------------------------------------------------------------------------------

// Move particle to global memory.
template <class TYPE>
__global__ void Kernel_ParticleRemove2( TYPE *d_ParticleProperty , float *d_BufferParticleProperty , int ParticleNum ){
	int	tidx = threadIdx.x , bidx = blockIdx.x , tdimx = blockDim.x , bdimx = gridDim.x , tdim = tdimx*bdimx , idx = bidx*tdimx + tidx ;

	for ( int i=idx ; i<ParticleNum ; i+=tdim )
		d_ParticleProperty[i]	=	d_BufferParticleProperty[i] ;
}

//------------------------------------------------------------------------------------------------------------

// Enter new particle.
void EnterParticle(	float *d_ParticleXCoord ,
			float *d_ParticleYCoord ,
			float *d_ParticleXVel ,
			float *d_ParticleYVel ,
			float *d_ParticleZVel ,
			int *d_ParticleCellNo ,
			int *d_ParticleLastCollide ,
			int *d_ParticleIn ,
			int *d_InletCellNo ,
			int *d_InletEdgeNo ,
			float *d_InletXVel ,
			float *d_InletYVel ,
			float *d_InletTemp ,
			float *d_InletNum ,
			float *d_InletNumR ,
			DSMC_DOMAIN *h_Domain ,
			DSMC_DOMAIN *d_Domain ,
			DSMC_NODE *d_Node ,
			DSMC_CELL *d_Cell ,
			RANDOM *d_Rand ){

	int		*d_EnterNum , *d_EnterNo , ParticleNum , NewParticleNum ;
	int		EnterNum , EnterNo ;


	ParticleNum	=	h_Domain->ParticleNum ;


	cudaMalloc((void**)&d_EnterNum , sizeof(int)*h_Domain->InletFaceNum) ;
	cudaMalloc((void**)&d_EnterNo , sizeof(int)*h_Domain->InletFaceNum) ;


	// Calculate Enter Particle Number.
	Kernel_CalEnterNum<<< 1024 , 32 >>>( d_EnterNum , d_InletNum , d_InletNumR , h_Domain->InletFaceNum ) ;
	cudaThreadSynchronize() ;



	// Scan number of enter particle.
	preallocBlockSums( h_Domain->InletFaceNum ) ;
	prescanArray( d_EnterNo , d_EnterNum , h_Domain->InletFaceNum ) ;
	cudaThreadSynchronize() ;
	deallocBlockSums() ;



	// Enter particle.
	Kernel_EnterParticle<<< 8192 , 32 >>>( d_ParticleXCoord , d_ParticleYCoord , d_ParticleXVel , d_ParticleYVel , d_ParticleZVel , 
					d_ParticleCellNo , d_ParticleLastCollide , d_ParticleIn , d_InletCellNo , d_InletEdgeNo , 
					d_InletXVel , d_InletYVel , d_InletTemp , d_EnterNum , d_EnterNo , h_Domain->InletFaceNum , 
					ParticleNum , d_Node , d_Cell , d_Rand ) ;


	cudaMemcpy( &EnterNum , &d_EnterNum[h_Domain->InletFaceNum-1] , sizeof(int) , cudaMemcpyDeviceToHost ) ;
	cudaMemcpy( &EnterNo , &d_EnterNo[h_Domain->InletFaceNum-1] , sizeof(int) , cudaMemcpyDeviceToHost ) ;


	// Updata Total Particle Number save at NewParticleNum.
	NewParticleNum	=	ParticleNum + EnterNum + EnterNo ;



	// Move New Particle.
	Kernel_NewParticleMove<<< 8192 , 32 >>>( d_ParticleXCoord , d_ParticleYCoord , d_ParticleXVel , d_ParticleYVel , d_ParticleCellNo , d_ParticleIn , 
						NewParticleNum , ParticleNum , d_Domain , d_Rand ) ;



	// Updata Total Particle Number.
	h_Domain->ParticleNum	=	NewParticleNum ;


	cudaFree(d_EnterNum) ;
	cudaFree(d_EnterNo) ;
}

//------------------------------------------------------------------------------------------------------------

// Calculate number of enter particle.
__global__ void Kernel_CalEnterNum( int *d_EnterNum , float *d_InletNum , float *d_InletNumR , int InletFaceNum ){

	int	tidx = threadIdx.x , bidx = blockIdx.x , tdimx = blockDim.x , bdimx = gridDim.x , tdim = tdimx*bdimx , idx = bidx*tdimx + tidx ;
	float	EnterNum ;

	for ( int i=idx ; i<InletFaceNum ; i+=tdim ){
		EnterNum	=	d_InletNum[i] + d_InletNumR[i] ;
		d_EnterNum[i]	=	EnterNum ;
		d_InletNumR[i]	=	EnterNum - d_EnterNum[i] ;
	}
}

//------------------------------------------------------------------------------------------------------------

// Enter new particle.
__global__ void Kernel_EnterParticle(	float *d_ParticleXCoord , 
					float *d_ParticleYCoord , 
					float *d_ParticleXVel , 
					float *d_ParticleYVel , 
					float *d_ParticleZVel , 
					int *d_ParticleCellNo , 
					int *d_ParticleLastCollide , 
					int *d_ParticleIn , 
					int *d_InletCellNo , 
					int *d_InletEdgeNo , 
					float *d_InletXVel , 
					float *d_InletYVel , 
					float *d_InletTemp , 
					int *d_EnterNum , 
					int *d_EnterNo , 
					int InletFaceNum , 
					int ParticleNum , 
					DSMC_NODE *d_Node ,
					DSMC_CELL *d_Cell ,
					RANDOM *d_Rand ){

	int	tidx = threadIdx.x , bidx = blockIdx.x , tdimx = blockDim.x , bdimx = gridDim.x , idx = bidx*tdimx + tidx ;

	float	MostProbableSpeed , Mass , NVel , PVel ;
	float	FS1 , FS2 , QA , VP , U , UN , A ;
	int	Node1 , Node2 ;


	Mass	=	MASS ;


	for ( int i=bidx ; i<InletFaceNum ; i+=bdimx ){
		MostProbableSpeed	=	sqrtf(2.*BOLTZ*d_InletTemp[i]/Mass) ;


		if ( d_InletEdgeNo[i] == 0 ){
			NVel	=	d_InletYVel[i]/MostProbableSpeed ;
			PVel	=	d_InletXVel[i]/MostProbableSpeed ;
		}else if ( d_InletEdgeNo[i] == 1 ){
			NVel	=	-d_InletXVel[i]/MostProbableSpeed ;
			PVel	=	d_InletYVel[i]/MostProbableSpeed ;
		}else if ( d_InletEdgeNo[i] == 2 ){
			NVel	=	-d_InletYVel[i]/MostProbableSpeed ;
			PVel	=	-d_InletXVel[i]/MostProbableSpeed ;
		}else if ( d_InletEdgeNo[i] == 3 ){
			NVel	=	d_InletXVel[i]/MostProbableSpeed ;
			PVel	=	-d_InletYVel[i]/MostProbableSpeed ;
		}


		FS1	=	NVel + sqrtf(NVel*NVel+2.) ;
		FS2	=	0.5 * (1.+NVel*(2.*NVel-FS1)) ;

		for ( int j=(ParticleNum+d_EnterNo[i]+tidx) ; j<(ParticleNum+d_EnterNo[i]+d_EnterNum[i]) ; j+=tdimx ){
			if ( fabsf(NVel*MostProbableSpeed) > 1.e-6 ){
				QA	=	3. ;
				if ( NVel < -3. ) QA = fabsf(NVel)+1. ;

				do{
					do{
						U	=	-QA + 2.*QA*Randn(&d_Rand[idx],idx) ;
						UN	=	U + NVel ;
					}while ( UN < 0. ) ;

					A	=	(2.*UN/FS1) * expf(FS2-U*U) ;
				}while ( A < Randn(&d_Rand[idx],idx) ) ;

				RandVel( &VP , &d_ParticleZVel[j] , MostProbableSpeed , &d_Rand[idx] , idx ) ;

				if ( d_InletEdgeNo[i] == 0 ){
					d_ParticleXVel[j]	=	VP + PVel*MostProbableSpeed ;
					d_ParticleYVel[j]	=	UN*MostProbableSpeed ;
				}else if ( d_InletEdgeNo[i] == 1 ){
					d_ParticleXVel[j]	=	-UN*MostProbableSpeed ;
					d_ParticleYVel[j]	=	VP + PVel*MostProbableSpeed ;
				}else if ( d_InletEdgeNo[i] == 2 ){
					d_ParticleXVel[j]	=	-(VP + PVel*MostProbableSpeed) ;
					d_ParticleYVel[j]	=	-UN*MostProbableSpeed ;
				}else if ( d_InletEdgeNo[i] == 3 ){
					d_ParticleXVel[j]	=	UN*MostProbableSpeed ;
					d_ParticleYVel[j]	=	-(VP + PVel*MostProbableSpeed) ;
				}

			}else{
				RandVel( &VP , &d_ParticleZVel[j] , MostProbableSpeed , &d_Rand[idx] , idx ) ;

				if ( d_InletEdgeNo[i] == 0 ){
					d_ParticleXVel[j]	=	VP + PVel*MostProbableSpeed ;
					d_ParticleYVel[j]	=	sqrtf(-logf(Randn(&d_Rand[idx],idx)))*MostProbableSpeed + NVel*MostProbableSpeed ;
				}else if ( d_InletEdgeNo[i] == 1 ){
					d_ParticleXVel[j]	=	-sqrtf(-logf(Randn(&d_Rand[idx],idx)))*MostProbableSpeed - NVel*MostProbableSpeed ;
					d_ParticleYVel[j]	=	VP + PVel*MostProbableSpeed ;
				}else if ( d_InletEdgeNo[i] == 2 ){
					d_ParticleXVel[j]	=	-(VP + PVel*MostProbableSpeed) ;
					d_ParticleYVel[j]	=	-sqrtf(-logf(Randn(&d_Rand[idx],idx)))*MostProbableSpeed - NVel*MostProbableSpeed ;
				}else if ( d_InletEdgeNo[i] == 3 ){
					d_ParticleXVel[j]	=	sqrtf(-logf(Randn(&d_Rand[idx],idx)))*MostProbableSpeed + NVel*MostProbableSpeed ;
					d_ParticleYVel[j]	=	-(VP + PVel*MostProbableSpeed) ;
				}
			}



			VP	=	Randn(&d_Rand[idx],idx) ;


			Node1	=	d_Cell[d_InletCellNo[i]].Node[d_InletEdgeNo[i]] ;
			if ( d_InletEdgeNo[i]==3 )
				Node2	=	d_Cell[d_InletCellNo[i]].Node[0] ;
			else
				Node2	=	d_Cell[d_InletCellNo[i]].Node[d_InletEdgeNo[i]+1] ;


			d_ParticleXCoord[j]	=	d_Node[Node1].XCoord + (d_Node[Node2].XCoord-d_Node[Node1].XCoord)*VP ;
			d_ParticleYCoord[j]	=	d_Node[Node1].YCoord + (d_Node[Node2].YCoord-d_Node[Node1].YCoord)*VP ;

			d_ParticleCellNo[j]	=	d_InletCellNo[i] ;
			d_ParticleLastCollide[j]=	-1 ;
			d_ParticleIn[j]		=	1 ;
		}
	}
}

//------------------------------------------------------------------------------------------------------------

// Move new particle.
__global__ void Kernel_NewParticleMove( float *d_ParticleXCoord , 
					float *d_ParticleYCoord , 
					float *d_ParticleXVel , 
					float *d_ParticleYVel , 
					int *d_ParticleCellNo , 
					int *d_ParticleIn , 
					int ParticleNum , 
					int Num ,
					DSMC_DOMAIN *d_Domain , 
					RANDOM *d_Rand ){


	int	tidx = threadIdx.x , bidx = blockIdx.x , tdimx = blockDim.x , bdimx = gridDim.x , tdim = tdimx*bdimx , idx = bidx*tdimx + tidx ;
	
	float	DomainXL , DomainXH , DomainYL , DomainYH , XCellSize , YCellSize , Time ;
	int	XCellNum , YCellNum ;


	DomainXL	=	d_Domain->XL ;
	DomainXH	=	d_Domain->XH ;
	DomainYL	=	d_Domain->YL ;
	DomainYH	=	d_Domain->YH ;
	XCellSize	=	d_Domain->XCellSize ;
	YCellSize	=	d_Domain->YCellSize ;
	XCellNum	=	d_Domain->XCellNum ;
	YCellNum	=	d_Domain->YCellNum ;



	for ( int i=(idx+Num) ; i<ParticleNum ; i+=tdim ){
		Time	=	d_Domain->Timestep * Randn(&d_Rand[idx],idx) ;
		
		Move( &d_ParticleXCoord[i] , &d_ParticleYCoord[i] , &d_ParticleXVel[i] , &d_ParticleYVel[i] , Time ) ;

		if ( d_ParticleXCoord[i] < DomainXL || d_ParticleXCoord[i] > DomainXH || d_ParticleYCoord[i] < DomainYL || d_ParticleYCoord[i] > DomainYH ){
			d_ParticleIn[i] = 0 ;
		}else{
			d_ParticleCellNo[i]	=	GetCellNo( d_ParticleXCoord[i] , d_ParticleYCoord[i] , DomainXL , DomainYL , XCellSize , YCellSize , XCellNum , 
								   YCellNum ) ;
		}
	}
}

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

// Index.
void Index( int *d_IndexParticle , int *d_IndexCell1 , int *d_IndexCell2 , int *d_ParticleCellNo , int ParticleNum , int CellNum , DSMC_TIME *Time ){

	cudaMemset( d_IndexCell2 , 0 , sizeof(int)*CellNum ) ;


	Kernel_Index1<<< 8192 , 32 >>>( d_IndexCell2 , d_ParticleCellNo , ParticleNum ) ;
	cudaThreadSynchronize() ;


	// CUDA SDK
	preallocBlockSums( CellNum ) ;
	prescanArray( d_IndexCell1, d_IndexCell2, CellNum ) ;
	cudaThreadSynchronize() ;
	deallocBlockSums() ;


	cudaMemset( d_IndexCell2 , 0 , sizeof(int)*CellNum ) ;
	Kernel_Index2<<< 8192 , 32 >>>( d_IndexParticle , d_IndexCell1 , d_IndexCell2 , d_ParticleCellNo , ParticleNum , CellNum ) ;
	cudaThreadSynchronize() ;
}

//------------------------------------------------------------------------------------------------------------

// Count number of particle in each cell.
__global__ void Kernel_Index1( int *d_IndexCell2 , int *d_ParticleCellNo , int ParticleNum ){
	int	tidx = threadIdx.x , bidx = blockIdx.x , tdimx = blockDim.x , bdimx = gridDim.x , tdim = tdimx*bdimx , idx = bidx*tdimx + tidx ;
	int	CellNo ;

	for ( int i=idx ; i<ParticleNum ; i+=tdim ){
		CellNo	=	d_ParticleCellNo[i] ;
		atomicAdd( &d_IndexCell2[CellNo] , 1 ) ;
	}
}

//------------------------------------------------------------------------------------------------------------

// Create relation between cells and particle.
__global__ void Kernel_Index2( int *d_IndexParticle , int *d_IndexCell1 , int *d_IndexCell2 , int *d_ParticleCellNo , int ParticleNum , int CellNum ){
	int	tidx = threadIdx.x , bidx = blockIdx.x , tdimx = blockDim.x , bdimx = gridDim.x , tdim = tdimx*bdimx , idx = bidx*tdimx + tidx ;
	int	CellNo , RParticleNo , BufferNo ;
	
	for ( int i=idx ; i<ParticleNum ; i+=tdim ){
		CellNo	=	d_ParticleCellNo[i] ;
		BufferNo=	atomicAdd( &d_IndexCell2[CellNo] , 1 ) ;

		RParticleNo	=	d_IndexCell1[CellNo] + BufferNo ;
		d_IndexParticle[RParticleNo]  =       i ;
	}
}

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

// Collision.
void Collision(	int *d_IndexParticle ,
		int *d_IndexCell1 ,
		int *d_IndexCell2 ,
		int *d_ParticleLastCollide ,
		float *d_ParticleXVel ,
		float *d_ParticleYVel ,
		float *d_ParticleZVel ,
		int *d_SampleParticleNum ,
		DSMC_DOMAIN *d_Domain ,
		DSMC_CELL *d_Cell ,
		RANDOM *d_Rand ){


	Kernel_Collision<<< 8192 , 32 >>>( d_IndexParticle , d_IndexCell1 , d_IndexCell2 , d_ParticleLastCollide , d_ParticleXVel , 
					d_ParticleYVel , d_ParticleZVel , d_SampleParticleNum , d_Domain , d_Cell , d_Rand ) ;
	cudaThreadSynchronize() ;
}

//------------------------------------------------------------------------------------------------------------

// Collision (kernel function).
__global__ void Kernel_Collision(	int *d_IndexParticle ,
					int *d_IndexCell1 ,
					int *d_IndexCell2 ,
					int *d_ParticleLastCollide ,
					float *d_ParticleXVel ,
					float *d_ParticleYVel ,
					float *d_ParticleZVel ,
					int *d_SampleParticleNum ,
					DSMC_DOMAIN *d_Domain ,
					DSMC_CELL *d_Cell ,
					RANDOM *d_Rand ){

	int		tidx = threadIdx.x , bidx = blockIdx.x , tdimx = blockDim.x , bdimx = gridDim.x , tdim = tdimx*bdimx , idx = bidx*tdimx + tidx ;

	int		SampleParticleNum , CollisionPairs , SelectParticleNo[2] ;
	float		RCollisionPairs , VelCrossSection=0. , MaxVelCrossSection , Relative_Vel[3] , RelativeVel , RelativeVelSq ;


	for ( int i=idx ; i<d_Domain->CellNum ; i+=tdim ){
		if ( d_Cell[i].Type == 2 ) continue ;

		SampleParticleNum	=	d_SampleParticleNum[i] ;

		if ( SampleParticleNum > 1 )
			SampleParticleNum = SampleParticleNum/(float)(d_Domain->SampleNum) ;
		else
			SampleParticleNum = d_IndexCell2[i] ;

		RCollisionPairs			=	0.5*d_IndexCell2[i]*SampleParticleNum*d_Domain->ParticleWeight*d_Cell[i].MaxCrossSectionSpeed*d_Domain->Timestep/d_Domain->CellVolume + d_Cell[i].RemainderCollisionPair ;
		CollisionPairs			=	RCollisionPairs ;
		d_Cell[i].RemainderCollisionPair=	RCollisionPairs - CollisionPairs ;

		if ( CollisionPairs > 0){
			if ( d_IndexCell2[i] < 2 ){
				d_Cell[i].RemainderCollisionPair	+=	CollisionPairs ;
			}else{
				MaxVelCrossSection	=	d_Cell[i].MaxCrossSectionSpeed ;
				for ( int j=0 ; j<CollisionPairs ; j++ ){
					VelCrossSection = SelectParticle( d_IndexParticle , &d_IndexCell1[i] , &d_IndexCell2[i] , d_ParticleLastCollide , 
									  d_ParticleXVel , d_ParticleYVel , d_ParticleZVel , SelectParticleNo , Relative_Vel , 
									  &RelativeVel , &RelativeVelSq , &d_Rand[idx] , idx ) ;

					if ( VelCrossSection > MaxVelCrossSection ) MaxVelCrossSection = VelCrossSection ;


					if ( Randn(&d_Rand[idx] , idx) < (VelCrossSection/d_Cell[i].MaxCrossSectionSpeed) ){
						d_ParticleLastCollide[SelectParticleNo[0]]	=	SelectParticleNo[1] ;
						d_ParticleLastCollide[SelectParticleNo[1]]	=	SelectParticleNo[0] ;

	
						// To calculate velocity of particle in post-collision.
						Elastic( &d_ParticleXVel[SelectParticleNo[0]] , &d_ParticleYVel[SelectParticleNo[0]] , &d_ParticleZVel[SelectParticleNo[0]] ,
							 &d_ParticleXVel[SelectParticleNo[1]] , &d_ParticleYVel[SelectParticleNo[1]] , &d_ParticleZVel[SelectParticleNo[1]] ,
							 &RelativeVel , &d_Rand[idx] , idx ) ;
					}
				}
				d_Cell[i].MaxCrossSectionSpeed = MaxVelCrossSection ;
			}
		}
	}
}

//------------------------------------------------------------------------------------------------------------

// Select two particles to collide.
__device__ float SelectParticle(int *d_IndexParticle , 
				int *_IndexCell1 , 
				int *_IndexCell2 ,
				int *d_ParticleLastCollide ,
				float *d_ParticleXVel ,
				float *d_ParticleYVel ,
				float *d_ParticleZVel ,
				int *SelectParticleNo , 
				float *Relative_Vel , 
				float *RelativeVel ,
				float *RelativeVelSq ,
				RANDOM *_Rand ,
				int idx ){

	float		VelCrossSection , MixCrossSection , MixRefTemp , ReduceMass , MixVisTempIdx ;
	int		Select , Select1 , Select2 ;
	bool		BSelect = true ;


	MixCrossSection	=	MIX_CROSS_SECTION ;
	MixRefTemp	=	MIX_RT ;
	ReduceMass	=	MIX_REDUCE_MASS ;
	MixVisTempIdx	=	MIX_VTI ;


	Select	=	(*_IndexCell2-0.00001) * Randn(_Rand , idx) ;
	Select1	=	d_IndexParticle[(*_IndexCell1)+Select] ;
	while ( BSelect ){
		Select	=	(*_IndexCell2-0.00001) * Randn(_Rand , idx) ;
		Select2 =	d_IndexParticle[(*_IndexCell1)+Select] ;
		
		if ( Select1 != Select2 ){
			if ( d_ParticleLastCollide[Select1] != Select2 || Select1 != d_ParticleLastCollide[Select2] || 
			     d_ParticleLastCollide[Select1] < 0 || d_ParticleLastCollide[Select2] < 0 ){
				break ;
			}else if ( (*_IndexCell2) == 2 ){
				break ;
			}
		}
	}


	SelectParticleNo[0] = Select1 ;
	SelectParticleNo[1] = Select2 ;
	

	Relative_Vel[0]	=	d_ParticleXVel[Select1] - d_ParticleXVel[Select2] ;
	Relative_Vel[1]	=	d_ParticleYVel[Select1] - d_ParticleYVel[Select2] ;
	Relative_Vel[2]	=	d_ParticleZVel[Select1] - d_ParticleZVel[Select2] ;


	*RelativeVelSq	=	Relative_Vel[0]*Relative_Vel[0] + Relative_Vel[1]*Relative_Vel[1] + Relative_Vel[2]*Relative_Vel[2] ;
	*RelativeVel	=	sqrtf((*RelativeVelSq)) ;

	VelCrossSection	=	(*RelativeVel) * MixCrossSection * (powf((2.*BOLTZ*MixRefTemp/(ReduceMass*(*RelativeVelSq))),(MixVisTempIdx-0.5))) / GAM(2.5-MIX_VTI) ;

	return	VelCrossSection ;
}

//------------------------------------------------------------------------------------------------------------

// Calculate particles properties after collision.
__device__ void Elastic(float *_1ParticleXVel ,
			float *_1ParticleYVel ,
			float *_1ParticleZVel ,
			float *_2ParticleXVel ,
			float *_2ParticleYVel ,
			float *_2ParticleZVel ,
			float *RelativeVel ,
			RANDOM *_Rand ,
			int idx ){

	float		ReduceMass = MIX_REDUCE_MASS , Mass = MASS ;
	float		MassRatio[2] ,  RelativeVelPost[3] , CenterMassVel[3] ;
	float		CosineElevationAngle , AzimuthAngle , Angle ;


	for ( int i=0 ; i<3 ; i++ ) RelativeVelPost[i] = 0. ;

	CosineElevationAngle	=	0. ;
	AzimuthAngle		=	0. ;
	Angle			=	0. ;


	MassRatio[0]	=	ReduceMass/Mass ;
	MassRatio[1]	=	ReduceMass/Mass ;


	// Calculate the components of the centre-of-mass velocity, eqn (2.1)
	CenterMassVel[0]	=	MassRatio[0]*(*_1ParticleXVel) + MassRatio[1]*(*_2ParticleXVel) ;
	CenterMassVel[1]	=	MassRatio[0]*(*_1ParticleYVel) + MassRatio[1]*(*_2ParticleYVel) ;
	CenterMassVel[2]	=	MassRatio[0]*(*_1ParticleZVel) + MassRatio[1]*(*_2ParticleZVel) ;


	// VHS
	CosineElevationAngle	=	2.*Randn(_Rand , idx) - 1. ;
	Angle			=	sqrtf(1.-(CosineElevationAngle*CosineElevationAngle)) ;

	RelativeVelPost[0]	=	CosineElevationAngle*(*RelativeVel) ;

	AzimuthAngle		=	2.*PI*Randn(_Rand , idx) ;

	RelativeVelPost[1]	=	Angle * cosf(AzimuthAngle) * (*RelativeVel) ;
	RelativeVelPost[2]	=	Angle * sinf(AzimuthAngle) * (*RelativeVel) ;


	(*_1ParticleXVel)	=	CenterMassVel[0] + RelativeVelPost[0]*MassRatio[1] ;
	(*_1ParticleYVel)	=	CenterMassVel[1] + RelativeVelPost[1]*MassRatio[1] ;
	(*_1ParticleZVel)	=	CenterMassVel[2] + RelativeVelPost[2]*MassRatio[1] ;

	(*_2ParticleXVel)	=	CenterMassVel[0] - RelativeVelPost[0]*MassRatio[0] ;
	(*_2ParticleYVel)	=	CenterMassVel[1] - RelativeVelPost[1]*MassRatio[0] ;
	(*_2ParticleZVel)	=	CenterMassVel[2] - RelativeVelPost[2]*MassRatio[0] ;
}

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

// Sample.
void Sample( 	int *d_IndexParticle ,
		int *d_IndexCell1 ,
		int *d_IndexCell2 ,
		int *d_SampleParticleNum ,
		float *d_SampleXVel , 
		float *d_SampleYVel , 
		float *d_SampleZVel , 
		float *d_SampleXVelSq , 
		float *d_SampleYVelSq ,
		float *d_SampleZVelSq ,
		float *d_ParticleXVel ,
		float *d_ParticleYVel ,
		float *d_ParticleZVel ,
		DSMC_CELL *d_Cell ,
		DSMC_DOMAIN *d_Domain ){

	Kernel_Sample<<< 8192 , 32 >>>( d_IndexParticle , d_IndexCell1 , d_IndexCell2 , d_SampleParticleNum , d_SampleXVel , 
					d_SampleYVel , d_SampleZVel , d_SampleXVelSq , d_SampleYVelSq , d_SampleZVelSq , d_ParticleXVel , 
					d_ParticleYVel , d_ParticleZVel , d_Cell , d_Domain ) ;
	cudaThreadSynchronize() ;
}

//------------------------------------------------------------------------------------------------------------

// Sample (kernel function).
__global__ void Kernel_Sample(	int *d_IndexParticle ,
				int *d_IndexCell1 ,
				int *d_IndexCell2 ,
				int *d_SampleParticleNum ,
				float *d_SampleXVel , 
				float *d_SampleYVel , 
				float *d_SampleZVel , 
				float *d_SampleXVelSq , 
				float *d_SampleYVelSq ,
				float *d_SampleZVelSq ,
				float *d_ParticleXVel ,
				float *d_ParticleYVel ,
				float *d_ParticleZVel ,
				DSMC_CELL *d_Cell ,
				DSMC_DOMAIN *d_Domain ){

	int	tidx = threadIdx.x , bidx = blockIdx.x , tdimx = blockDim.x , bdimx = gridDim.x , tdim = tdimx*bdimx , idx = bidx*tdimx + tidx ;
	int	ParticleNo , ParticleNum , IndexCell1 ;

	__shared__ float s_XVel[32] , s_YVel[32] , s_ZVel[32] , s_XVelSq[32] , s_YVelSq[32] , s_ZVelSq[32] ;


	s_XVel[tidx]	=	0. ;
	s_YVel[tidx]	=	0. ;
	s_ZVel[tidx]	=	0. ;
	s_XVelSq[tidx]	=	0. ;
	s_YVelSq[tidx]	=	0. ;
	s_ZVelSq[tidx]	=	0. ;


	if ( idx==0 ) (d_Domain->SampleNum)++ ;
	
	for ( int i=idx ; i<d_Domain->CellNum ; i+=tdim ){
		if ( d_Cell[i].Type == 2 ) continue ;

		s_XVel[tidx]	=	0. ;
		s_YVel[tidx]	=	0. ;
		s_ZVel[tidx]	=	0. ;
		s_XVelSq[tidx]	=	0. ;
		s_YVelSq[tidx]	=	0. ;
		s_ZVelSq[tidx]	=	0. ;

		ParticleNum	=	d_IndexCell2[i] ;
		IndexCell1	=	d_IndexCell1[i] ;

		for ( int j=0 ; j<ParticleNum ; j++ ){
			ParticleNo	=	d_IndexParticle[IndexCell1+j] ;

			s_XVel[tidx]	+=	d_ParticleXVel[ParticleNo] ;
			s_YVel[tidx]	+=	d_ParticleYVel[ParticleNo] ;
			s_ZVel[tidx]	+=	d_ParticleZVel[ParticleNo] ;
			s_XVelSq[tidx]	+=	(d_ParticleXVel[ParticleNo]*d_ParticleXVel[ParticleNo]) ;
			s_YVelSq[tidx]	+=	(d_ParticleYVel[ParticleNo]*d_ParticleYVel[ParticleNo]) ;
			s_ZVelSq[tidx]	+=	(d_ParticleZVel[ParticleNo]*d_ParticleZVel[ParticleNo]) ;
		}
	
		d_SampleParticleNum[i]	+=	ParticleNum ;

		d_SampleXVel[i]		+=	s_XVel[tidx] ;
		d_SampleYVel[i]		+=	s_YVel[tidx] ;
		d_SampleZVel[i]		+=	s_ZVel[tidx] ;
		d_SampleXVelSq[i]	+=	s_XVelSq[tidx] ;
		d_SampleYVelSq[i]	+=	s_YVelSq[tidx] ;
		d_SampleZVelSq[i]	+=	s_ZVelSq[tidx] ;
	}
}

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

// Calculate macroscopic properties.
void CalculateResult( int *h_SampleParticleNum , float *h_SampleXVel , float *h_SampleYVel , float *h_SampleZVel , float *h_SampleXVelSq , 
		      float *h_SampleYVelSq , float *h_SampleZVelSq , DSMC_OUTPUT *h_Result , DSMC_DOMAIN *h_Domain , DSMC_CELL *h_Cell ){

	float		Mass ;
	float		Weight , SumParticleNum , SumMass , SumMassVel[3] , SumMassVelSq ;
	float		NumDensity , Density , Vel[3] , VelSq , TotalTemp , TransTemp ; 
	float		AveParticleNum ;

	Mass	= 	MASS ;

	for ( int i=0 ; i<h_Domain->CellNum ; i++ ){
		if ( h_Cell[i].Type == 2 ) continue ;

		Weight		=	h_Domain->ParticleWeight/(h_Domain->CellVolume*h_Domain->SampleNum) ;

		SumParticleNum	=	0. ;
		SumMass		=	0. ;
		SumMassVel[0]	=	0. ;
		SumMassVel[1]	=	0. ;
		SumMassVel[2]	=	0. ;
		SumMassVelSq	=	0. ;

		SumParticleNum	+=	h_SampleParticleNum[i] ;
		SumMass		+=	(Mass*h_SampleParticleNum[i]) ;

		SumMassVel[0]	+=	Mass*h_SampleXVel[i] ;
		SumMassVel[1]	+=	Mass*h_SampleYVel[i] ;
		SumMassVel[2]	+=	Mass*h_SampleZVel[i] ;

		SumMassVelSq	+=	Mass * (h_SampleXVelSq[i]+h_SampleYVelSq[i]+h_SampleZVelSq[i]) ;


		// Calculate Results
		NumDensity	=	SumParticleNum * Weight ;
		Density		=	NumDensity*SumMass/SumParticleNum ;
		for ( int j=0 ; j<3 ; j++ )
			Vel[j]		=	SumMassVel[j]/SumMass ;
	
		VelSq		=	Vel[0]*Vel[0] + Vel[1]*Vel[1] + Vel[2]*Vel[2] ;
		TransTemp	=	(SumMassVelSq - SumMass*VelSq)/(3.*BOLTZ*SumParticleNum) ;
		
		TotalTemp	=	TransTemp ;

		AveParticleNum	=	SumParticleNum/h_Domain->SampleNum ;


		h_Result[i].NumDensity	=	NumDensity ;
		h_Result[i].Density	=	Density	;
		h_Result[i].XVel	=	Vel[0] ;
		h_Result[i].YVel	=	Vel[1] ;
		h_Result[i].ZVel	=	Vel[2] ;
		h_Result[i].Temp	=	TotalTemp ;
		h_Result[i].AveParticleNum=	AveParticleNum ;
	}
}

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

// Transfer data between host and device (direction: 0 is from host to device, and 1 is device to host).
void MemoryCopy(int direction , 
		int *d_IndexParticle		, int *h_IndexParticle , 
		int *d_IndexCell1		, int *h_IndexCell1 , 
		int *d_IndexCell2		, int *h_IndexCell2 , 
		float *d_ParticleXCoord		, float *h_ParticleXCoord , 
		float *d_ParticleYCoord		, float *h_ParticleYCoord , 
		float *d_ParticleXVel		, float *h_ParticleXVel , 
		float *d_ParticleYVel		, float *h_ParticleYVel , 
		float *d_ParticleZVel		, float *h_ParticleZVel , 
		int *d_ParticleCellNo		, int *h_ParticleCellNo , 
		int *d_ParticleLastCollide	, int *h_ParticleLastCollide , 
		int *d_ParticleIn		, int *h_ParticleIn , 
		int *d_SampleParticleNum	, int *h_SampleParticleNum , 
		float *d_SampleXVel		, float *h_SampleXVel , 
		float *d_SampleYVel		, float *h_SampleYVel , 
		float *d_SampleZVel		, float *h_SampleZVel , 
		float *d_SampleXVelSq		, float *h_SampleXVelSq , 
		float *d_SampleYVelSq		, float *h_SampleYVelSq , 
		float *d_SampleZVelSq		, float *h_SampleZVelSq , 
		int *d_InletCellNo 		, int *h_InletCellNo ,
		int *d_InletEdgeNo 		, int *h_InletEdgeNo ,
		float *d_InletNum 		, float *h_InletNum ,
		float *d_InletNumR 		, float *h_InletNumR ,
		float *d_InletXVel 		, float *h_InletXVel ,
		float *d_InletYVel 		, float *h_InletYVel ,
		float *d_InletZVel 		, float *h_InletZVel ,
		float *d_InletNumDen 		, float *h_InletNumDen ,
		float *d_InletTemp 		, float *h_InletTemp ,
		float *d_InletArea		, float *h_InletArea ,
		DSMC_DOMAIN *d_Domain		, DSMC_DOMAIN *h_Domain ,
		DSMC_NODE *d_Node		, DSMC_NODE *h_Node ,
		DSMC_CELL *d_Cell		, DSMC_CELL *h_Cell ,
		DSMC_BLOCK *d_Block		, DSMC_BLOCK *h_Block ){

	int	NodeNum , CellNum , BlockNum , MaxParticleNum , InletFaceNum ;

	NodeNum		=	h_Domain->NodeNum ;
	CellNum		=	h_Domain->CellNum ;
	BlockNum	=	h_Domain->BlockNum ;
	MaxParticleNum	=	h_Domain->MaxParticleNum ;
	InletFaceNum	=	h_Domain->InletFaceNum ;


	if ( direction == 0 ){
		cudaMemcpy(d_Domain		, h_Domain		, sizeof(DSMC_DOMAIN)		, cudaMemcpyHostToDevice) ;
		cudaMemcpy(d_Node		, h_Node		, sizeof(DSMC_NODE)*NodeNum	, cudaMemcpyHostToDevice) ;
		cudaMemcpy(d_Cell		, h_Cell		, sizeof(DSMC_CELL)*CellNum	, cudaMemcpyHostToDevice) ;
		cudaMemcpy(d_Block		, h_Block		, sizeof(DSMC_BLOCK)*BlockNum	, cudaMemcpyHostToDevice) ;
		
		cudaMemcpy(d_IndexParticle	, h_IndexParticle	, sizeof(int)*MaxParticleNum	, cudaMemcpyHostToDevice) ;
		cudaMemcpy(d_IndexCell1		, h_IndexCell1		, sizeof(int)*CellNum		, cudaMemcpyHostToDevice) ;
		cudaMemcpy(d_IndexCell2		, h_IndexCell2		, sizeof(int)*CellNum		, cudaMemcpyHostToDevice) ;
	
		cudaMemcpy(d_ParticleXCoord	, h_ParticleXCoord	, sizeof(float)*MaxParticleNum	, cudaMemcpyHostToDevice) ;
		cudaMemcpy(d_ParticleYCoord	, h_ParticleYCoord	, sizeof(float)*MaxParticleNum	, cudaMemcpyHostToDevice) ;
		cudaMemcpy(d_ParticleXVel	, h_ParticleXVel	, sizeof(float)*MaxParticleNum	, cudaMemcpyHostToDevice) ;
		cudaMemcpy(d_ParticleYVel	, h_ParticleYVel	, sizeof(float)*MaxParticleNum	, cudaMemcpyHostToDevice) ;
		cudaMemcpy(d_ParticleZVel	, h_ParticleZVel	, sizeof(float)*MaxParticleNum	, cudaMemcpyHostToDevice) ;
		cudaMemcpy(d_ParticleCellNo	, h_ParticleCellNo	, sizeof(int)*MaxParticleNum	, cudaMemcpyHostToDevice) ;
		cudaMemcpy(d_ParticleLastCollide, h_ParticleLastCollide	, sizeof(int)*MaxParticleNum	, cudaMemcpyHostToDevice) ;
		cudaMemcpy(d_ParticleIn		, h_ParticleIn		, sizeof(int)*MaxParticleNum	, cudaMemcpyHostToDevice) ;

		cudaMemcpy(d_SampleParticleNum	, h_SampleParticleNum	, sizeof(int)*CellNum		, cudaMemcpyHostToDevice) ;
		cudaMemcpy(d_SampleXVel		, h_SampleXVel		, sizeof(float)*CellNum		, cudaMemcpyHostToDevice) ;
		cudaMemcpy(d_SampleYVel		, h_SampleYVel		, sizeof(float)*CellNum		, cudaMemcpyHostToDevice) ;
		cudaMemcpy(d_SampleZVel		, h_SampleZVel		, sizeof(float)*CellNum		, cudaMemcpyHostToDevice) ;
		cudaMemcpy(d_SampleXVelSq	, h_SampleXVelSq	, sizeof(float)*CellNum		, cudaMemcpyHostToDevice) ;
		cudaMemcpy(d_SampleYVelSq	, h_SampleYVelSq	, sizeof(float)*CellNum		, cudaMemcpyHostToDevice) ;
		cudaMemcpy(d_SampleZVelSq	, h_SampleZVelSq	, sizeof(float)*CellNum		, cudaMemcpyHostToDevice) ;
		
		cudaMemcpy(d_InletCellNo 	, h_InletCellNo 	, sizeof(int)*InletFaceNum	, cudaMemcpyHostToDevice) ;
		cudaMemcpy(d_InletEdgeNo 	, h_InletEdgeNo 	, sizeof(int)*InletFaceNum	, cudaMemcpyHostToDevice) ;
		cudaMemcpy(d_InletNum 		, h_InletNum 		, sizeof(float)*InletFaceNum	, cudaMemcpyHostToDevice) ;
		cudaMemcpy(d_InletNumR 		, h_InletNumR 		, sizeof(float)*InletFaceNum	, cudaMemcpyHostToDevice) ;
		cudaMemcpy(d_InletXVel 		, h_InletXVel 		, sizeof(float)*InletFaceNum	, cudaMemcpyHostToDevice) ;
		cudaMemcpy(d_InletYVel 		, h_InletYVel 		, sizeof(float)*InletFaceNum	, cudaMemcpyHostToDevice) ;
		cudaMemcpy(d_InletZVel 		, h_InletZVel 		, sizeof(float)*InletFaceNum	, cudaMemcpyHostToDevice) ;
		cudaMemcpy(d_InletNumDen 	, h_InletNumDen 	, sizeof(float)*InletFaceNum	, cudaMemcpyHostToDevice) ;
		cudaMemcpy(d_InletTemp 		, h_InletTemp 		, sizeof(float)*InletFaceNum	, cudaMemcpyHostToDevice) ;
		cudaMemcpy(d_InletArea		, h_InletArea 		, sizeof(float)*InletFaceNum	, cudaMemcpyHostToDevice) ; 

	}else if ( direction == 1 ){
		cudaMemcpy(h_Domain		, d_Domain		, sizeof(DSMC_DOMAIN)		, cudaMemcpyDeviceToHost) ;
		cudaMemcpy(h_Node		, d_Node		, sizeof(DSMC_NODE)*NodeNum	, cudaMemcpyDeviceToHost) ;
		cudaMemcpy(h_Cell		, d_Cell		, sizeof(DSMC_CELL)*CellNum	, cudaMemcpyDeviceToHost) ;
		cudaMemcpy(h_Block		, d_Block		, sizeof(DSMC_BLOCK)*BlockNum	, cudaMemcpyDeviceToHost) ;

		cudaMemcpy(h_IndexParticle	, d_IndexParticle	, sizeof(int)*MaxParticleNum	, cudaMemcpyDeviceToHost) ;
		cudaMemcpy(h_IndexCell1		, d_IndexCell1		, sizeof(int)*CellNum		, cudaMemcpyDeviceToHost) ;
		cudaMemcpy(h_IndexCell2		, d_IndexCell2		, sizeof(int)*CellNum		, cudaMemcpyDeviceToHost) ;
	
		cudaMemcpy(h_ParticleXCoord	, d_ParticleXCoord	, sizeof(float)*MaxParticleNum	, cudaMemcpyDeviceToHost) ;
		cudaMemcpy(h_ParticleYCoord	, d_ParticleYCoord	, sizeof(float)*MaxParticleNum	, cudaMemcpyDeviceToHost) ;
		cudaMemcpy(h_ParticleXVel	, d_ParticleXVel	, sizeof(float)*MaxParticleNum	, cudaMemcpyDeviceToHost) ;
		cudaMemcpy(h_ParticleYVel	, d_ParticleYVel	, sizeof(float)*MaxParticleNum	, cudaMemcpyDeviceToHost) ;
		cudaMemcpy(h_ParticleZVel	, d_ParticleZVel	, sizeof(float)*MaxParticleNum	, cudaMemcpyDeviceToHost) ;
		cudaMemcpy(h_ParticleCellNo	, d_ParticleCellNo	, sizeof(int)*MaxParticleNum	, cudaMemcpyDeviceToHost) ;
		cudaMemcpy(h_ParticleLastCollide, d_ParticleLastCollide	, sizeof(int)*MaxParticleNum	, cudaMemcpyDeviceToHost) ;
		cudaMemcpy(h_ParticleIn		, d_ParticleIn		, sizeof(int)*MaxParticleNum	, cudaMemcpyDeviceToHost) ;

		cudaMemcpy(h_SampleParticleNum	, d_SampleParticleNum	, sizeof(int)*CellNum		, cudaMemcpyDeviceToHost) ;
		cudaMemcpy(h_SampleXVel		, d_SampleXVel		, sizeof(float)*CellNum		, cudaMemcpyDeviceToHost) ;
		cudaMemcpy(h_SampleYVel		, d_SampleYVel		, sizeof(float)*CellNum		, cudaMemcpyDeviceToHost) ;
		cudaMemcpy(h_SampleZVel		, d_SampleZVel		, sizeof(float)*CellNum		, cudaMemcpyDeviceToHost) ;
		cudaMemcpy(h_SampleXVelSq	, d_SampleXVelSq	, sizeof(float)*CellNum		, cudaMemcpyDeviceToHost) ;
		cudaMemcpy(h_SampleYVelSq	, d_SampleYVelSq	, sizeof(float)*CellNum		, cudaMemcpyDeviceToHost) ;
		cudaMemcpy(h_SampleZVelSq	, d_SampleZVelSq	, sizeof(float)*CellNum		, cudaMemcpyDeviceToHost) ;

		cudaMemcpy(h_InletCellNo 	, d_InletCellNo 	, sizeof(int)*InletFaceNum	, cudaMemcpyDeviceToHost) ;
		cudaMemcpy(h_InletEdgeNo 	, d_InletEdgeNo 	, sizeof(int)*InletFaceNum	, cudaMemcpyDeviceToHost) ;
		cudaMemcpy(h_InletNum 		, d_InletNum 		, sizeof(float)*InletFaceNum	, cudaMemcpyDeviceToHost) ;
		cudaMemcpy(h_InletNumR 		, d_InletNumR 		, sizeof(float)*InletFaceNum	, cudaMemcpyDeviceToHost) ;
		cudaMemcpy(h_InletXVel 		, d_InletXVel 		, sizeof(float)*InletFaceNum	, cudaMemcpyDeviceToHost) ;
		cudaMemcpy(h_InletYVel 		, d_InletYVel 		, sizeof(float)*InletFaceNum	, cudaMemcpyDeviceToHost) ;
		cudaMemcpy(h_InletZVel 		, d_InletZVel 		, sizeof(float)*InletFaceNum	, cudaMemcpyDeviceToHost) ;
		cudaMemcpy(h_InletNumDen 	, d_InletNumDen 	, sizeof(float)*InletFaceNum	, cudaMemcpyDeviceToHost) ;
		cudaMemcpy(h_InletTemp 		, d_InletTemp 		, sizeof(float)*InletFaceNum	, cudaMemcpyDeviceToHost) ;
		cudaMemcpy(h_InletArea		, d_InletArea 		, sizeof(float)*InletFaceNum	, cudaMemcpyDeviceToHost) ;
	}
}

//------------------------------------------------------------------------------------------------------------

// Transfer data from device to host for sampling data.
void MemoryCopy(int *d_SampleParticleNum	, int *h_SampleParticleNum , 
		float *d_SampleXVel		, float *h_SampleXVel , 
		float *d_SampleYVel		, float *h_SampleYVel , 
		float *d_SampleZVel		, float *h_SampleZVel , 
		float *d_SampleXVelSq		, float *h_SampleXVelSq , 
		float *d_SampleYVelSq		, float *h_SampleYVelSq , 
		float *d_SampleZVelSq		, float *h_SampleZVelSq , 
		DSMC_DOMAIN *d_Domain		, DSMC_DOMAIN *h_Domain ){

		int	CellNum = h_Domain->CellNum ;

		cudaMemcpy(h_Domain		, d_Domain		, sizeof(DSMC_DOMAIN)		, cudaMemcpyDeviceToHost) ;

		cudaMemcpy(h_SampleParticleNum	, d_SampleParticleNum	, sizeof(int)*CellNum		, cudaMemcpyDeviceToHost) ;
		cudaMemcpy(h_SampleXVel		, d_SampleXVel		, sizeof(float)*CellNum		, cudaMemcpyDeviceToHost) ;
		cudaMemcpy(h_SampleYVel		, d_SampleYVel		, sizeof(float)*CellNum		, cudaMemcpyDeviceToHost) ;
		cudaMemcpy(h_SampleZVel		, d_SampleZVel		, sizeof(float)*CellNum		, cudaMemcpyDeviceToHost) ;
		cudaMemcpy(h_SampleXVelSq	, d_SampleXVelSq	, sizeof(float)*CellNum		, cudaMemcpyDeviceToHost) ;
		cudaMemcpy(h_SampleYVelSq	, d_SampleYVelSq	, sizeof(float)*CellNum		, cudaMemcpyDeviceToHost) ;
		cudaMemcpy(h_SampleZVelSq	, d_SampleZVelSq	, sizeof(float)*CellNum		, cudaMemcpyDeviceToHost) ;
}

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

// Random number generator.
__device__ float Randn( RANDOM *Rand , int index_n ){
	int             i , j ;
	float           rf ;
	const int       MBIG    =       1000000000 ;
	const int       MSEED   =       161803398 ;
        const int       MZ      =       0 ;
        const float     FAC     =       1.E-9 ;
        int             MJ , MK , II ;

        if ( index_n < 0 || Rand->iff == 0 ){
                Rand->iff	=       1 ;
                MJ		=       MSEED - index_n ;
                MJ		=       MJ%MBIG ;
                Rand->ma[55]	=       MJ ;
                MK		=       1 ;

                for ( i=1 ; i<55 ; i++ ){
                        II		=       21*i%55 ;
                        Rand->ma[II]	=       MK ;
                        MK		=       MJ-MK ;

                        if ( MK < MZ ) MK = MK+MBIG ;
                        MJ	=	Rand->ma[II] ;
                }

                for ( j=0 ; j<4 ; j++ ){
                        for ( i=1 ; i<56 ; i++ ){
                                Rand->ma[i]	=	Rand->ma[i] - Rand->ma[1+(i+30)%55] ;
                                if ( Rand->ma[i] < MZ ) Rand->ma[i] = Rand->ma[i]+MBIG ;
                        }
                }

                Rand->inext   =       0 ;
		Rand->inextp  =       31 ;
        }


	do{
		Rand->inext   =       Rand->inext+1 ;
		if ( Rand->inext == 56 ) Rand->inext = 1 ;

                Rand->inextp  =       Rand->inextp+1 ;

                if ( Rand->inextp == 56 ) Rand->inextp = 1 ;

                MJ      =       Rand->ma[Rand->inext] - Rand->ma[Rand->inextp] ;
                if ( MJ < MZ ) MJ = MJ+MBIG ;

                Rand->ma[Rand->inext] = MJ ;

                rf = MJ*FAC ;
        }while ( !( rf > 1.E-8 && rf < 0.99999999 ) ) ;


        return rf ;
}

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

__device__ void RandVel( float *u , float *v , float most_probable_speed , RANDOM *_Rand , int idx ){
	float		a ;
	float		b ;

	(*u)	=	0. ;
	(*v)	=	0. ;
	
	a	=	sqrtf(-logf( Randn(_Rand,idx) )) ;
	b	=	6.283185308 * Randn(_Rand,idx) ;

	(*u)	=	a*sinf(b)*most_probable_speed ;
	(*v)	=	a*cosf(b)*most_probable_speed ;
}

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

__device__ float GAM( float x ){
	float          aa , y , gam ;

	aa      =       1. ;
	y       =       x ;

	if ( y < 1. ){
		aa      =       aa/y ;
	}else{
		y       =       y-1. ;
		while( y>= 1. ){
			aa      =       aa*y ;
			y       =       y-1. ;
		}
	}

	gam = aa*(1.-0.5748646*y+0.9512363*y*y-0.6998588*y*y*y+0.4245549*y*y*y*y-0.1010678*y*y*y*y*y) ;

	return  gam ;
}