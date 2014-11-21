#include "dsmcsg_class.h"

#if !defined(__DSMCSG_CUDAFUNCTION_H)
#define __DSMCSG_CUDAFUNCTION_H


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
			RANDOM *d_Rand ) ;


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
					RANDOM *d_Rand ) ;


// Update location of particle over a time.
__device__ void Move( float *XCoord , float *YCoord , float *XVel , float *YVel , float time ) ;


// Calculte cell no. of particle.
__device__ int GetCellNo( float XCoord , float YCoord , float XLower , float YLower , float XCellSize , float YCellSize , int XCellNum , int YCellNum ) ;


__device__ bool InCell( float XCoord , float YCoord , DSMC_NODE *d_Node , DSMC_CELL *_Cell ) ;


// Calculate time and edge number for particle collide edge of cell.
__device__ void GetOutEdge(	float BeforeXCoord , 
				float BeforeYCoord , 
				float XCoord , 
				float YCoord ,
				int FromCellNo ,
				int *EdgeNo , 
				float *MoveTime , 
				float RMoveTime , 
				DSMC_CELL *_Cell ) ;


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
					int idx ) ;


// Remove particle.
template <class TYPE>
void ParticleRemove( int *d_ParticleIn , int *d_ParticleWrite ,	TYPE *d_ParticleProperty , float *d_BufferParticleProperty , int ParticleNum , int RParticleNum ) ;


// Move particle to buffer.
template <class TYPE>
__global__ void Kernel_ParticleRemove1( int *d_ParticleIn , int *d_ParticleWrite , TYPE *d_ParticleProperty , float *d_BufferParticleProperty , int ParticleNum ) ;


// Move particle to global memory.
template <class TYPE>
__global__ void Kernel_ParticleRemove2(	TYPE *d_ParticleProperty , float *d_BufferParticleProperty , int ParticleNum ) ;


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
			RANDOM *d_Rand ) ;


// Calculate number of enter particle.
__global__ void Kernel_CalEnterNum( int *d_EnterNum , float *d_InletNum , float *d_InletNumR , int InletFaceNum ) ;


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
					RANDOM *d_Rand ) ;


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
					RANDOM *d_Rand ) ;


//------------------------------------------------------------------------------------------------------------


// Index.
void Index( int *d_IndexParticle , int *d_IndexCell1 , int *d_IndexCell2 , int *d_ParticleCellNo , int ParticleNum , int CellNum , DSMC_TIME *Time ) ;


// Count number of particle in each cell.
__global__ void Kernel_Index1( int *d_IndexCell2 , int *d_ParticleCellNo , int ParticleNum ) ;


// Create relation between cells and particle.
__global__ void Kernel_Index2( int *d_IndexParticle , int *d_IndexCell1 , int *d_IndexCell2 , int *d_ParticleCellNo , int ParticleNum , int CellNum ) ;


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
		RANDOM *d_Rand ) ;


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
					RANDOM *d_Rand ) ;


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
				int idx ) ;


// Calculate particles properties after collision.
__device__ void Elastic(float *_1ParticleXVel ,
			float *_1ParticleYVel ,
			float *_1ParticleZVel ,
			float *_2ParticleXVel ,
			float *_2ParticleYVel ,
			float *_2ParticleZVel ,
			float *RelativeVel ,
			RANDOM *_Rand ,
			int idx ) ;


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
		DSMC_DOMAIN *d_Domain ) ;


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
				DSMC_DOMAIN *d_Domain ) ;


// Calculate macroscopic properties.
void CalculateResult( int *h_SampleParticleNum , float *h_SampleXVel , float *h_SampleYVel , float *h_SampleZVel , float *h_SampleXVelSq , 
		      float *h_SampleYVelSq , float *h_SampleZVelSq , DSMC_OUTPUT *h_Result , DSMC_DOMAIN *h_Domain , DSMC_CELL *h_Cell ) ;


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
		DSMC_BLOCK *d_Block		, DSMC_BLOCK *h_Block ) ;


// Transfer data from device to host for sampling data.
void MemoryCopy(int *d_SampleParticleNum	, int *h_SampleParticleNum , 
		float *d_SampleXVel		, float *h_SampleXVel , 
		float *d_SampleYVel		, float *h_SampleYVel , 
		float *d_SampleZVel		, float *h_SampleZVel , 
		float *d_SampleXVelSq		, float *h_SampleXVelSq , 
		float *d_SampleYVelSq		, float *h_SampleYVelSq , 
		float *d_SampleZVelSq		, float *h_SampleZVelSq , 
		DSMC_DOMAIN *d_Domain		, DSMC_DOMAIN *h_Domain ) ;


//------------------------------------------------------------------------------------------------------------


// Random number generator.
__device__ float Randn( RANDOM *Rand , int index_n ) ;


__device__ void RandVel( float *u , float *v , float most_probable_speed , RANDOM *_Rand , int idx ) ;


__device__ float GAM( float x ) ;


extern "C"
void preallocBlockSums(unsigned int maxNumElements);


extern "C"
void prescanArray(int *outArray, int *inArray, int numElements);


extern "C"
void deallocBlockSums();


#endif
