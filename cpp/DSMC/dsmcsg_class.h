#include <sys/time.h>

#if !defined(__DSMCSG_CLASS_H)
#define __DSMCSG_CLASS_H


// Simulation conditions.
class DSMC_DOMAIN{
	public:
		int		TimestepNum , SampleNo ;
		int		ParticleNum , MaxParticleNum ;

		int		NodeNum , CellNum , XNodeNum , YNodeNum , XCellNum , YCellNum ; 
		int		BlockNum , BCBlockNum , InternalBlockNum , InletFaceNum ;

		int		SampleNum ;
		float		XCellSize , YCellSize , CellVolume , Timestep , ParticleWeight , WeightRatio ;
		float		XVel , YVel , ZVel , Temp , Density , NumDensity ;
		float		XL , XH , YL , YH ;


		DSMC_DOMAIN() ;
		
		void		Dump() ;
} ;

//------------------------------------------------------------------------------------------------------------

// Coordinate of node.
class DSMC_NODE{
	public:
		float		XCoord , YCoord ;

		DSMC_NODE() ;
} ;

//------------------------------------------------------------------------------------------------------------

// Cell data.
class DSMC_CELL{
	public:
		int		Type , Node[4] , Neighbor[4] ;        // Type: 0) general cell, 1) needed to tracking for particle, 2) solid cell.
		float		XCenter , YCenter , MaxCrossSectionSpeed , RemainderCollisionPair ;
		float		EdgeFA[4] , EdgeFB[4] , EdgeFC[4] ;

		DSMC_CELL() ;
} ;

//------------------------------------------------------------------------------------------------------------

// Block data for setup boundary conditions and solid cells inside domain.
class DSMC_BLOCK{
	public:
		int		Type ;
		float		XVel , YVel , NumDensity , Temp ;
		float		XL , XH , YL , YH ;

		int		NodeNum ;
		DSMC_NODE	Node[4] ;

		DSMC_BLOCK() ;
} ;

//------------------------------------------------------------------------------------------------------------

// Simulation results.
class DSMC_OUTPUT{
	public:
		float		NumDensity , Density , XVel , YVel , ZVel , Temp ;
		float		AveParticleNum ;

		DSMC_OUTPUT() ;	
} ;

//------------------------------------------------------------------------------------------------------------

// Simulation time.
class DSMC_TIME{
	private:
		struct timeval	start , end ;
		long double	time , tstart ;

	public:
		long double	Total , Init , Move , Index , Collision , Sample ;
		
		DSMC_TIME() ;

		void		Start() ;
		void		End() ;
		void		Time() ;
		void		Time( long double* ) ;
} ;

//------------------------------------------------------------------------------------------------------------

// Random number generator.
class RANDOM{
	public:
        	int		inext , inextp , iff , ma[56] ;
		
		RANDOM() ;
} ;

//------------------------------------------------------------------------------------------------------------

float randn() ;
void RandVel( float *u , float *v , float most_probable_speed ) ;


#endif
