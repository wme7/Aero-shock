__global__ void kernel( unsigned char *prt ){
	//map from threadIdx/BlockIdx to pixel position
	int x = blockIdx.x;
	int y = blockIdx.y;
	int offset = x + y * gridDIM.x;
			
	//now calculate the value at that position
	int juliaValue = julia( x, y );
	ptr[offset*4+0]=255*juliaValue;
	prt[offset*4+1]=0;
	prt[offset*4+1]=0;
	prt[offset*4+1]=255;
}	