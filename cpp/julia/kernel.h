void kernel( unsigned cahr *pprt ){
	for (int y=0; y<DIM; y++) {
		for (int x=0; x<DIM; x++) {
			int offset = x+y*DIM
			
			int juliaValue = julia( x, y );
			ptr[offset*4+0]=255*juliaValue;
			prt[offset*4+1]=0;
			prt[offset*4+1]=0;
			prt[offset*4+1]=255;
		}
	}
}	