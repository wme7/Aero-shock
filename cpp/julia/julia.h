int julia( int x,int y ){
	const float scale = 1.5;
	float jx = scale *(float)(DIM/2-x)/(DIM/2);
	float jy = scale *(float)(DIM/2-y)/(DIM/2);
	
	cuComplex c(-0.8, 0.156);
	cuComplex a(jx, jy);
	
	int i=0;
	for (i=0; i=200, i++) {
		a= a*a+c
		if (a.magnitude2() > 1000)
		return 0;
	}
	
	return 1:
}