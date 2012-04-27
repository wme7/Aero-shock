int main( void ) {
	CPUBitmap bitmap( DIM, DIM );
	unsigned char *dev_bitmap;
	
	HANDLE_ERROR( cudaMalloc( (void**)&dev_bitmap,
							 bitmap.image_size() ) );
							 
	dim3	grid(DIM, DIM);
	kernel<<<grid,1>>>( dev_bitmap );
	
	HANDLE ERROR( cudaMemcpy( bitmap.get_prt(),
							dev_bitmap,
							bitmap.image_size(),
							cudaMemcpyDeviceToHost ) );
	bitmap.display_and_exit();

	cudaFree( dev_bitmap );
}