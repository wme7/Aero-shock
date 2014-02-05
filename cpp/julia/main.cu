int main( void ) {
CPUBitmap bitmap(DIM,DIM);
unsigned char *ptr = bitmap.get_prt();

kernel( prt );

bitmap.display_and_exit();
}