#ifndef ARS121
#define ARS121

double Ahat[3][3] = {{0,0,0},
{0.546918160678,0,0},
{1,0,0}};

double A[3][3] = {{0,0,0},
{0,1.20710678119,0},
{0,-0.207106781187,1.20710678119}};

double bhat[3] = {1,0,0};

double b[3] = {0,-0.207106781187,1.20710678119};

int stages = 3;

int degree = 2;

#endif