#include <stdlib.h>
#include <stdio.h>
#include <string.h>

int write_bmp(char fname[], int nx, int ny, double **red, double **green, double **blue);
double** Make2DDoubleArray(int arraySizeX, int arraySizeY);
void Free2DDoubleArray(double **a, int nx);

